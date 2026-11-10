#!/usr/bin/env python3
"""
Conservation-guided miniaturisation of CasRx.

Builds a ladder of shortened CasRx variants by deleting regions the Cas13d family
treats as dispensable, then validates them by phylogenetic placement and
Needleman-Wunsch identity. No residues are generated: every variant is a strict
subsequence of wild-type CasRx.

Pipeline
  1. load 203 Cas13d sequences + the CasRx template
  2. cluster at 90% 5-mer containment -> per-sequence weights (deflates the
     dominant near-identical clique: 203 sequences are ~108 independent ones)
  3. MAFFT alignment, CasRx first
  4. per column: weighted gap frequency, weighted mean pairwise BLOSUM62
  5. dispensability = 0.6*gap + 0.4*(1 - chem_conservation)
  6. smooth, find runs, protect the catalytic sites, delete -> D1..D9
  7. FastTree placement of the variants against the natural distance band
  8. Needleman-Wunsch identity of each variant to all 203 naturals

Requires: mafft, FastTree on PATH; biopython, numpy, pandas.
  apt-get install -y mafft fasttree
  pip install biopython numpy pandas

Usage
  python casrx_deletion_pipeline.py \
      --dataset data/master_list_cas13d/master_list_cas13d_dataset.csv \
      --template data/evaluation_data/target_data.fasta \
      --template-name CasRx \
      --outdir results/v2
"""

import argparse
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Align import PairwiseAligner, substitution_matrices

AA20 = "ACDEFGHIKLMNPQRSTVWY"
AASET = set(AA20)
BLOSUM62 = substitution_matrices.load("BLOSUM62")

# ----------------------------------------------------------------------------
# defaults -- every one of these is a choice, not a derived constant
# ----------------------------------------------------------------------------
CLUSTER_K = 5            # k-mer size for containment clustering
CLUSTER_THR = 0.90       # containment above which two sequences share a cluster
W_GAP = 0.6              # weight on gap frequency in the dispensability score
W_CHEM = 0.4             # weight on chemical promiscuity
SMOOTH = 7               # smoothing window, residues
MIN_RUN = 6              # shortest deletable run, residues
PROTECT_UP = 15          # catalytic buffer upstream of each R-X4-H
PROTECT_DOWN = 20        # ...and downstream
THRESHOLDS = [round(x, 2) for x in np.arange(0.75, 0.34, -0.05)]
HEPN_RE = r"(?=(R.{4}H))"


# ----------------------------------------------------------------------------
# io
# ----------------------------------------------------------------------------
def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def write_fasta(path, seqs, width=60):
    with open(path, "w") as fh:
        for name, s in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(s), width):
                fh.write(s[i:i + width] + "\n")


def run(cmd, outfile=None):
    print(f"  $ {' '.join(cmd)}", file=sys.stderr)
    with open(outfile, "w") if outfile else subprocess.DEVNULL as out:
        p = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if p.returncode:
        raise RuntimeError(f"{cmd[0]} failed:\n{p.stderr[-2000:]}")


# ----------------------------------------------------------------------------
# step 2: redundancy weights
# ----------------------------------------------------------------------------
def cluster_weights(seqs, k=CLUSTER_K, thr=CLUSTER_THR):
    """Greedy clustering on k-mer containment. Returns {name: 1/cluster_size}.

    The corpus is dominated by near-identical strain variants; without this every
    column statistic measures the largest clique rather than the family.
    """
    names = list(seqs)
    kmers = {n: {seqs[n][i:i + k] for i in range(len(seqs[n]) - k + 1)} for n in names}
    reps, assign = [], {}
    for n in sorted(names, key=lambda x: -len(seqs[x])):     # longest first
        for r in reps:
            denom = min(len(kmers[n]), len(kmers[r]))
            if denom and len(kmers[n] & kmers[r]) / denom >= thr:
                assign[n] = r
                break
        else:
            reps.append(n)
            assign[n] = n
    size = Counter(assign.values())
    return {n: 1.0 / size[assign[n]] for n in names}, len(reps)


# ----------------------------------------------------------------------------
# steps 4-5: per-column scoring
# ----------------------------------------------------------------------------
def column_scores(aln, weights):
    """Weighted gap frequency and normalised BLOSUM62 conservation per column.

    gap[j]  - weighted fraction of sequences with no residue here.
              High => the family already lives without this position.
    chem[j] - weighted mean pairwise BLOSUM62 among observed residues, min-max
              normalised. High => the column varies only among chemically similar
              residues (D/E, I/L/V); low => anything goes.
    """
    names = list(aln)
    L = len(aln[names[0]])
    M = np.array([list(aln[n]) for n in names])
    w = np.array([weights[n] for n in names], dtype=float)
    w /= w.sum()

    gap = np.zeros(L)
    chem = np.zeros(L)
    for j in range(L):
        col = M[:, j]
        ok = np.fromiter((c in AASET for c in col), bool, len(col))
        gap[j] = w[~ok].sum()
        if ok.sum() < 2:
            chem[j] = -4.0                       # effectively unconstrained
            continue
        ww = w[ok] / w[ok].sum()
        freq = {}
        for a, x in zip(col[ok], ww):
            freq[a] = freq.get(a, 0.0) + x
        chem[j] = sum(pa * pb * BLOSUM62[a][b]
                      for a, pa in freq.items() for b, pb in freq.items())

    chem_n = (chem - chem.min()) / (chem.max() - chem.min())
    disp = W_GAP * gap + W_CHEM * (1.0 - chem_n)
    return gap, chem_n, disp


def project_to_template(aln, template_name, *arrays):
    """Map column-indexed arrays onto the template's own residue indices."""
    row = aln[template_name]
    cols = [j for j, ch in enumerate(row) if ch != "-"]
    return [a[cols] for a in arrays], cols


# ----------------------------------------------------------------------------
# step 6: deletion
# ----------------------------------------------------------------------------
def protected_positions(seq, up=PROTECT_UP, down=PROTECT_DOWN):
    """Residue indices around every R-X4-H that must never be deleted."""
    prot = set()
    for m in re.finditer(HEPN_RE, seq):
        for i in range(max(0, m.start() - up), min(len(seq), m.start() + 6 + down)):
            prot.add(i)
    return prot


def deletable_runs(disp, thr, protected, smooth=SMOOTH, min_run=MIN_RUN):
    sm = np.convolve(disp, np.ones(smooth) / smooth, mode="same")
    mask = sm >= thr
    runs, start = [], None
    for i, v in enumerate(mask):
        if v and start is None:
            start = i
        elif not v and start is not None:
            if i - start >= min_run:
                runs.append((start, i))
            start = None
    if start is not None and len(mask) - start >= min_run:
        runs.append((start, len(mask)))
    return [r for r in runs if not any(p in protected for p in range(*r))]


def build_ladder(template, disp, gap, chem_n, thresholds=THRESHOLDS):
    prot = protected_positions(template)
    variants, spec, segments = {}, [], []
    for i, thr in enumerate(thresholds, 1):
        runs = deletable_runs(disp, thr, prot)
        cut = {p for a, b in runs for p in range(a, b)}
        seq = "".join(c for j, c in enumerate(template) if j not in cut)
        key = f"D{i}"
        variants[key] = seq
        spec.append(dict(
            key=key, thr=thr, length=len(seq), cut=len(cut),
            pct_cut=round(len(cut) / len(template), 3), n_segments=len(runs),
            n_hepn=len(re.findall(HEPN_RE, seq)),
            mean_gap_support=round(float(np.mean([gap[p] for p in cut])), 3) if cut else 0.0,
            mean_chem_cons=round(float(np.mean([chem_n[p] for p in cut])), 3) if cut else 0.0,
        ))
        for a, b in runs:
            segments.append(dict(
                key=key, thr=thr, start=a, end=b, length=b - a,
                rel_start=round(a / len(template), 3),
                mean_gap=round(float(gap[a:b].mean()), 3),
                mean_chem_cons=round(float(chem_n[a:b].mean()), 3),
                seq=template[a:b]))
    return variants, pd.DataFrame(spec), pd.DataFrame(segments)


# ----------------------------------------------------------------------------
# step 7: phylogenetic placement
# ----------------------------------------------------------------------------
def place_variants(aln_path, variants, template, template_name, outdir, n_sample=40):
    """Add variants to the alignment, build a tree, compare patristic distances.

    The test is a necessary condition, not a discriminator: deletions target
    columns the family already varies in, so the conserved columns carrying the
    phylogenetic signal are untouched by construction. It catches a cut into
    conserved material; it says nothing about folding.
    """
    add_path = outdir / "variants_for_tree.fasta"
    write_fasta(add_path, {f"{template_name}_wt": template,
                           **{f"var_{k}": v for k, v in variants.items()}})
    merged = outdir / "ladder_aln.fa"
    run(["mafft", "--add", str(add_path), "--keeplength", "--thread", "4", str(aln_path)],
        outfile=merged)
    tre = outdir / "ladder.tre"
    run(["FastTree", "-quiet", "-nosupport", str(merged)], outfile=tre)

    tree = Phylo.read(tre, "newick")
    tips = {c.name: c for c in tree.get_terminals() if c.name}
    tree.root_at_midpoint()
    var_tips = {k: f"var_{k}" for k in variants if f"var_{k}" in tips}
    nat = [n for n in tips if not n.startswith("var_") and not n.endswith("_wt")]

    rng = np.random.default_rng(0)
    sample = [nat[i] for i in rng.choice(len(nat), min(n_sample, len(nat)), replace=False)]
    nn = np.array([min(tree.distance(tips[a], tips[b]) for b in nat if b != a)
                   for a in sample])
    band95 = float(np.percentile(nn, 95))

    wt = f"{template_name}_wt"
    rows = [dict(key="wild-type", d_nearest_natural=round(
        min(tree.distance(tips[wt], tips[b]) for b in nat), 4), d_to_wt=0.0, inside=True)]
    for k, tip in var_tips.items():
        dn = min(tree.distance(tips[tip], tips[b]) for b in nat)
        rows.append(dict(key=k, d_nearest_natural=round(float(dn), 4),
                         d_to_wt=round(float(tree.distance(tips[tip], tips[wt])), 4),
                         inside=bool(dn <= band95)))
    return pd.DataFrame(rows), dict(band_median=float(np.median(nn)),
                                    band_95th=band95, band_max=float(nn.max()))


# ----------------------------------------------------------------------------
# step 8: Needleman-Wunsch identity
# ----------------------------------------------------------------------------
def make_aligner(free_ends=True):
    a = PairwiseAligner()
    a.mode = "global"                      # Needleman-Wunsch
    a.substitution_matrix = BLOSUM62
    a.open_gap_score = -11
    a.extend_gap_score = -1
    if free_ends:                          # don't penalise terminal truncation
        try:
            a.end_insertion_score = 0.0    # biopython >= 1.85
            a.end_deletion_score = 0.0
        except AttributeError:
            a.target_end_gap_score = 0.0
            a.query_end_gap_score = 0.0
    return a


def nw_identity(aligner, a, b):
    """Three conventions, because the denominator matters when lengths differ."""
    x = aligner.align(a, b)[0]
    s, t = x[0], x[1]
    matches = sum(p == q for p, q in zip(s, t) if p != "-" and q != "-")
    aligned = sum(1 for p, q in zip(s, t) if p != "-" and q != "-")
    return dict(score=float(x.score), aln_len=len(s), matches=matches,
                id_aligned=matches / aligned if aligned else 0.0,
                id_alnlen=matches / len(s),
                id_shorter=matches / min(len(a), len(b)))


def identity_table(variants, naturals, template, template_name):
    al = make_aligner()
    rows = []
    for key, v in variants.items():
        r = nw_identity(al, v, template)
        rows.append(dict(key=key, length=len(v), against=template_name, **r))
        stats = [nw_identity(al, v, s) for s in naturals.values()]
        ids = np.array([s["id_aligned"] for s in stats])
        best = max(zip(ids, naturals), key=lambda t: t[0])
        rows.append(dict(key=key, length=len(v), against="203 naturals",
                         score=np.nan, aln_len=np.nan, matches=np.nan,
                         id_aligned=round(float(np.median(ids)), 3),
                         id_alnlen=np.nan, id_shorter=np.nan))
        rows[-1]["note"] = f"median; best non-self {best[1]} at {best[0]:.3f}"
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", required=True, help="CSV with accession,sequence[,length]")
    ap.add_argument("--template", required=True, help="FASTA holding the template")
    ap.add_argument("--template-name", default="CasRx")
    ap.add_argument("--outdir", default="results/v2")
    ap.add_argument("--min-len", type=int, default=600)
    ap.add_argument("--max-len", type=int, default=1200)
    ap.add_argument("--skip-tree", action="store_true", help="skip FastTree placement")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    # 1. load -----------------------------------------------------------------
    template = read_fasta(args.template)[args.template_name]
    df = pd.read_csv(args.dataset).drop_duplicates("sequence")
    if "length" not in df:
        df["length"] = df.sequence.str.len()
    df = df[(df.length.between(args.min_len, args.max_len)) &
            (~df.sequence.str.contains("X"))]
    naturals = {str(r.accession): r.sequence for r in df.itertuples()
                if r.sequence != template}
    print(f"[1] template {args.template_name} {len(template)} aa | "
          f"{len(naturals) + 1} sequences after filtering")
    hepn = [m.start() for m in re.finditer(HEPN_RE, template)]
    print(f"    catalytic R-X4-H at {hepn} "
          f"(relative {[round(h/len(template), 3) for h in hepn]})")

    # 2. weights --------------------------------------------------------------
    allseq = {args.template_name: template, **naturals}
    weights, n_clusters = cluster_weights(allseq)
    print(f"[2] {n_clusters} clusters at {CLUSTER_THR:.0%} containment | "
          f"effective N {sum(weights.values()):.1f} of {len(allseq)}")

    # 3. align ----------------------------------------------------------------
    raw = out / "input.fasta"
    write_fasta(raw, allseq)
    aln_path = out / "aln.fasta"
    run(["mafft", "--auto", "--thread", "4", str(raw)], outfile=aln_path)
    aln = read_fasta(aln_path)
    print(f"[3] alignment {len(aln)} x {len(aln[args.template_name])} columns")

    # 4-5. score --------------------------------------------------------------
    gap_c, chem_c, disp_c = column_scores(aln, weights)
    (gap, chem_n, disp), _ = project_to_template(
        aln, args.template_name, gap_c, chem_c, disp_c)
    print(f"[4] gap frequency median {np.median(gap):.3f} | "
          f"dispensability median {np.median(disp):.3f}, p90 {np.percentile(disp, 90):.3f}")
    pd.DataFrame(dict(index=np.arange(len(template)), residue=list(template),
                      gap_frac=np.round(gap, 4), chem_cons=np.round(chem_n, 4),
                      dispensability=np.round(disp, 4))
                 ).to_csv(out / "per_residue_scores.csv", index=False)

    # 6-7. ladder -------------------------------------------------------------
    variants, spec, segments = build_ladder(template, disp, gap, chem_n)
    print("[5] ladder:")
    print(spec.to_string(index=False))
    spec.to_csv(out / "scaffold_spec.csv", index=False)
    segments.to_csv(out / "deleted_segments.csv", index=False)
    write_fasta(out / "top_sequences.fasta",
                {f"{r.key}|rank0|thr{r.thr}|len{r.length}": variants[r.key]
                 for r in spec.itertuples()})

    bad = spec[spec.n_hepn < len(hepn)]
    if len(bad):
        print(f"    WARNING: {len(bad)} variant(s) lost a catalytic motif: "
              f"{list(bad.key)}", file=sys.stderr)

    # 8. placement ------------------------------------------------------------
    if not args.skip_tree:
        placement, band = place_variants(aln_path, variants, template,
                                         args.template_name, out)
        print(f"[6] natural nearest-neighbour band: median {band['band_median']:.3f}, "
              f"95th {band['band_95th']:.3f}, max {band['band_max']:.3f}")
        print(placement.to_string(index=False))
        placement.to_csv(out / "phylo_placement.csv", index=False)

    # 9. identity -------------------------------------------------------------
    ident = identity_table(variants, naturals, template, args.template_name)
    print("[7] Needleman-Wunsch identity:")
    print(ident.round(3).to_string(index=False))
    ident.to_csv(out / "identity.csv", index=False)

    print(f"\nwritten to {out}/")


if __name__ == "__main__":
    main()
