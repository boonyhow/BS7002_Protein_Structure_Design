#!/usr/bin/env python3
"""
Phylogenetic miniaturisation of CasRx.

A genuinely tree-based version of the conservation-guided deletion pipeline. The
earlier method scored alignment columns directly (weighted gap frequency plus
BLOSUM62 conservation), which is alignment-based rather than phylogenetic: it
approximates the correction for shared ancestry by down-weighting sequence
clusters. This version uses the tree itself.

Two tree-aware signals replace the two column statistics:

  evolutionary rate   rate4site, maximum-likelihood per-site rate under JTT with a
                      Gamma prior and an empirical-Bayes estimate. Accounts for
                      branch lengths, so a column that differs between two distant
                      lineages is not scored like one that differs between two
                      sister strains. Replaces BLOSUM62 column conservation.

  indel events        Fitch parsimony on the gap/residue character, giving the
                      minimum number of independent insertion or deletion events
                      on the tree. A region absent from 70 members of one clique
                      scores 1; a region lost independently in 8 lineages scores 8.
                      Replaces raw weighted gap frequency.

  dispensability      w_rate * norm(rate) + w_indel * norm(indel events)

Everything downstream (smoothing, run detection, catalytic protection, the
threshold ladder, phylogenetic placement, Needleman-Wunsch identity) is unchanged
from the alignment-based pipeline, so the two are directly comparable.

Requires: mafft, FastTree, rate4site on PATH; biopython, numpy, pandas.
  apt-get install -y mafft fasttree rate4site
  pip install biopython numpy pandas

Usage
  python casrx_phylo_pipeline.py \
      --dataset data/master_list_cas13d/master_list_cas13d_dataset.csv \
      --template data/evaluation_data/target_data.fasta \
      --template-name CasRx --outdir results/phylo

rate4site takes roughly 10-15 minutes on ~108 sequences of ~1900 columns. Its
output is cached in the output directory and reused on re-runs.
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

AASET = set("ACDEFGHIKLMNPQRSTVWY")
BLOSUM62 = substitution_matrices.load("BLOSUM62")

CLUSTER_K = 5
CLUSTER_THR = 0.90       # representatives for rate4site; also the weighting basis
W_RATE = 0.4             # weight on ML evolutionary rate
W_INDEL = 0.6            # weight on independent indel events
SMOOTH = 7
MIN_RUN = 6
PROTECT_UP = 15
PROTECT_DOWN = 20
THRESHOLDS = [round(x, 2) for x in np.arange(0.75, 0.34, -0.05)]
HEPN_RE = r"(?=(R.{4}H))"


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


def run(cmd, outfile=None, stdin_text=None):
    print(f"  $ {' '.join(cmd)}", file=sys.stderr)
    out = open(outfile, "w") if outfile else subprocess.DEVNULL
    try:
        p = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE,
                           input=stdin_text, text=True)
    finally:
        if outfile:
            out.close()
    if p.returncode:
        raise RuntimeError(f"{cmd[0]} failed:\n{p.stderr[-2000:]}")


# ----------------------------------------------------------------------------
def pick_representatives(seqs, template_name, k=CLUSTER_K, thr=CLUSTER_THR):
    """One sequence per 90%-containment cluster, template forced first.

    rate4site reports rates against the first sequence, so the template must lead.
    Using representatives rather than all sequences also keeps the ML optimisation
    tractable and stops one near-identical clique from dominating branch lengths.
    """
    raw = {n: s.replace("-", "") for n, s in seqs.items()}
    kmers = {n: {raw[n][i:i + k] for i in range(len(raw[n]) - k + 1)} for n in raw}
    reps = []
    for n in sorted(seqs, key=lambda x: -len(raw[x])):
        for r in reps:
            denom = min(len(kmers[n]), len(kmers[r]))
            if denom and len(kmers[n] & kmers[r]) / denom >= thr:
                break
        else:
            reps.append(n)
    if template_name in reps:
        reps.remove(template_name)
    return [template_name] + reps


# ----------------------------------------------------------------------------
def fitch_changes(tree_clades, states):
    """Minimum number of state changes on the tree (Fitch parsimony, down-pass).

    This is the phylogenetic correction that matters most here. Counting the
    fraction of sequences carrying a gap conflates one ancient deletion inherited
    by a large clade with many independent losses. Counting events on the tree
    separates them.
    """
    sets, changes = {}, 0
    for c in tree_clades:                      # postorder
        if c.is_terminal():
            sets[c] = states.get(c.name)
        else:
            kids = [sets[k] for k in c.clades if sets[k] is not None]
            if not kids:
                sets[c] = None
                continue
            inter = set.intersection(*kids) if len(kids) > 1 else kids[0]
            if inter:
                sets[c] = inter
            else:
                sets[c] = set.union(*kids)
                changes += 1
    return changes


def indel_events(aln, tree):
    """Independent gap/residue transitions per alignment column."""
    clades = list(tree.find_clades(order="postorder"))
    names = list(aln)
    L = len(aln[names[0]])
    out = np.zeros(L)
    for j in range(L):
        states = {n: ({"gap"} if aln[n][j] == "-" else {"res"}) for n in names}
        out[j] = fitch_changes(clades, states)
    return out


def run_rate4site(aln_path, tree_path, outdir):
    """ML per-site rates. Cached: rate4site is the slow step."""
    res = outdir / "rates.txt"
    if res.exists() and res.stat().st_size > 0:
        print("    (cached)", file=sys.stderr)
    else:
        run(["rate4site", "-s", str(aln_path), "-t", str(tree_path), "-o", str(res)],
            stdin_text="0\n")
    rates = {}
    for line in open(res):
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 3:
            try:
                rates[int(parts[0])] = float(parts[2])
            except ValueError:
                pass
    if not rates:
        raise RuntimeError("rate4site produced no rates")
    return np.array([rates[i] for i in sorted(rates)])


# ----------------------------------------------------------------------------
def norm(x):
    return (x - x.min()) / (x.max() - x.min()) if x.max() > x.min() else np.zeros_like(x)


def protected_positions(seq, up=PROTECT_UP, down=PROTECT_DOWN):
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


def build_ladder(template, disp, rate, indel, thresholds=THRESHOLDS):
    prot = protected_positions(template)
    variants, spec, segments = {}, [], []
    for i, thr in enumerate(thresholds, 1):
        runs = deletable_runs(disp, thr, prot)
        cut = {p for a, b in runs for p in range(a, b)}
        seq = "".join(c for j, c in enumerate(template) if j not in cut)
        key = f"P{i}"
        variants[key] = seq
        spec.append(dict(
            key=key, thr=thr, length=len(seq), cut=len(cut),
            pct_cut=round(len(cut) / len(template), 3), n_segments=len(runs),
            n_hepn=len(re.findall(HEPN_RE, seq)),
            mean_rate=round(float(np.mean([rate[p] for p in cut])), 3) if cut else 0.0,
            mean_indel_events=round(float(np.mean([indel[p] for p in cut])), 2) if cut else 0.0))
        for a, b in runs:
            segments.append(dict(key=key, thr=thr, start=a, end=b, length=b - a,
                                 rel_start=round(a / len(template), 3),
                                 mean_rate=round(float(rate[a:b].mean()), 3),
                                 mean_indel_events=round(float(indel[a:b].mean()), 2),
                                 seq=template[a:b]))
    return variants, pd.DataFrame(spec), pd.DataFrame(segments)


# ----------------------------------------------------------------------------
def make_aligner():
    a = PairwiseAligner()
    a.mode = "global"
    a.substitution_matrix = BLOSUM62
    a.open_gap_score = -11
    a.extend_gap_score = -1
    try:
        a.end_insertion_score = 0.0
        a.end_deletion_score = 0.0
    except AttributeError:
        a.target_end_gap_score = 0.0
        a.query_end_gap_score = 0.0
    return a


def nw_identity(aligner, a, b):
    x = aligner.align(a, b)[0]
    s, t = x[0], x[1]
    matches = sum(p == q for p, q in zip(s, t) if p != "-" and q != "-")
    aligned = sum(1 for p, q in zip(s, t) if p != "-" and q != "-")
    return dict(score=float(x.score), aln_len=len(s),
                id_aligned=matches / aligned if aligned else 0.0,
                id_alnlen=matches / len(s))


def place_variants(aln_path, variants, template, outdir):
    add = outdir / "variants_for_tree.fasta"
    write_fasta(add, {"TPL_wt": template, **{f"var_{k}": v for k, v in variants.items()}})
    merged = outdir / "ladder_aln.fa"
    run(["mafft", "--add", str(add), "--keeplength", "--thread", "4", str(aln_path)],
        outfile=merged)
    tre = outdir / "ladder.tre"
    run(["FastTree", "-quiet", "-nosupport", str(merged)], outfile=tre)
    tree = Phylo.read(tre, "newick")
    tips = {c.name: c for c in tree.get_terminals() if c.name}
    tree.root_at_midpoint()
    nat = [n for n in tips if not n.startswith("var_") and n != "TPL_wt"]
    rng = np.random.default_rng(0)
    sample = [nat[i] for i in rng.choice(len(nat), min(40, len(nat)), replace=False)]
    nn = np.array([min(tree.distance(tips[a], tips[b]) for b in nat if b != a)
                   for a in sample])
    band = float(np.percentile(nn, 95))
    rows = []
    for k in variants:
        tip = f"var_{k}"
        if tip not in tips:
            continue
        dn = min(tree.distance(tips[tip], tips[b]) for b in nat)
        rows.append(dict(key=k, d_nearest_natural=round(float(dn), 4),
                         d_to_wt=round(float(tree.distance(tips[tip], tips["TPL_wt"])), 4),
                         inside=bool(dn <= band)))
    return pd.DataFrame(rows), band


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--template", required=True)
    ap.add_argument("--template-name", default="CasRx")
    ap.add_argument("--outdir", default="results/phylo")
    ap.add_argument("--min-len", type=int, default=600)
    ap.add_argument("--max-len", type=int, default=1200)
    ap.add_argument("--skip-placement", action="store_true")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    template = read_fasta(args.template)[args.template_name]
    df = pd.read_csv(args.dataset).drop_duplicates("sequence")
    if "length" not in df:
        df["length"] = df.sequence.str.len()
    df = df[df.length.between(args.min_len, args.max_len) & ~df.sequence.str.contains("X")]
    allseq = {args.template_name: template,
              **{str(r.accession): r.sequence for r in df.itertuples()
                 if r.sequence != template}}
    hepn = [m.start() for m in re.finditer(HEPN_RE, template)]
    print(f"[1] template {len(template)} aa, R-X4-H at {hepn} | {len(allseq)} sequences")

    reps = pick_representatives(allseq, args.template_name)
    print(f"[2] {len(reps)} cluster representatives at {CLUSTER_THR:.0%} containment")

    sub = out / "reps.fasta"
    write_fasta(sub, {n: allseq[n] for n in reps})
    aln_path = out / "aln.fasta"
    run(["mafft", "--auto", "--thread", "4", str(sub)], outfile=aln_path)
    aln = read_fasta(aln_path)
    tree_path = out / "tree.tre"
    run(["FastTree", "-quiet", "-nosupport", str(aln_path)], outfile=tree_path)
    print(f"[3] alignment {len(aln)} x {len(aln[args.template_name])} columns, tree built")

    print("[4] rate4site (ML per-site rates, slow)...", file=sys.stderr)
    rate = run_rate4site(aln_path, tree_path, out)
    tree = Phylo.read(tree_path, "newick")
    indel_col = indel_events(aln, tree)
    cols = [j for j, ch in enumerate(aln[args.template_name]) if ch != "-"]
    indel = indel_col[cols]
    if len(rate) != len(template):
        raise RuntimeError(f"rate4site returned {len(rate)} positions, "
                           f"template is {len(template)}")
    print(f"[4] rate range {rate.min():.2f}..{rate.max():.2f} | "
          f"indel events median {np.median(indel):.0f}, max {indel.max():.0f} | "
          f"corr(rate, indel) {np.corrcoef(rate, indel)[0, 1]:.3f}")

    disp = W_RATE * norm(rate) + W_INDEL * norm(indel)
    pd.DataFrame(dict(index=np.arange(len(template)), residue=list(template),
                      ml_rate=np.round(rate, 4), indel_events=indel.astype(int),
                      dispensability=np.round(disp, 4))
                 ).to_csv(out / "per_residue_scores.csv", index=False)

    variants, spec, segments = build_ladder(template, disp, rate, indel)
    print("[5] ladder:")
    print(spec.to_string(index=False))
    spec.to_csv(out / "scaffold_spec.csv", index=False)
    segments.to_csv(out / "deleted_segments.csv", index=False)
    write_fasta(out / "top_sequences.fasta",
                {f"{r.key}|rank0|thr{r.thr}|len{r.length}": variants[r.key]
                 for r in spec.itertuples()})
    bad = spec[spec.n_hepn < len(hepn)]
    if len(bad):
        print(f"    WARNING: catalytic motif lost in {list(bad.key)}", file=sys.stderr)

    if not args.skip_placement:
        placement, band = place_variants(aln_path, variants, template, out)
        print(f"[6] natural 95th-percentile band {band:.3f}")
        print(placement.to_string(index=False))
        placement.to_csv(out / "phylo_placement.csv", index=False)

    al = make_aligner()
    rows = [dict(key=k, length=len(v), **nw_identity(al, v, template))
            for k, v in variants.items()]
    ident = pd.DataFrame(rows)
    print("[7] Needleman-Wunsch vs template:")
    print(ident.round(3).to_string(index=False))
    ident.to_csv(out / "identity.csv", index=False)
    print(f"\nwritten to {out}/")


if __name__ == "__main__":
    main()
