# Conservation-guided miniaturisation of CasRx

Code, data and results for a coursework project on shortening the CRISPR-Cas13d
effector CasRx by deleting regions that the Cas13d family itself varies in.

> **This is a module project.** The work has not been peer reviewed, and has not
> been subject to any external scrutiny. Nothing here has been tested in the lab:
> every result is computational, and the designed sequences are candidates for
> testing rather than validated enzymes. Treat the numbers accordingly.

## What this does

203 Cas13d sequences are aligned, weighted to discount near-identical database
entries, and scored per alignment column on two things: how often the family has no
residue there (gap frequency) and how little the family constrains the chemistry of
the residue it does have (BLOSUM62 column conservation). Regions scoring high on the
combination are deleted from CasRx, with both HEPN catalytic motifs protected. Nine
variants result, removing 56 to 311 of 967 residues.

No residues are generated. Every variant is a strict subsequence of wild-type CasRx.

Variants are then checked three ways: phylogenetic placement against the natural
family, Needleman-Wunsch identity, and structure prediction with the crRNA bound,
compared against the experimental CasRx-crRNA complex (PDB 9M31).

## Layout

```
analysis/     the pipeline, as runnable scripts
notebooks/    Colab notebooks for the GPU steps
pymol/        scripts that render the structure figures
data/         input sequences
results/      everything the pipeline produced
```

### `analysis/`

| file | what it does |
| --- | --- |
| `casrx_deletion_pipeline.py` | the main pipeline: cluster, align, score, delete, place, compare |
| `casrx_phylo_pipeline.py` | a tree-based variant using rate4site ML rates and Fitch indel counts instead of column statistics |

```bash
python analysis/casrx_deletion_pipeline.py \
    --dataset data/cas13d_203_sequences.csv \
    --template data/reference_targets.fasta \
    --template-name CasRx \
    --outdir results/rerun
```

Needs `mafft` and `FastTree` on PATH plus `biopython numpy pandas`. Add
`--skip-tree` if FastTree is unavailable. The phylogenetic variant additionally
needs `rate4site` and takes 10-15 minutes on this dataset.

### `notebooks/`

GPU steps, written for Colab.

- `Cas13d_compression_v3_esm2.ipynb` — the motif-scaffold route with ESM-2 infilling.
  Kept for the record. This route did not produce foldable sequences and is not what
  the deletion results come from.
- `Cas13d_boltz_structure_eval.ipynb` — folds variants with Boltz-2, apo and with
  crRNA, and computes TM-score, pLDDT, ipTM and catalytic geometry.

Two things that will waste your time if you miss them. Boltz must run with
`--use_msa_server`: in single-sequence mode wild-type CasRx folds at pLDDT 36 and
every comparison becomes meaningless. And the notebook asserts the wild-type control
clears pLDDT 80 before any design metric is read, which is the check that catches it.

### `pymol/`

```bash
pymol -cq pymol/render_figures.pml
```

Bare PyMOL commands silently rename short object names (`load f, x` creates `x_`),
after which selections match nothing and renders come out blank. These scripts drive
everything through `cmd.*` inside a `python` block for that reason.

### `data/`

| file | contents |
| --- | --- |
| `cas13d_203_sequences.csv` / `.fasta` | the 203 Cas13d sequences, NCBI Entrez free text on "Cas13d", filtered to 600-1200 aa |
| `reference_targets.fasta` | CasRx and three other Cas13d used as references |

### `results/`

| file | contents |
| --- | --- |
| `sequences/casrx_deletion_ladder.fasta` | wild-type CasRx and the nine variants D1-D9 |
| `alignment_203_mafft.fasta` | the MAFFT alignment the scores come from |
| `trees/ladder_fasttree.nwk` | FastTree tree of the 203 naturals plus the variants |
| `tables/deletion_ladder.csv` | one row per variant: length, residues removed, evidence strength |
| `tables/per_residue_scores.csv` | gap frequency, chemical conservation and dispensability for all 967 CasRx positions |
| `tables/structure_metrics_deletion_ladder.csv` | Boltz-2 output: TM, RMSD, pLDDT, ipTM, catalytic separation |
| `structures/*.pdb` | predicted complexes for wild type, D7 and D9 |

## Headline numbers

| | wild type | D7 | D9 |
| --- | --- | --- | --- |
| length (aa) | 967 | 799 | 656 |
| TM to PDB 9M31 | 0.912 | 0.901 | 0.720 |
| pLDDT | 87.1 | 85.2 | 74.5 |
| ipTM with crRNA | 0.849 | 0.754 | 0.689 |
| catalytic Cα R-R (Å) | 14.92 | 14.96 | 14.99 |

D7 removes 17% of the protein and scores within 0.011 TM units of wild-type CasRx
predicted under identical settings.

## Caveats worth reading before using any of this

Every structural number is a prediction. Boltz-2 places the crRNA roughly 58 Å from
its experimental position even for wild-type CasRx, so guide placement is only
interpretable relative to that control, not in absolute terms. The weighting between
the two scoring axes (0.6 / 0.4), the smoothing window, the minimum run length and
the catalytic buffer are all choices rather than derived values. The clustering is a
greedy k-mer method rather than MMseqs2, and cluster membership depends on input
order.

## Acknowledgements

Sequence retrieval uses the scraper from
[NTU-Singapore iGEM 2025](https://github.com/boonyhow/igem-2025).

## Licence

MIT for the code. The sequence data is from NCBI and PDB and carries their terms.
