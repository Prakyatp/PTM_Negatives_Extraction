# Acetylation Extraction Pipelines

Three datasets re-extracted from the DeepMVP acetylation subset with strict terminal-K bias correction (reduceK).

## Datasets

| Folder | Dataset ID | Name | Repository |
|--------|-----------|------|------------|
| `MSV000078509_kac/` | MSV000078509 | Kac from High_Kac + pSTY sheets | MassIVE |
| `PXD005793_SILAC_acetylome/` | PXD005793 | SILAC Acetylome (3 reps + 3 proteome) | PRIDE |
| `MSV000082644_medulloblastoma/` | MSV000082644 | Medulloblastoma (mmc3.xlsx) | MassIVE |

## Common pipeline philosophy

- Target residue: **K** (lysine). Trypsin cleaves after unmodified K, so every tryptic peptide ending in K is systematically under-represented as Kac.
- **reduceK** two-stage balancing: Stage 1 downsamples negatives until K_first% matches positives; Stage 2 downsamples positives until K_second_last% matches.
- Each dataset produces an **unfiltered**, **filtered_unbalanced**, and **filtered_balanced** variant.
- `filtered_unbalanced` is the primary training-ready output.
