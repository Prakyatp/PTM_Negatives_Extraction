# Phosphorylation Extraction Pipelines

Two datasets re-extracted from the DeepMVP phosphorylation subset with stricter tryptic filtering and explicit negative extraction.

## Datasets

| Folder | Dataset ID | Name | Repository |
|--------|-----------|------|------------|
| `PXD012174_human_phosphoproteome/` | PXD012174 | Human Phosphoproteome Map | PRIDE |
| `MSV000078509_serial_enrichment/` | MSV000078509 | Serial Enrichment – Supp Tables 1 & 2 | MassIVE |

## Common pipeline philosophy

- Every S/T/Y residue is labelled Positive (observed phosphorylated) or Negative (observed unmodified).
- MaxQuant output files (`evidence.txt` / XLSX) are the entry point.
- FASTA-based 4-check verification per site: AA match, FASTA position, start_pos agreement, tryptic filter.
- Both datasets produce an **unfiltered** and a **retention-time filtered** variant.
- Final outputs are split by amino acid (S / T / Y).
