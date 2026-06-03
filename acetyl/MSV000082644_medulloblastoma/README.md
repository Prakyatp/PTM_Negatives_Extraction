# MSV000082644 — Medulloblastoma (mmc3.xlsx)

**PTM type:** Acetylation  
**Repository:** MassIVE  
**Dataset:** `MSV000082644`  

## Workflow

### Step 1 — Extract Kac positives (`ac_extract.py`)

Reads Acetylomics-acK sheet from mmc3.xlsx in 4 stages: add_scans counts non-null sample columns; delete_acnu removes multi-accession peptides; label_peptides walks each modified sequence for K with prob 0.9 or higher (Positive) or prob 0.0 (Negative); map_to_fasta uses Aho-Corasick with strict tryptic rule.

**Key stats:**

- 10,988 rows loaded
- Positive: 10,574, Negative: 294
- 9,574 output rows

### Step 2 — Extract Kac negatives (`phospho_extract.py`)

Run on Phospho-pSTY (56,000 rows) and Phospho-pY (2,144 rows) sheets. Same add_scans + delete_acnu pipeline. extract_k_negatives finds the nearest internal K to the first phospho site. Terminal-K rescue is inverted vs ac_extract (K must be unmodified).

**Key stats:**

- pSTY: 17,094 mapped
- pY: 330 mapped
- Total: 17,424 K-negatives

### Step 3 — Merge all 3 CSVs (`merge.py`)

Simple concatenation of ac_mapped.csv (positives + Kac sheet negatives) and both phospho mapped files (pSTY and pY negatives). No deduplication at this stage.

**Key stats:**

- 3 input CSVs
- 26,998 total rows
- 9,497 pos / 17,501 neg

### Step 4 — Clean columns (`6_clean_columns.py`)

Standardises the merged file. Cleans the modified peptide sequence (strips probability annotations like (0.99) and underscores). Computes K_first and K_second_last flags. Renames columns to a clean schema.

**Key stats:**

- Strips (0.99), (1.0)
- Renames protein columns
- Adds K_first / K_second_last

### Step 5 — Dedup + reduceK (`7_reduceK.py`)

Dedup collapses 26,998 rows to 22,357 unique sites (unfiltered). Two-stage K-position balancing: Stage 1 downsamples Negatives until K_first matches Positives. Stage 2 downsamples Positives until K_second_last matches updated Negatives.

**Key stats:**

- 26,998 to 22,357 unique
- Unfiltered: 9,496 / 12,861
- Balanced: 8,453 / 11,919

### Step 6 — Scan filter + balance (`8_filter_balance.py`)

Collapses merged_clean.csv to unique sites and filters by number_of_scans greater than 1. For MSV8 only 3 sites removed - very high coverage dataset. We primarily use filtered_unbalanced as the scan filter has minimal impact.

**Key stats:**

- 22,357 to 22,354 unique (-3)
- filtered_unbalanced: 9,493 / 12,861
- filtered_balanced: 8,535 / 11,860

## Scripts (run in order)

1. `ac_extract.py`
2. `phospho_extract.py`
3. `merge.py`
4. `6_clean_columns.py`
5. `7_reduceK.py`
6. `8_filter_balance.py`
