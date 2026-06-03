# PXD005793 — SILAC Acetylome (3 reps + 3 proteome)

**PTM type:** Acetylation  
**Repository:** PRIDE  
**Dataset:** `PXD005793`  

## Workflow

### Step 1 — Extract (6 runs) (`1_extract.py`)

Run separately for 3 Kac enrichment folders (positives + negatives) and 3 Proteome folders (all negatives). Streams evidence.txt, filters quality, parses modified sequences for K(ac), applies terminal-K rescue, then verifies each site against FASTA with tryptic check.

**Key stats:**

- Kac 1+2+3: 8,958 rows
- Proteome 1+2+3: 36,539 rows
- Total per replicate

### Step 2 — Merge all 6 CSVs (`merge.py`)

Concatenates all 3 Kac and 3 Proteome extracted CSVs into one merged file without deduplication. Reports per-file and total label counts.

**Key stats:**

- 6 input CSVs
- 45,497 total rows
- 7,697 pos / 37,800 neg

### Step 3 — Add K position flags (`4_add_last_columns.py`)

Flag K_first (K at position 0) and K_second_last (K at penultimate position). These expose systematic trypsin cleavage bias: negatives are over-represented at terminal K positions because trypsin cleaves after K.

**Key stats:**

- K_first overall: 29.2%
- Pos K_sec_last: 47.45%
- Neg K_sec_last: 16.37%

### Step 4 — Dedup + reduceK (`5_reduceK.py`)

Collapse to unique sites. Then two-stage K-position balancing: Stage 1 downsamples negatives until K_first percent matches positives. Stage 2 downsamples positives until K_second_last percent matches updated negatives.

**Key stats:**

- 45,497 to 9,444 unique
- Unfiltered: 1,775 / 7,669
- Balanced: 1,377 / 5,472

### Step 5 — Scan filter + balance (`6_filter_balance.py`)

Collapse to unique sites, filter out those seen in only 1 scan, then apply the same two-stage K-position balancing. Produces filtered_unbalanced (our primary output) and filtered_balanced.

**Key stats:**

- 9,444 to 8,796 unique
- filtered_unbalanced: 1,679 / 7,117
- filtered_balanced: 1,323 / 5,088

## Scripts (run in order)

1. `1_extract.py`
2. `merge.py`
3. `4_add_last_columns.py`
4. `5_reduceK.py`
5. `6_filter_balance.py`
