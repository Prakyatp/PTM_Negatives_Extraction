# MSV000078509 — Kac from High_Kac + pSTY sheets

**PTM type:** Acetylation  
**Repository:** MassIVE  
**Dataset:** `MSV000078509`  

## Workflow

### Step 1 — Extract Table 1 (`11_extract.py`)

One script reads two sheets from Supp_Table1.xlsx. Step 1 (High_Kac): single-protein filter, K-only filter, max_loc_prob > 0.99 identifies Positives. Step 2 (High_pSTY): uses terminal-K prefix keys built from positives to extract K residues in phospho-peptides as Negatives.

**Key stats:**
- High_Kac: 3,159 rows
- 2,243 positives mapped
- 5,941 K-negatives mapped

**Run:**  
This script uses a hardcoded config block at the top. Edit the following variables before running:
```python
TABLE1_FILE  = "Supp_Table1.xlsx"   # path to your XLSX
FASTA_FILE   = "human.fasta.gz"     # path to your FASTA
OUT_STRICT   = "TP_TN_Kac_pSTY_STRICT.csv"
```
Then run:
```bash
python 11_extract.py
```

---

### Step 2 — Extract Table 2 (`22_extract_table2.py`)

Loads High_Proteome sheet. Isoform-aware single-protein filter: collapses isoform variants to canonical base IDs, accepts only peptides with exactly one canonical protein. Extracts ALL K residues from uniquely tryptic-mapped peptides as Negatives.

**Key stats:**
- 189,732 input rows
- 63,205 after protein filter
- 47,945 K-negatives

**Run:**  
Edit the config block at the top:
```python
IN_FILE    = "Supp_Table2.xlsx"    # path to your XLSX
FASTA_FILE = "human.fasta.gz"      # path to your FASTA
OUT_FILE   = "Final_Supp_Table2_extracted.csv"
```
Then run:
```bash
python 22_extract_table2.py
```

---

### Step 3 — Merge T1 + T2 (`merge.py`)

Concatenates TP_TN (positives + pSTY negatives) with Table 2 K-negatives. Applies the terminal-K rule and decorates Kac-positive sequences. K_first and K_second_last flags computed directly in this script.

**Key stats:**
- 47,945 neg + 8,184 TP_TN
- 56,129 merged rows
- 2,243 pos / 53,886 neg

**Run:**
```bash
python merge.py \
  --inputs TP_TN_Kac_pSTY_STRICT.csv Final_Supp_Table2_extracted.csv \
  --output merged.csv
```

---

### Step 4 — Dedup + reduceK (`8_reduceK.py`)

Dedup collapses 56,129 rows to 35,307 unique sites (unfiltered). Then two-stage K-position balancing: Stage 1 downsamples Negatives until K_first % matches Positives. Stage 2 downsamples Positives until K_second_last % matches.

**Key stats:**
- 56,129 to 35,307 unique
- Unfiltered: 2,213 / 33,094
- Balanced: 2,213 / 31,743

**Run:**
```bash
python 8_reduceK.py \
  --input merged.csv \
  --pre-out unfiltered_unique.csv \
  --out balanced_unique.csv
```

---

### Step 5 — Filter (Table 2 only) (`filter.py`)

Goes back to `Final_Supp_Table2_extracted.csv` and filters by distinct scan observations. A K site must appear with more than 1 distinct Number of Scans value to be kept. TP_TN Table 1 sites are pre-validated.

**Key stats:**
- 47,945 to 20,077 rows
- 8,135 sites kept
- 27,860 sites removed

**Run:**
```bash
python filter.py \
  --input Final_Supp_Table2_extracted.csv \
  --out filtered_table2.csv
```

---

### Step 6 — Merge + reduceK (filtered) (`8_reduceK.py`)

Re-run merge.py with filtered_table2.csv and the original TP_TN positives. Then run 8_reduceK.py on the filtered merged output. Produces filtered_unbalanced and filtered_balanced.

**Key stats:**
- 28,261 merged rows
- filtered_unbalanced: 2,213 / 11,628
- filtered_balanced: 2,213 / 11,259

**Run:**
```bash
python merge.py \
  --inputs TP_TN_Kac_pSTY_STRICT.csv filtered_table2.csv \
  --output filtered_merged.csv

python 8_reduceK.py \
  --input filtered_merged.csv \
  --pre-out filtered_unbalanced.csv \
  --out filtered_balanced.csv
```

---

## Scripts (run in order)

1. `11_extract.py`
2. `22_extract_table2.py`
3. `merge.py`
4. `8_reduceK.py`
5. `filter.py`
6. `merge.py` + `8_reduceK.py` (filtered pass)
