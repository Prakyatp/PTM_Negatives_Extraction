# MSV000078509 — Serial Enrichment – Supp Tables 1 & 2

**PTM type:** Phosphorylation  
**Repository:** MassIVE  
**Dataset:** `MSV000078509`  

## Workflow

### Step 1 — Extract Table 1 (`extract_table1.py`)

Parse the MaxQuant Phospho(STY)Sites XLSX. For each row, parse the probability string to get per-site probabilities. Apply threshold rules (prob ≥ 0.99 → Positive, prob ≤ 0.01 → Negative). Handle multi-phosphopeptide ambiguity. Protein_pos comes from the Position column.

**Key stats:**
- 27,392 input rows
- Outputs 74,625 site rows
- 3,486 multi-prot skipped

**Run:**
```bash
python extract_table1.py \
  --xlsx path/to/Supp_Table1.xlsx \
  --out t1_extracted.csv
```

---

### Step 2 — Extract Table 2 (`extract_table2.py`)

Parse the baseline proteome XLSX. Every S/T/Y in a uniquely-mapped, single-canonical-protein peptide becomes a Negative. Requires exactly one occurrence of the peptide in the protein (no multi-match ambiguity).

**Key stats:**
- 189,732 input rows
- 63,205 after protein filter
- Outputs 143,777 STY-neg

**Run:**
```bash
python extract_table2.py \
  --xlsx path/to/Supp_Table2.xlsx \
  --fasta path/to/human.fasta.gz \
  --out t2_extracted.csv
```

---

### Step 3 — Map (`new_map.py`)

Same FASTA verification as PXD. Pass 1 resolves IDs; Pass 2 runs 4-check verification (AA, position, tryptic) in parallel. Run once for Table 1 output and once for Table 2 output.

**Key stats:**
- T1: 74,625 to 41,182
- T2: 143,777 to 141,374
- 62,166 combos verified

**Run:**
```bash
# Map Table 1
python new_map.py \
  --csv t1_extracted.csv \
  --fasta path/to/human.fasta.gz \
  --out t1_mapped.csv \
  --workers 35

# Map Table 2
python new_map.py \
  --csv t2_extracted.csv \
  --fasta path/to/human.fasta.gz \
  --out t2_mapped.csv \
  --workers 35
```

---

### Step 4 — Merge T1 + T2 (`new_merge.py`)

Concatenate Table 1 mapped (positives) and Table 2 mapped (negatives) into one unified CSV. Both tables share the same column schema after mapping.

**Key stats:**
- 41,182 pos + 143,777 neg
- Merged 100,122 rows

**Run:**
```bash
python new_merge.py \
  --pos t1_mapped.csv \
  --neg t2_mapped.csv \
  --out merged.csv
```

---

### Step 5 — Dedup (Unfiltered) (`dedup.py`)

Deduplicate merged file on (protein, position, label). Flag ambiguous sites where the same position has both labels. Produces the unfiltered MSV phospho dataset of 58,824 unique sites.

**Key stats:**
- 100,122 to 58,824
- 12,826 pos / 45,998 neg
- 2,184 ambiguous

**Run:**
```bash
python dedup.py \
  --csv merged.csv \
  --out unfiltered_dedup.csv
```

---

### Step 6 — Filter + Dedup (`filter.py`)

Apply the retention-time filter (sequences seen at 2 or more RTs) to the merged file, then dedup again to produce the filtered MSV phospho dataset.

**Key stats:**
- Filtered unique: 58,824
- 12,826 pos / 45,998 neg

**Run:**
```bash
python filter.py \
  --input merged.csv \
  --out filtered.csv

python dedup.py \
  --csv filtered.csv \
  --out filtered_dedup.csv
```

---

### Step 7 — Split S / T / Y (`split_dataset.py`)

Split both unfiltered and filtered MSV datasets by amino acid. Six output CSVs total.

**Key stats:**
- 3 amino acids
- 6 output CSVs
- Unfiltered + Filtered

**Run:**
```bash
# Split unfiltered
python split_dataset.py \
  --inp unfiltered_dedup.csv \
  --out-prefix unfiltered

# Split filtered
python split_dataset.py \
  --inp filtered_dedup.csv \
  --out-prefix filtered
```
Produces `{prefix}_S.csv`, `{prefix}_T.csv`, `{prefix}_Y.csv` for each run.

---

## Scripts (run in order)

1. `extract_table1.py`
2. `extract_table2.py`
3. `new_map.py` (×2 — once per table)
4. `new_merge.py`
5. `dedup.py`
6. `filter.py` → `dedup.py`
7. `split_dataset.py`
