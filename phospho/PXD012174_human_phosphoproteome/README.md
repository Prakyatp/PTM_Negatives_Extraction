# PXD012174 — Human Phosphoproteome Map

**PTM type:** Phosphorylation  
**Repository:** PRIDE  
**Dataset:** `PXD012174`  

## Workflow

### Step 1 — Extract (`extract.py`)

Stream evidence.txt in 50k-row chunks. Filter reverse hits, contaminants, PEP ≥ 0.99, and multi-protein rows. Join peptides.txt for per-protein start positions. Label every S/T/Y site Positive (annotated ph) or Negative.

**Key stats:**
- evidence.txt + peptides.txt
- 24.4M raw input rows
- Outputs extracted.csv

**Run:**
```bash
python extract.py \
  --evidence path/to/evidence.txt \
  --peptides path/to/peptides.txt \
  --out 1_extracted.csv
```

---

### Step 2 — Map (`new_map.py`)

Pass 1 resolves protein IDs to canonical UniProt accessions (handles isoforms, _HUMAN suffixes). Pass 2 runs in parallel with 35 workers: 4-check verification per unique combo (AA match, FASTA position, start_pos agreement, tryptic filter).

**Key stats:**
- 535,231 unique combos
- 35 parallel workers
- 20.3M rows kept

**Run:**
```bash
python new_map.py \
  --csv 1_extracted.csv \
  --fasta path/to/human.fasta.gz \
  --out 2_mapped.csv \
  --workers 35
```

---

### Step 3 — Dedup (Unfiltered) (`dedup.py`)

Collapse 20.3M mapped rows to unique sites on key (leading_protein, Protein_pos, label). Flag ambiguous positions where the same site appears with both labels. Produces the unfiltered dataset of 361,789 unique sites.

**Key stats:**
- 20.3M to 361,789 sites
- Removed 19.9M duplicates
- 80,452 pos / 281,337 neg

**Run:**
```bash
python dedup.py \
  --csv 2_mapped.csv \
  --out 3_unfiltered_dedup.csv
```

---

### Step 4 — Filter (`filter.py`)

Applied to the full mapped file before the filtered dedup. Retains only peptides observed at 2 or more distinct retention times. Removes single-observation peptides that are likely spurious identifications.

**Key stats:**
- 20.3M to 20.2M
- Removed 106,251 rows
- 34,438 sequences removed

**Run:**
```bash
python filter.py \
  --input 2_mapped.csv \
  --out 4_filtered.csv
```

---

### Step 5 — Dedup (Filtered) (`dedup.py`)

Same dedup logic applied to the retention-time filtered mapped file. Produces the filtered dataset of 295,894 unique sites with higher confidence.

**Key stats:**
- 20.2M to 295,894 sites
- 67,847 pos / 228,047 neg

**Run:**
```bash
python dedup.py \
  --csv 4_filtered.csv \
  --out 5_filtered_dedup.csv
```

---

### Step 6 — Split S / T / Y (`split_dataset.py`)

Split both unfiltered and filtered datasets by amino acid. Each residue type has distinct phosphorylation biology, frequency, and class balance, so they are kept as separate downstream datasets.

**Key stats:**
- 3 amino acids
- 6 output CSVs total
- Unfiltered + Filtered

**Run:**
```bash
# Split unfiltered
python split_dataset.py \
  --inp 3_unfiltered_dedup.csv \
  --out-prefix unfiltered

# Split filtered
python split_dataset.py \
  --inp 5_filtered_dedup.csv \
  --out-prefix filtered
```
Produces `{prefix}_S.csv`, `{prefix}_T.csv`, `{prefix}_Y.csv` for each run.

---

## Scripts (run in order)

1. `extract.py`
2. `new_map.py`
3. `dedup.py`
4. `filter.py`
5. `dedup.py` (on filtered output)
6. `split_dataset.py`
