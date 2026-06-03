#!/usr/bin/env python3
import pandas as pd
import re
import gzip

# ========= CONFIG =========
IN_FILE    = "Supp_Table2.xlsx"
SHEET      = "High_Proteome"
FASTA_FILE = "/home/pp6448/human.fasta.gz"
OUT_FILE   = "Final_Supp_Table2_extracted.csv"
# ==========================


# ==========================
# FASTA loader
# ==========================
def load_fasta(path):
    seqs = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        current_id, current_seq = None, []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(current_seq)
                parts = line[1:].split("|")
                current_id = parts[1].strip() if len(parts) >= 2 else line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            seqs[current_id] = "".join(current_seq)
    return seqs


# ==========================
# Tryptic mapping helper
# ==========================
def find_tryptic_start(pep: str, prot_seq: str):
    """
    Find all occurrences of pep in prot_seq.
    Keep only tryptic ones: residue before start must be K, R, or start==0.
    Returns:
        (start, "unique") — exactly one tryptic occurrence
        (None,  "none")   — no tryptic occurrence
        (None,  "multi")  — multiple tryptic occurrences (ambiguous, drop)
    """
    tryptic_starts = []
    pos = 0
    while True:
        idx = prot_seq.find(pep, pos)
        if idx == -1:
            break
        if idx == 0 or prot_seq[idx - 1] in ("K", "R"):
            tryptic_starts.append(idx)
        pos = idx + 1

    if len(tryptic_starts) == 0:
        return None, "none"
    if len(tryptic_starts) == 1:
        return tryptic_starts[0], "unique"
    return None, "multi"


# ==========================
# Helpers
# ==========================
def clean_pepseq(seq: str) -> str:
    if pd.isna(seq):
        return ""
    s = str(seq).strip().strip("_")
    s = re.sub(r"\(.*?\)", "", s)
    return s.upper()

def k_indices_keep_all(pep: str):
    """Return 0-based indices of ALL Ks in the peptide."""
    if not pep:
        return []
    return [i for i, aa in enumerate(pep) if aa == "K"]

def is_single_base_protein(val) -> bool:
    if pd.isna(val):
        return False
    s = str(val).strip()
    if s == "" or s.lower() in {"na", "nan", "none", "null"}:
        return False
    tokens = [t.strip() for t in s.split(";") if t.strip()]
    if not tokens:
        return False
    base_ids = set(t.split("-")[0] for t in tokens)
    return len(base_ids) == 1

def get_base_protein_id(val) -> str:
    tokens = [t.strip() for t in str(val).split(";") if t.strip()]
    return tokens[0].split("-")[0]


# ==========================
# STEP 0: Load FASTA
# ==========================
print("=" * 50)
print("STEP 0: Loading FASTA")
print("=" * 50)
print(f"  File: {FASTA_FILE}")
seq_dict = load_fasta(FASTA_FILE)
print(f"  ✅ Loaded {len(seq_dict)} protein sequences")


# ==========================
# STEP 1: Load & Filter Table 2
# ==========================
print()
print("=" * 50)
print("STEP 1: Loading & Filtering Supp_Table2")
print("=" * 50)

print(f"  File: {IN_FILE}  Sheet: {SHEET}")
df = pd.read_excel(IN_FILE, sheet_name=SHEET)
print(f"  Total rows loaded:                         {len(df)}")

if "Peptide sequence" not in df.columns:
    raise KeyError(f"Missing required column 'Peptide sequence' in sheet '{SHEET}'")

before = len(df)
mask = df["Proteins"].apply(is_single_base_protein)
df = df[mask].copy()
print(f"  After isoform-aware single-protein filter: {len(df)} (removed {before - len(df)})")

df["sequence"]   = df["Peptide sequence"].apply(clean_pepseq)
df["protein_id"] = df["Proteins"].apply(get_base_protein_id)


# ==========================
# STEP 2: Extract Ks & FASTA map (tryptic)
# ==========================
print()
print("=" * 50)
print("STEP 2: Extracting Ks & Mapping to FASTA (tryptic)")
print("=" * 50)

rows             = []
missing_fasta    = 0
not_found        = 0
multi_tryptic    = 0

for _, row in df.iterrows():
    pep = row["sequence"]
    pid = row["protein_id"]
    k_idxs_0 = k_indices_keep_all(pep)

    if not k_idxs_0:
        continue

    prot_seq = seq_dict.get(pid)
    if prot_seq is None:
        missing_fasta += 1
        continue

    prot_seq = prot_seq.upper()
    start, status = find_tryptic_start(pep, prot_seq)

    if status == "none":
        not_found += 1
        continue
    if status == "multi":
        multi_tryptic += 1
        continue

    num_scans_val = row.get("Retention time", None)

    for ki0 in k_idxs_0:
        protein_pos = start + 1 + ki0
        k_index_1   = ki0 + 1
        rows.append({
            "sequence":      pep,
            "k_index":       k_index_1,
            "protein_id":    pid,
            "protein_pos":   protein_pos,
            "Proteins":      row.get("Proteins", None),
            "Protein Names": row.get("Protein Names", None),
            "Experiment":    row.get("Experiment", None),
            "Modifications": row.get("Modifications", None),
            "Charge":        row.get("Charge", None),
            "m/z":           row.get("m/z", None),
            "Number of Scans": num_scans_val,
            "Score":         row.get("Score", None),
            "Ratio H/L":     row.get("Ratio H/L", None),
            "label":         "Negative",
            "source":        "Supp_Table2_internalK",
        })

out_df = pd.DataFrame(rows)
print(f"  ✅ Successfully mapped:          {len(out_df)}")
print(f"  ❌ Protein missing in FASTA:    {missing_fasta}")
print(f"  ❌ No tryptic occurrence:       {not_found}")
print(f"  ❌ Multiple tryptic (dropped):  {multi_tryptic}")
print(f"  Unique peptides with K:         {out_df['sequence'].nunique()}")

out_df.to_csv(OUT_FILE, index=False)
print()
print(f"  💾 Saved to {OUT_FILE} with {len(out_df)} rows.")
