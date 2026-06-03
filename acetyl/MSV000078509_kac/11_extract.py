#!/usr/bin/env python3
import pandas as pd
import re
import gzip

# ==============================
# CONFIG
# ==============================
TABLE1_FILE  = "Supp_Table1.xlsx"
FASTA_FILE   = "/home/pp6448/human.fasta.gz"

SHEET_KAC    = "High_Kac"
SHEET_PSTY   = "High_pSTY"

OUT_STRICT   = "TP_TN_Kac_pSTY_STRICT.csv"

POS_THRESH   = 0.99


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
        (start, "unique")   — exactly one tryptic occurrence
        (None,  "none")     — no tryptic occurrence found
        (None,  "multi")    — multiple tryptic occurrences (ambiguous, drop)
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
def clean_from_prob_str(s):
    if pd.isna(s):
        return ""
    s = str(s).strip()
    s = re.sub(r"\(.*?\)", "", s)
    return s.upper()

def get_mod_index_from_prob_str(s):
    """Returns 0-based index of modified residue."""
    if pd.isna(s):
        return None
    s = str(s)
    pep_index = -1
    for ch in s:
        if ch.isalpha() and ch.isupper():
            pep_index += 1
        if ch == "(":
            return pep_index
    return None

def clean_modseq(seq):
    if pd.isna(seq):
        return ""
    s = str(seq).strip().strip("_")
    s = re.sub(r"\(.*?\)", "", s)
    return s.upper()

def internal_k_indices_only(pep: str):
    if not pep or len(pep) < 3:
        return []
    return [i for i, aa in enumerate(pep) if aa == "K" and 0 < i < len(pep) - 1]

def k_indices_with_terminal_rule(pep: str, terminal_prefix_keys: set):
    idxs = list(internal_k_indices_only(pep))
    if not pep:
        return idxs
    n = len(pep)
    if pep[0] == "K" and (pep, 0) in terminal_prefix_keys:
        idxs.append(0)
    if n > 1 and pep[-1] == "K" and (pep, n - 1) in terminal_prefix_keys:
        idxs.append(n - 1)
    return sorted(set(idxs))

def is_single_leading_protein(val) -> bool:
    if pd.isna(val):
        return False
    s = str(val).strip()
    if s == "" or s.lower() in {"na", "nan", "none", "null"}:
        return False
    return ";" not in s

def get_first_ph_index(modseq):
    """Returns 0-based peptide index of first (ph) residue."""
    if pd.isna(modseq):
        return None
    s = str(modseq).strip().strip("_")
    pep_index = -1
    first_ph_idx = None
    i = 0
    while i < len(s):
        ch = s[i]
        if ch.isalpha() and ch.isupper():
            pep_index += 1
            if first_ph_idx is None and s[i+1:].startswith("(ph)"):
                first_ph_idx = pep_index
            i += 1
        elif ch == "(":
            while i < len(s) and s[i] != ")":
                i += 1
            if i < len(s) and s[i] == ")":
                i += 1
        else:
            i += 1
    return first_ph_idx


# ==========================
# Load FASTA
# ==========================
print("=" * 50)
print("STEP 0: Loading FASTA")
print("=" * 50)
print(f"  File: {FASTA_FILE}")
seq_dict = load_fasta(FASTA_FILE)
print(f"  ✅ Loaded {len(seq_dict)} protein sequences")


# ==========================
# STEP 1: High_Kac
# ==========================
print()
print("=" * 50)
print("STEP 1: High_Kac")
print("=" * 50)

kac = pd.read_excel(TABLE1_FILE, sheet_name=SHEET_KAC)
print(f"  Total rows loaded:              {len(kac)}")

needed_kac_cols = [
    "Amino acid", "Acetyl (K) Probabilities",
    "Localization prob Exp1", "Localization prob Exp2", "Localization prob Exp3",
    "Leading proteins", "Protein names",
]
for c in needed_kac_cols:
    if c not in kac.columns:
        raise KeyError(f"Missing column in High_Kac: {c}")

before_lp = len(kac)
kac = kac[kac["Leading proteins"].apply(is_single_leading_protein)].copy()
print(f"  After single Leading proteins:  {len(kac)} (removed {before_lp - len(kac)})")

kac = kac[kac["Amino acid"] == "K"].copy()
print(f"  After K-only filter:            {len(kac)}")

loc_cols = ["Localization prob Exp1", "Localization prob Exp2", "Localization prob Exp3"]
kac["max_loc_prob"] = kac[loc_cols].max(axis=1, skipna=True).fillna(0.0)
kac["pep_clean"]    = kac["Acetyl (K) Probabilities"].apply(clean_from_prob_str)
kac["k_index_0"]    = kac["Acetyl (K) Probabilities"].apply(get_mod_index_from_prob_str)

df_pos = kac[kac["max_loc_prob"] > POS_THRESH].copy()
df_neg = kac[kac["max_loc_prob"] == 0].copy()
print(f"  Positives (max_loc_prob > {POS_THRESH}): {len(df_pos)}")
print(f"  Negatives (max_loc_prob == 0):  {len(df_neg)}")

# FASTA map positives using tryptic rule
print()
print("  --- Mapping High_Kac Positives to FASTA (tryptic) ---")
pos_rows          = []
pos_missing_fasta = 0
pos_not_found     = 0
pos_multi_tryptic = 0

for _, row in df_pos.iterrows():
    pid = str(row["Leading proteins"]).strip()
    pep = row["pep_clean"]
    ki0 = row["k_index_0"]
    if pd.isna(ki0) or not pep:
        continue

    prot_seq = seq_dict.get(pid)
    if prot_seq is None:
        pos_missing_fasta += 1
        continue

    prot_seq = prot_seq.upper()
    start, status = find_tryptic_start(pep, prot_seq)

    if status == "none":
        pos_not_found += 1
        continue
    if status == "multi":
        pos_multi_tryptic += 1
        continue

    ki0         = int(ki0)
    protein_pos = start + 1 + ki0
    k_index_1   = ki0 + 1

    pos_rows.append({
        "sequence":     pep,
        "k_index":      k_index_1,
        "protein_id":   pid,
        "protein_name": row.get("Protein names", None),
        "protein_pos":  protein_pos,
        "max_loc_prob": row["max_loc_prob"],
        "source":       "High_Kac_pos",
        "label":        "Positive",
    })

df_pos_out = pd.DataFrame(pos_rows)
print(f"  ✅ Successfully mapped:         {len(df_pos_out)}")
print(f"  ❌ Protein missing in FASTA:   {pos_missing_fasta}")
print(f"  ❌ No tryptic occurrence:      {pos_not_found}")
print(f"  ❌ Multiple tryptic (dropped): {pos_multi_tryptic}")

# Build prefix keys for terminal-K rule from positives
prefix_keys_for_terminal = set()
for _, row in df_pos.iterrows():
    pep = row["pep_clean"]
    ki0 = row["k_index_0"]
    if pd.isna(ki0):
        continue
    ki0 = int(ki0)
    if 0 <= ki0 < len(pep):
        prefix_keys_for_terminal.add((pep[: ki0 + 1], ki0))
print(f"  Prefix keys for terminal-K:    {len(prefix_keys_for_terminal)}")


# ==========================
# STEP 2: High_pSTY
# ==========================
print()
print("=" * 50)
print("STEP 2: High_pSTY (Negatives)")
print("=" * 50)

psty = pd.read_excel(TABLE1_FILE, sheet_name=SHEET_PSTY)
print(f"  Total rows loaded:              {len(psty)}")

needed_psty_cols = ["Modified sequence", "Leading proteins", "Protein names", "Gene names"]
for c in needed_psty_cols:
    if c not in psty.columns:
        raise KeyError(f"Missing '{c}' in High_pSTY")

before = len(psty)
psty = psty[psty["Leading proteins"].apply(is_single_leading_protein)].copy()
print(f"  After single Leading proteins:  {len(psty)} (removed {before - len(psty)})")

print()
print("  --- Mapping High_pSTY Negatives to FASTA (tryptic) ---")
neg_rows          = []
neg_missing_fasta = 0
neg_not_found     = 0
neg_multi_tryptic = 0
neg_no_k          = 0
neg_no_ph         = 0

for _, row in psty.iterrows():
    modseq = row["Modified sequence"]
    pep    = clean_modseq(modseq)
    pid    = str(row["Leading proteins"]).strip()

    ph_idx = get_first_ph_index(modseq)
    if ph_idx is None:
        neg_no_ph += 1
        continue

    k_idxs_0 = k_indices_with_terminal_rule(pep, prefix_keys_for_terminal)
    if not k_idxs_0:
        neg_no_k += 1
        continue

    prot_seq = seq_dict.get(pid)
    if prot_seq is None:
        neg_missing_fasta += 1
        continue

    prot_seq = prot_seq.upper()
    start, status = find_tryptic_start(pep, prot_seq)

    if status == "none":
        neg_not_found += 1
        continue
    if status == "multi":
        neg_multi_tryptic += 1
        continue

    for ki0 in k_idxs_0:
        protein_pos = start + 1 + ki0
        k_index_1   = ki0 + 1
        neg_rows.append({
            "sequence":     pep,
            "k_index":      k_index_1,
            "protein_id":   pid,
            "protein_name": row.get("Protein names", None),
            "protein_pos":  protein_pos,
            "max_loc_prob": None,
            "source":       "High_pSTY",
            "label":        "Negative",
        })

df_neg_out = pd.DataFrame(neg_rows)
print(f"  ✅ Successfully mapped:         {len(df_neg_out)}")
print(f"  ❌ Protein missing in FASTA:   {neg_missing_fasta}")
print(f"  ❌ No tryptic occurrence:      {neg_not_found}")
print(f"  ❌ Multiple tryptic (dropped): {neg_multi_tryptic}")
print(f"  ❌ No ph index found:          {neg_no_ph}")
print(f"  ❌ No valid K in peptide:      {neg_no_k}")


# ==========================
# STEP 3: Merge
# ==========================
print()
print("=" * 50)
print("STEP 3: Merging (no dedup)")
print("=" * 50)

all_df = pd.concat([df_pos_out, df_neg_out], ignore_index=True)
print(f"  Positives:                      {len(df_pos_out)}")
print(f"  Negatives:                      {len(df_neg_out)}")
print(f"  Total rows:                     {len(all_df)}")

site_label   = all_df[["protein_id", "protein_pos", "label"]].drop_duplicates()
n_pos_sites  = len(site_label[site_label["label"] == "Positive"])
n_neg_sites  = len(site_label[site_label["label"] == "Negative"])
label_counts = site_label.groupby(["protein_id", "protein_pos"])["label"].nunique()
n_ambiguous  = (label_counts > 1).sum()
print(f"  Unique positive protein sites:  {n_pos_sites}")
print(f"  Unique negative protein sites:  {n_neg_sites}")
print(f"  Ambiguous sites (both labels):  {n_ambiguous}")

all_df.to_csv(OUT_STRICT, index=False)
print()
print(f"  💾 Saved to {OUT_STRICT} with {len(all_df)} rows.")
