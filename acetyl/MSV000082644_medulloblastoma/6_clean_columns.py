#!/usr/bin/env python3
import pandas as pd
import re
import argparse

ap = argparse.ArgumentParser(description="Clean MSV/MOD dataset and compute K-first and K-second-last.")
ap.add_argument("--input", required=True, help="Input CSV")
ap.add_argument("--out", default="cleaned_final.csv", help="Output CSV")
args = ap.parse_args()

df = pd.read_csv(args.input)

# ---------------------------
# Clean peptide sequence
# ---------------------------
def clean_peptide(seq):
    if pd.isna(seq):
        return None
    s = str(seq)
    s = re.sub(r"\(\d+(\.\d+)?\)", "", s)   # remove (0.99), (1.0)
    s = s.replace("_", "")
    return s.strip().upper()

# ---------------------------
# K-first / K-second-last
# ---------------------------
def has_k_first(seq):
    if not isinstance(seq, str) or len(seq) == 0:
        return 0
    return 1 if seq[0] == "K" else 0

def has_k_second_last(seq):
    if not isinstance(seq, str) or len(seq) < 2:
        return 0
    return 1 if seq[-2] == "K" else 0

# ---------------------------
# Create cleaned output
# ---------------------------
clean = pd.DataFrame()

clean["sequence"] = df["sequenceVariablyModifiedLocations"].apply(clean_peptide)
clean["protein_id"] = df["mapped_uniprot_acc"]
clean["protein_pos"] = df["mapped_position"].astype("Int64")
clean["number_of_scans"] = df["Number_of_scans"].astype(int)
clean["label"] = df["label"].astype(str).str.capitalize()

# compute K-first and K-second-last
clean["K_first"] = clean["sequence"].apply(has_k_first)
clean["K_second_last"] = clean["sequence"].apply(has_k_second_last)

# ---------------------------
# Save
# ---------------------------
clean.to_csv(args.out, index=False)

print(f"Saved cleaned file → {args.out}")
print("Preview:")
print(clean.head(10))
