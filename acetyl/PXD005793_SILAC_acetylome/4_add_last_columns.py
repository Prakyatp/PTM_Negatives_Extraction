#!/usr/bin/env python3
import pandas as pd
import argparse

ap = argparse.ArgumentParser(description="Check how many peptides have K at first or second-last position.")
ap.add_argument('--input', required=True, help='Input CSV file with at least a sequence column')
ap.add_argument('--out', default=None, help='Optional CSV to save with K flags added')
args = ap.parse_args()

# -----------------------------
# Load
# -----------------------------
df = pd.read_csv(args.input)
if 'sequence' not in df.columns:
    raise SystemExit("❌ Missing 'sequence' column in input CSV.")

# -----------------------------
# Compute flags
# -----------------------------
def has_k_first(seq):
    if not isinstance(seq, str) or len(seq) == 0:
        return 0
    return 1 if seq[0].upper() == 'K' else 0

def has_k_second_last(seq):
    if not isinstance(seq, str) or len(seq) < 2:
        return 0
    return 1 if seq[-2].upper() == 'K' else 0

df['K_first'] = df['sequence'].apply(has_k_first)
df['K_second_last'] = df['sequence'].apply(has_k_second_last)

# -----------------------------
# Overall stats
# -----------------------------
total = len(df)
k_first_total = df['K_first'].sum()
k_second_total = df['K_second_last'].sum()

print("\n--- Trypsin Boundary K Analysis (Overall) ---")
print(f"Total sites analyzed: {total}")
print(f"K at first position:       {k_first_total}  ({k_first_total/total*100:.2f}%)")
print(f"K at second-last position: {k_second_total}  ({k_second_total/total*100:.2f}%)")

# -----------------------------
# By label (if present)
# -----------------------------
if 'label' in df.columns:
    print("\n--- K Occurrence by Label ---")
    stats = (
        df.groupby('label')[['K_first', 'K_second_last']]
        .sum()
        .assign(total=df.groupby('label').size())
    )
    stats['K_first_%'] = (stats['K_first'] / stats['total'] * 100).round(2)
    stats['K_second_last_%'] = (stats['K_second_last'] / stats['total'] * 100).round(2)
    print(stats)

# -----------------------------
# Optional save
# -----------------------------
if args.out:
    df.to_csv(args.out, index=False)
    print(f"\n✅ Wrote {args.out} (includes K_first and K_second_last columns)")
