#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse

ap = argparse.ArgumentParser(
    description="Two-stage downsampling: (1) Negatives on K_first, then (2) Positives on K_second_last, with UNIQUE-site summaries and unique-only CSVs."
)
ap.add_argument('--input', required=True,
                help='Input CSV with columns: label, protein_id, protein_pos, K_first, K_second_last (0/1)')
ap.add_argument('--pre-out', default='acetylation_sites_unique_prebalance.csv',
                help='Output CSV path for UNIQUE sites before balancing')
ap.add_argument('--out', default='acetylation_sites_downsampled_unique.csv',
                help='Output CSV path for UNIQUE sites after balancing')
ap.add_argument('--seed', type=int, default=42,
                help='Random seed for reproducibility')
args = ap.parse_args()

rng = np.random.default_rng(args.seed)
df = pd.read_csv(args.input)

# ---------------------------
# sanity checks
# ---------------------------
req = ['label', 'K_first', 'K_second_last', 'protein_id', 'protein_pos']
missing = [c for c in req if c not in df.columns]
if missing:
    raise SystemExit(
        f"❌ Missing required columns: {missing}. "
        f"Need label, protein_id, protein_pos, K_first, K_second_last."
    )

def pct(series_bool):
    return float(series_bool.mean()) if len(series_bool) else 0.0

def summarize(name, frame):
    """
    Summarize in terms of UNIQUE SITES (already unique by protein_id, protein_pos, label).
    K_first and K_second_last are taken directly per row (since rows are unique sites).
    """
    if frame.empty:
        print(f"{name}: 0 sites")
        return

    for lbl in ['Positive', 'Negative']:
        sub = frame[frame['label'] == lbl]
        if sub.empty:
            print(f"{name} — {lbl}: 0 sites")
            continue

        n_sites = len(sub)
        kf = pct(sub['K_first'] == 1) * 100
        ks = pct(sub['K_second_last'] == 1) * 100
        print(f"{name} — {lbl}: unique_sites={n_sites} | K_first={kf:.2f}% | K_second_last={ks:.2f}%")

# ---------------------------
# 1) Deduplicate to UNIQUE sites pre-balance
# ---------------------------
df_unique = df.drop_duplicates(subset=['protein_id', 'protein_pos', 'label']).copy()
print(f"Initial rows (input): {len(df)}")
print(f"Rows after unique-site dedup (protein_id, protein_pos, label): {len(df_unique)}")
print(f"Duplicates removed: {len(df) - len(df_unique)}")

# Save pre-balance unique CSV
df_unique.to_csv(args.pre_out, index=False)
print(f"💾 Saved pre-balance UNIQUE sites to: {args.pre_out}")

print("\n=== BEFORE (UNIQUE-SITE SUMMARY) ===")
summarize("Before", df_unique)

pos = df_unique[df_unique['label'] == 'Positive'].copy()
neg = df_unique[df_unique['label'] == 'Negative'].copy()

# ---------------------------
# 2) Stage 1: Downsample NEGATIVES on K_first to match POSITIVE K_first rate (row-based on unique sites)
# ---------------------------
target_kf_pct = pct(pos['K_first'] == 1)

def downsample_flag_to_target(frame, flag_col, target_pct, rng):
    """Drop rows with flag_col==1 to reach target percentage; keep others intact (row-based)."""
    if len(frame) == 0:
        return frame
    cur_pct = pct(frame[flag_col] == 1)
    if cur_pct <= target_pct or cur_pct == 0.0:
        return frame  # already at/below target (or none to drop)

    n1 = int((frame[flag_col] == 1).sum())
    n0 = int((frame[flag_col] == 0).sum())

    if target_pct >= 1.0:
        return frame

    # target_pct = x / (x + n0)  =>  x = target_pct*n0 / (1 - target_pct)
    x = int(np.floor((target_pct * n0) / (1.0 - target_pct)))
    x = max(0, min(x, n1))

    flag1_idx = frame.index[frame[flag_col] == 1].to_numpy()
    keep_idx1 = rng.choice(flag1_idx, size=x, replace=False) if x < n1 else flag1_idx

    keep_mask = frame[flag_col] == 0
    keep_mask.loc[keep_idx1] = True
    return frame[keep_mask]

neg_stage1 = downsample_flag_to_target(neg, 'K_first', target_kf_pct, rng)

# ---------------------------
# 3) Stage 2: Downsample POSITIVES on K_second_last to match UPDATED NEGATIVE K_second_last rate
# ---------------------------
updated_ks_pct_neg = pct(neg_stage1['K_second_last'] == 1)
pos_stage2 = downsample_flag_to_target(pos, 'K_second_last', updated_ks_pct_neg, rng)

# Recombine
downsampled = pd.concat([pos_stage2, neg_stage1], ignore_index=True)

# Ensure uniqueness again after downsampling (paranoia / safety)
downsampled_unique = downsampled.drop_duplicates(subset=['protein_id', 'protein_pos', 'label']).copy()

print("\n=== AFTER (UNIQUE-SITE SUMMARY, TWO-STAGE) ===")
summarize("After", downsampled_unique)

# Save post-balance unique CSV
downsampled_unique.to_csv(args.out, index=False)
print(f"\n✅ Wrote post-balance UNIQUE sites to: {args.out}")
