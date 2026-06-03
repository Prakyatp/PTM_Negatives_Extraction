#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse

ap = argparse.ArgumentParser(
    description=(
        "Merge acetylation datasets, collapse to unique sites, filter by number_of_scans>1, "
        "then perform two-stage downsampling (Negatives on K_first, Positives on K_second_last) "
        "with UNIQUE-site summaries."
    )
)

ap.add_argument(
    '--input', nargs='+', required=True,
    help='Input CSV(s) with columns: sequence, protein_id, protein_pos, number_of_scans, label, K_first, K_second_last'
)
ap.add_argument(
    '--pre-out', default='filtered_unique_prebalance.csv',
    help='Output CSV for filtered UNIQUE sites before balancing'
)
ap.add_argument(
    '--out', default='filtered_balanced_unique.csv',
    help='Output CSV for filtered UNIQUE sites after balancing'
)
ap.add_argument(
    '--seed', type=int, default=42,
    help='Random seed for reproducibility'
)
args = ap.parse_args()

rng = np.random.default_rng(args.seed)

# ---------------------------
# 0) Load and merge all inputs
# ---------------------------
frames = []
for path in args.input:
    print(f"Loading: {path}")
    df_part = pd.read_csv(path)
    df_part["source_file"] = path  # optional, can be kept or dropped later
    frames.append(df_part)

df = pd.concat(frames, ignore_index=True)
print(f"\nTotal rows after merging {len(args.input)} file(s): {len(df)}")

# ---------------------------
# sanity checks
# ---------------------------
req = [
    'sequence',
    'protein_id',
    'protein_pos',
    'number_of_scans',
    'label',
    'K_first',
    'K_second_last'
]
missing = [c for c in req if c not in df.columns]
if missing:
    raise SystemExit(
        f"❌ Missing required columns: {missing}. "
        f"Need sequence, protein_id, protein_pos, number_of_scans, label, K_first, K_second_last."
    )

def pct(series_bool):
    return float(series_bool.mean()) if len(series_bool) else 0.0

def summarize(name, frame):
    """
    Summarize in terms of UNIQUE SITES (already unique by protein_id, protein_pos, label).
    K_first and K_second_last are taken directly per row (site-level).
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
# 1) Collapse to UNIQUE sites (protein_id, protein_pos, label) and aggregate
# ---------------------------
print("\nCollapsing to unique (protein_id, protein_pos, label) sites...")

site_cols = ['protein_id', 'protein_pos', 'label']
agg = {
    'sequence': 'first',
    'number_of_scans': 'max',   # max scans across replicates/datasets
    'K_first': 'max',           # if any row at that site has 1 → site=1
    'K_second_last': 'max',
    'source_file': 'first'
}

sites = df.groupby(site_cols, as_index=False).agg(agg)

print(f"Total unique sites before scan filtering: {len(sites)}")

# ---------------------------
# 2) Filter: keep only sites with number_of_scans > 1
# ---------------------------
filtered = sites[sites['number_of_scans'] > 1].copy()
removed = sites[sites['number_of_scans'] <= 1].copy()

print(f"Sites removed (number_of_scans <= 1): {len(removed)}")
print(f"Sites kept   (number_of_scans > 1): {len(filtered)}")

print("\n=== UNIQUE-SITE SUMMARY AFTER SCAN FILTER (PRE-BALANCE) ===")
summarize("Filtered", filtered)

# Save pre-balance, filtered unique sites
filtered.to_csv(args.pre_out, index=False)
print(f"\n💾 Saved filtered UNIQUE sites (pre-balance) to: {args.pre_out}")

# ---------------------------
# 3) Two-stage balancing on filtered unique sites
# ---------------------------
pos = filtered[filtered['label'] == 'Positive'].copy()
neg = filtered[filtered['label'] == 'Negative'].copy()

def downsample_flag_to_target(frame, flag_col, target_pct, rng):
    """Drop rows with flag_col==1 to reach target percentage; keep others intact."""
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

# Stage 1: match Negative K_first rate to Positive K_first rate
target_kf_pct = pct(pos['K_first'] == 1)
neg_stage1 = downsample_flag_to_target(neg, 'K_first', target_kf_pct, rng)

# Stage 2: match Positive K_second_last rate to UPDATED Negative K_second_last rate
updated_ks_pct_neg = pct(neg_stage1['K_second_last'] == 1)
pos_stage2 = downsample_flag_to_target(pos, 'K_second_last', updated_ks_pct_neg, rng)

# Recombine
downsampled = pd.concat([pos_stage2, neg_stage1], ignore_index=True)

# Ensure uniqueness (just in case)
downsampled_unique = downsampled.drop_duplicates(subset=['protein_id', 'protein_pos', 'label']).copy()

print("\n=== AFTER BALANCING (UNIQUE-SITE SUMMARY) ===")
summarize("Balanced", downsampled_unique)

# ---------------------------
# 4) Save final filtered + balanced unique sites
# ---------------------------
downsampled_unique.to_csv(args.out, index=False)
print(f"\n✅ Saved filtered + balanced UNIQUE sites to: {args.out}")
