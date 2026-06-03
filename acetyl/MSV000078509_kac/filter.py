#!/usr/bin/env python3
import pandas as pd
import argparse

ap = argparse.ArgumentParser(
    description="Keep only sites seen with more than one distinct Number of Scans value."
)
ap.add_argument("--input",  required=True, help="Input CSV (Supp_Table2_extracted.csv)")
ap.add_argument("--out",    required=True, help="Output filtered CSV")
args = ap.parse_args()

df = pd.read_csv(args.input)
print(f"Total rows loaded: {len(df)}")

if "Number of Scans" not in df.columns:
    raise SystemExit("❌ Missing 'Number of Scans' column.")

# Define site key
site_cols = ["sequence", "k_index", "protein_id", "protein_pos"]

# Count distinct Number of Scans values per site
distinct_scans = (
    df.groupby(site_cols)["Number of Scans"]
    .nunique()
    .reset_index()
    .rename(columns={"Number of Scans": "distinct_scan_count"})
)

# Keep sites with more than 1 distinct scan value
keep_sites = distinct_scans[distinct_scans["distinct_scan_count"] > 1]
print(f"Sites with > 1 distinct scan observation: {len(keep_sites)}")
print(f"Sites with <= 1 distinct scan observation (removed): {len(distinct_scans) - len(keep_sites)}")

# Merge back to keep all rows belonging to valid sites
filtered = df.merge(keep_sites[site_cols], on=site_cols, how="inner")
print(f"Rows after filter: {len(filtered)}")

filtered.to_csv(args.out, index=False)
print(f"💾 Saved to {args.out}")
