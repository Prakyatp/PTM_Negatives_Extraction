#!/usr/bin/env python3
import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser(
        description="Merge labeled CSVs without deduplication."
    )
    ap.add_argument("--inputs", nargs="+", required=True,
                    help="List of labeled CSV files to merge.")
    ap.add_argument("--output", default="merged_labeled_sites.csv",
                    help="Output CSV file.")
    args = ap.parse_args()

    dfs = []
    for path in args.inputs:
        df_i = pd.read_csv(path)
        print(f"  {path}: {len(df_i):,} rows")
        dfs.append(df_i)

    df = pd.concat(dfs, ignore_index=True)

    pos_count = (df["label"] == "Positive").sum()
    neg_count = (df["label"] == "Negative").sum()

    print(f"\n=== Merge Summary ===")
    print(f"Total rows  : {len(df):,}")
    print(f"Positive    : {pos_count:,}")
    print(f"Negative    : {neg_count:,}")

    df.to_csv(args.output, index=False)
    print(f"✅ Saved → {args.output}")


if __name__ == "__main__":
    main()

