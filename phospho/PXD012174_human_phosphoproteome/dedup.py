#!/usr/bin/env python3
import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Remove duplicate phospho sites and summarize counts.\n"
            "Duplicates removed based on (leading_protein, Protein_pos, label). "
            "Ambiguous = positions with both Positive and Negative."
        )
    )
    ap.add_argument("--csv", required=True, help="Input CSV from previous script")
    ap.add_argument(
        "--out",
        default=None,
        help="Output CSV path (default: <input>_dedup.csv)",
    )
    args = ap.parse_args()

    in_path = args.csv
    out_path = args.out
    if out_path is None:
        if in_path.lower().endswith(".csv"):
            out_path = in_path[:-4] + "_dedup.csv"
        else:
            out_path = in_path + "_dedup.csv"

    # Load CSV
    df = pd.read_csv(in_path)

    # Required columns
    required_cols = ["leading_protein", "Protein_pos", "label"]
    for col in required_cols:
        if col not in df.columns:
            raise KeyError(f"Required column '{col}' not found in input CSV.")

    # Normalize label case
    df["label"] = df["label"].astype(str).str.capitalize()

    # 1) Remove duplicates by (leading_protein, Protein_pos, label)
    before = len(df)
    df = df.drop_duplicates(subset=["leading_protein", "Protein_pos", "label"])
    after = len(df)
    n_deleted_duplicates = before - after

    # 2) Count positives and negatives (ambiguous still included)
    n_pos = (df["label"] == "Positive").sum()
    n_neg = (df["label"] == "Negative").sum()

    # 3) Ambiguous positions: same (leading_protein, Protein_pos) having both pos & neg
    grouped = df.groupby(["leading_protein", "Protein_pos"])
    ambiguous = grouped["label"].agg(lambda x: {"Positive", "Negative"}.issubset(set(x)))
    n_ambiguous = ambiguous.sum()

    # 4) Print summary only
    print("=== Summary after deduplication ===")
    print(f"Deleted duplicates: {n_deleted_duplicates}")
    print(f"Positive sites:      {n_pos}")
    print(f"Negative sites:      {n_neg}")
    print(f"Ambiguous positions: {n_ambiguous}")

    # 5) Save deduplicated dataframe (including ambiguous rows)
    df.to_csv(out_path, index=False)


if __name__ == "__main__":
    main()
