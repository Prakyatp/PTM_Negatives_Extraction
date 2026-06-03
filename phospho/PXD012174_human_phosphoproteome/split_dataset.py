#!/usr/bin/env python3
import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser(
        description="Split dataset by Amino_acid (S/T/Y) and report label counts"
    )
    ap.add_argument("--inp", required=True, help="Input CSV")
    ap.add_argument("--out-prefix", default="dataset", help="Output file prefix")
    args = ap.parse_args()

    df = pd.read_csv(args.inp)

    required_cols = {"Amino_acid", "label"}
    missing = required_cols - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    # Normalize label just in case
    df["label"] = df["label"].astype(str).str.capitalize()

    for aa in ["S", "T", "Y"]:
        sub = df[df["Amino_acid"] == aa].copy()

        out_file = f"{args.out_prefix}_{aa}.csv"
        sub.to_csv(out_file, index=False)

        n_total = len(sub)
        n_pos = (sub["label"] == "Positive").sum()
        n_neg = (sub["label"] == "Negative").sum()

        print(f"\n=== Dataset for Amino Acid {aa} ===")
        print(f"Total rows:   {n_total}")
        print(f"Positive:     {n_pos}")
        print(f"Negative:     {n_neg}")
        print(f"Saved to:     {out_file}")


if __name__ == "__main__":
    main()
