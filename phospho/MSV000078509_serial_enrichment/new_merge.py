#!/usr/bin/env python3
import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Merge positive + negative site CSVs.\n"
            "Final output columns must be: leading_protein, Protein_pos, Sequence, Amino_acid, probability, label."
        )
    )
    ap.add_argument("--pos", required=True, help="Positive sites CSV")
    ap.add_argument("--neg", required=True, help="Negative sites CSV")
    ap.add_argument("--out", default="merged_sites.csv", help="Output merged CSV")
    args = ap.parse_args()

    # Load
    df_pos = pd.read_csv(args.pos)
    df_neg = pd.read_csv(args.neg)

    print(f"Loaded {len(df_pos)} positive rows")
    print(f"Loaded {len(df_neg)} negative rows")

    # Merge
    df = pd.concat([df_pos, df_neg], ignore_index=True)

    # Count pos/neg
    df["label"] = df["label"].astype(str).str.capitalize()
    n_pos = (df["label"] == "Positive").sum()
    n_neg = (df["label"] == "Negative").sum()

    # Save final
    df.to_csv(args.out, index=False)

    print("\n=== Merge Summary ===")
    print(f"Total rows:             {len(df)}")
    print(f"Positive sites:         {n_pos}")
    print(f"Negative sites:         {n_neg}")
    print(f"Final output saved to:  {args.out}")


if __name__ == "__main__":
    main()
