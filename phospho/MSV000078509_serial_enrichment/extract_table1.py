#!/usr/bin/env python3
"""
Extract phospho STY sites from MaxQuant Phospho(STY)Sites XLSX.

Outputs per-site rows with:
    Protein_pos — taken directly from the Position column (no calculation)
    pep_pos     — 1-based position of the S/T/Y within the peptide sequence

map.py will:
    1. Find the peptide in the FASTA → start_idx
    2. Verify: start_idx + pep_pos == Protein_pos
    3. Verify: pep_seq[pep_pos - 1] == Amino_acid
    Both must pass → KEEP

Output columns:
    leading_protein | Protein_pos | pep_pos | Sequence
    Amino_acid | probability | label
"""

import argparse
import pandas as pd


def colfind(df, *cands, required=True):
    """Case-insensitive column finder."""
    low = {c.lower(): c for c in df.columns}
    for name in cands:
        if name and name.lower() in low:
            return low[name.lower()]
    if required:
        raise KeyError(f"Missing required column. Tried: {cands}")
    return None


def parse_prob_string(s):
    """
    Parse the 'Phospho (STY) Probabilities' string.

    Example:
        'INS(1)APS(0.251)S(0.749)PIKTNK'
            seq   = 'INSAPSSPIKTNK'
            probs = {3: 1.0, 6: 0.251, 7: 0.749}
    """
    if pd.isna(s):
        return None, {}

    s     = str(s).strip()
    seq   = []
    probs = {}
    i     = 0
    pos   = 0   # 1-based peptide position

    while i < len(s):
        ch = s[i]
        if ch.isalpha():
            pos += 1
            seq.append(ch)
            i += 1
            if i < len(s) and s[i] == "(":
                j = s.find(")", i + 1)
                if j == -1:
                    break
                try:
                    probs[pos] = float(s[i + 1 : j])
                except ValueError:
                    pass
                i = j + 1
        else:
            i += 1

    return "".join(seq), probs


def count_ph(mod_seq):
    """Count '(ph)' occurrences in the Modified sequence."""
    if pd.isna(mod_seq):
        return 0
    return str(mod_seq).lower().count("(ph)")


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Extract per-site phospho labels from MaxQuant Phospho(STY)Sites XLSX.\n"
            "Outputs Protein_pos (from Position column, unchanged) and\n"
            "pep_pos (1-based offset of the STY within the peptide).\n"
            "map.py verifies both against the FASTA."
        )
    )
    ap.add_argument("--xlsx",  required=True)
    ap.add_argument("--sheet", default=None)
    ap.add_argument("--out",   default="1_extracted.csv")
    args = ap.parse_args()

    df = pd.read_excel(args.xlsx, sheet_name=args.sheet)
    print(f"Loaded {len(df):,} rows from {args.xlsx}")

    c_leading = colfind(df, "Leading proteins")
    c_pos     = colfind(df, "Position")
    c_mod     = colfind(df, "Modified sequence")
    c_prob    = colfind(df, "Phospho (STY) Probabilities")

    records = []
    n_skipped_multi   = 0
    n_skipped_noph    = 0
    n_skipped_noseq   = 0
    n_skipped_nopos   = 0
    n_skipped_multiph = 0

    for _, row in df.iterrows():
        leading = str(row[c_leading]).strip()

        # Only single leading protein
        if ";" in leading:
            n_skipped_multi += 1
            continue

        # Must have at least one phospho
        n_ph = count_ph(row[c_mod])
        if n_ph == 0:
            n_skipped_noph += 1
            continue

        # Parse probability string → clean sequence + per-position probs
        seq, pos_probs = parse_prob_string(row[c_prob])
        if not seq:
            n_skipped_noseq += 1
            continue

        # All STY positions in the peptide (1-based)
        sty_positions = [i for i, aa in enumerate(seq, start=1) if aa in "STY"]
        if not sty_positions:
            continue

        site_probs     = {p: pos_probs.get(p, 0.0) for p in sty_positions}
        positive_sites  = [p for p in sty_positions if site_probs[p] >= 0.99]
        ambiguous_sites = [p for p in sty_positions if 0.01 < site_probs[p] < 0.99]

        # Multi-ph rule 1: all ambiguous, no positives → drop entire peptide
        if n_ph >= 2 and ambiguous_sites and not positive_sites:
            n_skipped_multiph += 1
            continue

        # Multi-ph rule 2: both positives AND ambiguous → suppress negatives
        suppress_negatives = (n_ph >= 2 and ambiguous_sites and positive_sites)

        # Get Protein_pos directly from the Position column — no calculation
        try:
            protein_pos = int(row[c_pos])
        except (TypeError, ValueError):
            n_skipped_nopos += 1
            continue

        # The Position column gives the protein position of the anchor site.
        # Anchor = first positive if any, else first STY.
        pep_pos_anchor = positive_sites[0] if positive_sites else sty_positions[0]

        # Emit one row per Positive or Negative site
        for pep_pos in sty_positions:
            aa = seq[pep_pos - 1]
            p  = site_probs[pep_pos]

            if p >= 0.99:
                label = "Positive"
            elif p <= 0.01:
                if suppress_negatives:
                    continue
                label = "Negative"
            else:
                # Ambiguous → skip
                continue

            # Protein_pos for THIS site, calculated from the anchor
            # anchor_protein_pos is for pep_pos_anchor
            # so this site's protein pos = protein_pos + (pep_pos - pep_pos_anchor)
            site_protein_pos = protein_pos + (pep_pos - pep_pos_anchor)

            records.append({
                "leading_protein" : leading,
                "Protein_pos"     : site_protein_pos,  # from Position column, unchanged
                "pep_pos"         : pep_pos,            # 1-based offset within peptide
                "Sequence"        : seq,
                "Amino_acid"      : aa,
                "probability"     : p,
                "label"           : label,
            })

    out_df = pd.DataFrame.from_records(records)
    out_df.to_csv(args.out, index=False)

    print(f"\n=== Extraction summary ===")
    print(f"Input rows                        : {len(df):,}")
    print(f"Skipped — multiple proteins       : {n_skipped_multi:,}")
    print(f"Skipped — no phospho              : {n_skipped_noph:,}")
    print(f"Skipped — no sequence parsed      : {n_skipped_noseq:,}")
    print(f"Skipped — no valid Position       : {n_skipped_nopos:,}")
    print(f"Skipped — multi-ph all ambiguous  : {n_skipped_multiph:,}")
    print(f"Output site rows                  : {len(out_df):,}")
    print(f"✅ Wrote: {args.out}")


if __name__ == "__main__":
    main()
