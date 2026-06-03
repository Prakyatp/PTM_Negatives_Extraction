#!/usr/bin/env python3
"""
Extract phospho site labels from MaxQuant evidence.txt + peptides.txt.

Outputs per-site rows with:
    pep_pos   — 1-based position of the S/T/Y within the peptide
    start_pos — 1-based start position of the peptide in the protein
                (from peptides.txt, used as confirmation in map.py)

The map script will:
    1. Find the peptide in the FASTA → get start_idx (0-based)
    2. Calculate Protein_pos = start_idx + pep_pos
    3. Verify start_idx + 1 == start_pos  (FASTA agrees with peptides.txt)
    4. Verify pep_seq[pep_pos - 1] == Amino_acid

Output columns:
    leading_protein | pep_pos | start_pos | Sequence | Amino_acid
    probability | label | Retention_time
"""

import argparse
import math
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


def parse_start_single(x):
    """Parse a single Start position value e.g. '123' or '123.0'."""
    try:
        return int(float(x))
    except Exception:
        return math.nan


def pick_start_for_protein(prot_str, start_str, lead_prot):
    """
    From peptides.txt:
        Proteins       : 'P12345;Q99999;O60341'
        Start position : '109;57;1'
    From evidence.txt:
        Leading protein: 'P12345'

    Finds the start position for the specific leading_protein.

    Handles isoforms:
        leading_protein 'P12345' matches 'P12345-2' in Proteins list
        by stripping isoform suffix before comparing.
    """
    if pd.isna(prot_str) or pd.isna(start_str) or pd.isna(lead_prot):
        return math.nan

    lead_prot  = str(lead_prot).strip()
    lead_base  = lead_prot.split("-")[0]   # strip isoform suffix

    pro_list   = [p.strip() for p in str(prot_str).split(";") if p.strip()]
    start_list = [s.strip() for s in str(start_str).split(";") if s.strip()]

    if not pro_list or not start_list:
        return math.nan

    # If only one start value given for multiple proteins, replicate it
    if len(start_list) == 1 and len(pro_list) > 1:
        start_list = start_list * len(pro_list)

    if len(pro_list) != len(start_list):
        return math.nan

    # Try exact match first, then base accession match (handles isoforms)
    for i, prot in enumerate(pro_list):
        if prot == lead_prot:
            return parse_start_single(start_list[i])

    for i, prot in enumerate(pro_list):
        if prot.split("-")[0] == lead_base:
            return parse_start_single(start_list[i])

    return math.nan


def parse_phospho_probs(prob_str, n_sites):
    """
    Parse 'Phospho (STY) Probabilities' as ';'-separated floats.
    Returns a list of length n_sites.

    Example:
        "0.991;0.009;0.251"  →  [0.991, 0.009, 0.251]
    """
    probs = [0.0] * n_sites
    if isinstance(prob_str, str) and prob_str.strip():
        parts = [p.strip() for p in prob_str.split(";")]
        vals  = []
        for p in parts:
            if p in ("", "NaN", "nan"):
                continue
            try:
                vals.append(float(p))
            except ValueError:
                vals.append(0.0)
        for i in range(min(len(vals), n_sites)):
            probs[i] = vals[i]
    return probs


def get_phospho_sites_from_modified(mod_seq):
    """
    Return a set of 1-based peptide positions annotated with '(ph)'.

    Handles MaxQuant Modified sequence formats:
        _AAAAAAS(ph)AASAGGK_
        R.AAAAAAS(ph)AASAGGK.K
        AAAAAAS(ph)AASAGGK
    """
    if pd.isna(mod_seq):
        return set()

    s = str(mod_seq)

    # Strip flanking protein context  'K.PEPTIDE.R' → 'PEPTIDE'
    if s.count(".") == 2:
        s = s.split(".")[1]
    s = s.strip("_")

    pos          = 0
    ph_positions = set()
    i            = 0

    while i < len(s):
        ch = s[i]
        if ch.isalpha() and ch.isupper():
            pos += 1
            j = i + 1
            if j < len(s) and s[j] == "(":
                k = s.find(")", j + 1)
                if k != -1:
                    annot = s[j + 1 : k]
                    if "ph" in annot.lower():
                        ph_positions.add(pos)
                    i = k + 1
                    continue
            i += 1
        else:
            i += 1

    return ph_positions


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Extract per-site phospho labels from MaxQuant evidence.txt.\n"
            "Outputs pep_pos (1-based offset within peptide) and\n"
            "start_pos (1-based peptide start in protein from peptides.txt).\n"
            "Both are used in map.py for double verification."
        )
    )
    ap.add_argument("--evidence",  required=True, help="evidence.txt from MaxQuant")
    ap.add_argument("--peptides",  required=True, help="peptides.txt from MaxQuant")
    ap.add_argument("--out",       default="1_extracted.csv")
    ap.add_argument("--chunksize", type=int, default=50000)
    args = ap.parse_args()

    # ── Load peptides.txt (only the 3 columns we need) ────────────────────────
    print(f"Reading peptides: {args.peptides}")
    pep_head    = pd.read_csv(args.peptides, sep="\t", nrows=0)
    c_seq_pep   = colfind(pep_head, "Sequence")
    c_prot_pep  = colfind(pep_head, "Proteins")
    c_start_pep = colfind(pep_head, "Start position")

    pep = pd.read_csv(
        args.peptides,
        sep="\t",
        usecols=[c_seq_pep, c_prot_pep, c_start_pep],
        low_memory=False,
    )
    pep = pep.dropna(subset=[c_seq_pep, c_prot_pep, c_start_pep])
    print(f"  Peptide rows loaded: {len(pep):,}")

    # ── Inspect evidence header ───────────────────────────────────────────────
    print(f"Reading evidence header: {args.evidence}")
    ev_head    = pd.read_csv(args.evidence, sep="\t", nrows=0)
    c_seq      = colfind(ev_head, "Sequence")
    c_modseq   = colfind(ev_head, "Modified sequence")
    c_ph_probs = colfind(ev_head, "Phospho (STY) Probabilities")
    c_leadprot = colfind(ev_head, "Leading proteins")
    c_rt       = colfind(ev_head, "Retention time")
    c_pep      = colfind(ev_head, "PEP")
    c_rev      = colfind(ev_head, "Reverse")
    c_cont     = colfind(ev_head, "Potential contaminant")

    ev_usecols = [
        c_seq, c_modseq, c_ph_probs,
        c_leadprot, c_rt, c_pep, c_rev, c_cont,
    ]

    # ── Stream evidence in chunks ─────────────────────────────────────────────
    first_out       = True
    total_site_rows = 0
    chunk_index     = 0

    print(f"\nStreaming evidence in chunks of {args.chunksize:,} ...")

    for chunk in pd.read_csv(
        args.evidence,
        sep="\t",
        usecols=ev_usecols,
        chunksize=args.chunksize,
        low_memory=False,
    ):
        chunk_index += 1
        print(f"\n--- Chunk {chunk_index} | rows before filters: {len(chunk):,} ---")

        # Basic quality filters
        chunk = chunk[chunk[c_rev].fillna("") != "+"]
        chunk = chunk[chunk[c_cont].fillna("") != "+"]
        chunk = chunk[chunk[c_pep].fillna(1.0) < 0.99]
        chunk = chunk[~chunk[c_leadprot].fillna("").astype(str).str.contains(";")]
        print(f"  After filters: {len(chunk):,}")

        if len(chunk) == 0:
            continue

        # Merge with peptides on Sequence to get Proteins + Start position
        merged = chunk.merge(
            pep,
            left_on=c_seq,
            right_on=c_seq_pep,
            how="left",
            suffixes=("", "_pep"),
        )

        # Compute start_pos for this specific leading_protein
        # This now handles isoforms correctly via pick_start_for_protein
        merged["start_pos"] = merged.apply(
            lambda r: pick_start_for_protein(
                r[c_prot_pep],
                r[c_start_pep],
                r[c_leadprot],
            ),
            axis=1,
        )

        # Drop rows where start_pos couldn't be determined
        merged = merged.dropna(subset=["start_pos"])
        merged["start_pos"] = merged["start_pos"].astype(int)
        print(f"  After start_pos join: {len(merged):,}")

        if len(merged) == 0:
            continue

        # Expand to site level
        out_rows = []

        for _, row in merged.iterrows():
            seq        = str(row[c_seq]).strip()
            mod_seq    = row[c_modseq]
            lead_prot  = str(row[c_leadprot]).strip()
            rt         = row[c_rt]
            start_pos  = int(row["start_pos"])
            ph_prob_str = row.get(c_ph_probs, None)

            sty_positions = []
            sty_aas       = []
            for i, aa in enumerate(seq, start=1):
                if aa in "STY":
                    sty_positions.append(i)
                    sty_aas.append(aa)

            if not sty_positions:
                continue

            site_probs = parse_phospho_probs(ph_prob_str, len(sty_positions))
            ph_sites   = get_phospho_sites_from_modified(mod_seq)

            for site_idx, (pep_pos, aa) in enumerate(zip(sty_positions, sty_aas)):

                if pep_pos in ph_sites:
                    label = "Positive"
                    prob  = site_probs[site_idx]
                else:
                    label = "Negative"
                    prob  = 0.0

                out_rows.append({
                    "leading_protein" : lead_prot,
                    "pep_pos"         : pep_pos,    # 1-based offset in peptide
                    "start_pos"       : start_pos,  # 1-based start in protein (from peptides.txt)
                    "Sequence"        : seq,
                    "Amino_acid"      : aa,
                    "probability"     : prob,
                    "label"           : label,
                    "Retention_time"  : rt,
                })

        if not out_rows:
            continue

        out_df = pd.DataFrame(out_rows)
        total_site_rows += len(out_df)
        print(f"  Site rows this chunk: {len(out_df):,}  (total: {total_site_rows:,})")

        if first_out:
            out_df.to_csv(args.out, index=False, mode="w")
            first_out = False
        else:
            out_df.to_csv(args.out, index=False, mode="a", header=False)

    print(f"\n✅ Done.")
    print(f"   Total site rows : {total_site_rows:,}")
    print(f"   Output file     : {args.out}")


if __name__ == "__main__":
    main()
