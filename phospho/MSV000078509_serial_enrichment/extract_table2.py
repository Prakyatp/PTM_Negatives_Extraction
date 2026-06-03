#!/usr/bin/env python3
import argparse
import gzip
import pandas as pd


# =========================
# FASTA loader
# =========================
def load_uniprot_fasta(path):
    """
    Load UniProt FASTA into a dict: {uniprot_id: sequence}
    Handles gzipped or plain FASTA.
    """
    seqs = {}
    opener = gzip.open if path.endswith(".gz") else open

    with opener(path, "rt") as fh:
        current_id = None
        current_seq = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(current_seq)
                parts = line[1:].split("|")
                if len(parts) >= 2 and parts[0] in ("sp", "tr"):
                    current_id = parts[1].strip()
                else:
                    current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            seqs[current_id] = "".join(current_seq)
    return seqs


# =========================
# Helpers
# =========================
def normalize_protein_id(pid):
    pid = str(pid).strip()
    if "-" in pid:
        pid = pid.split("-")[0]
    return pid


def pick_single_canonical_protein(proteins_str):
    if pd.isna(proteins_str):
        return None
    parts = [p.strip() for p in str(proteins_str).split(";") if p.strip()]
    canon_ids = {normalize_protein_id(p) for p in parts}
    if len(canon_ids) == 1:
        return list(canon_ids)[0]
    return None


def normalize_pepseq(seq):
    if pd.isna(seq):
        return None
    s = str(seq).strip().strip("_")
    letters = []
    in_paren = False
    for ch in s:
        if ch == "(":
            in_paren = True
            continue
        if ch == ")":
            in_paren = False
            continue
        if in_paren:
            continue
        if ch.isalpha():
            letters.append(ch.upper())
    return "".join(letters)


def find_all_occurrences(haystack, needle):
    starts = []
    start = 0
    while True:
        idx = haystack.find(needle, start)
        if idx == -1:
            break
        starts.append(idx)
        start = idx + 1
    return starts


# =========================
# Main
# =========================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--sheet", default=None)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", default="negative_STY_from_acetylomics.csv")
    args = ap.parse_args()

    print(f"Loading UniProt FASTA from: {args.fasta}")
    uniprot = load_uniprot_fasta(args.fasta)
    print(f"Loaded {len(uniprot)} UniProt entries")

    print(f"Loading XLSX from: {args.xlsx}")
    df = pd.read_excel(args.xlsx, sheet_name=args.sheet)
    print(f"Total input peptide rows: {len(df)}")

    col_pep = "Peptide sequence"
    col_prot = "Proteins"
    col_rt = "Retention time"

    for col in [col_pep, col_prot, col_rt]:
        if col not in df.columns:
            raise KeyError(f"Missing required column: {col}")

    # Protein filter
    df["canonical_protein"] = df[col_prot].apply(pick_single_canonical_protein)
    df = df[df["canonical_protein"].notna()].copy()
    print(f"After protein filter: {len(df)} rows")

    # Normalize peptide
    df["Sequence_norm"] = df[col_pep].apply(normalize_pepseq)
    df = df[df["Sequence_norm"].notna()].copy()
    print(f"After sequence normalization: {len(df)} rows")

    records = []
    n_missing_protein = n_no_match = n_multi_match = n_mapped = 0

    for _, row in df.iterrows():
        prot_id = row["canonical_protein"]
        pep_seq = row["Sequence_norm"]
        rt = row[col_rt]

        prot_seq = uniprot.get(prot_id)
        if prot_seq is None:
            n_missing_protein += 1
            continue

        starts = find_all_occurrences(prot_seq.upper(), pep_seq.upper())
        if not starts:
            n_no_match += 1
            continue
        if len(starts) > 1:
            n_multi_match += 1
            continue

        start_idx = starts[0]
        n_mapped += 1

        for pep_pos, aa in enumerate(pep_seq, start=1):
            if aa not in "STY":
                continue
            records.append(
                {
                    "leading_protein": prot_id,
                    "Protein_pos": start_idx + pep_pos,
                    "Sequence": pep_seq,
                    "Amino_acid": aa,
                    "Retention_time": rt,
                    "probability": 0.0,
                    "label": "Negative",
                }
            )

    out_df = pd.DataFrame(records)
    out_df = out_df[
        [
            "leading_protein",
            "Protein_pos",
            "Sequence",
            "Amino_acid",
            "Retention_time",
            "probability",
            "label",
        ]
    ]

    out_df.to_csv(args.out, index=False)

    print("\n=== Mapping summary ===")
    print(f"Mapped peptides uniquely: {n_mapped}")
    print(f"Missing proteins:         {n_missing_protein}")
    print(f"No match in protein:      {n_no_match}")
    print(f"Multiple matches:         {n_multi_match}")
    print(f"Total STY negatives:      {len(out_df)}")
    print(f"Output written to:        {args.out}")


if __name__ == "__main__":
    main()
