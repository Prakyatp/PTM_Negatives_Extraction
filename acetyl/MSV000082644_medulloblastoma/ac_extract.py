#!/usr/bin/env python3
"""
ac_extract.py
Pipeline:
    1. Read Acetylomics-acK sheet from mmc3.xlsx (add_scans)
    2. Remove peptides with >1 accession_number (delete_acnu)
    3. Label peptides Positive / Negative, save _pep_kidx (K offset in peptide)
    4. Map to UniProt FASTA using strict tryptic verification + _pep_kidx
       - peptide must have exactly 1 tryptic hit in exactly 1 protein
       - if peptide appears in ANY other protein (even non-tryptically) → dropped
"""

import re, gzip, argparse
import pandas as pd
import ahocorasick

# ── Sample columns ─────────────────────────────────────────────────────────────
SAMPLE_COLS = [
    "MB278","MB282","MB275","MB247","MB091","MB234","MB199","MB136","MB288",
    "Burdenko_1360","MB226","MB118","MB166","MB239","MB102","MB271","MB177",
    "MB037","MB265","MB170","MB244","MB260","MB274","MB264","MB268","MB095",
    "MB227","MB269","MB088","MB281","MB284","MB174","MB277","MB164","x5M15",
    "MB287","MB099","MB266","MB104","MB248","x1M6","MB106","MB270","MB018",
    "MB206",
]

# ── Step 1: add_scans ──────────────────────────────────────────────────────────

def colfind(df, *cands):
    low = {c.lower(): c for c in df.columns}
    for cand in cands:
        if cand and cand.lower() in low:
            return low[cand.lower()]
    raise KeyError(f"Could not find any of {cands} in columns: {list(df.columns)}")

def add_scans(xlsx_path, sheet):
    df = pd.read_excel(xlsx_path, sheet_name=sheet, header=1)
    seq_col = colfind(df, "sequenceVariablyModifiedLocations", "sequenceVML")
    meta_actual = [
        seq_col,
        colfind(df, "variableSites"),
        colfind(df, "accession_number"),
        colfind(df, "geneSymbol", "GeneSymbol"),
        colfind(df, "entry_name"),
        colfind(df, "accessionNumber_VMsites_numVMsitesPresent_numVMsitesLocalizedBest_earliestVMsiteAA_latestVMsiteAA"),
    ]
    sample_actual = [c for c in SAMPLE_COLS if c in df.columns]
    df["Number_of_scans"] = df[sample_actual].notna().sum(axis=1)
    out = df[meta_actual + ["Number_of_scans"]].copy()
    out = out.rename(columns={seq_col: "sequenceVariablyModifiedLocations"})
    print(f"[add_scans]   Rows loaded: {len(out):,}")
    return out

# ── Step 2: delete_acnu ────────────────────────────────────────────────────────

def delete_acnu(df):
    pep_col = "sequenceVariablyModifiedLocations"
    acc_col = "accession_number"
    df = df.copy()
    df["_acc_list"] = (
        df[acc_col].astype(str).str.split(";")
        .apply(lambda xs: [x.strip() for x in xs if x.strip()])
    )
    pep_to_accs = df.explode("_acc_list").groupby(pep_col)["_acc_list"].nunique()
    multi_acc   = pep_to_accs[pep_to_accs > 1].index
    bad_mask    = df[pep_col].isin(multi_acc)
    df_clean    = df[~bad_mask].drop(columns=["_acc_list"]).copy()
    print(f"[delete_acnu] Peptides with >1 accession: {len(multi_acc):,}")
    print(f"[delete_acnu] Rows removed: {bad_mask.sum():,}  |  Remaining: {len(df_clean):,}")
    return df_clean

# ── Step 3: label + save _pep_kidx ────────────────────────────────────────────

def parse_peptide(seq):
    if pd.isna(seq):
        return "", []
    s = str(seq).strip()
    plain, ks = [], []
    i = pos = 0
    while i < len(s):
        ch = s[i]
        if ch.isalpha():
            aa = ch
            plain.append(aa)
            prob = None
            i += 1
            if i < len(s) and s[i] == "(":
                j = s.find(")", i + 1)
                if j != -1:
                    try: prob = float(s[i+1:j])
                    except: prob = None
                    i = j + 1
            if aa == "K":
                ks.append({"pos": pos, "prob": prob})
            pos += 1
        else:
            i += 1
    return "".join(plain), ks

def rescued_terminal_K(short_plain, k_pos, all_plain, all_ks):
    Ls = len(short_plain)
    for big_plain, ks_big in zip(all_plain, all_ks):
        if len(big_plain) > Ls and big_plain.startswith(short_plain):
            for site in ks_big:
                if site["pos"] == k_pos and site["prob"] is not None and site["prob"] >= 0.9:
                    if site["pos"] < len(big_plain) - 1:
                        return True
    return False

def label_peptides(df):
    """
    Label each row Positive/Negative and save _pep_kidx:
      Positive: 0-based index of the acetylated K (prob >= 0.9, non-terminal)
      Negative: 0-based index of the unmodified K (prob == 0.0, non-terminal)
    """
    seq_col = "sequenceVariablyModifiedLocations"
    all_plain, all_ks = [], []
    for s in df[seq_col]:
        p, ks = parse_peptide(s)
        all_plain.append(p)
        all_ks.append(ks)

    labels, pep_kidxs = [], []
    pos_count = neg_count = 0

    for i, (idx, row) in enumerate(df.iterrows()):
        plain = all_plain[i]
        ks    = all_ks[i]
        L     = len(plain)

        found_pos_kidx = None
        found_neg_kidx = None

        for site in ks:
            prob = site["prob"]
            if prob is None:
                continue
            pos0        = site["pos"]
            is_terminal = (pos0 == L - 1)

            if prob >= 0.9 and not is_terminal:
                if found_pos_kidx is None:
                    found_pos_kidx = pos0
            elif prob == 0.0:
                if not is_terminal:
                    if found_neg_kidx is None:
                        found_neg_kidx = pos0
                else:
                    if rescued_terminal_K(plain, pos0, all_plain, all_ks):
                        if found_neg_kidx is None:
                            found_neg_kidx = pos0

        if found_pos_kidx is not None:
            labels.append("Positive")
            pep_kidxs.append(found_pos_kidx)
            pos_count += 1
        else:
            labels.append("Negative")
            pep_kidxs.append(found_neg_kidx)
            neg_count += 1

    df = df.copy()
    df["label"]     = labels
    df["_pep_kidx"] = pep_kidxs
    print(f"[label]       Positive: {pos_count:,}  |  Negative: {neg_count:,}")
    return df

# ── Step 4: FASTA mapping with strict tryptic verification ────────────────────

def clean_peptide(seq):
    if pd.isna(seq): return None
    s = re.sub(r'\(\d+(\.\d+)?\)', '', str(seq))
    s = s.replace('_', '').strip().upper()
    return s if s else None

def build_automaton(peptides):
    A = ahocorasick.Automaton()
    for idx, pep in enumerate(peptides):
        if pep: A.add_word(pep, (idx, pep))
    A.make_automaton()
    return A

def scan_fasta(fasta_path, automaton):
    matches, proteins = {}, {}
    opener = gzip.open if fasta_path.endswith(".gz") else open
    with opener(fasta_path, "rt") as fh:
        acc, seq = None, []
        def _flush():
            if acc is None or not seq: return
            protein = "".join(seq)
            proteins[acc] = protein
            for end_idx, (idx, pep) in automaton.iter(protein):
                start = end_idx - len(pep) + 1
                if pep not in matches: matches[pep] = []
                matches[pep].append((acc, start))
        for line in fh:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                _flush()
                parts = line[1:].split("|")
                acc = parts[1] if len(parts) >= 2 and parts[0] in ("sp", "tr") else line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        _flush()
    return matches, proteins

def get_tryptic_hits(pep, matches, proteins):
    """
    Strict tryptic filter:
      1. Find all FASTA hits (tryptic or not)
      2. Find tryptic hits (start == 0 OR protein[start-1] in K/R)
      3. If peptide appears in ANY protein beyond the single tryptic one
         (even non-tryptically) → return empty → dropped as ambiguous
      4. If exactly 1 tryptic hit and no other protein hits → return it
    """
    hits = matches.get(pep, [])

    # All proteins this peptide appears in (tryptic or not)
    all_proteins = {acc for acc, _ in hits}

    # Tryptic hits only
    tryptic = [
        (acc, start) for acc, start in hits
        if start == 0 or proteins[acc][start - 1] in ("K", "R")
    ]
    tryptic_proteins = {acc for acc, _ in tryptic}

    # If peptide exists in any protein beyond the single tryptic one → ambiguous
    if len(tryptic_proteins) == 1 and len(all_proteins) > 1:
        return []   # drop — non-tryptic hit in another protein makes it ambiguous

    return tryptic

def map_to_fasta(df, fasta_path, mismatch_out=None):
    df = df.copy()
    df["_clean"] = df["sequenceVariablyModifiedLocations"].apply(clean_peptide)

    # Drop rows with no valid K offset
    no_kidx = df["_pep_kidx"].isna().sum()
    if no_kidx > 0:
        print(f"[map]         Dropped (no valid K offset): {no_kidx:,}")
        df = df[df["_pep_kidx"].notna()].copy()
    df["_pep_kidx"] = df["_pep_kidx"].astype(int)

    unique_peps = [p for p in df["_clean"].unique() if p]
    print(f"[map]         Unique peptides: {len(unique_peps):,}")

    A = build_automaton(unique_peps)
    print(f"[map]         Scanning FASTA ...")
    matches, proteins = scan_fasta(fasta_path, A)

    # Pre-compute strict tryptic hits per unique peptide
    tryptic_map = {pep: get_tryptic_hits(pep, matches, proteins) for pep in unique_peps}

    # Peptide-level stats
    n_absent     = sum(1 for p in unique_peps if p not in matches)
    n_no_tryptic = sum(1 for p in unique_peps if p in matches and len(tryptic_map[p]) == 0)
    n_multi_prot = sum(1 for p in unique_peps if len({a for a, _ in tryptic_map[p]}) > 1)
    n_multi_site = sum(1 for p in unique_peps
                       if len({a for a, _ in tryptic_map[p]}) == 1 and len(tryptic_map[p]) > 1)
    n_good       = sum(1 for p in unique_peps if len(tryptic_map[p]) == 1)

    print(f"\n[map] Peptide-level tryptic stats:")
    print(f"  Absent                         : {n_absent:,}")
    print(f"  No tryptic / ambiguous in FASTA: {n_no_tryptic:,}  → dropped")
    print(f"  Multi-protein tryptic          : {n_multi_prot:,}  → dropped")
    print(f"  Multi-site (1 protein)         : {n_multi_site:,}  → dropped")
    print(f"  Unique tryptic hit             : {n_good:,}  → kept")

    verdicts, final_pos, final_acc = [], [], []

    for _, row in df.iterrows():
        pep  = row["_clean"]
        kidx = int(row["_pep_kidx"])

        if pep is None or pep not in matches:
            verdicts.append("absent");        final_pos.append(None); final_acc.append(None); continue

        tryptic_hits = tryptic_map[pep]

        if len(tryptic_hits) == 0:
            verdicts.append("no_tryptic");    final_pos.append(None); final_acc.append(None); continue

        distinct_proteins = {acc for acc, _ in tryptic_hits}
        if len(distinct_proteins) > 1:
            verdicts.append("multi_protein"); final_pos.append(None); final_acc.append(None); continue

        if len(tryptic_hits) > 1:
            verdicts.append("multi_site");    final_pos.append(None); final_acc.append(None); continue

        # Exactly one tryptic hit in exactly one protein, no other hits anywhere
        acc, start = tryptic_hits[0]
        prot_pos   = start + kidx + 1    # 1-based protein position of K
        verdicts.append("confirmed")
        final_pos.append(prot_pos)
        final_acc.append(acc)

    df["_verdict"]   = verdicts
    df["_final_pos"] = final_pos
    df["_final_acc"] = final_acc

    print(f"\n[map] Row-level Verification Summary:")
    for v in ["confirmed", "absent", "no_tryptic", "multi_protein", "multi_site"]:
        print(f"  {v:20s}: {(df['_verdict'] == v).sum():,}")

    # Save all dropped rows
    if mismatch_out:
        dropped_df = df[df["_verdict"] != "confirmed"].copy()
        dropped_df = dropped_df.drop(columns=["_clean", "_pep_kidx", "_verdict",
                                               "_final_pos", "_final_acc"])
        dropped_df.to_csv(mismatch_out, index=False)
        print(f"[dropped]     {len(dropped_df):,} rows saved to: {mismatch_out}")

    out_df = df[df["_verdict"] == "confirmed"].copy()
    out_df["variableSites"]      = out_df["_final_pos"].apply(lambda p: f"K{int(p)}k")
    out_df["mapped_position"]    = out_df["_final_pos"].astype(int)
    out_df["mapped_uniprot_acc"] = out_df["_final_acc"]
    out_df = out_df.drop(columns=["_clean", "_pep_kidx", "_verdict", "_final_pos", "_final_acc"])
    return out_df

# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="Acetylomics full pipeline: scans → delete_acnu → label → map")
    ap.add_argument("--input",        required=True,             help="mmc3.xlsx")
    ap.add_argument("--sheet",        default="Acetylomics-acK", help="Sheet name")
    ap.add_argument("--fasta",        required=True,             help="UniProt FASTA (.gz ok)")
    ap.add_argument("--output",       default="ac_mapped.csv")
    ap.add_argument("--mismatch-out", default="ac_dropped.csv",  help="CSV for all dropped rows")
    args = ap.parse_args()

    df = add_scans(args.input, args.sheet)
    df = delete_acnu(df)
    df = label_peptides(df)
    df = map_to_fasta(df, args.fasta, args.mismatch_out)

    df.to_csv(args.output, index=False)
    print(f"\nTotal output rows: {len(df):,}")
    print(f"Written to: {args.output}")

if __name__ == "__main__":
    main()
