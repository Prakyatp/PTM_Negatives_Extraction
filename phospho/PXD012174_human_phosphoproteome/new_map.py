#!/usr/bin/env python3
"""
Map phospho sites to canonical UniProt sequences.

Input CSV must have:
    leading_protein | Protein_pos | pep_pos | Sequence | Amino_acid | ...

Verification logic per unique (canonical_acc, Sequence, Amino_acid, pep_pos, Protein_pos):

    1. Find all occurrences of peptide in the canonical protein → start_idx (0-based)
    2. Check AA:      pep_seq[pep_pos - 1] == Amino_acid
    3. Check pos:     start_idx + pep_pos   == Protein_pos
    4. Check tryptic: start_idx == 0  OR  prot_seq[start_idx - 1] in {K, R}
    All must pass → KEEP (Protein_pos unchanged)

Drop reasons:
    absent       — peptide not found in protein
    aa_mismatch  — letter at pep_pos doesn't match Amino_acid
    pos_mismatch — FASTA-derived position doesn't match given Protein_pos
    not_tryptic  — position-verified match exists but is not preceded by K/R
    ambiguous    — multiple occurrences pass all checks
"""

import argparse
import gzip
import math
import pandas as pd
from multiprocessing import Pool, cpu_count


# ── FASTA loading ─────────────────────────────────────────────────────────────

def load_uniprot_fasta(path):
    acc_to_seq  = {}
    name_to_seq = {}
    name_to_acc = {}

    opener = gzip.open if path.endswith(".gz") else open

    with opener(path, "rt") as fh:
        cur_acc = cur_name = None
        cur_seq = []

        def _save():
            if not cur_seq:
                return
            seq = "".join(cur_seq)
            if cur_acc:
                acc_to_seq[cur_acc] = seq
            if cur_name:
                name_to_seq[cur_name] = seq
                if cur_acc:
                    name_to_acc[cur_name] = cur_acc

        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                _save()
                parts = line[1:].split()[0].split("|")
                if len(parts) >= 3:
                    cur_acc, cur_name = parts[1], parts[2]
                elif len(parts) == 2:
                    cur_acc, cur_name = parts[1], None
                else:
                    cur_acc, cur_name = None, parts[0]
                cur_seq = []
            else:
                cur_seq.append(line)
        _save()

    return acc_to_seq, name_to_seq, name_to_acc


# ── ID normalisation ──────────────────────────────────────────────────────────

def to_canonical_acc(pid, acc_to_seq, name_to_seq, name_to_acc):
    if pid is None or (isinstance(pid, float) and math.isnan(pid)):
        return None, None
    pid = str(pid).strip()
    if not pid:
        return None, None

    base = pid.split("-")[0] if "-" in pid else None

    for cand in dict.fromkeys(filter(None, [
        pid,
        base,
        None if pid.endswith("_HUMAN") else pid + "_HUMAN",
        None if (not base or base.endswith("_HUMAN")) else base + "_HUMAN",
    ])):
        if cand in acc_to_seq:
            return cand, acc_to_seq[cand]
        if cand in name_to_seq:
            acc = name_to_acc.get(cand, cand)
            return acc, name_to_seq[cand]

    return None, None


# ── Tryptic check ─────────────────────────────────────────────────────────────

def is_tryptic(start_idx, prot_seq):
    """
    Return True if the peptide starting at start_idx (0-based) is tryptic:
      - starts at the very beginning of the protein (N-terminal peptide), OR
      - the residue immediately before it is K or R.
    """
    return start_idx == 0 or prot_seq[start_idx - 1] in ("K", "R")


# ── Core verification ─────────────────────────────────────────────────────────

def verify_combo(acc, pep_seq, aa, pep_pos, protein_pos, acc_seq_map):
    """
    Verify a unique (canonical_acc, Sequence, Amino_acid, pep_pos, Protein_pos) combo.

    Steps:
        1. Check pep_seq[pep_pos - 1] == aa
           (the letter at pep_pos in the peptide must match Amino_acid)

        2. Find all occurrences of pep_seq in the protein.

        3. For each occurrence:
               derived_pos = start_idx + pep_pos
               Check A (pos):     derived_pos == Protein_pos
               Check B (tryptic): start_idx == 0 OR prot_seq[start_idx-1] in {K, R}

        4. Exactly 1 occurrence passes all checks → KEEP
           0 occurrences found in protein        → absent
           AA check fails                        → aa_mismatch
           Position check fails for all hits     → pos_mismatch
           Position check passes but none tryptic→ not_tryptic
           Multiple pass all checks              → ambiguous
    """
    prot_seq = acc_seq_map.get(acc, "")
    if not prot_seq:
        return "absent"

    pep_len = len(pep_seq)

    # Check 1 — AA at pep_pos in the peptide must match Amino_acid
    if pep_pos < 1 or pep_pos > pep_len:
        return "aa_mismatch"
    if pep_seq[pep_pos - 1] != aa:
        return "aa_mismatch"

    # Find all occurrences of the peptide in the protein
    pos_valid     = 0   # occurrences passing position check
    tryptic_valid = 0   # occurrences passing both position + tryptic check
    start         = 0
    found_any     = False

    while True:
        idx = prot_seq.find(pep_seq, start)
        if idx == -1:
            break
        found_any = True

        # Check 2A — does the FASTA-derived position match given Protein_pos?
        derived_pos = idx + pep_pos    # 0-based start + 1-based pep_pos = 1-based protein pos
        if derived_pos == protein_pos:
            pos_valid += 1

            # Check 2B — tryptic: preceded by K/R or at protein N-terminus
            if is_tryptic(idx, prot_seq):
                tryptic_valid += 1
                if tryptic_valid > 1:
                    return "ambiguous"

        start = idx + 1

    if not found_any:
        return "absent"
    if pos_valid == 0:
        return "pos_mismatch"
    if tryptic_valid == 0:
        return "not_tryptic"
    # tryptic_valid == 1 (ambiguous already returned above)
    return "keep"


# ── Multiprocessing worker ────────────────────────────────────────────────────

_worker_acc_seq_map = {}

def _init_worker(acc_seq_map):
    global _worker_acc_seq_map
    _worker_acc_seq_map = acc_seq_map


def _worker_task(args):
    acc, pep_seq, aa, pep_pos, protein_pos = args
    return verify_combo(acc, pep_seq, aa, pep_pos, protein_pos, _worker_acc_seq_map)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Verify phospho sites against canonical UniProt FASTA.\n"
            "Uses pep_pos and Protein_pos from extraction for double verification,\n"
            "plus tryptic check (preceding residue must be K or R)."
        )
    )
    ap.add_argument("--csv",       required=True)
    ap.add_argument("--fasta",     required=True)
    ap.add_argument("--out",       default="2_mapped.csv")
    ap.add_argument("--workers",   type=int,
                    default=max(1, cpu_count() - 1),
                    help="Parallel workers (default: all CPUs - 1)")
    ap.add_argument("--chunksize", type=int, default=2000)
    args = ap.parse_args()

    # ── Load ──────────────────────────────────────────────────────────────────
    print(f"Loading site table : {args.csv}")
    df = pd.read_csv(args.csv, low_memory=False)

    for c in ["leading_protein", "Protein_pos", "pep_pos", "Amino_acid", "Sequence"]:
        if c not in df.columns:
            raise KeyError(
                f"Missing required column: '{c}'\n"
                f"Available: {list(df.columns)}"
            )

    print(f"  Rows loaded      : {len(df):,}")

    print(f"Loading FASTA      : {args.fasta}")
    acc_to_seq, name_to_seq, name_to_acc = load_uniprot_fasta(args.fasta)
    print(f"  Accessions : {len(acc_to_seq):,}")
    print(f"  Entry names: {len(name_to_seq):,}")

    # ═════════════════════════════════════════════════════════════════════════
    # PASS 1 — Resolve protein IDs to canonical accessions
    # ═════════════════════════════════════════════════════════════════════════
    print("\nPass 1: resolving protein IDs ...")

    unique_ids = df["leading_protein"].unique()
    print(f"  Unique protein IDs   : {len(unique_ids):,}")

    id_map = {
        pid: to_canonical_acc(pid, acc_to_seq, name_to_seq, name_to_acc)
        for pid in unique_ids
    }

    acc_seq_map = {
        acc: seq for acc, seq in id_map.values() if acc is not None
    }

    n_unmapped = sum(1 for v in id_map.values() if v[0] is None)
    n_isoforms = sum(
        1 for pid in unique_ids
        if "-" in str(pid) and id_map[pid][0] is not None
    )
    print(f"  Not found in FASTA   : {n_unmapped:,}")
    print(f"  Isoforms → canonical : {n_isoforms:,}")

    df["_was_isoform"]   = df["leading_protein"].str.contains("-", na=False)
    df["_canonical_acc"] = df["leading_protein"].map(lambda p: id_map[p][0])

    n_before     = len(df)
    df           = df[df["_canonical_acc"].notna()].copy()
    n_dropped_p1 = n_before - len(df)
    print(f"  Rows dropped (not in FASTA): {n_dropped_p1:,}")
    print(f"  Rows carried forward       : {len(df):,}")

    df["leading_protein"] = df["_canonical_acc"]

    # Validate Protein_pos and pep_pos are numeric
    df["_protein_pos_int"] = pd.to_numeric(df["Protein_pos"], errors="coerce")
    df["_pep_pos_int"]     = pd.to_numeric(df["pep_pos"],     errors="coerce")

    n_bad = int((df["_protein_pos_int"].isna() | df["_pep_pos_int"].isna()).sum())
    if n_bad:
        print(f"  Rows with invalid Protein_pos/pep_pos dropped: {n_bad:,}")
        df = df[df["_protein_pos_int"].notna() & df["_pep_pos_int"].notna()].copy()

    df["_protein_pos_int"] = df["_protein_pos_int"].astype(int)
    df["_pep_pos_int"]     = df["_pep_pos_int"].astype(int)

    # ═════════════════════════════════════════════════════════════════════════
    # PASS 2 — Verification + tryptic filter
    #
    # Deduplicate on (acc, Sequence, Amino_acid, pep_pos, Protein_pos).
    # Run verify_combo() once per unique combo, merge back to all rows.
    # ═════════════════════════════════════════════════════════════════════════
    print("\nPass 2: verifying sites ...")

    combo_cols = [
        "_canonical_acc", "Sequence", "Amino_acid",
        "_pep_pos_int", "_protein_pos_int"
    ]
    combos = (
        df[combo_cols]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    print(f"  Unique combos to verify : {len(combos):,}")
    print(f"  Workers                 : {args.workers}")

    tasks = [
        (row["_canonical_acc"],
         str(row["Sequence"]),
         str(row["Amino_acid"]),
         int(row["_pep_pos_int"]),
         int(row["_protein_pos_int"]))
        for _, row in combos.iterrows()
    ]

    print("  Running verification ...")
    with Pool(
        processes=args.workers,
        initializer=_init_worker,
        initargs=(acc_seq_map,)
    ) as pool:
        results = pool.map(_worker_task, tasks, chunksize=args.chunksize)

    combos["_verdict"] = results
    df = df.merge(combos, on=combo_cols, how="left")

    # ── Counts ────────────────────────────────────────────────────────────────
    n_kept         = int((df["_verdict"] == "keep").sum())
    n_absent       = int((df["_verdict"] == "absent").sum())
    n_aa_mismatch  = int((df["_verdict"] == "aa_mismatch").sum())
    n_pos_mismatch = int((df["_verdict"] == "pos_mismatch").sum())
    n_not_tryptic  = int((df["_verdict"] == "not_tryptic").sum())
    n_ambiguous    = int((df["_verdict"] == "ambiguous").sum())
    n_iso_kept     = int(((df["_verdict"] == "keep") & df["_was_isoform"]).sum())

    # ── Write output ──────────────────────────────────────────────────────────
    drop_cols = [
        "_canonical_acc", "_was_isoform",
        "_protein_pos_int", "_pep_pos_int", "_verdict",
        "pep_pos",    # internal, not needed downstream
    ]

    out_df = df[df["_verdict"] == "keep"].copy()
    out_df["Protein_pos"] = out_df["_protein_pos_int"]

    fixed_cols = ["leading_protein", "Protein_pos", "Sequence", "Amino_acid"]
    remaining  = [
        c for c in out_df.columns
        if c not in fixed_cols and c not in drop_cols
    ]
    out_df = out_df[fixed_cols + remaining]
    out_df = out_df.drop(columns=[c for c in drop_cols if c in out_df.columns])
    out_df.to_csv(args.out, index=False)

    print("\n=== Summary ===")
    print(f"Total input rows                             : {n_before:,}")
    print(f"--- Pass 1 (ID normalisation) ---")
    print(f"  Dropped — protein not in FASTA             : {n_dropped_p1:,}")
    print(f"--- Pass 2 (verification + tryptic filter) ---")
    print(f"  Kept    — all checks passed (tryptic)      : {n_kept:,}")
    print(f"    of which were isoforms                   :   {n_iso_kept:,}")
    print(f"  Dropped — peptide not in protein           : {n_absent:,}")
    print(f"  Dropped — AA mismatch at pep_pos           : {n_aa_mismatch:,}")
    print(f"  Dropped — pos mismatch (FASTA vs extracted): {n_pos_mismatch:,}")
    print(f"  Dropped — not tryptic (no K/R before pep)  : {n_not_tryptic:,}")
    print(f"  Dropped — ambiguous (multiple tryptic hits): {n_ambiguous:,}")
    print(f"Total output rows                            : {len(out_df):,}")
    print(f"✅ Wrote: {args.out}")


if __name__ == "__main__":
    main()
