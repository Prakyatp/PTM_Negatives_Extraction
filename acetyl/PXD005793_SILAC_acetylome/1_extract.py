#!/usr/bin/env python
import pandas as pd
import numpy as np
import re
import gzip
import argparse
from pathlib import Path
from multiprocessing import Pool, cpu_count


# ---------------------------------------------------------
# STEP 0 — Helper utilities
# ---------------------------------------------------------
def colfind(df, *cands, required=True):
    low = {c.lower(): c for c in df.columns}
    for name in cands:
        if name and name.lower() in low:
            return low[name.lower()]
    if required:
        raise KeyError("Missing required column. Tried: %s" % (cands,))
    return None

def is_plus(x):
    return str(x).strip() == '+' if pd.notna(x) else False

def is_missing_or_blank(x):
    if pd.isna(x):
        return True
    s = str(x).strip()
    return s == "" or s.lower() in {"na", "nan", "none", "null"}

def normalize_pepseq(seq):
    if pd.isna(seq):
        return None
    s = str(seq).strip()
    if s.count('.') == 2:
        s = s.split('.')[1]
    return s.strip('_').upper()

def parse_modified_sequence(modseq):
    if pd.isna(modseq):
        return {'all_k_idx': [], 'k_ac_idx': [], 'k_anymod_idx': [], 'residues': [], 'pep_len': 0}
    s = str(modseq)
    if s.count('.') == 2:
        s = s.split('.')[1]
    s = s.strip('_')

    all_k_idx    = []
    k_ac_idx     = []
    k_anymod_idx = []
    residues     = []
    i            = 0
    res_idx      = -1
    last_res     = None

    while i < len(s):
        ch = s[i]
        if ch.isalpha():
            res_idx += 1
            last_res = ch.upper()
            residues.append(last_res)
            if last_res == 'K':
                all_k_idx.append(res_idx)
            i += 1
        elif ch == '(':
            j = s.find(')', i + 1)
            if j == -1:
                break
            modtext = s[i+1:j].lower()
            if last_res == 'K':
                if modtext.strip():
                    k_anymod_idx.append(all_k_idx[-1])
                if 'acetyl' in modtext or re.search(r'\bac\b', modtext):
                    k_ac_idx.append(all_k_idx[-1])
            i = j + 1
        else:
            i += 1

    return {
        'all_k_idx'    : all_k_idx,
        'k_ac_idx'     : k_ac_idx,
        'k_anymod_idx' : k_anymod_idx,
        'residues'     : residues,
        'pep_len'      : len(residues),
    }

def derive_prefix_from_evidence(evidence_path):
    p = Path(evidence_path).resolve().parts
    if len(p) >= 3:
        if p[-2].lower() == 'txt' and len(p) >= 4:
            return "%s_%s" % (p[-4], p[-3])
        return "%s_%s" % (p[-3], p[-2])
    return Path(evidence_path).stem

def filter_single_leading_protein(ev_raw, c_leadprot):
    s            = ev_raw[c_leadprot].astype(str)
    missing_mask = ev_raw[c_leadprot].apply(is_missing_or_blank)
    multi_mask   = s.str.contains(";", na=False)
    keep_mask    = (~missing_mask) & (~multi_mask)
    out          = ev_raw[keep_mask].copy()
    return out, len(out), int((~keep_mask).sum())


# ---------------------------------------------------------
# STEP 1 — FASTA loading
# ---------------------------------------------------------
def load_uniprot_fasta(path):
    seqs   = {}
    opener = gzip.open if path.endswith(".gz") else open

    with opener(path, "rt") as fh:
        current_id  = None
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
                    current_id = parts[1]
                else:
                    current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            seqs[current_id] = "".join(current_seq)
    return seqs


# ---------------------------------------------------------
# STEP 2 — FASTA verification + tryptic check
# ---------------------------------------------------------
_worker_seqs = {}

def _init_worker(seqs):
    global _worker_seqs
    _worker_seqs = seqs


def is_tryptic(prot_seq: str, start_idx: int) -> bool:
    """
    Returns True if the peptide at 0-based start_idx is tryptic:
      - start_idx == 0  → N-terminus of protein, always valid
      - prot_seq[start_idx - 1] in (K, R) → preceded by trypsin cleavage residue
    """
    if start_idx == 0:
        return True
    return prot_seq[start_idx - 1].upper() in ('K', 'R')


def verify_site(args):
    """
    Verify K position, enforce tryptic rule, and remap start_pos if needed.

    Inputs:
        protein_id — canonical UniProt accession (isoform suffix already stripped)
        pep_seq    — plain peptide sequence (no mod annotations)
        start_pos  — 1-based peptide start from peptides.txt
        kidx       — 0-based index of K within the peptide

    Returns (verdict, fasta_pos, fasta_start):

      Kept verdicts:
        "keep"             — single locus, tryptic, K confirmed, start_pos matches
        "remapped"         — single locus, tryptic, K confirmed, start_pos corrected
        "rescued_tryptic"  — multiple loci, exactly one is tryptic, K confirmed

      Dropped verdicts:
        "absent"           — peptide not found in protein at all
        "not_tryptic"      — peptide found but no locus has K/R before it
        "multi_tryptic"    — 2+ tryptic loci, cannot resolve which is correct
        "mismatch_bad_kidx"   — pep_seq[kidx] is not K
        "mismatch_bad_k_prot" — protein residue at computed position is not K
    """
    protein_id, pep_seq, start_pos, kidx = args
    prot_seq = _worker_seqs.get(protein_id, "")

    if not prot_seq:
        return ("absent", None, None)

    # Find every position the peptide occurs in this protein
    starts = []
    pos    = 0
    while True:
        idx = prot_seq.find(pep_seq, pos)
        if idx == -1:
            break
        starts.append(idx)
        pos = idx + 1

    if not starts:
        return ("absent", None, None)

    # ---- Single locus ----
    if len(starts) == 1:
        start_idx = starts[0]

        # Tryptic check: residue at start_pos - 1 (0-based: start_idx - 1) must be K or R
        if not is_tryptic(prot_seq, start_idx):
            return ("not_tryptic", None, None)

        # K integrity checks — unrecoverable if these fail
        if pep_seq[kidx] != 'K':
            return ("mismatch_bad_kidx", None, None)
        if prot_seq[start_idx + kidx] != 'K':
            return ("mismatch_bad_k_prot", None, None)

        protein_pos = start_idx + kidx + 1   # 1-based
        fasta_start = start_idx + 1           # 1-based

        # start_pos from peptides.txt is wrong — remap to FASTA value
        if fasta_start != start_pos:
            return ("remapped", protein_pos, fasta_start)

        return ("keep", protein_pos, fasta_start)

    # ---- Multiple loci: keep only tryptic ones ----
    tryptic_starts = [s for s in starts if is_tryptic(prot_seq, s)]

    if len(tryptic_starts) == 0:
        return ("not_tryptic", None, None)

    if len(tryptic_starts) > 1:
        # Cannot resolve which tryptic locus the modification came from
        return ("multi_tryptic", None, None)

    # Exactly one tryptic locus — use it
    start_idx = tryptic_starts[0]

    if pep_seq[kidx] != 'K':
        return ("mismatch_bad_kidx", None, None)
    if prot_seq[start_idx + kidx] != 'K':
        return ("mismatch_bad_k_prot", None, None)

    protein_pos = start_idx + kidx + 1
    fasta_start = start_idx + 1
    return ("rescued_tryptic", protein_pos, fasta_start)


# ---------------------------------------------------------
# CLI
# ---------------------------------------------------------
ap = argparse.ArgumentParser(
    description="Extract Kac sites from MaxQuant evidence+peptides with "
                "FASTA position verification and tryptic filtering."
)
ap.add_argument('--evidence',   required=True)
ap.add_argument('--peptides',   required=True)
ap.add_argument('--fasta',      required=True)
ap.add_argument('--out-prefix', default=None)
ap.add_argument('--out-dir',    default='.')
ap.add_argument('--workers',    type=int, default=max(1, cpu_count() - 1))
args = ap.parse_args()

out_dir = Path(args.out_dir)
out_dir.mkdir(parents=True, exist_ok=True)
prefix  = args.out_prefix or derive_prefix_from_evidence(args.evidence)


# ---------------------------------------------------------
# STEP 3 — Load FASTA
# ---------------------------------------------------------
print(f"Loading FASTA: {args.fasta}")
fasta_seqs = load_uniprot_fasta(args.fasta)
print(f"  Loaded {len(fasta_seqs):,} proteins")


# ---------------------------------------------------------
# STEP 4 — Load and filter evidence.txt
# ---------------------------------------------------------
ev_raw = pd.read_csv(args.evidence, sep='\t', low_memory=False)
pep    = pd.read_csv(args.peptides, sep='\t', low_memory=False)

c_PEP      = colfind(ev_raw, 'PEP')
c_ModSeq   = colfind(ev_raw, 'Modified sequence')
c_Seq      = colfind(ev_raw, 'Sequence')
c_Proteins = colfind(ev_raw, 'Proteins')
c_UniqProt = colfind(ev_raw, 'Unique (Proteins)', required=False)
c_Reverse  = colfind(ev_raw, 'Reverse', required=False)
c_Contam   = colfind(ev_raw, 'Potential contaminant', 'Potential contaminant ')
c_NumScans = colfind(ev_raw, 'Number of scans')
c_Gene     = colfind(ev_raw, 'Gene names', required=False)
c_LeadProt = colfind(ev_raw, 'Leading proteins', 'Leading Proteins', 'Leading protein')

print("Total evidence rows: %d" % len(ev_raw))

ev_raw, kept_lp, removed_lp = filter_single_leading_protein(ev_raw, c_LeadProt)
print("After Leading proteins single-entry filter: %d (removed %d)" % (kept_lp, removed_lp))

ev_raw['protein_id'] = ev_raw[c_LeadProt].astype(str).str.strip()

ev = ev_raw[ev_raw[c_PEP] < 0.01].copy()
if c_Reverse:
    ev = ev[~ev[c_Reverse].apply(is_plus)]
ev = ev[~ev[c_Contam].apply(is_plus)]

def is_proteotypic(row):
    if c_UniqProt and pd.notna(row[c_UniqProt]):
        try:
            return int(row[c_UniqProt]) == 1
        except Exception:
            pass
    prots = [p for p in str(row[c_Proteins]).split(';') if p.strip()]
    return len(prots) == 1

ev = ev[ev.apply(is_proteotypic, axis=1)].copy()
ev['seq_norm'] = ev[c_Seq].map(normalize_pepseq)
ev['modseq']   = ev[c_ModSeq].astype(str)
print("After evidence filters (PEP, Reverse, contaminant, proteotypic): %d" % len(ev))


# ---------------------------------------------------------
# STEP 5 — Build peptide start_pos map from peptides.txt
# ---------------------------------------------------------
p_Seq      = colfind(pep, 'Sequence')
p_Start    = colfind(pep, 'Start position')
p_LRP      = colfind(pep, 'Leading razor protein', required=False)
p_LeadProt = colfind(pep, 'Leading proteins', 'Leading Proteins', 'Leading protein', required=False)
p_Proteins = colfind(pep, 'Proteins', required=False)

pep['seq_norm'] = pep[p_Seq].map(normalize_pepseq)

def pick_peptide_protein(row):
    if p_LeadProt and pd.notna(row[p_LeadProt]) and str(row[p_LeadProt]).strip():
        v = str(row[p_LeadProt]).strip()
        if ";" in v:
            return None
        return v
    if p_LRP and pd.notna(row[p_LRP]) and str(row[p_LRP]).strip():
        return str(row[p_LRP]).strip()
    if p_Proteins and pd.notna(row[p_Proteins]):
        plist = [p.strip() for p in str(row[p_Proteins]).split(';') if p.strip()]
        return plist[0] if len(plist) == 1 else None
    return None

pep['protein_id'] = pep.apply(pick_peptide_protein, axis=1)
pep['_start']     = pep[p_Start].apply(lambda x: int(x) if pd.notna(x) else np.nan)
pep_map           = pep[['seq_norm', 'protein_id', '_start']].dropna().drop_duplicates()

ev = ev.merge(pep_map, on=['seq_norm', 'protein_id'], how='inner')
ev = ev.rename(columns={'_start': 'start_pos'})
ev['start_pos'] = ev['start_pos'].astype(int)


# ---------------------------------------------------------
# STEP 6 — Parse modified sequences and label sites
# ---------------------------------------------------------
parsed = ev['modseq'].apply(parse_modified_sequence)
ev     = pd.concat([ev, parsed.apply(pd.Series)], axis=1)

ev['internal_k_idx'] = ev.apply(
    lambda r: [i for i in r['all_k_idx'] if i < r['pep_len'] - 1], axis=1
)
ev = ev[(ev['pep_len'] > 0) & (ev['internal_k_idx'].map(len) > 0)]

site_rows         = ev.explode('internal_k_idx').rename(columns={'internal_k_idx': 'kidx'})
site_rows['kidx'] = site_rows['kidx'].astype(int)
site_rows['protein_pos'] = site_rows['start_pos'].astype(int) + site_rows['kidx'].astype(int)

site_rows['next_res'] = site_rows.apply(
    lambda r: r['residues'][r['kidx'] + 1] if (r['kidx'] + 1) < r['pep_len'] else None, axis=1
)

site_rows['is_pos'] = site_rows.apply(lambda r: r['kidx'] in r['k_ac_idx'], axis=1)
site_rows['is_neg'] = site_rows.apply(
    lambda r: (r['kidx'] not in r['k_anymod_idx']) and (r['next_res'] != 'P'), axis=1
)

labeled                    = site_rows[(site_rows['is_pos']) | (site_rows['is_neg'])].copy()
labeled['label']           = np.where(labeled['is_pos'], 'Positive', 'Negative')
labeled['is_terminal']     = False
labeled['number_of_scans'] = labeled[c_NumScans].fillna(0).astype(int)


# ---------------------------------------------------------
# STEP 7 — Terminal K rescue
# ---------------------------------------------------------
ev['terminal_k_idx'] = ev.apply(
    lambda r: [r['pep_len'] - 1] if (r['pep_len'] > 0 and r['residues'][-1] == 'K') else [], axis=1
)
term_rows                    = ev.explode('terminal_k_idx').dropna(subset=['terminal_k_idx'])
term_rows                    = term_rows.rename(columns={'terminal_k_idx': 'kidx'})
term_rows['kidx']            = term_rows['kidx'].astype(int)
term_rows['protein_pos']     = term_rows['start_pos'].astype(int) + term_rows['kidx'].astype(int)
term_rows['is_neg_terminal'] = term_rows.apply(lambda r: r['kidx'] not in r['k_anymod_idx'], axis=1)
term_neg                     = term_rows[term_rows['is_neg_terminal']].copy()

pos_sites = labeled[labeled['label'] == 'Positive'][['protein_id', 'protein_pos']].drop_duplicates()
rescued   = term_neg.merge(pos_sites, on=['protein_id', 'protein_pos'], how='inner')
rescued['label']           = 'Negative'
rescued['is_terminal']     = True
rescued['number_of_scans'] = rescued[c_NumScans].fillna(0).astype(int)


# ---------------------------------------------------------
# STEP 8 — Gene name mapping
# ---------------------------------------------------------
if c_Gene and c_Gene in ev.columns:
    gene_map = ev[['protein_id', c_Gene]].drop_duplicates().rename(columns={c_Gene: 'gene_name'})
else:
    gene_map = None


# ---------------------------------------------------------
# STEP 9 — Combine labeled and rescued rows
# ---------------------------------------------------------
rows_all = pd.concat([labeled, rescued], ignore_index=True)
rows_all = rows_all[[
    'protein_id', 'protein_pos', 'seq_norm', c_PEP,
    'label', 'start_pos', 'kidx', 'pep_len', 'is_terminal', 'number_of_scans'
]].rename(columns={'seq_norm': 'sequence', c_PEP: 'pep'})

if gene_map is not None:
    rows_all = rows_all.merge(gene_map, on='protein_id', how='left')
else:
    rows_all['gene_name'] = ''

rows_all = rows_all.sort_values(['protein_id', 'protein_pos', 'label'])
print(f"\nBefore FASTA verification: {len(rows_all):,} rows")


# ---------------------------------------------------------
# STEP 10 — FASTA verification + tryptic check (parallel)
# ---------------------------------------------------------
rows_all['_canon_id'] = rows_all['protein_id'].apply(
    lambda x: str(x).split('-')[0] if pd.notna(x) else None
)

combo_cols = ['_canon_id', 'sequence', 'start_pos', 'kidx']
combos     = rows_all[combo_cols].drop_duplicates().reset_index(drop=True)
print(f"Unique combos to verify : {len(combos):,}")
print(f"Workers                 : {args.workers}")

tasks = [
    (row['_canon_id'],
     str(row['sequence']),
     int(row['start_pos']),
     int(row['kidx']))
    for _, row in combos.iterrows()
]

canon_seqs = {}
for acc, seq in fasta_seqs.items():
    canon_seqs[acc] = seq
    base = acc.split('-')[0]
    if base not in canon_seqs:
        canon_seqs[base] = seq

print("Running verification + tryptic check ...")
with Pool(
    processes=args.workers,
    initializer=_init_worker,
    initargs=(canon_seqs,)
) as pool:
    results = pool.map(verify_site, tasks)

combos['_verdict']     = [r[0] for r in results]
combos['_fasta_pos']   = [r[1] for r in results]
combos['_fasta_start'] = [r[2] for r in results]

rows_all = rows_all.merge(combos, on=combo_cols, how='left')

# Overwrite positions with FASTA-derived values for all kept verdicts
kept_verdicts = ['keep', 'remapped', 'rescued_tryptic']
keep_mask     = rows_all['_verdict'].isin(kept_verdicts)
rows_all.loc[keep_mask, 'protein_pos'] = rows_all.loc[keep_mask, '_fasta_pos'].astype(int)
rows_all.loc[keep_mask, 'start_pos']   = rows_all.loc[keep_mask, '_fasta_start'].astype(int)

n_keep            = int((rows_all['_verdict'] == 'keep').sum())
n_remapped        = int((rows_all['_verdict'] == 'remapped').sum())
n_rescued_tryptic = int((rows_all['_verdict'] == 'rescued_tryptic').sum())
n_not_tryptic     = int((rows_all['_verdict'] == 'not_tryptic').sum())
n_multi_tryptic   = int((rows_all['_verdict'] == 'multi_tryptic').sum())
n_absent          = int((rows_all['_verdict'] == 'absent').sum())
n_bad_kidx        = int((rows_all['_verdict'] == 'mismatch_bad_kidx').sum())
n_bad_k_prot      = int((rows_all['_verdict'] == 'mismatch_bad_k_prot').sum())

print("\n=== FASTA Verification + Tryptic Summary ===")
print(f"  keep             — confirmed, tryptic, position exact     : {n_keep:,}")
print(f"  remapped         — tryptic, start_pos corrected via FASTA : {n_remapped:,}")
print(f"  rescued_tryptic  — multi-locus, one tryptic locus kept    : {n_rescued_tryptic:,}")
print(f"  ---")
print(f"  not_tryptic      — no K/R before peptide start            : {n_not_tryptic:,}")
print(f"  multi_tryptic    — 2+ tryptic loci, ambiguous             : {n_multi_tryptic:,}")
print(f"  absent           — peptide not found in FASTA             : {n_absent:,}")
print(f"  bad_kidx         — kidx not K in peptide                  : {n_bad_kidx:,}")
print(f"  bad_k_prot       — kidx not K in protein                  : {n_bad_k_prot:,}")


# ---------------------------------------------------------
# STEP 11 — Build final output
# ---------------------------------------------------------
final = rows_all[keep_mask].copy()

final['position_source'] = final['_verdict'].map({
    'keep'            : 'fasta_verified',
    'remapped'        : 'fasta_remapped',
    'rescued_tryptic' : 'fasta_rescued_tryptic',
})

final = final.drop(columns=['_canon_id', '_verdict', '_fasta_pos', '_fasta_start'])
final = final.sort_values(['protein_id', 'protein_pos', 'label'])

pos_count = int((final['label'] == 'Positive').sum())
neg_count = int((final['label'] == 'Negative').sum())

print(f"\nFinal output rows : {len(final):,}")
print(f"  Positive        : {pos_count:,}")
print(f"  Negative        : {neg_count:,}")
print(f"  number_of_scans stored for downstream filtering")

out_file = out_dir / ("%s_acetylation_sites.csv" % prefix)
final.to_csv(out_file, index=False)
print(f"✅ Wrote {out_file}")


# ---------------------------------------------------------
# STEP 12 — Diagnostic CSVs
# ---------------------------------------------------------
drop_cols = ['_canon_id', '_verdict', '_fasta_pos', '_fasta_start']

def write_diag(df_in, verdict_val, suffix, label):
    sub = df_in[df_in['_verdict'] == verdict_val].copy()
    sub = sub.drop(columns=drop_cols)
    if len(sub) > 0:
        fpath = out_dir / ("%s_%s.csv" % (prefix, suffix))
        sub.to_csv(fpath, index=False)
        print(f"⚠️  Wrote {len(sub):,} {label} rows → {fpath}")

write_diag(rows_all, 'remapped',            'remapped_sites',       'remapped')
write_diag(rows_all, 'rescued_tryptic',     'rescued_tryptic',      'rescued_tryptic')
write_diag(rows_all, 'not_tryptic',         'dropped_not_tryptic',  'not_tryptic dropped')
write_diag(rows_all, 'multi_tryptic',       'dropped_multi_tryptic','multi_tryptic dropped')
write_diag(rows_all, 'absent',              'dropped_absent',       'absent dropped')
write_diag(rows_all, 'mismatch_bad_kidx',   'dropped_bad_kidx',     'bad_kidx dropped')
write_diag(rows_all, 'mismatch_bad_k_prot', 'dropped_bad_k_prot',   'bad_k_prot dropped')
