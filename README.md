# PTM Extraction Pipelines

Re-extraction of the phosphorylation and acetylation slices of the [DeepMVP PTMAtlas](https://github.com/DeepMVP) dataset with stricter quality controls.

## Structure

```
ptm-extraction-pipelines/
├── phospho/
│   ├── PXD012174_human_phosphoproteome/   # PRIDE — Human Phosphoproteome Map
│   └── MSV000078509_serial_enrichment/    # MassIVE — Serial Enrichment
└── acetyl/
    ├── MSV000078509_kac/                  # MassIVE — Kac from High_Kac + pSTY
    ├── PXD005793_SILAC_acetylome/         # PRIDE — SILAC Acetylome
    └── MSV000082644_medulloblastoma/      # MassIVE — Medulloblastoma
```

See the `README.md` inside each folder for the step-by-step workflow.

## PTM types

- **[Phosphorylation →](phospho/)** — 2 datasets, labels on S/T/Y residues, split by amino acid.
- **[Acetylation →](acetyl/)** — 3 datasets, labels on K residues, terminal-K bias corrected via reduceK.

## Dependencies

```
pandas
numpy
biopython       # for FASTA parsing
pyteomics       # optional, for tryptic rules
ahocorasick     # MSV000082644 pipeline only
```

## Reference

Wen et al. *DeepMVP: Multi-View Pooling Recovers Missing Annotations from Protein Post-Translational Modification Databases.* Nature Methods, 2025.
