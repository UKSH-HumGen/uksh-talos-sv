# uksh-talos-sv

**Talos Structural Variant Extension — UKSH Fork**

> **Status:** This repository is a fork of
> [populationgenomics/talos](https://github.com/populationgenomics/talos)
> at upstream version **v8.2.0**, extended with structural-variant (SV)
> support. It is **not** rebased onto newer upstream releases. The
> changes here were developed for, and used to produce the results in,
> the Kaschta *et al.* preprint on automated rare disease reanalysis at
> UKSH (see [Citing this fork](#citing-this-fork) below).

[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

---

## Overview

[Talos](https://github.com/populationgenomics/talos) is an open-source
variant prioritisation framework built by the Centre for Population
Genomics (CPG) for automated reanalysis of rare-disease cohorts. It
combines static annotations (gnomAD, predicted consequence, in-silico
scores) with dynamic evidence (ClinVar, PanelApp) and applies
ACMG/AMP-aligned, MOI-aware rules to flag variants in known disease
genes.

`uksh-talos-sv` extends Talos with **structural-variant support**, so
the same reanalysis logic can run end-to-end on a joint VCF that
contains both small variants (SNVs/indels) and SVs (CNVs, deletions,
duplications, inversions, complex events). The fork was developed at
the Institute of Human Genetics, Universitaetsklinikum
Schleswig-Holstein (UKSH).

For the unchanged upstream documentation of the small-variant pipeline,
the configuration, and the philosophy of the tool, please refer to the
[upstream README](https://github.com/populationgenomics/talos).

## Key features added by this fork

- **`CategorizeSVs.py`** — GFF3-based gene/exon annotation for SVs.
  Each SV is intersected with Ensembl gene/transcript/exon intervals to
  attach gene IDs, affected exon counts, and predicted consequence
  classes used downstream by `ValidateMOI`.
- **Genotype-aware multi-sample SV merging.** SVs in joint VCFs are
  assigned only to samples whose genotype confirms carriage (`0/1`,
  `1/1`), instead of being broadcast to every sample in the joint
  matrix. This fixes false-positive assignments seen with DRAGEN
  joint-genotyped CNVs.
- **MOI validation for SVs.** `utils.py` and `ValidateMOI.py` recognise
  structural variants alongside small variants, so AR/AD/XL inheritance
  reasoning works for SV/SNV+SV combinations.
- **SV coordinate preservation** during VCF merging. The
  `MergeVcfsWithBcftools` step is modified so that SV ALT alleles and
  END/SVLEN/SVTYPE fields are not flattened.
- **DRAGEN CNV display fix** in the HTML report.

A more detailed walk-through of each modification (with code anchors)
is in
[`docs/TALOS_SV_MODIFICATIONS.md`](docs/TALOS_SV_MODIFICATIONS.md).

## Differences from upstream

| Aspect | Upstream Talos v8.2.0 | uksh-talos-sv 1.0.0 |
|---|---|---|
| Small variants | Supported | Supported |
| Structural variants | Limited (SV labels propagated, no MOI) | Fully integrated through MOI validation |
| Multi-sample SV assignment | Broadcast to all samples | Genotype-restricted per sample |
| `cpg-utils` requirement | Hard dependency | Optional extra (`pip install ".[cpg]"`) |
| CPG-internal scripts | `src/talos/cpg_internal_scripts/` | Removed (metamist/seqr, cpg-flow stages) |
| Docker image | `talos:8.2.0`, GCR-pushed | `uksh-talos-sv:1.0.0`, local build only |
| Experimental modules | Mixed into main workflow | Quarantined under `experimental/` |

The fork stays on the v8.2.0 base. Upstream features added in v8.3+ are
**not** ported.

## Installation

The package is published as `uksh-talos-sv` but installs as the Python
module `talos` (drop-in name-compatible with upstream).

```bash
git clone https://github.com/UKSH-HumGen/uksh-talos-sv.git
cd uksh-talos-sv

# Standalone install (no Google Cloud helpers)
pip install .

# Or, with cpg-utils for Hail/Batch initialisation against GCP:
pip install ".[cpg]"

# Development install with the test suite:
pip install -e ".[test]"
```

A Docker image can be built locally with:

```bash
docker build -t uksh-talos-sv:1.0.0 .
```

## Quick start

The Nextflow pipeline under `nextflow/` orchestrates the full reanalysis
workflow. A synthetic example pedigree is shipped so you can test the
plumbing without supplying real patient data.

```bash
# 1. Either use the bundled synthetic pedigree directly:
PED=nextflow/inputs/pedigree_synthetic_example.ped

# 2. ...or generate one from a Cases.txt sheet:
python nextflow/inputs/generate_pedigree.py \
    --input nextflow/inputs/Cases_example.txt \
    --output nextflow/inputs/pedigree_synthetic_example.ped

# 3. Validate it:
python nextflow/inputs/validate_pedigree_simple.py \
    nextflow/inputs/pedigree_synthetic_example.ped

# 4. Run the workflow (configure paths in nextflow/talos.config first):
nextflow run nextflow/talos.nf -c nextflow/talos.config --pedigree $PED
```

Real patient PEDs and VCFs are **deliberately excluded from this
repository** (and from `.gitignore`'s whitelist) — keep them outside
the working tree.

## Repository layout

```
uksh-talos-sv/
├── src/talos/                — Python package (same import name as upstream)
├── nextflow/                 — active Nextflow workflow + modules
│   └── inputs/               — config and the synthetic example pedigree
├── experimental/             — not integrated; see experimental/README.md
├── docs/                     — SV modifications and known issues
├── test/                     — pytest test suite
└── pyproject.toml, Dockerfile, CITATION.cff
```

## Citing this fork

Until the preprint is published, please cite both the upstream tool and
this fork. A machine-readable `CITATION.cff` is in the repository root.

```bibtex
@software{uksh_talos_sv_2026,
  author       = {Kaschta, Daniel},
  title        = {uksh-talos-sv: Talos with structural variant support},
  version      = {1.0.0-uksh},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.XXXX/PLACEHOLDER-PREPRINT-DOI},
  url          = {https://github.com/UKSH-HumGen/uksh-talos-sv}
}
```

> **TODO Daniel:** replace `10.XXXX/PLACEHOLDER-PREPRINT-DOI` with the
> real DOI once the Zenodo release and the medRxiv preprint are live.
> The same placeholder appears in `CITATION.cff` — update both.

Please also cite the upstream tool:

> Centre for Population Genomics. **Talos**.
> https://github.com/populationgenomics/talos

### Generating a Zenodo DOI

Zenodo issues a citable DOI for each tagged GitHub release. Setup is a
one-time, ~5-minute job; afterwards every release auto-archives.

1. Sign in to Zenodo with the same GitHub account that owns the
   repository, then visit your Zenodo **GitHub** settings page and
   flip the toggle for `UKSH-HumGen/uksh-talos-sv` to **ON**.
2. In GitHub, create a release tagged `v1.0.0-uksh` (release title can
   match). Zenodo picks up the webhook automatically.
3. Wait a minute, then open the Zenodo deposit Zenodo created. Note the
   DOI it minted (looks like `10.5281/zenodo.XXXXXXX`).
4. Replace `10.XXXX/PLACEHOLDER-PREPRINT-DOI` in `CITATION.cff` and in
   the BibTeX block above with that real DOI. Commit and push.
5. Optional: add the Zenodo "DOI" badge to the top of this README.

## Acknowledgements

Upstream Talos was conceived and is maintained by the
[Centre for Population Genomics](https://github.com/populationgenomics)
team. The structural-variant extensions in this fork build directly on
their work and would not exist without it.

## License

MIT, with dual copyright. See [`LICENSE`](LICENSE).

```
Copyright (c) 2022 Centre for Population Genomics
Copyright (c) 2025-2026 Daniel Kaschta, UKSH (modifications and additions)
```
