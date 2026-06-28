# LUMEN

AI-driven de novo design of small molecules that restore the BMP signaling axis, for pulmonary arterial hypertension (PAH).

**Target:** FKBP12 (UniProt P62942) — see [`docs/decisions/ADR-001-target-redefinition.md`](docs/decisions/ADR-001-target-redefinition.md).
**Goal:** a non-immunosuppressive, oral small molecule that activates BMPR2-axis signaling (like FK506) while avoiding calcineurin inhibition.

> The original proposal (`proposal.html`) targeted mutant BMPR2 directly. Phase 0 validation found this premise scientifically flawed (BMPR2 haploinsufficiency; sotatercept already approved 2024) and the target was redefined. See [`docs/research-log/2026-06-28-premise-validation.md`](docs/research-log/2026-06-28-premise-validation.md).

## Layout
```
docs/decisions/      Architecture decision records (ADRs)
docs/research-log/   Dated research notes
docs/target-dossier-FKBP12.md   Single source of truth for the target
src/structure/       Structure fetch, pocket & interface analysis
src/generation/      (GPU) de novo generation — interfaces only, run on cluster
src/scoring/         ADMET / affinity / selectivity scoring
data/raw/            Downloaded structures (gitignored; provenance in manifests)
results/             Analysis outputs
```

## Environment
Conda (conda-forge + bioconda). The compchem stack needs conda, not pip wheels.
```bash
conda env create -f environment.yml        # main env: rdkit, biopython, fpocket
conda env create -f environment-dock.yml   # docking env: smina, openbabel (separate, avoids icu conflict)
```

## Reproduce Phase 1
```bash
python src/structure/fetch_structures.py     # download FKBP12 PDB structures
python src/structure/map_interfaces.py       # map ON/OFF-target interfaces
python src/data/collect_chembl_ligands.py    # collect FKBP12 ligands from ChEMBL
python src/scoring/dock_benchmark.py          # validate docking against known binders
```

## Compute note
This repo is developed on macOS without a GPU. AlphaFold3 inference, diffusion-based
generation (DiffSBDD), and RL optimization (REINVENT4) require a GPU cluster; those
steps are prepared here to the input/interface stage and offloaded.
