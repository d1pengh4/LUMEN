# Docking Benchmark — Scoring Pipeline Validation

Date: 2026-06-28 · Related: target dossier §4–5.2
Script: `src/scoring/dock_benchmark.py` · Output: `results/phase1/docking/docking_scores.csv`
Engine: smina (AutoDock Vina 1.1.2 scoring), env `lumen-dock`

## Setup
- Receptor: 1FKF with FK506 and waters removed → PDBQT (rigid, pH 7.4).
- Box: centered on native FK506 centroid (22.3, 30.2, 23.7), 24 Å cube.
- Ligands: native FK506 (positive control, redock) + top 6 drug-like ChEMBL anchors.

## Results
| Ligand | pKi/pChEMBL | Dock score (kcal/mol) |
|---|---|---|
| FK506 (native, redock) | high | **−10.3** |
| CHEMBL4090599 | 9.35 | −7.52 |
| CHEMBL4072643 | 8.15 | −7.57 |
| CHEMBL4075704 | 8.15 | −7.55 |
| CHEMBL4101268 | 8.15 | −6.93 |
| CHEMBL33094 | 9.00 | −6.85 |
| CHEMBL4102121 | 8.15 | −6.73 |

## Interpretation
- **Pipeline validated:** all known binders score strongly negative; the native high-affinity macrocycle FK506 scores best (−10.3). Receptor prep + box + docking are wired correctly.
- **Caveat (honest):** within the anchor set, dock score does not track pKi (e.g., the most potent, pKi 9.35, is not the best-scoring). Docking scores are weak affinity predictors and biased toward molecular size/contacts.
- **Consequence for Phase 2:** docking is fit for coarse pose generation / gross filtering only. Affinity ranking must use an ML model (GNN) and/or MM-GBSA rescoring, and scores should be size-normalized (ligand efficiency). This matches the multi-property scoring plan in ADR-001.

## Next
- Decoy enrichment test (property-matched decoys) to quantify discrimination.
- Build/port the GNN affinity model (Phase 2; GPU offload).
