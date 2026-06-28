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

---

## UPDATE — Enrichment result (2026-06-28): docking score has NO discriminative power
Scripts: `src/scoring/build_enrichment_set.py`, `src/scoring/run_enrichment.py`
Set: 20 known drug-like actives (pChEMBL≥6) vs 39 property-matched random-ChEMBL decoys (assumed non-binders). Docked 59/60.

| Metric | Value | Read |
|---|---|---|
| ROC-AUC | **0.455** | below random (0.5) |
| EF@10% | 1.18 | negligible early enrichment |
| EF@20% | 1.07 | negligible |
| active mean score | −7.28 | — |
| decoy mean score | −7.33 | decoys score *as well or better* |

**Conclusion:** smina/Vina docking score does NOT separate FKBP12 binders from random drug-like molecules. The redock RMSD (0.81 Å) shows pose geometry is correct, but the *scoring function* is uninformative for affinity ranking on this shallow, near-apolar pocket (druggability 0.517, polar SASA ≈ 0) — Vina score tracks hydrophobic size, which most drug-like molecules share here.

**Honest caveats:** decoys are presumed (not confirmed) non-binders; actives are the small drug-like subset (the most potent FKBP12 binders are macrocycles); rigid receptor, single conformation, SMILES-derived 3D. Flexible-receptor docking or MM-GBSA/FEP rescoring might help marginally, but the strategic conclusion is robust.

**Consequence (Phase 2 strategy change):**
1. **Do NOT use docking score for ranking/selection.** Use docking only for pose generation / gross steric filtering.
2. **Affinity ranking must come from an ML model** trained on the 725 ChEMBL ligands (or MM-GBSA/FEP for top candidates).
3. Re-run this enrichment protocol as the acceptance test for any future scorer (require ROC-AUC ≥ 0.7 before trusting it).
