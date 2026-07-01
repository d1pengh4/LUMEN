# Developability panel + continuous verification (Track A, step 1)

Date: 2026-07-01
Scripts: `src/scoring/developability_panel.py`
Outputs: `results/developability/*.csv`, `results/phase2b/B_integrated_candidates.csv`

## Developability (local, rule-based proxies — honest)
RDKit-only ADMET-proxy + SA score + drug-likeness (NOT a trained ADMET model).

| Set | n | developable* | SA median | QED median | PAINS |
|---|---|---|---|---|---|
| B promising (dock-selected) | 8 | 4 | 4.88 | 0.72 | 0 |
| B all docked | 49 | 16 | 4.75 | 0.57 | 0 |
| A top-50 (BRICS) | 50 | 3 | 4.64 | 0.60 | 0 |
*developable = SA≤4.5 & QED≥0.5 & Lipinski & Veber & no PAINS

Honest finding: DiffSBDD (B) molecules have somewhat high SA scores (median 4.75–4.88);
synthetic accessibility is a real weak point of de novo molecules — only half of the
dock-promising B set is also easily synthesizable. Documented as a paper limitation.

## Integrated B candidates (TIER-1)
Combine independent filters: dock ≤ −7 AND developable AND novel (Tanimoto < 0.3).
→ 4 TIER-1 candidates (`B_integrated_candidates.csv`).

## Continuous verification — docking stability across seeds
Because docking SCORE ranks weakly (enrichment ROC-AUC 0.455), a −7 threshold is
only meaningful if the number is reproducible. Re-docked the 4 TIER-1 molecules
with seeds 1/7/99 (orig seed 42):

| orig | seed1 | seed7 | seed99 | std |
|---|---|---|---|---|
| −8.54 | −7.87 | −7.35 | −7.37 | 0.48 |
| −7.58 | −7.44 | −7.56 | −7.53 | 0.05 |
| −7.57 | −7.29 | −7.48 | −7.27 | 0.13 |
| −7.34 | −7.39 | −7.35 | −7.35 | 0.02 |

Findings:
- All 4 stably dock ≤ −7.3 across seeds → the filter is self-consistent, not noise. ✅
- **Correction:** candidate #1's −8.54 (seed 42) was seed-optimistic; true ≈ −7.5.
  Report the seed-averaged score, not the best single seed.
- **Caveat (unchanged):** stability ≠ discrimination. Random decoys also dock ~−7.3,
  so TIER-1 = "passes multiple independent weak filters (novel + developable +
  stably in binder-range)", a prioritization, NOT validated affinity. Real ranking
  awaits MM-GBSA (next) and ultimately experiment.
