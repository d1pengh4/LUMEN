# Selectivity Panel — and why ligand-based ML cannot solve selectivity here

Date: 2026-06-28 · Related: ml-scorer log, dossier §5
Script: `src/scoring/selectivity_panel.py` · Data: `data/processed/fkbp_offtargets.csv`

## Question
Selectivity (FKBP12 over calcineurin, FKBP12.6, other FKBPs) is the *core differentiator* of this project (non-immunosuppressive, FKBP12-specific). Can we score it computationally to drive the generative reward?

## Data availability (ChEMBL, pChEMBL)
| Counter-target | ligands | matched pairs w/ FKBP12 | usable for ML? |
|---|---|---|---|
| Calcineurin (PPP3CA) | 38 | — | No (sparse + direct-inhibitor mechanism, not the immunophilin ternary complex) |
| FKBP12.6 (FKBP1B) | 44 | 40 | No (83% identical to FKBP12) |
| FKBP4 | 101 | 62 | Marginal |
| FKBP5 | 246 | 142 | Yes (data-wise) |

FKBP12 vs FKBP12.6 sequence identity = **83%** (108 aa, BLOSUM62 global) — computed, confirms intrinsic hardness.

## Leakage-free selectivity test (matched pairs held out of training)
Predicted ΔpChEMBL (reg_FKBP12 − reg_offtarget) vs measured ΔpChEMBL:

| Off-target | measured mean Δ | FKBP12-selective (Δ>0.5) | predicted-vs-measured Spearman |
|---|---|---|---|
| FKBP5 | +0.33 | 60/142 | **−0.01 (no power)** |
| FKBP4 | +1.53 | 58/62 | 0.31 (weak) |
| FKBP12.6 | −0.14 | **1/40** | n/a (no model) |

## Honest conclusions
1. **Ligand-based ML cannot predict FKBP-family selectivity** from current data. Selectivity is a small difference between two noisy affinity predictions; errors compound to ≈0 signal (FKBP5 ρ=−0.01).
2. **FKBP12.6 is the critical, hardest case:** 83% identity, and known ligands do not discriminate it (only 1/40 are FKBP12-selective). Since FKBP12.6–RyR2 is the cardiac-safety liability (verification audit D3), an **experimental FKBP12.6 counter-screen is mandatory and non-substitutable by computation**.
3. **FKBP4 is empirically not a major problem** (ligands are naturally FKBP12-selective, mean Δ +1.53), though we can't predict it well.
4. **Calcineurin selectivity stays structural** (40s-loop avoidance from 1TCO), not ML.

## Impact on strategy (important)
- Do NOT build the generative reward around an ML selectivity score — it would optimize noise.
- Selectivity enforcement = (a) structural geometric constraints in generation (calcineurin 40s loop), (b) potency + drug-likeness + applicability-domain from ML, (c) **experimental counter-screens (esp. FKBP12.6) as the real gate** in Phase 3.
- The potency regressor (ρ 0.76) remains valid and useful; the per-off-target affinity regressors are saved but selectivity *differencing* is unreliable — use only as coarse "binds off-target at all?" flags, not for ranking.
