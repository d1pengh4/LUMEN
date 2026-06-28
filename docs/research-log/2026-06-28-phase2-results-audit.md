# Phase 2 Results Audit — skeptical analysis of the candidate set

Date: 2026-06-28 · Related: phase2-poc log
Scope: full re-verification of repo state + adversarial analysis of the 400 generated candidates.

## Repo integrity — PASS
11 commits, clean tree, pushed. All 10 key artifacts exist on disk (structures, ligand sets, enrichment, off-targets, interfaces, candidates, models).

## Candidate analysis

### 1. Chemical sanity — mostly OK, one real flaw
- 400/400 parsable; PAINS alerts ~0% (1 molecule). Good.
- **Lipophilicity skew (flaw):** logP median 5.3; 51/400 have logP>6 or >12 rotatable bonds. The shallow apolar FKBP12 pocket pulls generation toward greasy molecules → poor oral-drug properties. **Action: tighten filters (logP ≤ 5, rotatable ≤ 10) in future runs.**

### 2. Negative control — what binder_prob really means
BRICS molecules built from RANDOM ChEMBL fragments vs from FKBP12-binder fragments:
| Source | mean binder_prob | in-domain |
|---|---|---|
| random-ChEMBL fragments | **0.027** | 0% |
| FKBP12-binder fragments (our candidates) | 0.591 | 100% |

Reading: binder_prob is NOT "any BRICS molecule scores high" (control is 0.027) — so the signal is real in the triage sense. BUT our candidates score high largely because they descend from binder fragments. **binder_prob = "lives in FKBP12-binder chemical space," a valid triage signal, NOT a proof of affinity.** Consistent with the documented classifier limitation (FKBP-family triage, not FKBP12-specific affinity).

### 3. Applicability-domain gate — works
Glucose, a triglyceride, polyamine, benzene, metformin → all flagged in_domain=False (ad_tanimoto 0.08–0.29). Note glucose fools the regressor (pred 6.87) but is correctly caught by the AD gate + binder_prob 0.02 → multilayer defense holds.

### 4. Potency prediction — conservative, no over-extrapolation
Candidate pred_pchembl (median 6.29, max 6.93) sits within the known-ligand range (median 6.71, max 11.0) and never exceeds it. Honest uncertainty: regressor CV MAE ~0.66 log → ~5× in Ki. So "predicted ~0.3 µM" means "plausibly sub-to-low-µM," not a precise number.

## Verdict
The Phase-2 PoC output is **internally valid and honestly characterized**: real chemistry, working AD gate, conservative potency, controls behave correctly. Its genuine weaknesses are (a) low scaffold diversity (11%, from phase2-poc log), (b) lipophilicity skew, (c) binder_prob is triage not affinity, (d) no selectivity scoring. None of these are hidden; all motivate Phase 2-B (a real pocket-conditioned generative model) + tighter physchem filters.

**Bottom line (unchanged, reaffirmed):** these are testable hypotheses produced by a validated, non-self-deceiving pipeline — not drug candidates.
