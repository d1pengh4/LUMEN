# G2: Goal-directed optimization — attempt log (detailed)

Date: 2026-07-02 · Roadmap G2
Goal: push candidates beyond the brute-force plateau (CNNaffinity ~6.5) toward the
STOP target (≥7.0) via directed optimization, on a frugal (CPU, no GPU) budget.

## Design decision (why GA, not REINVENT4/RL)
- RL (REINVENT4) with a GNINA reward needs thousands of dockings → infeasible on
  free GPU (violates frugal budget).
- Chosen: **BRICS fragment-recombination genetic algorithm** optimizing the FAST
  validated ML scorer (potency regressor `pred_pchembl`), with hard constraints
  (in-domain ad≥0.30, QED≥0.55, MW 250–500). Optimized output then GNINA-validated.
- Script: `src/generation/ga_optimize.py` (local CPU).

## Attempt 1 — FLAWED (seed dominance). Recorded for the record.
Seeds: generated developable candidates + drug-like ChEMBL anchors (as fragment donors).
Result: best `pred_pchembl` = 8.39, **flat across all 8 generations** (no improvement).
Top "optimized" molecules had **ad_tanimoto = 1.0** — i.e. they ARE known ChEMBL
binders (identical to training molecules), e.g.
`COC[C@H]1CN(Cc2ccccn2)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1cc(Cl)cc(Cl)c1`.

**Why it failed (honest):**
1. Seeding the *population* with real known binders (anchors) put ad=1.0 molecules
   at the top. The regressor scores them highest (it was TRAINED on them — memorized).
2. The GA maximizing `pred_pchembl` just keeps those known binders on top; offspring
   rarely beat a memorized training molecule's score.
3. So "8.39" = "a known binder's regressor score," NOT a novel optimized candidate.
   This is seed-dominance / reward-hacking, not real optimization.

**Evidence it's an artifact:** ad_tanimoto = 1.0 (identical to training set); best
score never moved (8.39 → 8.39) → the GA never actually improved anything novel.

## Fix (attempt 2)
- **Do NOT seed the population with anchors.** Seed only with NOVEL generated
  developable candidates. Use anchors ONLY as fragment donors during breeding.
- **Add a novelty ceiling:** require 0.30 ≤ ad_tanimoto ≤ 0.65 (in-domain but NOT a
  near-duplicate of a known molecule). This forbids the GA from "winning" by
  regurgitating training molecules.
- Objective unchanged (maximize pred_pchembl) within the novelty+developability box.
- Re-run, then GNINA-validate the top novel optimized molecules vs the 6.56 baseline.

## Attempt 2 — novelty-constrained (fixed)
Changes: population seeded ONLY with novel generated developable candidates;
anchors used as fragment DONORS only; hard novelty box 0.30 ≤ ad_tanimoto ≤ 0.65.

Result (`results/phase2b/ga_optimized.csv`):
- best `pred_pchembl` trajectory: 5.16 → 7.34 → (flat to gen 8).
- Top molecules are NOVEL (ad_tanimoto 0.44–0.52, i.e. in-domain but not duplicates).
  Example: `C=C[C@H]1CN(CC(=O)O)C(=O)[C@@H]2CCC[C@H]1N2Cc1cccc(CC(N)=O)c1`
  (bicyclic diketopiperazine core of the FKBP12 chemotype with a novel
  benzyl-acetamide arm replacing the dichlorophenyl-sulfonyl).
- pred_pchembl 7.34 > developable-candidate levels and above the brute-force plateau.

**Honest caveats:**
1. Converged fast (best flat after gen 2) — limited exploration; the GA found one
   good region and stopped improving. Not a dramatic optimization curve.
2. **Model disagreement (yellow flag):** top molecules have pred_pchembl 7.34 but
   binder_prob only ~0.37. The potency REGRESSOR likes them; the binder CLASSIFIER
   is lukewarm. Likely the regressor extrapolating optimistically on a novel
   chemotype. So "7.34" must NOT be trusted at face value.
3. `pred_pchembl` is the FAST proxy, not the validated ranker. **Decisive test =
   GNINA** on these optimized molecules vs the 6.56 baseline / 7.0 STOP target.

## Next (decisive)
GNINA-score the top ~15 GA-optimized novel molecules. Only if GNINA confirms
improvement over 6.56 (ideally toward ≥7.0) is G2 a real success; otherwise the
ML-space gain is an artifact and we document the ceiling honestly.
