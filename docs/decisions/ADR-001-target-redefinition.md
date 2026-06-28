# ADR-001 — Target Redefinition: Mutant-BMPR2 Pocket Binder → FKBP12-Based BMP Signaling-Axis Restoration

Status: **Accepted** · 2026-06-28
Related: `docs/research-log/2026-06-28-premise-validation.md`

## Context
The original hypothesis in `proposal.html` was: "predict the 3D structure of mutant BMPR2 and de novo-design a small molecule for its binding pocket." Phase 0 validation found two problems:
1. ~68% of BMPR2 PAH variants are truncating → NMD → **haploinsufficiency (quantitative loss of protein)**. A pocket binder on the mutant protein cannot restore lost function; variants like W9X have no foldable structure at all.
2. The "zero remodeling-reversal drugs" premise is already broken by sotatercept (approved 2024).

## Decision
Redefine the target to **FKBP12 (UniProt P62942, gene FKBP1A)**.
- Design goal: a small molecule that disrupts the FKBP12–type I BMP receptor interaction (ALK1/ALK2/ALK3, GS domain) to activate downstream BMPR2 signaling (SMAD1/5/8, ID1).
- Key constraint (= the AI design differentiator): **avoid calcineurin inhibition** → eliminate FK506's immunosuppressive toxicity. Oral-compatible physicochemistry.
  - **Precision note (verification audit 2026-06-28):** this is an efficacy–toxicity *trade-off*, not a pure win. FK506 activates BMP signaling via BOTH calcineurin-dependent and calcineurin-independent pathways (JCI 2013). A calcineurin-sparing ligand captures the calcineurin-independent (FKBP12-displacement) arm and removes immunosuppression, but may deliver partial efficacy. Mitigation: engage all three ALK receptors (like FK506, unlike rapamycin) and avoid both calcineurin and mTOR liabilities. Efficacy sufficiency is a Phase 3 question.
- Differentiation: a different molecular node (FKBP12) + oral small molecule, vs sotatercept (ActRIIA ligand trap, injectable).

## Rationale
- FK506 (tacrolimus) activates BMP signaling at this same node → reverses rodent PAH, entered clinical trials (Spiekerkoetter, JCI 2013 / ERJ 2017). Mechanism and clinical precedent exist.
- FKBP12 has abundant experimental structures (many PDB entries; well-defined FK506/rapamycin pocket). Low AF3 dependence → Phase 1 can start immediately.
- Selectivity (engage FKBP12 to activate signaling, avoid calcineurin inhibition, avoid other FKBP family members) is a well-defined multi-objective problem → clear value for the de novo + RL pipeline.

## Alternatives considered but deferred
- **Wild-type BMPR2 stabilizer/activator**: target structure/pocket less well defined. Backup option.
- **Downstream nodes (e.g., SOX17)**: transcription factor, low small-molecule druggability.
- **Direct mutant-BMPR2 targeting (original)**: rejected for the reasons above.

## Consequences (what this decision changes)
- Phase 1: "BMPR2 variant structure atlas" → "**FKBP12–type I receptor interface & pocket definition + selectivity counter-target (calcineurin / FKBP family) panel**".
- Pipeline (AF3 → diffusion → RL) is preserved; only the input target is swapped.
- Positioning / paper framing: "world-first remodeling drug" → "non-immunosuppressive BMP-axis activator (complementary to sotatercept, oral)."

## Reversal
A backup target (wild-type BMPR2 stabilizer) is on record; if FKBP12 selectivity proves impossible, switch via ADR-002.
