# An honest de novo design pipeline for non-immunosuppressive FKBP12 ligands as a pulmonary-arterial-hypertension strategy

*Working draft (G8). Every number traces to `docs/PROJECT_LOG.md` and the dated
research logs. Candidates herein are computational hypotheses to test, not drugs.*

## Abstract
Pulmonary arterial hypertension (PAH) is driven in part by loss of BMPR2 signaling.
Because ~68% of BMPR2 PAH variants are truncating (haploinsufficiency), the receptor
itself is a poor small-molecule target; and sotatercept (approved 2024-03) already
addresses the ligand-trap axis. We instead target **FKBP12** (UniProt P62942), whose
FK506-bound state activates BMP signaling — aiming for a **non-immunosuppressive**
(calcineurin-sparing) binder. We report a fully documented de novo pipeline and, as
much as the honest negatives, its **validated limits**: docking reproduces poses
(redock RMSD 0.81 Å) but cannot rank affinity here (enrichment ROC-AUC 0.455); a
ligand-based ML model ranks potency (scaffold-CV Spearman 0.75) but **cannot** predict
FKBP12 vs paralog selectivity (ρ ≈ 0); GNINA CNNaffinity *can* rank affinity (enrichment
ROC-AUC 0.799) and was adopted. Pocket-conditioned generation (DiffSBDD) plus
goal-directed GA optimization produced three developable, novel leads (GNINA
CNNaffinity 6.95–7.21) that **all remain bound in 1 ns MD** (ligand RMSD ≤ 2.9 Å).
Retrosynthesis (AiZynthFinder) shows their novel cores require custom synthesis, and a
structural analysis shows FKBP12 and FKBP12.6 differ at only **two conservative pocket
residues (Trp59, Ile90)** — which the leads nonetheless contact, making selectivity
structurally addressable but narrow. We present this as a reproducible methodology and a
concrete, testable hypothesis set rather than a drug claim.

## 1. Introduction
- PAH and the BMPR2 axis; why direct BMPR2 pocket-binding is misguided under
  haploinsufficiency; sotatercept as the already-approved axis modulator (Phase 0
  premise correction, ADR-001).
- FKBP12 rationale: FK506 activates BMP/ALK signaling via FKBP12 (JCI 2013); GPI-1046 is
  a non-immunosuppressive FKBP12 ligand precedent. Goal: pocket-binding that spares the
  calcineurin-recruiting effector surface.
- Contribution: not a drug, but (i) a reproducible de novo pipeline and (ii) an honest
  map of what each computational method can and cannot deliver on this target.

## 2. Methods
- **Structures/pocket.** RCSB 1FKF (FKBP12–FK506), 1B6C (ALK5 proxy, on-target), 1TCO
  (calcineurin, off-target), 2PPN (apo), 1F40 (GPI-1046). fpocket on 1FKF (FK506
  removed): druggability 0.517, volume 277 Å³, hydrophobicity 61.2, polar SASA ≈ 0.
  Interfaces (BioPython, 5 Å): on-target 22 residues overlapping the pocket; off-target
  adds the 40s-loop patch (residues 32–45) → design rule "fill pocket, avoid 40s loop."
- **Ligand data.** ChEMBL CHEMBL1902: 821 activities → 725 unique binders (pChEMBL
  3.1–11.0, median 6.71); 68 drug-like anchors (MW ≤ 500, QED ≥ 0.5, pChEMBL ≥ 6).
- **Scoring stack.** smina docking (pose control + enrichment); RandomForest on ECFP4
  (binder triage + pChEMBL potency + a selectivity probe); GNINA CNNaffinity rescoring.
- **Generation.** BRICS baseline; DiffSBDD pocket-conditioned generation on Colab GPU;
  goal-directed genetic-algorithm optimization under an applicability-domain box.
- **MD.** 1 ns explicit-solvent (Amber ff14SB + GAFF-2.11 + TIP3P, 0.15 M NaCl,
  ~20.5k atoms), OpenMM, ligand heavy-atom RMSD vs the docked start (protein-Cα aligned).
- **Retrosynthesis.** AiZynthFinder 4.4.1, public USPTO expansion/filter models + ZINC
  stock (17.4 M). **Selectivity.** UniProt P62942 vs P68106 alignment + pocket mapping +
  local re-docking contact analysis.

## 3. Results
### 3.1 What the scoring methods can and cannot do (the core methodology)
- **Docking = poses only.** FK506 redock RMSD 0.81 Å, but 20-active/40-decoy enrichment
  ROC-AUC **0.455** (below random). Retained for geometry, not ranking.
- **ML potency, not selectivity.** Binder triage scaffold-split ROC-AUC 0.993; pChEMBL
  regressor scaffold-split Spearman 0.748 (R² 0.56). But FKBP12−paralog selectivity
  prediction ρ ≈ **−0.01**: ligand-based selectivity is **not** learnable here. The
  classifier is an FKBP-*family* detector, not FKBP12-specific.
- **GNINA works.** 60-molecule enrichment ROC-AUC **0.799** (actives 6.34 ± 0.36 vs
  decoys 5.80 ± 0.53). Adopted as the affinity signal with a Tanimoto ≥ 0.30
  applicability-domain gate.

### 3.2 Generation and optimization
BRICS gave drug-like but low-diversity output (scaffold diversity 11%). DiffSBDD gave
diverse, novel, pocket-conditioned scaffolds (diversity 89%, novelty 0.18, PAINS 0%) but
best GNINA only 6.56, and scale-up plateaued (~6.5). Goal-directed GA broke the plateau
(GNINA 6.56 → 7.06 R1; high-quality exhaustiveness-8 rescoring 7.06 → 7.21), yielding
three developable, novel leads: **L1 7.21, L2 7.04, L3 6.95** (QED 0.62–0.72, SA
4.05–4.39, Tanimoto < 0.4 to known).

### 3.3 MD pose-stability (G3) — all three leads stable
1 ns MD ligand RMSD (mean/max Å): **L1 1.15/1.70, L2 1.70/2.91, L3 1.05/1.60** — all stay
in the pocket (≤ 2.9 Å max), separating real binding from lucky poses.

### 3.4 Synthesizability (G4) — an affinity↔synthesizability trade-off
AiZynthFinder solves no complete purchasable route for L1–L3: peripheral chemistry is
routine (4/5–5/6 precursors catalog) but the **novel core is the irreducible precursor**.
Screening all 37 GNINA-ranked candidates: **7/37 fully synthesizable**, but capped at
GNINA ~6.34; and the best synthesizable candidate (`O=C(O)CCCCN(CCO)c1cccc(F)c1F`, GNINA
6.34, SA 2.28, 1-step) **fails MD** (RMSD final 5.45 Å). No single candidate is
simultaneously synthesizable, MD-stable, and high-affinity.

### 3.5 Selectivity feasibility (G7 precursor) — narrow but real
FKBP12 vs FKBP12.6: 83.3% identical; the binding pocket differs at only **two
conservative residues, Trp59→Phe and Ile90→Val** (12.6 smaller/less aromatic at both) —
explaining the ML failure. Re-docking shows L1–L3 contact **Trp59** (all) and **Ile90**
(L2, L3), so selectivity is structurally addressable for the leads via those two handles.
Magnitude (ΔΔG) requires relative FEP (resource-gated) or experiment.

## 4. Discussion
The honest negatives are the contribution: on a shallow hydrophobic pocket with a
93%-identical-pocket paralog, docking-score ranking and ligand-based selectivity both
fail, while CNN rescoring + physics (MD) + structural analysis together still yield a
defensible, testable shortlist. Two lead classes emerge — (A) high-affinity, MD-stable,
novel-core leads needing custom synthesis; (B) readily-synthesizable but weaker/less
stable binders — framing a real medicinal-chemistry choice. Selectivity is a specific,
falsifiable design objective (exploit Trp59/Ile90), not an assumption.

**Limitations.** GNINA is a validated-but-imperfect proxy (ROC 0.80, not measured
binding); 1 ns single-replica MD is a stability screen, not free energy; the ALK
interface is a proxy (1B6C); FEP and wet-lab (binding + FKBP12.6/cardiac counter-screen)
are not yet done. All candidates remain hypotheses.

## 5. Conclusion
A reproducible pipeline and three MD-stable, novel, developable FKBP12 lead hypotheses
with a specific selectivity rationale — plus an honest account of where standard methods
fail on this target. The next decisive experiments are relative FEP for the
FKBP12/FKBP12.6 window and wet-lab validation of binding and selectivity.

## Data & code availability
All scripts (`src/`), gate logs (`docs/research-log/`), the master record
(`docs/PROJECT_LOG.md`), and result tables (`results/`) are in the project repository.
Large model/stock files are gitignored with provenance documented.
