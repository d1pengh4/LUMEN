# Target Dossier — FKBP12

Status: draft v0.1 · 2026-06-28 · Related: ADR-001, Phase 1 plan
Purpose: single source of truth for the target definition / pocket / selectivity criteria used as conditioning input for structure-based de novo molecule design.

## 1. Protein identity
| Field | Value |
|---|---|
| Name | Peptidyl-prolyl cis-trans isomerase FKBP1A (FKBP12) |
| UniProt | P62942 (human) |
| Gene | FKBP1A |
| Length | 108 aa (small, single-domain → structurally very well defined) |
| Function | PPIase; FK506/rapamycin binding protein; binds the GS domain of type I TGF-β/BMP receptors and suppresses signaling |

## 2. Therapeutic hypothesis (why this protein)
FKBP12 binds BMP/TGF-β type I receptors (ALK1/2/3, etc.) and **keeps them locked in an inactive state**. FK506 binds FKBP12 and releases this lock → activates SMAD1/5/8 signaling → reverses rodent PAH (JCI 2013). However, the **FKBP12–FK506 complex also inhibits calcineurin**, causing immunosuppressive toxicity.

→ **Design goal:** a non-immunosuppressive, oral small molecule that binds FKBP12 and displaces it from type I receptors (ON-target → turns BMP signaling on) **without creating the calcineurin-inhibitory surface (OFF-target avoidance)**.

## 3. Structural basis (PDB structures collected)
`src/structure/fetch_structures.py` → `data/raw/structures/`

| PDB | Role | Design meaning |
|---|---|---|
| 1FKF | drug pocket | FKBP12–FK506. Defines the **pocket the molecule binds** |
| 1B6C | ON-target interface | FKBP12–TGF-β type I receptor. The **signaling interface to disrupt** (homologous to BMP type I) |
| 1TCO | OFF-target interface | FKBP12–FK506–calcineurin. The **surface that must NOT be created** = selectivity criterion |
| 2PPN | apo | Ligand-free reference |

## 4. Binding pocket (fpocket quantification complete — 2026-06-28)
Ran fpocket on 1FKF (FK506 complex) after removing FK506 → **Pocket 1 = the FK506 binding site**.
`results/phase1/pocket_analysis/1FKF_out/`

| Metric | Value | Interpretation |
|---|---|---|
| Druggability score | **0.517** | Moderate. Shallow hydrophobic pocket → affinity will be non-trivial (Phase 2 risk) |
| Volume | 277 Å³ | Fits a single small molecule |
| Hydrophobicity score | 61.2 | Strongly hydrophobic |
| Polar SASA | 0.000 | **Almost entirely apolar** → strategy: hydrophobic scaffold + a few H-bond anchors |
| Charge score | -2 | Weakly negative |

Literature key residues (human FKBP12 numbering): Tyr26, Phe36, Asp37, Arg42, Phe46, Val55, Ile56, Trp59, Tyr82, His87, Ile90, Ile91, Phe99.
- **Design implication:** Polar SASA ≈ 0 → affinity depends on hydrophobic contacts. Using Asp37/Arg42 as H-bond/salt-bridge anchors is the key to achieving affinity and selectivity together.
- apo (2PPN) has druggability 0.007 (open pocket) → use the **holo (1FKF) pocket geometry** as the generative-model conditioning.

## 5. Selectivity counter-target panel (things to avoid)
| Counter-target | UniProt | Reason |
|---|---|---|
| Calcineurin (PPP3CA / PPP3R1) | Q08209 / P63098 | Source of FK506's immunosuppressive toxicity. The complex must not create this surface |
| FKBP1B (FKBP12.6) | P68106 | Regulates cardiac RyR2. Cross-binding risks cardiac side effects |
| Other FKBP family (FKBP51/52, etc.) | - | Selectivity profiling panel |

Selectivity metrics enter the Phase 2 Pareto objective (ADR-001).

### 5.1 Interface mapping results (2026-06-28) — structural selectivity rule
`src/structure/map_interfaces.py` → `results/phase1/interfaces/interface_residues.json`

- **ON-target interface (1B6C, FKBP12↔type I receptor):** 22 residues. Largely overlaps the FK506-pocket residues (Tyr26, Phe36, Asp37, Arg42, Phe46, Val55, Ile56, Trp59, Tyr82, His87, Ile90/91, Phe99) → **occupying the pocket disrupts receptor binding** (mechanism confirmed).
- **OFF-target interface (1TCO, FKBP12↔calcineurin):** 23 residues. Beyond the pocket, it additionally uses the **"40s loop" patch (Thr27, Asp32, Lys34/35, Ser38, Arg40, Asn43, Lys44, Pro45)**. This patch is **absent** from the receptor interface.
- **➡ Design rule (non-immunosuppressive):** fill the hydrophobic pocket (Pocket 1) but use a **small molecule that does not touch or reorganize the calcineurin-recruiting 40s loop (residues 32–45)**. FK506 is immunosuppressive because its out-of-pocket effector portion recruits calcineurin → removing that portion is the key to selectivity. Impose "confined to pocket interior + avoid 40s loop" as a geometric constraint on the generative model.

## 6. Differentiation (vs competing drugs)
- **Sotatercept (approved 2024):** ActRIIA ligand trap, injectable. → This project: **different node (FKBP12) + oral small molecule**.
- **FK506 (clinical):** non-selective, immunosuppressive. → This project: **calcineurin-avoidance selectivity** is the key advance.

## Open / next
- [x] Quantify pocket with fpocket (1FKF)
- [x] Map ON-target interface (1B6C) and OFF-target surface (1TCO)
- [ ] Collect FKBP12-binding ligands from ChEMBL (training set for the affinity model)
- [ ] Convert pocket to generative-model input format (DiffSBDD conditioning)
