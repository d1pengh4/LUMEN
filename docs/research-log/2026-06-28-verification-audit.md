# Verification Audit — Phase 0–1

Date: 2026-06-28
Purpose: independent re-verification of all claims and computations to date (medical-grade rigor). Each item checked against primary evidence, not assumption.

## A. Technical / computational checks — ALL PASS

| Check | Method | Result |
|---|---|---|
| Structure identity (1FKF/1B6C/1TCO/2PPN) | PDB TITLE/COMPND headers | ✓ exactly as labeled |
| Chain assignment 1B6C | COMPND | ✓ A,C,E,G = FKBP12; B,D,F,H = TGF-β type I receptor |
| Chain assignment 1TCO | COMPND | ✓ A = calcineurin A (catalytic), B = calcineurin B, C = FKBP12 |
| fpocket Pocket 1 = FK506 site | centroid distance pocket1 ↔ native FK506 | ✓ 4.35 Å (same site) |
| ChEMBL completeness | API total_count vs fetched | ✓ 821 = 821 (no pagination loss); 725 unique after cleaning |
| Docking setup validity | redock FK506, RMSD vs crystal pose (obrms) | ✓ **0.81 Å** (success criterion <2.0 Å) |

Conclusion: the structural pipeline, data collection, and docking are geometrically and numerically sound.

## B. Scientific claim verification

### B1. BMPR2 haploinsufficiency (original premise rejection) — CONFIRMED
~68% of BMPR2 PAH variants are nonsense/frameshift → NMD → haploinsufficiency. Two independent sources. The rejection of the "mutant-BMPR2 pocket binder" premise stands.

### B2. Sotatercept already approved — CONFIRMED
WINREVAIR (sotatercept), FDA 2024-03, first-in-class disease-modifying PAH drug (ActRIIA ligand trap). The "zero remodeling drugs" premise is correctly rejected.

### B3. Core target mechanism — CONFIRMED with an important correction
Verified (Spiekerkoetter et al., JCI 2013, PMC3726153):
- FK506 releases FKBP12 from ALK1/ALK2/ALK3 → activates SMAD1/5, MAPK, ID1. FKBP12 is a *repressor* of BMP signaling. ✓ Our mechanism is correct.
- Receptor coverage: BMP4 signals mainly via ALK3 (some ALK1); **FK506 engages all three ALKs**, whereas **rapamycin only engages FKBP12 bound to ALK2** → FK506 superior.

**⚠ CORRECTION (error found in our own docs):** FK506 activates BMP signaling via **BOTH calcineurin-dependent AND calcineurin-independent pathways** (dual mechanism). Our ADR-001 / dossier framed "avoid calcineurin inhibition" as a pure benefit (toxicity removal only). That is **inaccurate** — calcineurin inhibition contributes to part of the efficacy. The honest framing is an **efficacy–toxicity trade-off**:
- A calcineurin-sparing FKBP12 ligand captures the **calcineurin-independent (FKBP12-displacement) arm** of the effect and removes immunosuppression.
- Risk: it may deliver **partial efficacy** vs FK506, losing the calcineurin-dependent component.
- Mitigation built into hypothesis: maximize the displacement arm by engaging **all three ALK receptors** (like FK506, unlike rapamycin) and avoid BOTH calcineurin and mTOR liabilities (rapamycin's). Precedent (rapamycin) shows pocket-binding-only ligands DO activate BMP signaling, supporting feasibility.

### B4. Structural proxy nuance (precision note)
1B6C is FKBP12–**TGF-β type I receptor (ALK5)**, not a BMP type I receptor (ALK1/2/3). We use it as a **homologous proxy** for the FKBP12–type I receptor interface, justified because FKBP12 binds the conserved GS domain shared across type I receptors. No experimental FKBP12–ALK1/2/3 complex was used. Docs updated to state this explicitly.

## C. Net assessment
- No computational errors found; all geometry/data checks pass.
- One scientific oversimplification found and corrected (B3): the calcineurin-avoidance benefit is a trade-off, not free. Target hypothesis remains viable and is now stated precisely.
- Open scientific risk (carried to risk register): whether calcineurin-sparing efficacy is sufficient — answerable only by Phase 3 cellular assays (PASMC SMAD1/5 readout).

## Sources
- JCI 2013 FK506/BMPR2: jci.org/articles/view/65592, PMC3726153, PubMed 23867624
- BMPR2 mutation spectrum: PMC2762975; AHA CIRCRESAHA.122.321554
- Sotatercept: Merck/FDA 2024-03; PubMed 39466552
