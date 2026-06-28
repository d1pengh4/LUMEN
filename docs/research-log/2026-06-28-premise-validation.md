# Phase 0 — Premise Validation

Date: 2026-06-28
Author: Pre-research literature check of the core claims in `proposal.html`
Conclusion: **Major scientific correction required before adopting the proposal.** Two headline premises are factually wrong. However, the project is salvageable with a redirected target.

---

## Check 1 — "Zero remodeling-reversal drugs exist" → **FALSE**

Proposal claim (Hero, problem-card 01, impact "0→1"):
> "There is not a single drug in the world that reverses the remodeling itself."

Literature:
- **Sotatercept (WINREVAIR, Merck)** — FDA approved **March 2024** (based on STELLAR Phase 3). The first **disease-modifying / first-in-class** drug for PAH.
- Mechanism: activin signaling inhibitor (ActRIIA ligand trap) — restores BMP/activin signaling balance to regulate vascular cell proliferation.
- Preclinical data show **reduced vessel-wall thickness, partial reversal of RV remodeling, improved hemodynamics** — i.e., exactly the effect the proposal claims to deliver "for the first time in the world."

Impact: The project's core differentiation claim ("world-first remodeling inhibitor") is **already preempted by an approved drug**. Positioning must shift from "0→1" to "1→N (better/oral/targeted)."

## Check 2 — "De novo small molecule binding the BMPR2 pocket" → **Flawed target strategy**

Proposal claim (Hypothesis, Phase 1·2):
> "Predict the 3D structure of BMPR2 variants and de novo-generate small molecules optimized for the binding pocket."
> Example variants: R899X, W9X, 1750delA

Literature:
- **~68% of BMPR2 PAH variants are nonsense/frameshift (truncating)** → **nonsense-mediated decay (NMD)** → **quantitative loss of protein = haploinsufficiency (loss of function)**.
- W9X = stop at residue 9 → essentially no protein. R899X = loss of C-terminus → NMD target. **A "3D structure atlas" or "binding pocket" of such variants is meaningless** (the protein is absent).
- Core contradiction: the problem is "too little BMPR2," yet designing a small molecule that binds the BMPR2 pocket (inhibitor-like) **does not restore function**. A lost protein cannot be rescued by a pocket binder.

Impact: "Directly target the mutant BMPR2 protein" is a biological dead end. The viable direction is **signaling-axis restoration/activation** (activate the residual wild-type allele, modulate upstream type I receptors, restore downstream SMAD1/5/8) — not a de novo ligand for the mutant receptor itself.

## Check 3 — Constructive alternative: a validated small-molecule precedent exists

Literature:
- **FK506 (tacrolimus)** — identified as a strong BMP-signaling activator in a screen of ~3,500 FDA-approved drugs (Spiekerkoetter et al., JCI 2013).
- Mechanism: **releases FKBP12 from type I receptors (ALK1/ALK2/ALK3), activating downstream BMPR2 signaling (SMAD1/5, MAPK, ID1)** → reverses occlusive vasculopathy in rodent PAH models.
- **TransformPAH** trial (Stanford, Spiekerkoetter) — low-dose FK506 safe, increased BMPR2 expression in a subset of patients (ERJ 2017).

Implication: The concept of "small-molecule restoration of the BMPR2 signaling axis" is sound and clinically precedented. The target, however, is not the mutant BMPR2 protein but the **FKBP12–type I receptor interaction** (or another signaling node). The real white space for de novo design is a **more selective BMPR2-axis activator than FK506**.

---

## Synthesis

| Proposal premise | Verdict |
|---|---|
| Zero remodeling-reversal drugs | ✗ Sotatercept (2024) already approved |
| De novo ligand for mutant BMPR2 | ✗ Most variants cause protein loss (haploinsufficiency) |
| AI molecule-generation tools mature | ✓ (DiffSBDD/AF3 etc. real — licensing/validation separate) |
| Small-molecule intervention on BMPR2 axis feasible | ✓ FK506 precedent, in clinical trials |

**Recommendation:** Do not execute the proposal as written. Redefine the target hypothesis from "mutant-BMPR2 pocket binder" to **"small molecule that restores the BMPR2/BMP signaling axis (e.g., FKBP12–ALK interface, or wild-type BMPR2 stabilizer)"**, and reframe positioning against sotatercept (oral, small-molecule, target selectivity) before starting Phase 1.

## Sources
- Merck/FDA: Winrevair (sotatercept) approval, 2024-03 — merck.com news; PubMed 39466552
- BMPR2 mutation spectrum: Springer Respir Res 10:87 (PMC2762975); AHA CircRes (CIRCRESAHA.122.321554)
- FK506/BMPR2: JCI 2013 (jci.org/articles/view/65592, PubMed 23867624); TransformPAH ERJ 2017 (PubMed 28893866); NCT01647945
