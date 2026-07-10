# LUMEN — Path to experimental validation & the startup case

Purpose: turn the computational results into an actionable, fundable plan for a wet-lab
program (goal = real drug / company, not just a paper). Grounded in project results
(`docs/PROJECT_LOG.md`); costs/timelines are rough order-of-magnitude estimates a CRO
would refine. Candidates are validated *hypotheses*, not confirmed binders.
Last updated: 2026-07-10.

## 0. Thesis in one paragraph
FK506 activates the BMP/BMPR2 axis (relevant to PAH) via FKBP12 but is immunosuppressive
because its effector arm recruits calcineurin. We designed de novo small molecules that
**fill only the FKBP12 pocket and avoid the calcineurin-recruiting 40s loop** — i.e.
BMP-axis engagement without immunosuppression. Three leads are computationally de-risked
(validated affinity signal + MD-stable + novel + contact the FKBP12-vs-FKBP12.6
selectivity residues). The value of the whole program hinges on one experiment class:
**does a lead bind FKBP12, spare calcineurin (non-immunosuppressive), and discriminate
FKBP12.6 (cardiac safety)?**

## 1. What is already de-risked (computational)
- **Target rationale** corrected and evidenced (BMPR2 haploinsufficiency → FKBP12 pivot;
  FK506/GPI-1046 precedent).
- **Affinity**: leads score above the known-active mean on GNINA (the one scorer we
  validated, enrichment ROC-AUC 0.799).
- **Real binding vs lucky pose**: all three leads stay bound in 1 ns MD (RMSD ≤ 2.9 Å).
- **Selectivity is structurally addressable**: FKBP12 vs FKBP12.6 pockets differ at only
  two conservative residues (Trp59, Ile90); the leads contact them.
- **Novelty** (IP-relevant): Tanimoto < 0.4 to known FKBP12 ligands; novel core scaffolds.

## 2. Lead triage — what to make and test first
| Lead | GNINA | MD (mean/max Å) | Selectivity handles contacted | Synthesis (AiZynth route) | QED / MW / SA |
|------|-------|-----------------|-------------------------------|---------------------------|---------------|
| **L3** | 6.95 | **1.05 / 1.60** (best) | **Trp59 + Ile90** | 4 steps (shortest), 4/5 catalog + custom core | 0.67 / 454 / 4.39 |
| **L2** | 7.04 | 1.70 / 2.91 | **Trp59 + Ile90** | 6 steps, 4/5 + custom oligo-pyrrolidine core (hardest) | **0.72 / 392 (smallest) / 4.05** |
| L1 | **7.21** (best) | 1.15 / 1.70 | Trp59 only | 5 steps, 5/6 + custom core | 0.67 / 454 / 4.34 |

**Recommendation:** synthesize **L3 first** (best MD stability, both selectivity handles,
shortest route), then **L2** (most drug-like/smallest, both handles). Keep **L1** as the
high-affinity comparator. Rationale for goal (2): because selectivity is the thesis,
prioritize the leads that engage *both* differing residues (L2, L3) over raw affinity (L1).

## 3. Synthesis (to get compound in hand)
AiZynthFinder confirms the **peripheral chemistry is routine** (sulfonylation, amide
coupling, reductive amination, Suzuki — 4/5–5/6 building blocks are catalog). The one
custom step per lead is the **novel bicyclic proline-derived core** (L1/L3) or the
oligo-pyrrolidine chain (L2). Action: request a medicinal-chemistry CRO quote to make
L3 + L2 (~10–50 mg each), flagging the core as custom asymmetric synthesis.
Rough estimate: **~$8–25k and 2–4 months** for 2 leads incl. custom core route dev.

## 4. Experimental validation cascade (each tier is a go/no-go)
**Tier 1 — Binding & potency (does it bind FKBP12 at all?)**
FK506-competition fluorescence-polarization assay and/or the classic FKBP12 PPIase
(rotamase) inhibition assay; SPR/ITC for Kd on a hit. *De-risks:* the core "does it
bind" question. *Est.:* ~$5–15k, 4–8 wks. **GO if measurable FKBP12 binding (µM or better).**

**Tier 2 — THE make-or-break: selectivity (two counter-screens)**
(a) **Non-immunosuppression:** calcineurin phosphatase / NFAT-luciferase (or T-cell IL-2)
assay with compound+FKBP12 — the compound must **not** recruit calcineurin (this is the
whole differentiator vs FK506; our design avoids the 40s loop, so this is the direct test
of the design rule). (b) **FKBP12.6 counter-screen:** repeat Tier-1 binding on FKBP12.6 —
tests the Trp59/Ile90 selectivity hypothesis; want FKBP12 ≫ FKBP12.6 (cardiac safety).
*De-risks:* the entire thesis. *Est.:* ~$10–25k, 6–10 wks. **GO if non-immunosuppressive
AND some FKBP12/12.6 window.**

**Tier 3 — Mechanism/efficacy (does it do the therapeutic thing?)**
BMP/SMAD1/5 phosphorylation or BRE-luciferase reporter in pulmonary-artery smooth-muscle
cells (PASMC) or C2C12; show BMP-axis activation like FK506/GPI-1046 without
immunosuppression. *Est.:* ~$15–40k, 2–3 mo. **GO → lead optimization + IND-enabling.**

**Tier 4 (later):** ADMET/PK, then in-vivo PAH models (monocrotaline / SuHx), GLP tox.

## 5. IP posture (for the company case)
Leads are novel (Tanimoto < 0.4; the bicyclic proline-derived sulfonamide and
oligo-pyrrolidine cores are not in the known FKBP12 ligand set). A **composition-of-matter
provisional** on the scaffold genus is plausible *before* public disclosure (note: the
methodology paper is a disclosure — file first if the startup track is chosen). What a
patent attorney needs from us: the scaffold Markush, the synthesized exemplars + any
assay data, and a prior-art list. **Sequence the paper vs the filing deliberately.**

## 6. The ask (minimum to reach the decisive data)
The first truly decisive, investable milestone is **Tier 1 + Tier 2** (binding +
selectivity). Minimum viable program to get there: **synthesis of L3 & L2 + Tier 1/2
assays ≈ $25–65k and ~4–6 months** via a CRO or an academic collaborator (FKBP/PPIase or
PAH lab). That data set — "binds FKBP12, non-immunosuppressive, FKBP12.6-selective" —
is the gate to seed investment; until then everything is computational.

## 7. Honest risk register
- **Selectivity magnitude unproven** — only 2 conservative pocket residues differ; the
  FKBP12/12.6 window may be small. Highest-impact unknown. (FEP could pre-estimate it but
  is resource-gated and still computational; the assay is more decisive.)
- **No measured binding yet** — GNINA/MD are validated proxies, not affinity.
- **Synthesis of the novel core** — needs custom route dev; a cost/time risk, not a
  scientific blocker (peripheral chemistry is routine).
- **Shallow pocket** (druggability 0.517) → potency may be hard to push.
- **Affinity↔synthesizability trade-off** (Figure 1): the makeable analogues are weaker;
  the potent leads need custom synthesis.

## 8. Decision this enables
If Tier 1–2 pass → **startup/IP track** (provisional patent on the scaffold, seed raise on
the safety-differentiated FKBP12 binder). If they fail → the honest methodology paper
(G8 draft ready) is the fallback deliverable. The next concrete action is a **CRO/collab
quote for synthesis of L3+L2 and the Tier-1/2 assay panel.**
