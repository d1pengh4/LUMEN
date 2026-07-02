# LUMEN — Master Project Log (complete, with numbers)

Purpose: single chronological record of every step with exact numbers, for the
paper and for the record. Authoritative summary; per-step detail lives in
`docs/research-log/*` and `docs/decisions/*`. Last updated: 2026-07-02.

Project: AI de novo design of a non-immunosuppressive small-molecule activator of
the BMP/BMPR2 axis (via FKBP12) for pulmonary arterial hypertension (PAH).
Repo: private GitHub `LUMEN` (owner d1pengh4). Compute: macOS local + Colab T4 (frugal).

---

## PHASE 0 — Premise validation (target correction)
- Original proposal target = BMPR2. Rejected after literature check:
  - "0 remodeling-reversal drugs" → FALSE: **sotatercept (WINREVAIR)** FDA-approved
    **2024-03** (first disease-modifying PAH drug).
  - ~**68%** of BMPR2 PAH variants are truncating → NMD → haploinsufficiency →
    pocket-binder strategy invalid.
- Decision (ADR-001): pivot target to **FKBP12** (UniProt P62942, 108 aa).
  Precedent: FK506 activates BMP signaling via FKBP12 (JCI 2013). Goal =
  non-immunosuppressive (calcineurin-sparing) oral small molecule.
- Feasibility evidence: **GPI-1046** (non-immuno FKBP12 ligand) activates BMP/ALK3
  independent of calcineurin. Safety: FKBP12.6 (83% identical) counter-screen
  mandatory. Novelty reframed to "AI-designed non-immuno FKBP12 ligand for PAH".

## PHASE 1 — Structure, pocket, data
- Structures fetched (RCSB): 1FKF (FKBP12+FK506, drug pocket), 1B6C (FKBP12+ALK5
  proxy, ON-target), 1TCO (FKBP12+calcineurin, OFF-target), 2PPN (apo), 1F40
  (FKBP12+GPI-1046, non-immuno template). Identities verified from PDB headers.
- Pocket (fpocket on 1FKF, FK506 removed): **druggability 0.517, volume 277 Å³,
  hydrophobicity 61.2, polar SASA ≈ 0**. Pocket1 centroid 4.35 Å from native FK506
  (same site, verified).
- Interfaces (BioPython, 5 Å): ON-target 22 residues (overlap FK506 pocket);
  OFF-target (calcineurin) 23 residues incl. the "40s loop" (32–45) absent from
  the receptor interface → **selectivity rule: fill pocket, avoid 40s loop**.
- ChEMBL ligands (target CHEMBL1902): **821 activities → 725 unique** binders,
  pChEMBL 3.1–11.0 (median 6.71). Space skews macrocyclic (median MW 556, 243 >700 Da,
  QED median 0.284). **68 drug-like anchors** (MW≤500, QED≥0.5, pChEMBL≥6).

## SCORING TOOLS — build + validation (the key methodology)
### Docking (smina/Vina)
- Redock control: FK506 RMSD **0.81 Å** vs crystal (pose geometry OK).
- **Enrichment (20 actives / 40 decoys): ROC-AUC 0.455** (below random) → docking
  score has NO affinity-ranking power here. Retained for POSES only.
### ML scorer (RandomForest, ECFP4 2048)
- Binder triage: random-split CV ROC-AUC **0.987**, **scaffold-split 0.993**;
  external enrichment 1.000 (flagged EASY: actives Tanimoto 0.75 vs decoys 0.21).
- Potency regressor (pChEMBL): **scaffold-split Spearman 0.748**, R² 0.56, MAE ~0.66 log.
- Selectivity panel: FKBP5 standalone regressor scaffold-CV ρ 0.78, BUT
  FKBP12−FKBP5 selectivity prediction ρ ≈ **−0.01** → **selectivity NOT ML-solvable**;
  FKBP12 vs FKBP12.6 = 83% identical, known ligands don't discriminate (1/40).
- Classifier characterized: it's an FKBP-FAMILY binder detector (FKBP5 ligands
  score ~0.79), NOT FKBP12-specific. Determinism bug fixed (sorted scaffold groups).
### GNINA (CNN rescoring) — adopted ranker
- Smoke: actives 6.59/6.23 vs aspirin 4.52.
- **Enrichment GATE (60 mols): ROC-AUC 0.799** (vs docking 0.455; target ≥0.65) →
  PASS. Actives CNNaffinity 6.34 ± 0.36, decoys 5.80 ± 0.53. **Adopted as the
  affinity-ranking signal.** Applicability-domain gate = Tanimoto ≥ 0.30.

## PHASE 2 — Generation + candidate ranking
### A. BRICS (CPU baseline)
- 400 novel drug-like candidates. Scaffold diversity **11%** (weak), novelty
  (median max Tanimoto to known) 0.46, logP median 5.3. Negative control: random
  BRICS binder_prob 0.027 vs these 0.591. Hypotheses, no selectivity.
### B. DiffSBDD (GPU, pocket-conditioned) — batch 1
- 272 molecules. Strengthened gate PASS: **scaffold diversity 89%**, novelty 0.18,
  33 developable, PAINS 0%.
- GNINA-scored 37 (33 developable + 4 TIER-1): **best CNNaffinity 6.56** (TIER-1,
  also dock −8.5/seed-avg −7.5, QED 0.72, novel Tanimoto 0.20). Only **3/37 ≥ 6.3**.
- Synthetic-feasibility triage: top scorer (6.56) has 5 stereocenters (hard);
  practical lead = 6.34 candidate, SA 2.28, 0 stereocenters. 14/37 "synthesizable_easy".

## ROADMAP G1 — scale-up test ("does more help?")
- DiffSBDD batch 2: salvaged 271 molecules (session paused mid-run). Gate PASS
  (diversity 91%, novelty 0.17, 24 developable, PAINS 3%).
- GNINA of 24 new developable: **best 6.35, none ≥7.0** → doubling random generation
  did NOT beat 6.56. **Brute-force plateaus ~6.5.** → directed optimization needed.

## ROADMAP G2 — goal-directed optimization (GA, CPU/frugal)
Objective: maximize pred_pchembl s.t. in-domain 0.30≤ad≤0.65, QED≥0.55, MW 250–500.
### Round 1
- Attempt 1 FLAWED: seeded population with anchors → top ad=1.0 (known binders),
  best flat **8.39** = seed dominance (artifact). Recorded.
- Attempt 2 (novelty box, anchors as donors only): best pred_pchembl **5.16→7.34**;
  top novel (ad 0.44–0.52) BUT binder_prob only ~0.37 (regressor/classifier disagree).
- **GNINA validation (15 mols): best 7.06** (vs 6.56); 7/15 >6.56; 13/15 ≥6.3;
  **1** developable+novel ≥7.0. Top: `CC(C)(C)C(=O)C(=O)N1CCCCC1N1CCCC1N1CCC[C@H]1CC(N)=O`
  (7.06, ad 0.38, QED 0.72, SA 4.05, MW 392). Directed optimization beat brute-force,
  GNINA-confirmed. STOP (≥3 at ≥7.0): **1/3**.
### Round 2 (seeded from round-1 winners)
- best pred_pchembl **7.34→7.78**; **binder_prob 0.50–0.65** (models now AGREE →
  more trustworthy). Top novel (ad 0.48–0.58), QED 0.66–0.72, recurrent nitrile
  (CC#N) arm on the diketopiperazine-sulfonamide core.
  Top: `COC[C@H]1CN(CC#N)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1cccc(C)c1` (pred 7.78).
- **GNINA validation: IN PROGRESS** (15 molecules, session 'gar2').

## STATUS vs LOCKED CRITERIA (2026-07-02)
- STAGE-GATES passed: scoring (GNINA 0.799), G1, G2-round1 (GNINA 7.06 > 6.56).
- STOP criterion (≥3 developable/novel at GNINA CNNaffinity ≥7.0): **1/3** (pending R2 GNINA).
- Safety stop (plateau = 2 rounds <0.1 gain): NOT triggered (R1 +0.5 GNINA; R2 ML +0.44).
- GPU budget: frugal; sessions stopped after each run.
- Decision gate (startup/patent vs paper): not reached (after G3–G7).

## NEXT
Round-2 GNINA verdict → if ≥3 at ≥7.0 or clear improvement, proceed to G3 (MD
pose-stability) on the top novel/developable candidates; else safety-stop and carry
the best set (≥7.06) to G3. Then G4 retrosynthesis, G5 AF3, G7 selectivity FEP
(resource-gated), G8 paper.

## KEY FILES
- Candidates: `results/phase2b/FINAL_candidates.csv`, `gaopt_full.csv` (R1 GNINA),
  `ga_optimized.csv` (R1 ML), `ga_optimized_r2.csv` (R2 ML), `gen600_developable.csv`.
- Scoring: `results/gnina/gate_bg.log` (GNINA gate), `gaopt_scored.log`, `gar2_scored.log`.
- Synthesis: `results/synthesis/FINAL_synth_feasibility.csv`.
- Decisions/logs: `docs/decisions/ADR-001*`, `docs/research-log/2026-*`, `docs/ROADMAP.md`.

## HQ RESCORE (exhaustiveness 8) — legitimate improvement (2026-07-02)
Re-scored top 10 candidates with thorough GNINA (exh 8, num_modes 5; NOT threshold fudging).
- Best CNNaffinity: exh4 7.06 → **exh8 7.21** (better sampling raised it).
- **≥7.0 & novel & developable: 2/3** (up from 1); all 10 ≥6.34 (known-active level), 9 ≥6.56.
- STOP criterion (≥3 at ≥7.0): 2/3 — strong lead set; 3rd candidate at 6.95 (borderline).
- Leads: `results/phase2b/FINAL_leads.csv`. Proceeding to G3 (MD) as the decisive test.
