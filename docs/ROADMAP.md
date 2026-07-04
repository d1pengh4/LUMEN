# LUMEN Roadmap (DEFINITIVE, v2)

Last updated: 2026-07-02
Target: FKBP12 (P62942) — non-immunosuppressive small-molecule activator of the
BMP/BMPR2 signaling axis for PAH. See ADR-001, target dossier.

This is the authoritative plan. It encodes three decisions the user locked in:
1. **Stop criterion** = target score reached.
2. **Validation scope** = ALL computational gates (MD, selectivity FEP, retrosynthesis, AF3).
3. **GPU budget** = frugal (free/cheap Colab T4 first).

An honest tension: "all gates incl. FEP" vs "frugal GPU" — FEP is the most
GPU-expensive step and clashes with a free-tier budget. Resolution is written
into step G7 (FEP only on the final 1–2 candidates; if free-tier can't do it,
it becomes an explicit resource-gated step, not silently skipped).

---

## GOAL / DEFINITION OF DONE
Produce a **defensible, honestly-validated shortlist of de novo FKBP12 candidates**
+ a reproducible pipeline, sufficient for a **methodology paper** and to attract a
wet-lab collaborator. "Done" = STOP criterion met AND all gates run (or explicitly
resource-gated) AND paper draft written.

## GATE STATUS (2026-07-04) — see docs/PROJECT_LOG.md for numbers
- G1 scale-up: DONE (brute-force plateaus ~6.5). G2 GA optimization: DONE (GNINA 7.06→7.21).
- **G3 MD pose-stability: PASS** — L1/L2/L3 all stay bound over 1 ns (RMSD max ≤ 2.9 Å).
- **G4 retrosynthesis: CONDITIONAL** — leads' novel cores need custom synthesis; pool
  screen: 7/37 synthesizable but only ≤ GNINA 6.34, and the best of those fails MD.
- **G7 selectivity feasibility: DONE (structural)** — FKBP12 vs FKBP12.6 pockets differ
  at only 2 conservative residues (Trp59, Ile90); L1–L3 contact them → selectivity is
  addressable for our leads. Full FEP still **resource-gated** (free-tier can't run it).
- Pending: G5 AF3 (supporting, interactive), G7 FEP (paid/borrowed GPU), G8 paper.
- STOP criterion (≥3 developable/novel at GNINA ≥7.0): 2/3 strictly (L3 at 6.95); all 3
  MD-stable. Decision gate (startup/patent vs paper): approaching — hinges on FEP/wet-lab.

## STAGE-GATE PRINCIPLE (unchanged, mandatory)
No step starts until the previous PASSES its concrete numeric/artifact gate.
A weak foundation collapses the stack. Every gate outcome recorded in a research log.

---

## STOP CRITERION (the generate↔optimize loop)
Primary: **≥ 3 candidates with GNINA CNNaffinity ≥ 7.0** that are also developable
(SA ≤ 4, QED ≥ 0.6, logP 1–4, no PAINS) and novel (Tanimoto < 0.4 to known).
- Rationale: gate actives averaged 6.34 ± 0.36; ≥ 7.0 is clearly top-tier ("better
  than a typical known binder" by our validated scorer).
Safety stop (prevents infinite loop): if **2 consecutive RL rounds do not raise the
best CNNaffinity by ≥ 0.1**, declare a plateau, stop, and carry the best achieved
(documented honestly). Frugal budget also caps total generation at ~3–4k molecules.

---

## PIPELINE (build order, with gates)

**Done so far:** ① target pivot ② structure/pocket/ligands ③ scoring built +
GNINA gate PASS (ROC 0.799) ④ generation (272) + first candidate ranking
(top CNNaffinity 6.56). Now scaling generation (600, in progress).

### G1 — Scale generation + GNINA rescore  [IN PROGRESS]
Generate more (→600 now; more if needed), filter (strengthened developability),
GNINA-rescore, rank. GPU: T4, batched, detached (frugal-OK).
**GATE:** pool passes strengthened diversity/novelty/developability gate; produce
ranked CNNaffinity table. Check: did the best score improve over 6.56?

### G2 — RL / multi-objective optimization (REINVENT4 or DiffSBDD-guided)  [NEXT]
"Smart" loop: learn from high-scoring molecules, steer generation toward higher
CNNaffinity + QED + logD, keeping the applicability-domain gate.
**GATE:** best CNNaffinity rises toward the STOP target (≥7.0), diversity kept,
stays in-domain. Loop G1↔G2 until STOP criterion or plateau.
GPU: T4 detached. Frugal: cap rounds; stop on plateau.

### G3 — MD pose-stability (OpenMM)  [after STOP criterion met]
Short explicit-solvent MD (a few ns) of the top ~3–5 candidate–FKBP12 complexes.
**GATE:** ligand stays bound (RMSD stable, key pocket contacts retained). Unstable
poses dropped — separates real binding from lucky docking poses.
GPU: T4, short runs, only top candidates (frugal-OK). Honest risk: ligand
parameterization (GAFF/OpenFF) setup may need iteration.

### G4 — Synthetic accessibility + retrosynthesis (AiZynthFinder)  [mostly CPU]
Confirm a chemist could make the finalists.
**GATE:** each finalist has SA ≤ ~4 AND a plausible retrosynthetic route. Cheap.

### G5 — AlphaFold3 of FKBP12–ALK complexes (AF Server / ColabFold)  [supporting]
Model the real ON-target interface (not just the ALK5 proxy) to strengthen the
mechanism narrative. **GATE:** none (supporting evidence). Free (AF Server web).

### G6 — MM-GBSA cross-check (optional, if budget allows)
Second physics-based opinion on the finalists to complement GNINA.
**GATE:** correlates with GNINA ranking on finalists (or documented divergence).

### G7 — Selectivity by relative FEP: FKBP12 vs FKBP12.6  [EXPENSIVE / resource-gated]
The ML-unsolvable limit. Precedent: MD+FEP resolves FKBP51/52 (83%-identical).
**GATE:** FEP reproduces the FKBP12-vs-12.6 selectivity sign on known reference
ligands (retrospective control) BEFORE trusting it on candidates.
**Frugal reconciliation:** run FEP ONLY on the final 1–2 top candidates. If the
free/cheap tier cannot complete FEP (likely), this step is **explicitly
resource-gated**: document FKBP12.6 selectivity as an open risk requiring paid
GPU or a collaborator — NOT silently skipped. This honors "all gates" honestly.

### DECISION GATE (after G3–G7, before G8) — data decides the endgame
Purpose stays OPEN until here; the accumulated data picks the path (user decision
2026-07-02). Concrete branch:
- **STRONG result** (STOP criterion met: ≥3 developable/novel candidates at
  CNNaffinity ≥ 7.0, surviving MD stability, synthesizable, selectivity at least
  not-ruled-out) → pursue **startup / patent** track: prioritize IP (provisional
  patent on the lead scaffolds), seek a wet-lab collaborator/CRO for experimental
  validation, then paper. Experimental data is the gate to investment.
- **LIMITED result** (plateaus below target, or candidates fail MD/selectivity, or
  scoring-ceiling reached) → wrap up as a **methodology paper**: the honest pipeline
  + negative findings are themselves the contribution. No overclaiming.
Either branch is a legitimate outcome by design — the roadmap makes "it didn't
reach a drug" a publishable finding, not a failure.

### G8 — Methodology paper (both branches converge here)
Write the paper: honest pipeline + negative findings (docking fails, ML selectivity
fails, GNINA works) + validated candidate shortlist + all gate results (incl.
whichever were resource-gated). In the STRONG branch, this accompanies IP/collab;
in the LIMITED branch, it is the final deliverable. Attracts wet-lab collaboration
either way.

---

## Track B — beyond compute (unchanged; the real validation)
Synthesis + binding assay + FKBP12.6 counter-screen → cell (PASMC SMAD1/5) →
animal (MCT/SuHx) → GLP tox → IND. Requires a wet lab / CRO / funding. The paper
(G8) is the lever to obtain these.

## Standing principles
- Stage-gate: pass the gate before advancing; no skipping.
- Commit every work unit; English; private repo; reply to user in Korean.
- Verify before believing; surface limits; no "world-first" inflation.
- Every candidate carries its applicability-domain + "hypothesis-not-drug" caveat.
- Frugal GPU: T4 first, batched, detached; escalate cost only where a gate demands it.
