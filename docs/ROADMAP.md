# LUMEN Roadmap

Last updated: 2026-06-28
Target: FKBP12 (P62942) — non-immunosuppressive small-molecule activator of the
BMP/BMPR2 signaling axis for PAH. See ADR-001 and the FKBP12 target dossier.

This roadmap is honest about what is done, what is computationally achievable
next, and what genuinely requires resources beyond a laptop + Colab (and how to
get them). Mapped against the original `proposal.html` 4-phase plan.

---

## Where we are (2026-06-28)

| Proposal phase | Status | Evidence |
|---|---|---|
| Premise validation (added) | ✅ done | premise-validation + verification-audit logs |
| P1 — target structure & pocket | ✅ done (experimental PDB, not AF3; no MD) | dossier §3–4, interfaces |
| P2 — generative design + scoring | 🟡 core done | A (BRICS) + B (DiffSBDD) ran; ML scorer validated |
| P2 — RL / Pareto / full ADMET | ⬜ not started | — |
| P3 — synthesis + cell + animal | ⬜ needs wet lab | — |
| P4 — GLP tox + IND | ⬜ needs CRO + funding | — |

Deliverable so far: 8 structurally-promising, novel de novo FKBP12 candidates
(`results/phase2b/diffsbdd_promising.csv`) — testable hypotheses, not drugs.

Key validated methodology facts (the spine of credibility):
- docking = poses only (enrichment ROC-AUC 0.455, no affinity ranking);
- ML potency scorer works (triage scaffold-CV 0.993, potency ρ 0.76);
- selectivity is NOT ligand-ML-predictable → structural + experimental;
- calcineurin-sparing is an efficacy/toxicity trade-off (GPI-1046 supports feasibility);
- FKBP12.6 counter-screen is safety-critical.

---

## Track A — Computational, achievable now (laptop + Colab GPU)

Ordered by value. All produce real, publishable methodology.

1. **MM-GBSA / FEP rescoring of B candidates** (Colab GPU)
   Better-than-docking physics-based affinity estimate for the 8 promising B
   molecules. Goal: replace weak docking score with a defensible ranking.
   Verify: rescoring separates known binders from decoys better than ROC-AUC 0.455.

2. **RL / multi-objective optimization of B scaffolds** (Colab GPU, REINVENT4)
   Take B's novel scaffolds and optimize for predicted potency + QED + logD,
   keeping the applicability-domain gate. Fills the proposal's P2 RL/Pareto gap.
   Verify: optimized set beats current B on potency proxy without leaving domain.

3. **Full ADMET profiling** (local / API)
   ADMET (absorption, hERG, CYP, clearance) for promising candidates.
   Verify: candidates pass basic developability filters.

4. **Pocket-conditioned generation at scale + tightened filters** (Colab GPU)
   Larger DiffSBDD run (1000s), filter logP ≤ 4, then ML + docking + AD gate.
   Verify: a clean shortlist of ~20 diverse, drug-like, pocket-fitting molecules.

5. **AlphaFold3 / MD where it adds value** (Colab GPU)
   AF3 of FKBP12–ALK complexes (not just experimental analogs); short MD to check
   pose stability of top candidates. Closes the proposal's P1 MD gap.

6. **Synthetic accessibility + retrosynthesis** (local, AiZynthFinder/SA score)
   Filter candidates a chemist could actually make. Verify: SA score / routes exist.

---

## Track B — Requires resources beyond compute (and how to get them)

These cannot be done on a computer. The honest path to each:

1. **Wet-lab synthesis + binding assay (the real validation).**
   Needs a medicinal-chemistry lab. Path: (a) collaborate with an academic group
   (e.g., PAH/BMP labs — Spiekerkoetter/Stanford, Imperial College, domestic
   PAH registries); (b) contract research org (CRO) for synthesis + FKBP12 SPR/ITC
   binding — costs $ per compound but is buyable; (c) FKBP12 + FKBP12.6 counter-screen.

2. **Cell efficacy (PASMC SMAD1/5 readout).** Academic collaborator or CRO.
3. **Animal models (MCT/SuHx).** Animal facility + ethics approval.
4. **GLP tox / IND.** CRO + substantial funding (grant or investor).

How determination actually translates here: not "do it on the laptop," but
**produce computational results strong enough to win a collaborator or grant.**
That is the lever. A tight, honest dossier + novel validated-in-silico candidates
is what gets a wet lab to say yes.

---

## Immediate next step (recommended)
Track A #1 (MM-GBSA rescoring) — it directly attacks our weakest link (no
trustworthy affinity ranking) and strengthens the case for the 8 B candidates.

## Standing principles (do not drop)
- Commit every work unit; English; private repo.
- Verify before believing; surface limitations; no "world-first" inflation.
- Every candidate carries its applicability-domain caveat. In-silico = hypothesis.
