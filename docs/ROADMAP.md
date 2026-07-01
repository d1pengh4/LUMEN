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

### STAGE-GATE PRINCIPLE (mandatory — "a weak foundation collapses")
No step begins until the previous step PASSES its gate. Each gate is a concrete,
verified numeric/artifact criterion, not a feeling. If a step fails its gate, we
fix or redesign it (or honestly declare the limit) before advancing. This is why
we already caught: docking (ROC 0.455 → rejected), ML selectivity (ρ≈0 → rejected).
Record every gate outcome in a research-log entry.

Order below is the build order. Steps done so far: developability panel (step done),
scaled generation (in progress). Physics-based steps (1,5,7) are the core of
overcoming the "computer ≠ reality" and "selectivity" limits — literature confirms
ABFE/FEP improves enrichment over docking, and MD+FEP resolves selectivity between
closely related FKBP paralogs (FKBP51/52 precedent).

1. **MM-GBSA rescoring of candidates** (Colab GPU)
   Physics-based affinity estimate to replace the weak docking score.
   **GATE:** on our enrichment set (20 actives/40 decoys), MM-GBSA ROC-AUC must
   beat docking's 0.455 (target ≥ 0.65). If it fails, do NOT use it for ranking —
   escalate to ABFE/FEP or declare the limit. Only after PASS, rescore candidates.

2. **RL / multi-objective optimization of B scaffolds** (Colab GPU, REINVENT4)
   Optimize B's novel scaffolds for potency + QED + logD, keeping the AD gate.
   **GATE:** optimized set beats the current pool on the (validated) scorer without
   leaving the applicability domain; diversity preserved.

3. **Full ADMET profiling** (local / API)
   **GATE:** a defined developability filter (SA≤4.5, QED≥0.5, Lipinski+Veber,
   no PAINS, hERG/CYP flags) applied; shortlist passes.

4. **Pocket-conditioned generation at scale + tightened filters** (Colab GPU)
   **GATE (strengthened):** scaffold diversity > 70%, novelty (median Tanimoto to
   known) < 0.35, ≥ 15 developable (SA ≤ 4.0, QED ≥ 0.6, 1 ≤ logP ≤ 4,
   Lipinski+Veber, no PAINS), PAINS fraction < 5%. Efficacy is NOT gated here
   (unreliable for out-of-domain novel mols) — verified later by MM-GBSA.

5. **MD pose-stability of top candidates** (Colab GPU, OpenMM)
   Short explicit-solvent MD of the top complexes.
   **GATE:** candidate stays bound (ligand RMSD stable, key contacts retained)
   over the run; unstable poses are dropped. This is the real test that a novel
   molecule actually occupies the pocket (vs a lucky docking pose).

6. **Synthetic accessibility + retrosynthesis** (local, AiZynthFinder/SA score)
   **GATE:** each final candidate has SA ≤ ~4 AND a plausible retrosynthetic route.

7. **[NEW] Selectivity by relative free energy (FKBP12 vs FKBP12.6)** (Colab GPU, FEP)
   The limit ML could not solve. Literature precedent: MD+FEP resolves FKBP51/52
   selectivity (83%-identical paralogs), the same problem class as FKBP12/12.6.
   **GATE:** FEP reproduces the FKBP12-vs-12.6 selectivity sign/magnitude on a few
   known reference ligands (retrospective control) BEFORE trusting it on candidates.
   If the control fails, selectivity stays a wet-lab counter-screen (honest limit).

8. **AlphaFold3 of FKBP12–ALK complexes** (Colab GPU / AF Server)
   Model the actual ON-target interface (not just the ALK5 proxy). Supporting, not gating.

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

## Immediate next step
Finish scaled generation (step in progress) and lock its gate (a clean,
diverse candidate pool). ONLY THEN start Track A #1 (MM-GBSA), which must pass
its enrichment gate before we trust it on candidates.

## Standing principles (do not drop)
- **Stage-gate:** each step must PASS its numeric/artifact gate before the next
  begins. A weak foundation collapses the whole stack. No skipping ahead.
- Commit every work unit; English; private repo.
- Verify before believing; surface limitations; no "world-first" inflation.
- Every candidate carries its applicability-domain caveat. In-silico = hypothesis.
