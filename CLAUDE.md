# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

LUMEN is a computational drug-discovery pipeline, not a conventional software app. Goal: AI de novo design of a **non-immunosuppressive small-molecule activator of the BMP/BMPR2 axis via FKBP12** (UniProt P62942) for pulmonary arterial hypertension. The deliverable is a defensible, honestly-validated candidate shortlist plus a reproducible pipeline — sufficient for a methodology paper and a wet-lab collaborator.

Read these first to get current: `docs/PROJECT_LOG.md` (master numeric record, single source of truth for status), `docs/ROADMAP.md` (locked plan with stage-gates G1–G8), `docs/decisions/ADR-001-target-redefinition.md` (why the target is FKBP12, not BMPR2 as in `proposal.html`), and `docs/target-dossier-FKBP12.md`. The write-up itself lives in `docs/MANUSCRIPT_DRAFT.md` and `docs/paper/manuscript_{EN,KO}.md` (+ `docs/paper/figures/`) — every number in the paper must trace back to `docs/PROJECT_LOG.md`.

## Environment & commands

Two separate conda envs (the compchem stack needs conda, not pip wheels; docking is split out to avoid an icu conflict):
```bash
conda env create -f environment.yml        # main `lumen`: rdkit, biopython, fpocket, scikit-learn
conda env create -f environment-dock.yml   # `lumen-dock`: smina, openbabel
```
Reproduce Phase 1 (structure → data → scoring), from README:
```bash
python src/structure/fetch_structures.py      # download FKBP12 PDBs (1FKF, 1B6C, 1TCO, 2PPN, 1F40)
python src/structure/map_interfaces.py        # ON/OFF-target interfaces -> selectivity rule
python src/data/collect_chembl_ligands.py     # CHEMBL1902 ligands
python src/scoring/dock_benchmark.py          # smina redock RMSD gate (needs lumen-dock env)
python src/scoring/train_affinity_model.py    # trains models/*.joblib (gitignored, reproducible)
python src/scoring/run_enrichment.py          # docking enrichment ROC-AUC (documents the smina negative)
python src/scoring/selectivity_panel.py       # FKBP12 vs paralog ML selectivity check (documents the ML negative)
python src/scoring/score_molecules.py cands.csv > scored.csv   # unified scoring interface
```
There is no test suite, linter, or build step. "Verification" here means a numeric/artifact **gate** (e.g. redock RMSD, enrichment ROC-AUC), recorded in a research log — see the stage-gate principle below.

## Local vs GPU split (important)

This repo is developed on **macOS without a GPU**. The GPU steps — DiffSBDD generation, GNINA rescoring, OpenMM MD — do NOT run locally. They are written as **self-contained `colab_*.py` scripts** (in `src/scoring/`, `src/generation/`, `src/md/`) that are pasted into a Colab T4 session: each installs its own deps, downloads its binary (gnina) or structure (1FKF from RCSB), embeds molecules with RDKit, runs, and prints `STEP`/marker lines that get salvaged back into `results/` logs. Molecule lists are frequently **hardcoded inline** as JSON in these scripts (e.g. the `DATA=[...]` in `colab_gnina_hq.py`) — that is the established pattern, not a smell to fix. GPU budget is deliberately **frugal**: T4, batched, detached, sessions stopped after each run.

## Generation pipeline (the G1–G4 flow across `src/generation/`, `src/scoring/`, `src/md/`)

Two parallel generation routes feed the same downstream optimize/validate loop:
- **Route A — BRICS baseline** (`src/generation/brics_generate.py`, local CPU): fragments known FKBP12 binders and recombines them with RDKit. The honest CPU stand-in for a diffusion model.
- **Route B — DiffSBDD** (`src/generation/prep_diffsbdd_input.py` builds the pocket from 1FKF, `colab_diffsbdd_run.py` runs pocket-conditioned 3D diffusion on a Colab GPU, `score_colab_output.py` parses/scores the returned SMILES). Both routes are scored/gated the same way (`score_and_gate.py` applies the diversity/novelty/developability gate from `docs/ROADMAP.md`).
- **G2 optimization** (`ga_optimize.py`, `_r2.py`, `_synth.py`): a BRICS-recombination genetic algorithm that hill-climbs the fast ML potency scorer under hard developability/applicability-domain constraints (resists reward-hacking), since RL against a GNINA reward is too GPU-expensive.
- **GNINA rescoring** (`src/scoring/colab_gnina_*.py` — `gate`, `smoke`, `candidates`, `gaopt`, `gaopt_r2`, `gen600`, `hq`): near-identical self-contained Colab scripts, one per generation/round; they differ mainly in which molecule list is hardcoded in and exhaustiveness settings. This is the expensive, validated affinity check that GA-optimized/generated molecules must clear.
- **G3 MD validation** (`src/md/colab_md_setup.py` bootstraps micromamba on the Colab VM, `colab_md.py` runs short explicit-solvent OpenMM MD): separates real binding (stable ligand RMSD) from lucky docking poses, run only on the top GNINA survivors (frugal-GPU discipline).
- **G4 synthesis** (`src/scoring/synth_feasibility.py` for SA/screening, AiZynthFinder externally for retrosynthesis): confirms a chemist could plausibly make the finalists; results here can override an otherwise-attractive candidate (see the affinity↔synthesizability trade-off in `docs/PROJECT_LOG.md`).

## Scoring stack — what works, and what deliberately does not

This is the core methodology; do not "fix" the honest negatives:
- **Docking (smina)**: correct poses (redock RMSD 0.81 Å) but **no affinity ranking** (enrichment ROC-AUC 0.455, below random). Used for POSES ONLY.
- **ML potency** (`score_molecules.py`, RandomForest/ECFP4): triage scaffold-CV 0.993, potency Spearman 0.76. Use for ranking, gated by a Tanimoto applicability-domain threshold (`AD_THRESHOLD = 0.30`). The classifier is an **FKBP-family** binder detector, not FKBP12-specific.
- **GNINA CNNaffinity**: the adopted affinity ranker (enrichment ROC-AUC 0.799, gate PASS). STOP-criterion currency is CNNaffinity ≥ 7.0.
- **Selectivity is NOT ML-solvable here** (`selectivity_panel.py`: FKBP12-vs-FKBP5 ρ≈0; FKBP12.6 is 83% identical). `score_molecules.py` deliberately emits no selectivity score. Selectivity is enforced structurally (avoid the calcineurin "40s loop", residues 32–45) and experimentally (FKBP12.6 counter-screen) — never claim a predicted selectivity number.

## Working conventions specific to LUMEN

- **Stage-gate, mandatory**: no step starts until the previous one PASSES its concrete numeric/artifact gate, and every gate outcome is recorded in `docs/research-log/YYYY-MM-DD-*.md`. A weak foundation collapses the stack.
- **Log everything with exact numbers**: after any work unit, update `docs/PROJECT_LOG.md` and add/extend the dated research log. Commit every work unit. Convert relative dates to absolute.
- **Honesty over hype**: every candidate carries its applicability-domain caveat and a "hypothesis, not a drug" framing. No "world-first" inflation; negative findings are results, not failures.
- `results/` and `data/` are gitignored (provenance lives in docs/manifests); `models/*.joblib` are gitignored and regenerated by `train_affinity_model.py`. Key candidate files: `results/phase2b/FINAL_leads.csv`, `FINAL_candidates.csv`.

---

## General coding guidelines

Behavioral guidelines to reduce common LLM coding mistakes. Merge with project-specific instructions as needed.

**Tradeoff:** These guidelines bias toward caution over speed. For trivial tasks, use judgment.

## 1. Think Before Coding

**Don't assume. Don't hide confusion. Surface tradeoffs.**

Before implementing:
- State your assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them - don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

## 2. Simplicity First

**Minimum code that solves the problem. Nothing speculative.**

- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.

Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.

## 3. Surgical Changes

**Touch only what you must. Clean up only your own mess.**

When editing existing code:
- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it - don't delete it.

When your changes create orphans:
- Remove imports/variables/functions that YOUR changes made unused.
- Don't remove pre-existing dead code unless asked.

The test: Every changed line should trace directly to the user's request.

## 4. Goal-Driven Execution

**Define success criteria. Loop until verified.**

Transform tasks into verifiable goals:
- "Add validation" → "Write tests for invalid inputs, then make them pass"
- "Fix the bug" → "Write a test that reproduces it, then make it pass"
- "Refactor X" → "Ensure tests pass before and after"

For multi-step tasks, state a brief plan:
```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

Strong success criteria let you loop independently. Weak criteria ("make it work") require constant clarification.

---

**These guidelines are working if:** fewer unnecessary changes in diffs, fewer rewrites due to overcomplication, and clarifying questions come before implementation rather than after mistakes.