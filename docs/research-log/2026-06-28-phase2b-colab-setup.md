# Phase 2-B — GPU generation via Colab CLI (DiffSBDD)

Date: 2026-06-28 · Related: phase2-poc (A baseline), target-dossier §4 (pocket)
Goal: replace the BRICS CPU baseline (A) with a real pocket-conditioned 3D
diffusion model (DiffSBDD) on a Colab GPU, then score with the SAME validated
pipeline so B is comparable to A head-to-head.

## Tooling
- `google-colab-cli` (installed via `uv tool install google-colab-cli`).
  `colab run --gpu L4 --timeout <s> script.py` runs a local .py on a fresh GPU
  VM and returns stdout; the VM is torn down automatically.

## Division of labor
- USER (only they can): `colab auth` (interactive Google OAuth), and `colab pay`
  if GPU compute units are needed.
- ME / automatable: the generation script, running it, parsing + scoring output.

## Files
- `src/generation/colab_diffsbdd_run.py` — self-contained GPU job: fetches 1FKF,
  builds the FKBP12 pocket (chain A + FK506 A:108), clones DiffSBDD, downloads
  the full-atom conditional checkpoint (Zenodo 8183747), generates N molecules,
  prints SMILES between `===LUMEN_SMILES_BEGIN/END===` markers.
- `src/generation/score_colab_output.py` — parses that log, scores with
  FKBP12Scorer, and prints the A-vs-B audit (diversity / novelty / logP).
- `src/generation/prep_diffsbdd_input.py` + `data/processed/diffsbdd_input/` —
  local copy of the pocket input (for reference / manual runs).

## Run sequence
```
colab auth                                   # user, one time
colab run --gpu L4 --timeout 2400 \
    src/generation/colab_diffsbdd_run.py 200 | tee results/phase2b/colab.log
python src/generation/score_colab_output.py results/phase2b/colab.log
```

## Honest risks (smoke-test run #1 first, n=50-100)
1. **Dependency/CUDA mismatch (most likely failure):** DiffSBDD pins torch 1.12.1
   / CUDA 10.2; Colab ships CUDA 12 / torch 2.x. The script installs torch-scatter
   matched to the VM torch, best-effort. If `generate_ligands.py` imports fail,
   that's expected on run #1 → iterate on pinned versions.
2. **GPU availability / compute units:** free Colab may not grant a GPU; may need
   Pro / units (`colab pay`).
3. **Checkpoint download:** Zenodo URL must resolve (~hundreds of MB).

## Acceptance criteria for B (vs A baseline)
- generation completes, SMILES returned;
- scaffold diversity clearly > 11% (A) and novelty ≥ ~0.46;
- then apply tightened physchem filters (logP ≤ 5, rot ≤ 10) and the FKBP12Scorer.
B output, like A, remains HYPOTHESES TO TEST, not drug candidates.
