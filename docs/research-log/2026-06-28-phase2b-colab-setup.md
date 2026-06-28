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

---

## RESULT — B works (2026-06-28). Dependency fixes + first generation.

### Fixes required (each found by reading the real error; honest record)
1. **macOS SSL** — Colab kernel WebSocket failed `CERTIFICATE_VERIFY_FAILED`.
   Fix: `SSL_CERT_FILE=$(certifi.where())` → wrapper `src/generation/colabx.sh`.
2. **ipykernel doesn't forward subprocess stdout/stderr** — errors were invisible.
   Fix: `sh()` captures and prints child output.
3. `ModuleNotFoundError: openbabel` → `pip install openbabel-wheel`.
4. `three_to_one` removed in Biopython ≥1.80 (1.79 has no py3.12 wheel) → source
   shim via `protein_letters_3to1` (no downgrade).
5. torch ≥2.6 `weights_only=True` rejects the checkpoint's argparse.Namespace →
   `sitecustomize.py` forcing `torch.load(weights_only=False)` (trusted Zenodo ckpt).

VM: Tesla T4, Python 3.12, torch 2.11+cu128. Checkpoint: crossdocked_fullatom_cond (18 MB).

### A (BRICS) vs B (DiffSBDD), N≈50, scored with the same FKBP12Scorer
| Metric | A (BRICS) | B (DiffSBDD) |
|---|---|---|
| scaffold diversity | 11% | **96%** |
| median logP | 5.3 | **1.7** |
| novelty (median max Tanimoto to known) | 0.46 | **0.17** |
| QED median | ~0.6 | 0.57 |
| MW median | ~450 | 296 |
| in-domain & binder_prob≥0.5 | 346 | **0** |

### Honest interpretation (key scientific point)
B generates **diverse, novel, drug-like, pocket-conditioned** scaffolds — exactly
what A could not. BUT all 49 are **out of the ligand-ML scorer's applicability
domain** (Tanimoto ≤ 0.29), so the scorer cannot vouch for them. This is NOT
evidence they don't bind — it means the ligand-based scorer is the wrong tool for
novel scaffolds. The exploitation (A) vs exploration (B) tension is real:
- A: high predicted-binder rate, low novelty (rearranged known binders).
- B: high novelty/diversity, unverifiable by ligand-ML.

**Correct validation for B = structure-based** (docking pose/score vs known
binders; ultimately experiment). Docking of B follows in `dock_generated.py`.
B remains HYPOTHESES TO TEST.
