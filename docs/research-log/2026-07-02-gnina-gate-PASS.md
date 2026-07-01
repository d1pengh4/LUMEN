# GNINA Rescoring GATE — PASS (Track A step 1)

Date: 2026-07-02 · Related: docking-benchmark (docking enrichment 0.455), ml-scorer
Result file: `results/gnina/gate_bg.log` · Script: `src/scoring/colab_gnina_gate.py`

## Purpose
Find a scoring method that actually ranks FKBP12 binders (docking's Vina score
failed: enrichment ROC-AUC 0.455, no discrimination). Test GNINA (CNN rescoring)
on the same 60-molecule enrichment set (20 actives / 40 decoys). GATE: must beat
0.455; target ≥ 0.65.

## Method (robust)
- RDKit ETKDG 3D per molecule (no obabel hang), one multi-mol SDF.
- Single gnina call, `--cnn_scoring rescore --exhaustiveness 4`, on a Colab T4.
- Run detached (nohup) on a persistent session -> survives client disconnects
  (earlier synchronous runs died on websocket drops / obabel hangs).

## Result — PASS
| Metric | GNINA | Docking (Vina) |
|---|---|---|
| **ROC-AUC** | **0.799** | 0.455 |
| EF@10% | 1.50 | ~1.0 |
| actives CNNaffinity | 6.34 ± 0.36 | — |
| decoys CNNaffinity | 5.80 ± 0.53 | — |

60/60 scored. GNINA clearly separates actives from decoys where Vina docking
could not. **GATE PASSED (0.799 ≥ 0.65).** GNINA CNNaffinity is adopted as the
affinity-ranking signal for candidates.

## Honest caveats
- 0.799 is good but not perfect; CNNaffinity is still a computational proxy, not
  measured binding. Decoys are presumed non-binders.
- This replaces docking for RANKING; it does not prove any candidate binds.
  Physical validation (MD stability, ultimately experiment) still required.

## Engineering lessons (documented so we don't repeat)
- gnina binary is 1.28 GB -> "setup hangs" were just the download, not stalls.
- obabel --gen3d hangs on some SMILES -> use RDKit ETKDG instead.
- long synchronous colab runs drop on websocket -> run detached + poll + download.
- reused-output-filename caused log corruption -> unique filenames per run.

## Next (gate passed -> advance)
Rescore the 33 developable generated candidates + 4 TIER-1 with GNINA; rank.
