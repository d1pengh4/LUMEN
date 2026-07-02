# G1 result: brute-force scale-up plateaus — motivates RL (G2)

Date: 2026-07-02 · Roadmap G1
Question (user): does generating more molecules yield better candidates?

## Experiment
Generated a second DiffSBDD batch (~271 more molecules; ~540 total with the first
272). Same strengthened developability gate (PASS: diversity 91%, novelty 0.17,
24 developable, PAINS 3%). GNINA-scored the 24 new developable candidates.

## Result — NO improvement from brute-force scaling
| | CNNaffinity best | ≥7.0 (STOP target) | ≥6.3 |
|---|---|---|---|
| Batch 1 (272) | 6.56 | 0 | 3/33 |
| Batch 2 (+271) | **6.35** | 0 | 1/24 |

Doubling the random sample did NOT beat 6.56; new best is 6.35. Random DiffSBDD
sampling plateaus around CNNaffinity ~6.5. This confirms (with data) the expected
diminishing returns of unguided scaling.

## Interpretation & decision
- **Brute-force "generate more" has hit its ceiling** for this pocket/model.
- Reaching the STOP target (≥7.0) requires **directed optimization (RL / G2)** —
  learn from high-scorers and steer generation, not just sample more randomly.
- This is itself a publishable methodology finding: unguided pocket-conditioned
  generation saturates; goal-directed optimization is necessary to push affinity.
- Best candidate overall remains batch-1's 6.56 (but it has 5 stereocenters / hard
  synthesis; the practical lead is the SA-2.28 candidate at 6.34 — see synth log).

## Next
Proceed to G2 (RL / goal-directed optimization). GPU session stopped (frugal).
