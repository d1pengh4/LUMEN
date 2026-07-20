# LUMEN

**An honest, reproducible de novo design pipeline for non-immunosuppressive FKBP12
ligands, targeting the BMP/BMPR2 signaling axis in pulmonary arterial hypertension
(PAH).**

Full write-up: [`docs/paper/manuscript_EN.md`](docs/paper/manuscript_EN.md) ([한국어](docs/paper/manuscript_KO.md)) ·
Supplementary data: [`docs/paper/manuscript_SI_EN.md`](docs/paper/manuscript_SI_EN.md) ·
Master numeric record: [`docs/PROJECT_LOG.md`](docs/PROJECT_LOG.md) ·
Blog series (한국어): [velog.io/@eurochoic/series/LUMEN](https://velog.io/@eurochoic/series/LUMEN)

> Candidate molecules reported here are **computational hypotheses for experimental
> testing** — not validated binders or drugs. Every number below traces to a dated
> research log in [`docs/research-log/`](docs/research-log/).

---

## The premise

Loss of BMPR2 signaling drives pulmonary arterial hypertension. Because most disease
BMPR2 variants are loss-of-function (haploinsufficiency), the receptor itself is a poor
direct drug target. Tacrolimus (FK506) activates the same axis through a different node —
FKBP12 — but its immunosuppression (via calcineurin) makes it unusable as a PAH therapy.
LUMEN designs FKBP12-pocket ligands that keep the BMP-activating effect while avoiding the
calcineurin-recruiting surface. The original proposal targeted BMPR2 directly; that
premise was rejected after literature review (see
[ADR-001](docs/decisions/ADR-001-target-redefinition.md)) and the project pivoted to
FKBP12.

## What actually works, and what doesn't

The core result of this project is not just candidate molecules — it's a validated map of
which computational tools can be trusted on this target class.

| Method | Task | Result | Verdict |
|---|---|---|---|
| Docking (smina) | reproduce pose | redock RMSD **0.81 Å** | trustworthy (geometry) |
| Docking (smina) | rank affinity | enrichment ROC-AUC **0.455** | **fails** — below random |
| Ligand-based ML (RF/ECFP4) | rank potency | scaffold-split ρ **0.748** | usable, domain-gated |
| Ligand-based ML | FKBP12-vs-paralog selectivity | ρ ≈ **−0.01** | **fails** — not learnable |
| GNINA CNNaffinity | rank affinity | enrichment ROC-AUC **0.799** | **adopted** |

Docking fails to rank affinity on this shallow, near-apolar pocket; a ligand-based
classifier detects the FKBP *family* but can't discriminate FKBP12 from its 83%-identical
paralog FKBP12.6. GNINA's CNN rescoring is the one tool that works, and selectivity is
handled structurally (avoid the calcineurin "40s loop"), not by a learned score.

## Result: three leads

Pocket-conditioned generation (DiffSBDD) + goal-directed genetic-algorithm optimization
produced three developable, novel leads (Tanimoto < 0.4 to any known FKBP12 ligand):

| Lead | GNINA CNNaffinity | MD stability (1 ns, max RMSD) | Selectivity handles | Synthesis |
|---|---|---|---|---|
| L1 | **7.21** | 1.70 Å — stable | Trp59 | custom core (5 steps, 5/6 catalog) |
| L2 | 7.04 | 2.91 Å — stable | Trp59 + Ile90 | custom core (6 steps, 4/5 catalog) |
| L3 | 6.95 | 1.60 Å — stable | Trp59 + Ile90 | custom core (4 steps, 4/5 catalog) |

All three stay bound over 1 ns of explicit-solvent MD (separating real binding from a
lucky docking pose), and structurally contact the only two pocket residues that
distinguish FKBP12 from the cardiac-safety-relevant paralog FKBP12.6. None has a fully
catalog-purchasable synthesis route — the novel cores need custom synthesis, exposing a
robust **affinity↔synthesizability trade-off** (Fig. 1 in the paper) that directed
optimization narrowed but did not eliminate.

Two documented failure modes are part of the result, not swept under the rug: a genetic
algorithm run that turned out to be reward-hacking the training set (seed-dominance
artifact), and a second optimization round that improved every proxy metric but was
rejected by the validated GNINA check and a pre-registered safety-stop rule. Both are
in the paper (§2.3) as the negative results the field usually doesn't publish.

## Status

| Gate | Result |
|---|---|
| G1 — scale-up generation | brute-force plateaus at CNNaffinity ≈ 6.5 |
| G2 — goal-directed optimization | **PASS** — 6.56 → 7.21, beating the plateau |
| G3 — MD pose-stability | **PASS** — all 3 leads stay bound |
| G4 — retrosynthesis | **CONDITIONAL** — periphery routine, novel cores need custom synthesis |
| G7 (precursor) — selectivity feasibility | **DONE** — narrow but structurally addressable |
| G5 (AF3), G7 (FEP), G8 (paper) | supporting / resource-gated / this repository |

Full stage-gate log: [`docs/ROADMAP.md`](docs/ROADMAP.md). What's still missing is
experimental: no compound in this pipeline has been synthesized or assayed. The
[path-to-validation plan](docs/PATH_TO_VALIDATION.md) lays out the specific wet-lab
cascade (binding → calcineurin/NFAT counter-screen → FKBP12.6 counter-screen) that would
turn these into confirmed leads.

## Repository layout

```
docs/paper/              The manuscript (EN/KO), supplementary data, figures
docs/PROJECT_LOG.md       Master chronological record with exact numbers
docs/ROADMAP.md           Stage-gate plan (G1–G8) and current status
docs/decisions/           ADRs — why the target and major method choices were made
docs/research-log/        Dated, granular logs — one per work unit, including failures
docs/target-dossier-FKBP12.md   Target/pocket/selectivity reference
src/structure/            Structure fetch, pocket & interface analysis
src/data/                 ChEMBL ligand collection
src/scoring/              Docking, ML scorers, GNINA rescoring, selectivity/developability
src/generation/           BRICS baseline, DiffSBDD (GPU), GA optimization
src/md/                   Molecular-dynamics pose-stability pipeline
```

`results/` and `data/` are gitignored (large/regeneratable); provenance for every number
lives in the docs above. The GPU steps (DiffSBDD, GNINA, OpenMM MD) run as self-contained
scripts on transient cloud GPUs — see the Supplementary Methods in the paper for the exact
setup and the engineering issues encountered along the way.

## License / status

Private research project by [@d1pengh4](https://github.com/d1pengh4). Not affiliated with
any institution. Contact: eurochoic@gmail.com.
