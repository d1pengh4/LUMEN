# Supplementary Information

*Accompanies: "What de novo design can and cannot deliver on a hard target: an honest,
reproducible pipeline for non-immunosuppressive FKBP12 ligands against pulmonary arterial
hypertension." All tables below are transcribed directly from the result files named in
each caption, which live in `results/` in the project repository (gitignored; provenance
and generating scripts are in `docs/PROJECT_LOG.md` and `docs/research-log/`).*

---

## Table S1. Docking-benchmark validation set (redock control + drug-like anchors)

Source: `results/phase1/docking/docking_scores.csv`; §2.2, §4.4. FK506 is the native
redock control; the six anchors are the most potent drug-like (MW ≤ 500) FKBP12 ligands in
the curated ChEMBL set.

| Ligand | pChEMBL | Dock score (kcal/mol) |
|---|---|---|
| FK506 (native, redock) | high | **−10.30** |
| CHEMBL4090599 | 9.35 | −7.52 |
| CHEMBL4072643 | 8.15 | **−7.57** (best of the six anchors) |
| CHEMBL4075704 | 8.15 | −7.55 |
| CHEMBL4101268 | 8.15 | −6.93 |
| CHEMBL33094 | 9.00 | −6.85 |
| CHEMBL4102121 | 8.15 | −6.73 |

Note the ranking inversion discussed in §2.2: CHEMBL4090599 (pChEMBL 9.35, the most potent
anchor) does not score best; CHEMBL4072643 (pChEMBL 8.15) does. This six-molecule set
motivated the formal 60-molecule enrichment test (Table 1).

---

## Table S2. Full 37-candidate GNINA-ranked pool (Route B, batch 1)

Source: `results/phase2b/FINAL_candidates.csv`; §2.3, §2.5, Fig. 1, Table 4. The 33
developable + 4 TIER-1 DiffSBDD candidates from generation batch 1, ranked by GNINA
CNNaffinity (exhaustiveness 4 screening pass). `ad_tanimoto` is populated only for the
`developable` source rows (the applicability-domain gate used for the ML potency model,
§2.2); `dock_score` is populated only for `tier1` rows (independently flagged by docking +
developability ahead of GNINA scoring, §2.3). Blank cells reflect which upstream filter a
molecule passed through, not missing data.

| CNNaffinity | source | SMILES | binder_like | QED | SA | logP | MW | AD Tanimoto | dock (kcal/mol) |
|---|---|---|---|---|---|---|---|---|---|
| 6.560 | tier1 | `CCC[C@@H]1CC(=O)N[C@@]2(C)[C@H](C1)[C@@H](Cc1ccc(F)cc1N)CC[C@@H]2O` | True | 0.719 | 4.39 | 3.42 | 362.2 | — | −8.54 |
| 6.336 | developable | `O=C(O)CCCCN(CCO)c1cccc(F)c1F` | True | 0.711 | 2.28 | 2.02 | 273.1 | 0.224 | — |
| 6.336 | tier1 | `CC1([C@H](Cc2cccc(Cl)c2)OCCNC(=O)O)CCC1` | True | 0.754 | 3.10 | 3.73 | 311.1 | — | −7.58 |
| 6.273 | developable | `CCSC1=NCC[C@@H](CN(CF)C(=O)NCCc2ccccc2)N1` | False | 0.742 | 3.45 | 2.64 | 352.2 | 0.338 | — |
| 6.164 | tier1 | `OCCCC[C@H]1CCCC[C@H](CC2c3ccccc32)N1` | False | 0.776 | 3.43 | 3.59 | 273.2 | — | −7.57 |
| 6.155 | developable | `C=C(C)NNS(=O)(=O)c1ccccc1Cc1ncccn1` | False | 0.788 | 2.63 | 1.38 | 304.1 | 0.197 | — |
| 6.067 | developable | `O=C1[C@@H](Cc2ccccc2)C[C@H]2C[C@H]2N1CC1CCOCC1` | False | 0.855 | 3.46 | 2.89 | 299.2 | 0.291 | — |
| 5.941 | developable | `Cc1ccc(Cn2ccnn2)cc1NC=C1CC1` | False | 0.893 | 2.75 | 2.72 | 240.1 | 0.213 | — |
| 5.820 | developable | `CCCCCN[C@@H]1[C@@H](c2ccccc2)CCCC[C@@H]1O` | False | 0.610 | 3.07 | 3.85 | 275.2 | 0.333 | — |
| 5.803 | developable | `NCCCCS[C@H]1[C@H](N2CCCCC2)CCCS1(=O)=O` | False | 0.757 | 3.78 | 1.85 | 320.2 | 0.169 | — |
| 5.761 | developable | `Cc1ccc(CN2CCC[C@@H](CCCC(=O)O)C2)cc1` | False | 0.864 | 2.27 | 3.46 | 275.2 | 0.328 | — |
| 5.647 | tier1 | `C[C@H]1CCCN(c2ccc(F)cc2)[C@@H]1CCCCC(=O)O` | False | 0.804 | 2.92 | 4.08 | 293.2 | — | −7.34 |
| 5.625 | developable | `C[C@]1(c2cccc(O)c2)CCOC[C@H](CC#N)N1` | False | 0.836 | 3.92 | 1.90 | 246.1 | 0.191 | — |
| 5.620 | developable | `CNC(=O)[C@H]1CCO[C@H](C[C@H]2CCCO2)[C@@H]1c1ccc(C#N)cc1` | False | 0.921 | 3.86 | 2.36 | 328.2 | 0.222 | — |
| 5.583 | developable | `C[C@H]1CC[C@@H]([C@@]2(O)CCC[C@H](N)C2)CC1` | False | 0.700 | 3.47 | 2.45 | 211.2 | 0.159 | — |
| 5.580 | developable | `CCCOC[C@@H]1CC[C@H](O)C[C@@H]1NCC(=O)C1CCCC1` | False | 0.675 | 3.71 | 2.29 | 297.2 | 0.222 | — |
| 5.510 | developable | `Nc1cc(F)ccc1[C@@H]1CCCC[C@@H](C=CC(=O)O)CO1` | False | 0.662 | 3.47 | 3.30 | 293.1 | 0.208 | — |
| 5.476 | developable | `C[C@H]1CCC[C@@](Oc2ccccc2)(C(=O)NCC(=O)O)[C@@H]1O` | False | 0.759 | 3.55 | 1.19 | 307.1 | 0.286 | — |
| 5.464 | developable | `O=C(O)CCCN1Cc2cc(F)ccc2[C@H](F)C1` | False | 0.898 | 2.90 | 2.52 | 255.1 | 0.228 | — |
| 5.183 | developable | `C1CCC(C[C@H]2CCCS2)CC1` | False | 0.628 | 3.20 | 3.85 | 184.1 | 0.173 | — |
| 5.166 | developable | `NC(=O)CCO[C@@H]1CC[C@H]2C[C@@H]2CC(=O)C1` | False | 0.776 | 3.59 | 1.03 | 225.1 | 0.206 | — |
| 5.137 | developable | `NC(=O)[C@@H]1CCCCCCc2c(OOF)cccc2O1` | False | 0.683 | 3.69 | 2.62 | 283.1 | 0.219 | — |
| 5.116 | developable | `CC[C@@H]1C[C@@]1(C)c1ccc(F)cc1N` | False | 0.717 | 3.43 | 3.10 | 193.1 | 0.176 | — |
| 5.060 | developable | `Nc1ccc(Cl)c(CC2NOCS2)c1` | False | 0.764 | 3.97 | 2.02 | 230.0 | 0.150 | — |
| 5.043 | developable | `CCCC[C@H]1CC=C[C@H]2CC[C@@H](C(=O)O)C[C@@H]12` | False | 0.752 | 3.86 | 3.87 | 236.2 | 0.212 | — |
| 4.878 | developable | `NCCc1cccc2c1C(=O)SC2=O` | False | 0.792 | 2.82 | 1.22 | 207.0 | 0.137 | — |
| 4.753 | developable | `COC(=O)N1CCC[C@@H]1C1=CCC=CC1=O` | False | 0.676 | 3.42 | 1.67 | 221.1 | 0.274 | — |
| 4.657 | developable | `CCCCCN1OC(=O)CCC12CC2` | False | 0.647 | 3.67 | 2.26 | 197.1 | 0.172 | — |
| 4.597 | developable | `CC[C@@H](CCC(=O)O)N(C(=O)O)[C@@H]1C(=O)[C@H]1C` | False | 0.730 | 3.77 | 1.20 | 243.1 | 0.188 | — |
| 4.459 | developable | `C[C@@H]1SCC[C@@H]1C(=O)O` | False | 0.604 | 3.88 | 1.21 | 146.0 | 0.170 | — |
| 4.307 | developable | `CCC(=O)CC[C@H](C)C(=O)O` | False | 0.659 | 2.73 | 1.47 | 158.1 | 0.191 | — |
| 4.289 | developable | `C[C@@H]1C(=O)COC(F)(F)[C@H]1c1ccccc1` | False | 0.735 | 3.56 | 2.60 | 226.1 | 0.226 | — |
| 4.202 | developable | `FC(F)(F)O[C@@H](c1ccccc1)C1CC1` | False | 0.748 | 2.79 | 3.67 | 216.1 | 0.230 | — |
| 4.027 | developable | `O=C(O)CNc1ccccc1` | False | 0.682 | 1.45 | 1.18 | 151.1 | 0.277 | — |
| 3.690 | developable | `C=C(C)N(CCC)C(=O)O` | False | 0.655 | 2.99 | 1.91 | 143.1 | 0.155 | — |
| 3.434 | developable | `O=C(O)Cc1ccccc1` | False | 0.665 | 1.21 | 1.31 | 136.1 | 0.343 | — |
| 3.094 | developable | `CN1COc2ccc(O)cc21` | False | 0.604 | 2.56 | 1.18 | 151.1 | 0.216 | — |

Only the top 3 rows clear the active-like threshold (≥ 6.3); the remaining 34 score at or
below decoy level (5.80 ± 0.53, Table 1) — the honest majority outcome of pocket-conditioned
generation referenced in §2.3.

---

## Table S3. Molecular-dynamics raw trajectory summary (leads + the failed synthesizable candidate)

Source: `results/md/g3_md_rmsd.csv`; §2.4, §2.5, Fig. 2. All runs: 1 ns, OpenCL platform
(§4.6, §S3 of this document below), 100 saved frames.

| Lead | CNNaffinity | Heavy atoms (ligand) | RMSD mean (Å) | RMSD max (Å) | RMSD final (Å) | Verdict |
|---|---|---|---|---|---|---|
| L1 | 7.21 | 32 | 1.15 | 1.70 | 1.67 | PASS |
| L2 | 7.04 | 28 | 1.70 | 2.91 | 1.52 | PASS |
| L3 | 6.95 | 32 | 1.05 | 1.60 | 0.59 | PASS |
| S1 (synthesizable backup, §2.5) | 6.34 | 19 | 2.93 | 6.04 | 5.45 | **FAIL** |

L3's final-frame RMSD (0.59 Å) is notably lower than its own mean (1.05 Å), i.e. it
re-approaches the initial pose by the end of the 1 ns window rather than drifting away —
consistent with a genuinely stable, not merely slowly-drifting-but-not-yet-escaped, pose.
S1's final RMSD (5.45 Å) exceeds its own mean (2.93 Å), i.e. it is still moving away from
the pocket when the simulation ends — consistent with ongoing unbinding rather than a
transient excursion.

---

## Table S4. Selectivity-handle contact geometry (raw re-docking output)

Source: `results/selectivity/leads_handle_contacts.csv`; §2.6, Table 3, Fig. 3. Re-docking
used smina at exhaustiveness 8 (§4.8); contact threshold ≤ 4.5 Å.

| Lead | Dock score (kcal/mol) | Trp59 min. distance (Å) | Ile90 min. distance (Å) | Contacts Trp59 | Contacts Ile90 |
|---|---|---|---|---|---|
| L1 | −7.70 | 3.44 | 5.16 | yes | no |
| L2 | −6.90 | 3.87 | 3.83 | yes | yes |
| L3 | **−9.00** | 3.67 | 3.59 | yes | yes |

L3 has the most favorable re-docking score of the three leads (−9.00 kcal/mol) in addition
to the best MD stability (Table S3) and both selectivity handles — the combined basis for
prioritizing it first in the experimental roadmap (§3).

---

## Table S5. Structural/synthesizability triage of the three leads

Source: `results/synthesis/leads_synth_feasibility.csv`; §2.5, §4.7. `synthesizable_easy`
requires SA ≤ 3.5 among other criteria (§4.4); all three leads fail this strict triage even
before retrosynthesis is attempted, foreshadowing the AiZynthFinder result (Table 3).

| Lead | SA score | Stereocenters | Total rings | Aromatic rings | Max ring size | Strained ring | Spiro centers | "Easy" synthesizability |
|---|---|---|---|---|---|---|---|---|
| L1 | 4.34 | 3 | 4 | 2 | 7 | No | 0 | No |
| L2 | 4.05 | 3 | 3 | 0 | 6 | No | 0 | No |
| L3 | 4.39 | 3 | 4 | 2 | 7 | No | 0 | No |

---

## Table S6. Developability panel across the three generation pools

Source: `docs/research-log/2026-07-01-developability-and-verification.md`
(`results/developability/*.csv`); §2.3, §2.8. RDKit-only rule-based developability proxy
(QED, SA, Lipinski/Veber, PAINS — §4.4), applied before any docking or GNINA filtering, on
an early pilot scale distinct from the full gated batches in Table 4.

| Pool | n | developable (SA≤4.5, QED≥0.5, Lipinski+Veber, no PAINS) | SA (median) | QED (median) | PAINS alerts |
|---|---|---|---|---|---|
| Route A, top 50 (BRICS) | 50 | 3 | 4.64 | 0.60 | 0 |
| Route B, all docked (DiffSBDD pilot) | 49 | 16 | 4.75 | 0.57 | 0 |
| Route B, dock-promising subset | 8 | 4 | 4.88 | 0.72 | 0 |

Only half of even the dock-promising DiffSBDD subset is simultaneously easily
synthesizable — an early, smaller-scale confirmation of the affinity↔synthesizability
tension formalized later at full scale in §2.5 and Fig. 1.

---

## Table S7. Extended reproducibility notes (see also §S1–S4 in the main text's Supplementary Methods)

| Issue | Symptom | Resolution |
|---|---|---|
| gnina binary size | Setup appeared to hang | 1.28 GB download, not a stall; added progress logging |
| `obabel --gen3d` | Hangs on certain SMILES | Switched to RDKit ETKDG embedding |
| Long synchronous Colab execs | Dropped on WebSocket disconnect | Run detached (`nohup`), poll, download on completion |
| Reused output filenames | Silent log corruption across runs | Unique filenames per run/session |
| arm64 Colab CLI wrapper | Forced x86_64 Python.framework; `cryptography` wheel failed to `dlopen` | Switched to native uv-installed `colab` binary |
| Mid-run VM reclaim | ~40 min synchronous `colab exec` lost the whole VM (L2's first MD attempt) | Prefer harness-tracked background execution that keeps the client connected |

---

*End of Supplementary Information. All source CSVs/logs referenced above are in the
project repository's `results/` and `docs/research-log/` directories.*
