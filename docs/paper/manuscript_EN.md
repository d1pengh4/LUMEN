# What de novo design can and cannot deliver on a hard target: an honest, reproducible pipeline for non-immunosuppressive FKBP12 ligands against pulmonary arterial hypertension

**Authors:** [Author One]¹, [Author Two]², [Corresponding Author]¹*
**Affiliations:** ¹[Institution]; ²[Institution]
**\*Correspondence:** [email]

*Preprint draft — LUMEN project. All quantitative results trace to the project record
(`docs/PROJECT_LOG.md`) and dated research logs; code and result tables are in the
project repository. Candidate molecules reported here are computational hypotheses for
experimental testing, not validated binders or drugs.*

---

## Abstract

**Background.** Loss of bone morphogenetic protein receptor type-2 (BMPR2) signaling is
central to pulmonary arterial hypertension (PAH). Because most disease-associated *BMPR2*
variants are loss-of-function (haploinsufficiency), the receptor is a poor direct
small-molecule target, and the immunophilin FKBP12 — through which tacrolimus (FK506)
activates BMP signaling — is a mechanistically attractive alternative, provided the
calcineurin-mediated immunosuppression of FK506 is avoided.

**Results.** We built and, critically, *validated the limits of*, a de novo design
pipeline for FKBP12 pocket ligands that spare the calcineurin-recruiting "40s loop." On
this shallow hydrophobic pocket we find that (i) molecular docking reproduces the
crystallographic pose (redock RMSD 0.81 Å) but has **no affinity-ranking power**
(enrichment ROC-AUC 0.455, below random); (ii) a ligand-based machine-learning model
ranks potency (scaffold-split Spearman ρ = 0.75) but **cannot** predict FKBP12-versus-
paralog selectivity (ρ ≈ 0); and (iii) a convolutional-neural-network rescorer (GNINA
CNNaffinity) *does* rank affinity (enrichment ROC-AUC 0.799) and was adopted. Pocket-
conditioned generative design (DiffSBDD) followed by goal-directed genetic-algorithm
optimization yielded three developable, novel leads (GNINA CNNaffinity 6.95–7.21) that
**all remain bound over 1 ns of explicit-solvent molecular dynamics** (ligand RMSD ≤ 2.9
Å). Retrosynthetic analysis (AiZynthFinder) shows their novel cores require custom
synthesis, exposing a robust **affinity↔synthesizability trade-off** that directed
optimization did not overcome; and a structural analysis shows the FKBP12 and FKBP12.6
pockets differ at only **two conservative residues (Trp59, Ile90)**, which the leads
nonetheless contact — making selectivity structurally addressable but narrow.

**Conclusions.** We present a reproducible pipeline, a small honestly-validated shortlist
of testable FKBP12 lead hypotheses with a concrete selectivity rationale, and a candid map
of where standard computational methods succeed and fail on this target class. The honest
negatives — not the candidates alone — are the contribution.

**Keywords:** pulmonary arterial hypertension; FKBP12; BMPR2; de novo drug design;
molecular dynamics; retrosynthesis; selectivity; GNINA; honest benchmarking.

---

## 1. Introduction

Pulmonary arterial hypertension (PAH) is a progressive vasculopathy in which reduced
BMPR2 signaling drives pulmonary vascular remodeling [1,2]. Approximately two-thirds of
heritable *BMPR2* variants are truncating and undergo nonsense-mediated decay, producing
haploinsufficiency rather than a druggable altered protein [3]; direct pocket-binding of
the mutant receptor is therefore ill-founded. The 2024 approval of sotatercept, an
activin-signaling ligand trap, established that rebalancing the BMP/TGF-β axis is
disease-modifying [4], but leaves room for orthogonal, orally available mechanisms.

Tacrolimus (FK506) activates BMP signaling and reverses experimental pulmonary
hypertension by binding FKBP12 and relieving its inhibition of BMP type-I receptors [5].
Its therapeutic use is precluded by immunosuppression: the FK506–FKBP12 complex recruits
calcineurin through an effector surface adjacent to the FKBP12 "40s loop." Non-
immunosuppressive FKBP ligands (e.g., the neuroimmunophilin class) demonstrate that pocket
occupancy can be decoupled from calcineurin recruitment [6]. This motivates a specific
design objective: **fill the FKBP12 FK506 pocket while avoiding the calcineurin-recruiting
40s loop (residues 32–45)**, and discriminate the 83%-identical paralog FKBP12.6 (FKBP1B)
to avoid cardiac ryanodine-receptor liabilities.

De novo generative design is often reported through its successes. Here we instead report
a fully documented pipeline **and its validated failure modes** on FKBP12 — a shallow,
hydrophobic, low-druggability pocket with a near-identical paralog — because knowing which
computational instruments are trustworthy on such a target is itself a transferable
result. Our contributions are (i) a reproducible generation-to-validation pipeline; (ii)
an explicit capability map of docking, ligand-ML, CNN rescoring, MD, retrosynthesis, and
structural selectivity analysis on this target; and (iii) three MD-stable, novel lead
hypotheses with a concrete, falsifiable selectivity rationale.

## 2. Results

### 2.1 Target, pocket, and the structural selectivity rule

After a premise correction (from BMPR2 to FKBP12; §4.1), we defined the pocket on the
FKBP12–FK506 complex (PDB 1FKF). With FK506 removed, fpocket identified the FK506 site as
the top pocket (druggability 0.517, volume 277 Å³, hydrophobicity 61.2, polar SASA ≈ 0;
pocket centroid 4.35 Å from native FK506) — a shallow, hydrophobic, moderately druggable
site. Interface mapping (BioPython, 5 Å) against an on-target proxy (FKBP12–ALK5, 1B6C; 22
residues, overlapping the pocket) and the off-target calcineurin complex (1TCO; 23
residues) showed the calcineurin interface uniquely uses the 40s-loop patch (Thr27, Asp32,
Lys34/35, Ser38, Arg40, Asn43, Lys44, Pro45), absent from the receptor interface. This
yields the design rule used throughout: **occupy the pocket, avoid the 40s loop.**

From ChEMBL (CHEMBL1902) we assembled 725 unique FKBP12 binders (821 activities, pChEMBL
3.1–11.0, median 6.71); the space skews macrocyclic (median MW 556), leaving 68 drug-like
anchors (MW ≤ 500, QED ≥ 0.5, pChEMBL ≥ 6) as fragment donors.

### 2.2 A capability map: which scoring methods can be trusted (Table 1)

We benchmarked each scoring instrument before using it.

**Table 1. Validated capabilities and limits of the scoring stack.**

| Method | Task | Metric | Verdict |
|--------|------|--------|---------|
| smina docking | reproduce pose | FK506 redock RMSD **0.81 Å** | trustworthy (geometry) |
| smina docking | rank affinity | enrichment ROC-AUC **0.455** | **fails** (poses only) |
| RF/ECFP4 classifier | binder triage | scaffold-split ROC-AUC **0.993** | trustworthy (vs random) |
| RF/ECFP4 regressor | rank potency | scaffold-split ρ **0.748** (R² 0.56) | usable, domain-gated |
| RF/ECFP4 | FKBP12 vs paralog selectivity | ρ ≈ **−0.01** | **fails** (not learnable) |
| GNINA CNNaffinity | rank affinity | enrichment ROC-AUC **0.799** | **adopted** |

Three conclusions are load-bearing. Docking score is uninformative for affinity here (ROC
0.455) and was restricted to pose generation. Ligand-based ML predicts potency (ρ = 0.75)
but is a family-level FKBP binder detector — it scores FKBP5 ligands ~0.79 and cannot
predict FKBP12/paralog selectivity (ρ ≈ 0). GNINA CNNaffinity separates actives (6.34 ±
0.36) from decoys (5.80 ± 0.53) with ROC-AUC 0.799 and was adopted as the affinity signal,
gated by an applicability-domain Tanimoto threshold ≥ 0.30. Selectivity was therefore
addressed structurally and (prospectively) experimentally, never by a learned reward.

### 2.3 Generation and goal-directed optimization

A BRICS baseline produced drug-like but low-diversity molecules (scaffold diversity 11%).
Pocket-conditioned diffusion (DiffSBDD) produced diverse, novel scaffolds (scaffold
diversity 89%, novelty 0.18, PAINS 0%) but a best GNINA CNNaffinity of only 6.56, and
doubling the sample did not improve it (brute-force plateau ≈ 6.5). Goal-directed genetic-
algorithm optimization of the validated potency model, under hard developability and
applicability-domain constraints, broke the plateau (GNINA 6.56 → 7.06; high-quality
exhaustiveness-8 rescoring 7.06 → 7.21). Three developable, novel leads emerged (Table 2).

**Table 2. The three lead candidates (structures and properties).**

| Lead | SMILES | GNINA CNNaff. | QED | MW | SA | AD Tanimoto | novel |
|------|--------|--------------|-----|----|----|-------------|-------|
| L1 | `C=C[C@H]1CN(c2cccc(NCC)c2)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1ccc(C)cc1` | 7.21 | 0.67 | 454 | 4.34 | 0.46 | yes |
| L2 | `CC(C)(C)C(=O)C(=O)N1CCCCC1N1CCCC1N1CCC[C@H]1CC(N)=O` | 7.04 | 0.72 | 392 | 4.05 | 0.38 | yes |
| L3 | `C=C[C@H]1CN(c2cccc(NCC)c2)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1cccc(C)c1` | 6.95 | 0.67 | 454 | 4.39 | 0.45 | yes |

The synthesizable candidate discussed in §2.5 is `O=C(O)CCCCN(CCO)c1cccc(F)c1F` (denoted
S1; GNINA 6.34, SA 2.28).

### 2.4 Molecular-dynamics pose stability: all three leads stay bound (Fig. 2)

For each lead we built the FKBP12–ligand complex (Amber ff14SB + GAFF-2.11 + TIP3P, 0.15 M
NaCl, ≈ 20,500 atoms) and ran 1 ns of explicit-solvent MD (OpenMM), tracking ligand
heavy-atom RMSD from the docked pose after Cα alignment. All three leads remained bound
(mean/max RMSD: L1 1.15/1.70; L2 1.70/2.91; L3 1.05/1.60 Å) — well below the ~3–4 Å
unbinding range (Fig. 2) — separating genuine binding from lucky docking poses.

### 2.5 Synthesizability and the affinity↔synthesizability trade-off (Fig. 1)

Retrosynthetic analysis (AiZynthFinder; USPTO models, ZINC purchasable stock of 17.4 M
compounds) found **no complete purchasable route** for L1–L3, even under deep search: the
peripheral chemistry is routine (sulfonylation, amide coupling, reductive amination,
Suzuki; 4/5–5/6 precursors catalog), but the irreducible precursor is in every case the
**novel core scaffold** (a fused proline-derived diaza-bicyclic sulfonamide for L1/L3; an
oligo-pyrrolidine for L2). Screening all 37 GNINA-ranked candidates, only 7/37 were fully
solvable, and these cap at GNINA ≈ 6.34 (known-active level); the best synthesizable
candidate (`O=C(O)CCCCN(CCO)c1cccc(F)c1F`; GNINA 6.34, SA 2.28, 1-step) **failed MD** (RMSD
final 5.45 Å). No single molecule was simultaneously high-affinity, MD-stable, and
synthesizable (Fig. 1). A genetic algorithm constrained to synthesizable-easy chemistry
(SA ≤ 3.5) reached predicted potency 7.09, yet retrosynthesis solved only 1/10 of its top
candidates and the one solvable molecule had disagreeing potency/binder models — the
trade-off is robust, and SA score is a weak surrogate for a real catalog route.

### 2.6 Selectivity feasibility: narrow but structurally addressable (Fig. 3)

Aligning FKBP12 (UniProt P62942) and FKBP12.6 (P68106; both 108 aa) shows 83.3% identity
(18 differences). Only **two of ~13 pocket-lining residues differ, both conservatively:
Trp59→Phe and Ile90→Val** (FKBP12.6 is smaller and less aromatic at each). This explains
the ligand-ML selectivity failure (§2.2) — there is almost no chemical difference to learn
— and simultaneously localizes the selectivity handle. Re-docking the leads and measuring
per-residue contacts shows all three contact Trp59 (3.4–3.9 Å) and L2/L3 additionally
contact Ile90 (3.6–3.8 Å) (Fig. 3): selectivity is addressable for the leads, not merely
for the target. Contact does not equal a favorable ΔΔG — that requires relative
free-energy perturbation (FEP) or experiment — but the prerequisite chain (the pockets
differ; the leads touch the differing residues) holds end-to-end.

### 2.7 End-to-end validation summary (Table 3)

Table 3 consolidates the four validation axes per lead. All three are affinity-plausible,
MD-stable, and engage the Trp59 selectivity residue; L2 and L3 additionally engage Ile90.
None is catalog-synthesizable (each requires custom core synthesis). This is the honest
lead package the pipeline delivers.

**Table 3. End-to-end validation of the three leads.**

| Lead | Affinity (GNINA) | MD max RMSD (Å) | MD verdict | Trp59 (Å) | Ile90 (Å) | selectivity handles | AiZynth route |
|------|------------------|-----------------|------------|-----------|-----------|---------------------|---------------|
| L1 | 7.21 | 1.70 | stable | 3.44 | 5.16 | Trp59 only | not solved (5 steps, 5/6 catalog) |
| L2 | 7.04 | 2.91 | stable | 3.87 | 3.83 | Trp59 + Ile90 | not solved (6 steps, 4/5 catalog) |
| L3 | 6.95 | 1.60 | stable | 3.67 | 3.59 | Trp59 + Ile90 | not solved (4 steps, 4/5 catalog) |

## 3. Discussion

On a shallow hydrophobic pocket with a near-identical paralog, two workhorse methods fail
in instructive ways: docking-*score* ranking (ROC 0.455) and ligand-based selectivity
prediction (ρ ≈ 0). What remains effective is a layered strategy — CNN rescoring for
affinity, explicit-solvent MD for pose reality, and structural analysis for selectivity —
which together still yield a defensible, testable shortlist. We argue that publishing this
capability map, including the failures, is more useful to the field than another
success-only generative case study, because the failures are target-class properties, not
project-specific accidents.

Two trade-offs structure the medicinal-chemistry outlook. First, **affinity↔
synthesizability** (Fig. 1): the highest-affinity, MD-stable leads are novel-core scaffolds
that require custom synthesis, while the readily-synthesizable candidates are weaker and
less stable; directed optimization narrowed but did not eliminate this tension, and — a
methodological point — SA score is not a reliable proxy for retrosynthetic solvability
(retrosynthesis must be the in-loop objective). Second, **potency↔selectivity feasibility**:
the same near-identical pocket that defeats ligand-ML selectivity also implies a narrow
selectivity window, yet it is a *specific* one — the Trp59 indole and the Ile90 methyl,
present in FKBP12 but not FKBP12.6, which our leads already engage. Selectivity thus becomes
a falsifiable design objective rather than an assumption.

**Limitations.** GNINA CNNaffinity is a validated-but-imperfect proxy (ROC 0.80, not
measured binding); 1 ns single-replica MD is a stability screen, not a free-energy
calculation; the on-target interface is modeled by an ALK5 proxy; and neither relative FEP
for the FKBP12/FKBP12.6 window nor any wet-lab binding, immunosuppression, or efficacy
assay has been performed. Every candidate is a hypothesis for experimental testing. The
decisive next experiments are (i) FKBP12 binding, (ii) a calcineurin/NFAT counter-screen to
confirm non-immunosuppression, and (iii) an FKBP12.6 counter-screen to test the Trp59/Ile90
selectivity hypothesis.

## 4. Methods

**4.1 Target definition.** The original BMPR2-pocket premise was rejected after a
literature check (BMPR2 haploinsufficiency; sotatercept approval) and the target redefined
to FKBP12 (recorded as ADR-001). Structures 1FKF, 1B6C, 1TCO, 2PPN, 1F40 were retrieved
from the RCSB and identities verified from headers.

**4.2 Pocket and interfaces.** fpocket [7] on 1FKF (FK506 removed) quantified the pocket;
interfaces were mapped with BioPython at a 5 Å cutoff.

**4.3 Ligand data.** FKBP12 activities were collected from ChEMBL (CHEMBL1902) [8];
drug-like anchors were filtered at MW ≤ 500, QED ≥ 0.5, pChEMBL ≥ 6.

**4.4 Scoring.** Docking used smina [9] (AutoDock Vina scoring [10]); redock RMSD and a
20-active/40-decoy enrichment set assessed pose and ranking. A RandomForest on 2048-bit
ECFP4 provided binder triage, a pChEMBL regressor, and a selectivity probe, evaluated by
scaffold-split cross-validation. GNINA CNNaffinity [11] rescoring was validated on a
60-molecule enrichment set and adopted with a Tanimoto ≥ 0.30 applicability-domain gate.

**4.5 Generation and optimization.** BRICS recombination provided a CPU baseline; DiffSBDD
[12] performed pocket-conditioned diffusion generation on a GPU; a BRICS-fragment genetic
algorithm optimized the validated potency model under hard developability and
applicability-domain constraints. GNINA rescoring at exhaustiveness 8 ranked the winners.

**4.6 Molecular dynamics.** Complexes were prepared with PDBFixer; ligands were docked into
the pocket (autobox on FK506). Systems (Amber ff14SB [13] + GAFF-2.11 [14] + TIP3P [15],
0.15 M NaCl) were built with OpenMM/openmmforcefields [16], energy-minimized, and run for
1 ns (LangevinMiddle integrator, 4 fs with hydrogen-mass repartitioning, Monte Carlo
barostat) using OpenCL (a CUDA/driver PTX incompatibility on the cloud GPU forced the
OpenCL platform). Ligand heavy-atom RMSD versus the initial pose was computed with MDTraj
after protein-Cα alignment.

**4.7 Retrosynthesis.** AiZynthFinder 4.4.1 [17] with public USPTO expansion/ring-breaker/
filter models and the ZINC purchasable stock (17.4 M) was run per molecule (default and a
deep 500-iteration search); a molecule was "solved" only if a complete route to catalog
building blocks was found.

**4.8 Selectivity analysis.** UniProt sequences P62942 and P68106 were aligned positionally
(no indels); differences were intersected with the pocket residue set. Leads were re-docked
into FKBP12 (smina, exhaustiveness 8) and per-residue minimum lead–protein distances were
computed; contact was defined at ≤ 4.5 Å.

**4.9 Reproducibility.** All scripts, per-step research logs, the master numeric record,
and result tables are in the project repository; large model/stock files are excluded from
version control with provenance documented. GPU steps run as self-contained scripts on
transient cloud GPUs.

## Data and code availability
All code (`src/`), gate logs (`docs/research-log/`), the master record
(`docs/PROJECT_LOG.md`), result tables (`results/`), and figures (`docs/paper/figures/`)
are available in the project repository. Structures are from the RCSB PDB; ligand data from
ChEMBL; retrosynthesis models/stock from the AiZynthFinder public data release.

## Author contributions
[To be completed.] Conceptualization, methodology, software, analysis, and writing.

## Competing interests
The authors declare no competing interests.

## Acknowledgements
[To be completed.]

## References
*(References 2, 4, 5, 11, 12, 17, 19 verified against the primary sources with DOIs below;
the remaining well-established method/database citations should have their page/DOI details
confirmed at copy-editing.)*

1. Morrell NW, et al. Cellular and molecular basis of pulmonary arterial hypertension.
   *J Am Coll Cardiol* / *Nat Rev Cardiol* (reviews).
2. Lane KB, Machado RD, Pauciulo MW, et al. Heterozygous germline mutations in BMPR2,
   encoding a TGF-β receptor, cause familial primary pulmonary hypertension. *Nat Genet*
   2000;26:81–84. doi:10.1038/79226.
3. Machado RD, et al. Mutations of the TGF-β type II receptor BMPR2 in pulmonary arterial
   hypertension. *Hum Mutat* 2006 (NMD/haploinsufficiency).
4. Hoeper MM, Badesch DB, Ghofrani HA, et al. Phase 3 trial of sotatercept for treatment
   of pulmonary arterial hypertension (STELLAR). *N Engl J Med* 2023;388(16):1478–1490.
   doi:10.1056/NEJMoa2213558.
5. Spiekerkoetter E, et al. FK506 activates BMPR2, rescues endothelial dysfunction, and
   reverses pulmonary hypertension. *J Clin Invest* 2013;123(8):3600–3613. doi:10.1172/JCI65592.
6. Gold BG, Steiner JP, et al. Neuroimmunophilin (non-immunosuppressive FKBP) ligands.
   *(neuroimmunophilin literature; e.g., GPI-1046).*
7. Le Guilloux V, Schmidtke P, Tuffery P. Fpocket: an open source platform for ligand
   pocket detection. *BMC Bioinformatics* 2009;10:168.
8. Mendez D, et al. ChEMBL: towards direct deposition of bioassay data. *Nucleic Acids Res*
   2019;47:D930–D940.
9. Koes DR, Baumgartner MP, Camacho CJ. Lessons learned in empirical scoring with smina.
   *J Chem Inf Model* 2013;53:1893–1904.
10. Trott O, Olson AJ. AutoDock Vina. *J Comput Chem* 2010;31:455–461.
11. McNutt AT, Francoeur P, Aggarwal R, et al. GNINA 1.0: molecular docking with deep
    learning. *J Cheminform* 2021;13:43. doi:10.1186/s13321-021-00522-2.
12. Schneuing A, Harris C, Du Y, et al. Structure-based drug design with equivariant
    diffusion models (DiffSBDD). *Nat Comput Sci* 2024;4(12):899–909.
    doi:10.1038/s43588-024-00737-x (preprint arXiv:2210.13695).
13. Maier JA, et al. ff14SB: improving the accuracy of protein side chain and backbone
    parameters. *J Chem Theory Comput* 2015;11:3696–3713.
14. Wang J, et al. Development and testing of a general Amber force field (GAFF).
    *J Comput Chem* 2004;25:1157–1174.
15. Jorgensen WL, et al. Comparison of simple potential functions for simulating liquid
    water (TIP3P). *J Chem Phys* 1983;79:926–935.
16. Eastman P, et al. OpenMM 8: molecular dynamics simulation with machine-learning
    potentials. *J Phys Chem B* 2024 (and OpenMM 7, *PLoS Comput Biol* 2017).
17. Genheden S, Thakkar A, Chadimová V, et al. AiZynthFinder: a fast, robust and flexible
    open-source software for retrosynthetic planning. *J Cheminform* 2020;12:70.
    doi:10.1186/s13321-020-00472-1.
18. Ertl P, Schuffenhauer A. Estimation of synthetic accessibility score. *J Cheminform*
    2009;1:8.
19. Van Duyne GD, Standaert RF, Karplus PA, Schreiber SL, Clardy J. Atomic structure of
    FKBP–FK506, an immunophilin–immunosuppressant complex. *Science* 1991;252:839–842.
    doi:10.1126/science.1709302.

---

## Figures

**Figure 1.** The affinity↔synthesizability trade-off. GNINA CNNaffinity vs SA score for
the 37 ranked candidates, colored by AiZynthFinder solvability; MD-stable leads (stars) and
the MD-failing best synthesizable candidate (red ring) are marked. The desirable
high-affinity/low-SA/synthesizable region is empty.
(`docs/paper/figures/fig1_affinity_synth_tradeoff.png`)

**Figure 2.** Pose stability over 1 ns MD. Ligand heavy-atom RMSD (mean and max) for L1–L3
(stable) and the synthesizable candidate S1* (drifts out).
(`docs/paper/figures/fig2_md_stability.png`)

**Figure 3.** Leads contact the FKBP12/FKBP12.6 selectivity residues. Minimum lead–residue
distances to Trp59 and Ile90; all leads engage Trp59, L2/L3 also engage Ile90.
(`docs/paper/figures/fig3_selectivity_contacts.png`)
