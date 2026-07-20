# What de novo design can and cannot deliver on a hard target: an honest, reproducible pipeline for non-immunosuppressive FKBP12 ligands against pulmonary arterial hypertension

**Author:** Yulho Choi¹*
**Affiliation:** ¹Independent Researcher
**\*Correspondence:** eurochoic@gmail.com

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
Mechanistically, FKBP12 acts as a constitutive repressor bound to the GS domain shared by
the type-I BMP/TGF-β receptors (ALK1, ALK2, ALK3); FK506 displaces it, permitting receptor
autophosphorylation and downstream SMAD1/5 (and MAPK/ID1) signaling. FK506 engages all
three ALK receptors, whereas rapamycin — which acts through a distinct FKBP12–mTOR
interaction rather than receptor displacement — engages this axis only via FKBP12 bound to
ALK2, a narrower and mechanistically distinct effect [5].

Its therapeutic use is precluded by immunosuppression: the FK506–FKBP12 complex recruits
calcineurin through an effector surface adjacent to the FKBP12 "40s loop." A precise point
matters here and is easy to over-simplify: FK506 activates BMP signaling through **both** a
calcineurin-dependent and a calcineurin-independent (FKBP12-displacement) arm [5], so
avoiding calcineurin is not a free benefit — a calcineurin-sparing ligand captures only the
displacement arm and may therefore deliver **partial** efficacy relative to FK506 itself.
This is an honest efficacy–toxicity trade-off, not a pure win, and we return to it as an
open risk in §3. The trade-off is de-risked, though not eliminated, by direct precedent:
GPI-1046, a non-immunosuppressive FKBP12 ligand from the neuroimmunophilin class, has been
reported to induce ALK3 transcription and BMP-responsive signaling to a comparable degree
independent of calcineurin modulation [6] — evidence that the displacement arm alone can be
sufficient, at least in the systems it was tested in. Non-immunosuppressive FKBP ligands are
therefore not a novel chemical concept (the neuroimmunophilin program dates to the late
1990s and included compounds that later had mixed clinical results and reports of
FKBP12-independent activity); what is novel here is the combination of this known design
principle with de novo generative design purpose-built for the PAH BMP-axis, multi-ALK
engagement, and FKBP12.6 selectivity — a narrower and more defensible claim than "first
FKBP12 ligand," and the one we make throughout.

This motivates a specific
design objective: **fill the FKBP12 FK506 pocket while avoiding the calcineurin-recruiting
40s loop (residues 32–45)**, and discriminate the 83%-identical paralog FKBP12.6 (FKBP1B).
The latter is not a generic selectivity nicety: FKBP12.6 physically tethers the cardiac
ryanodine receptor RyR2, and its displacement has been mechanistically linked to
arrhythmia and pulmonary-hypertension-related cardiac dysfunction, making FKBP12.6
cross-reactivity a specific safety liability to design against, not merely an
off-target-affinity concern. On the competitive landscape: sotatercept is an injectable
ligand trap; FK506 is the immunosuppressive mechanistic comparator; and at least one other
oral small-molecule BMP-axis enhancer with a distinct mechanism — an inhibitor of the
E3 ubiquitin ligase SMURF1, which normally targets BMP receptors and SMADs for degradation
— is already in PAH clinical trials (LTP001, ClinicalTrials.gov NCT05135000,
NCT05764265), underscoring that oral BMP-axis restoration is an active area of competing
mechanisms, not a hypothetical niche.

De novo generative design is often reported through its successes. Here we instead report
a fully documented pipeline **and its validated failure modes** on FKBP12 — a shallow,
hydrophobic, low-druggability pocket with a near-identical paralog — because knowing which
computational instruments are trustworthy on such a target is itself a transferable
result. Our contributions are (i) a reproducible generation-to-validation pipeline; (ii)
an explicit capability map of docking, ligand-ML, CNN rescoring, MD, retrosynthesis, and
structural selectivity analysis on this target; and (iii) three MD-stable, novel lead
hypotheses with a concrete, falsifiable selectivity rationale.

**Related work.** Two families of methods dominate directed molecular optimization:
policy-gradient reinforcement learning over a generative prior, exemplified by REINVENT
[20], and evolutionary search directly over molecular graphs, exemplified by graph-based
genetic algorithms (GB-GA) [21]. We adopt the latter (§2.3, §4.5) because it needs orders
of magnitude fewer reward evaluations to converge — a deliberate fit to a frugal,
GPU-constrained compute budget rather than a claim of general superiority. Our
developability filters follow established medicinal-chemistry heuristics — quantitative
drug-likeness (QED) [22], the Lipinski rule-of-five [23], and the Veber rules [24] — and
our fragment-recombination baseline and structural-alert screen use, respectively, the
BRICS fragmentation scheme [25] and PAINS filters [26], both applied on ECFP4/Morgan
fingerprints [27] via RDKit [28]. Unlike most generative-design case studies, which report
a single successful scoring/generation stack, we report the stack that *did not* work
(docking-score ranking, ligand-based selectivity) alongside the one that did, in the belief
that negative capability results generalize better across targets than positive ones.

## 2. Results

### 2.1 Target, pocket, and the structural selectivity rule

After a premise correction (from BMPR2 to FKBP12; §4.1), we assembled a five-structure
reference set from the RCSB, each serving a distinct role: **1FKF** (FKBP12–FK506 [19]) defines
the drug pocket itself; **1B6C** (FKBP12–ALK5 GS-domain peptide) is the on-target proxy
interface, since no experimental FKBP12–ALK1/2/3 structure exists; **1TCO**
(FKBP12–FK506–calcineurin) is the off-target interface whose avoidance defines the
selectivity rule; **2PPN** is the apo structure (pocket geometry without a bound ligand);
and **1F40** (FKBP12–GPI-1046) is the non-immunosuppressive template precedent showing
that pocket occupancy without calcineurin recruitment is structurally achievable. Chain and
residue identities were verified against each PDB header before use.

With FK506 removed from 1FKF, fpocket identified the FK506 site as
the top pocket (druggability 0.517, volume 277 Å³, hydrophobicity 61.2, polar SASA ≈ 0;
pocket centroid 4.35 Å from native FK506) — a shallow, hydrophobic, moderately druggable
site. Interface mapping (BioPython, 5 Å) against the on-target proxy (1B6C; 22
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

Three conclusions are load-bearing, and each is worth unpacking mechanistically rather than
taken at face value.

**Why docking score fails here.** The redock control shows the *pose* is right (FK506, the
native high-affinity macrocycle, scores best at −10.3 kcal/mol and reproduces the
crystallographic geometry to 0.81 Å), so receptor preparation and the search itself are not
the problem. What fails is *ranking*: among the six most potent drug-like anchors, the
compound with the highest reported potency (pChEMBL 9.35) does **not** score best by Vina;
scores cluster tightly (−6.7 to −7.6 kcal/mol) regardless of measured affinity. On the
formal 20-active/40-decoy enrichment set this compresses to ROC-AUC 0.455 with negligible
early enrichment (EF@10% = 1.18, EF@20% = 1.07 — both close to the chance value of 1.0);
active and decoy mean scores are statistically indistinguishable (−7.28 vs. −7.33
kcal/mol). The mechanistic reading: FKBP12's pocket is shallow and almost entirely
hydrophobic (polar SASA ≈ 0), so the Vina/smina scoring function — which rewards
hydrophobic contact area — cannot distinguish a true binder from any other molecule that
merely fills the same greasy cavity. This is a target-class property (shallow apolar
pockets defeat contact-count scoring functions generically), not an implementation error,
which is why we report it rather than tune around it.

**Why the ligand-based ML model is a family detector, not an FKBP12-specific one.**
The classifier and regressor were trained on 725 ChEMBL FKBP12 binders against 1,200
random drug-like ChEMBL decoys, with the 60-molecule enrichment set held out entirely from
training to prevent leakage. Scaffold-split cross-validation is the headline number
because a naive external-test ROC-AUC of 1.000 is misleadingly easy: a similarity audit
showed the enrichment actives sit at median Tanimoto 0.75 to the training positives versus
0.21 for decoys, so that number mostly reflects "near-FKBP12 chemical space vs. random
space," not fine-grained discrimination — hence we quote the harder scaffold-split numbers
(0.993 triage, ρ = 0.748 potency) throughout. The same classifier scores known FKBP5
ligands at mean P(binder) ≈ 0.79 (vs. 0.99 for FKBP12 actives and 0.02 for random
molecules) and an independent FKBP5 potency regressor is itself accurate (scaffold-CV ρ =
0.78) — yet the *difference* between the two paralogs' predicted potency is unlearnable
(ρ ≈ −0.01 for FKBP12-vs-FKBP5 selectivity, ρ ≈ 0 for FKBP12-vs-FKBP12.6 in §2.6). Two
individually good regressors whose difference carries no signal is the operational
definition of "not ligand-ML-solvable" used throughout this paper. (A weaker, low-n
FKBP4 regressor, ρ = 0.30 on 101 ligands, was also built as a secondary off-target check
and is not load-bearing.)

**Why GNINA works where the other two do not.** GNINA CNNaffinity separates actives (6.34
± 0.36) from decoys (5.80 ± 0.53) with ROC-AUC 0.799 and was adopted as the affinity
signal, gated by an applicability-domain Tanimoto threshold ≥ 0.30. Unlike Vina/smina's
hand-crafted contact terms, GNINA's convolutional scoring function was trained end-to-end
on 3D structural data to predict affinity directly, which evidently captures pocket-filling
quality beyond raw hydrophobic contact area on this target — though, as emphasized
throughout, ROC-AUC 0.799 is a validated proxy, not measured binding. Selectivity was
therefore addressed structurally and (prospectively) experimentally, never by a learned
reward.

### 2.3 Generation and goal-directed optimization

A BRICS baseline produced drug-like but low-diversity molecules (scaffold diversity 11%).
An early head-to-head pilot (N ≈ 50 per route, before the full gated batches) made the
exploitation-vs-exploration tension between the two routes concrete: BRICS candidates had
higher median logP (5.3 vs. 1.7) and larger median MW (450 vs. 296 Da) than DiffSBDD
candidates, and — most tellingly — 346 of the BRICS molecules fell inside the ligand-ML
scorer's applicability domain with binder\_prob ≥ 0.5, versus **zero** of the DiffSBDD
molecules. This is not evidence the DiffSBDD molecules do not bind; it is evidence that a
ligand-based scorer trained on known chemotypes is structurally blind to genuinely novel
ones, which is exactly why §2.2 restricts that scorer to in-domain candidates and why
structure-based tools (docking pose, then GNINA, then MD) carry the validation burden for
novel scaffolds. Structure-based screening of this pilot batch (41/49 successfully docked)
placed 16/41 at or better than known drug-like binders' typical range (≤ −7.0 kcal/mol) and
2/41 better than any known drug-like binder in our set (≤ −8.0), while explicitly *not*
treating this as validation, given docking score's known lack of ranking power (§2.2) —
only as a prioritization hint later confirmed or rejected by GNINA. As a rigor check on
that hint, we re-docked the four resulting TIER-1 candidates (novel, developable, and
dock ≤ −7 kcal/mol) at three additional random seeds: scores were stable (SD 0.02–0.48
kcal/mol across seeds) for three of four, but one candidate's initial best-seed score
(−8.54) proved seed-optimistic (true, seed-averaged ≈ −7.5) — establishing the practice,
used from then on, of reporting seed-averaged rather than best-single-seed docking scores.

Pocket-conditioned diffusion (DiffSBDD) at full scale produced diverse, novel scaffolds (scaffold
diversity 89%, novelty 0.18, PAINS 0%) but a best GNINA CNNaffinity of only 6.56. A second,
independent DiffSBDD batch (271 molecules, salvaged from a session that paused mid-run)
passed the same gate (diversity 91%, novelty 0.17, 24 developable, PAINS 3%) yet reached
only 6.35, with none ≥ 7.0 — confirming, rather than a single unlucky run, that brute-force
sampling plateaus at ≈ 6.5 regardless of batch and motivating directed optimization.

**A documented failure mode: seed dominance.** Our first optimization attempt seeded the
genetic-algorithm population with the drug-like ChEMBL anchors themselves (as a source of
"good" starting material) and maximized the potency regressor. The reported best score
(pred\_pchembl 8.39) never moved across 8 generations, and the top "optimized" molecules
had applicability-domain Tanimoto = 1.0 — i.e., they were verbatim training-set binders,
such as
`COC[C@H]1CN(Cc2ccccn2)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1cc(Cl)cc(Cl)c1`. The regressor was
simply recognizing molecules it had memorized during training; the GA had optimized
nothing. We record this as a general caution for ligand-based optimization loops: seeding
a population with training-set molecules under a training-set-derived reward is a latent
reward-hacking mode that a flat/high score alone will not reveal — only an
applicability-domain audit of the winners does. The fix used from then on: anchors serve
only as fragment **donors** during breeding, never as population seeds, and a hard novelty
box (0.30 ≤ AD Tanimoto ≤ 0.65) forbids near-duplicates of known binders from winning.

With that fix, round 1 of the genetic algorithm (8 generations, §4.5) raised the potency
proxy from 5.16 to 7.34 among genuinely novel molecules (AD Tanimoto 0.44–0.52), and GNINA
validation of the top 15 confirmed a real gain (best 7.06 vs. the 6.56 baseline; 7/15 above
6.56; 1 developable/novel candidate ≥ 7.0) — directed optimization beating brute-force
scaling, on the validated ranker, not just the fast proxy.

**A second, equally instructive failure: over-optimizing the proxy.** Because round 1's
gain (+0.5 GNINA) was well clear of our pre-registered safety-stop threshold (< 0.1 gain
over two rounds), we ran a round 2, reseeded from the round-1 winners. On the fast ML proxy
alone, round 2 looked like unambiguous further progress: the potency proxy rose again
(7.34 → 7.78), and — more suggestively — the binder-classifier and potency-regressor
**began to agree** (binder probability rose from ≈ 0.37 in round 1, itself flagged as a
yellow flag for over-optimistic extrapolation, to 0.50–0.65 in round 2, on molecules
sharing a recurring nitrile-diketopiperazine-sulfonamide motif). Every proxy-space signal
pointed toward a genuine improvement. GNINA validation of the round-2 top 15 said
otherwise: best CNNaffinity **6.64 — lower than round 1's 7.06**, not higher. The GA had
overfit to the fast regressor, drifting into a region the proxy rewards but the validated
tool does not — textbook reward-hacking, and a second concrete demonstration (after the
seed-dominance failure) that rising proxy-space metrics, even ones that look
mutually corroborating, are not evidence of real improvement without the validated check.
Our pre-registered safety-stop rule (stop if GNINA does not improve ≥ 0.1 over two rounds
and the STOP criterion is still unmet) triggered exactly as designed, and G2 optimization
halted with round 1's 7.06 candidate as the best validated result; round 2's molecules were
not carried forward.

Independently of round 2, we then applied a **higher-quality re-scoring pass**
(exhaustiveness 8 instead of 4, five poses instead of one — a resampling change, not a
threshold adjustment) to the existing round-1-era shortlist (the 7.06 candidate plus the
broader ≥ 6.56 group). More thorough conformational sampling alone raised the confirmed
best from 7.06 to 7.21 and moved the STOP-criterion count from 1/3 to 2/3 candidates ≥ 7.0.
Three developable, novel leads emerged (Table 2), all traceable to round-1-era
optimization plus this resampling step — none to round 2.

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
candidate (`O=C(O)CCCCN(CCO)c1cccc(F)c1F`; GNINA 6.34, SA [18] 2.28, 1-step) **failed MD** (RMSD
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

### 2.8 Pipeline scale: the full generation-to-lead funnel (Table 4)

Table 4 gives the scale of every stage, end to end, to make explicit how much winnowing
separates the initial pools from the three reported leads — a generation-and-optimization
pipeline is only as credible as its funnel is transparent.

**Table 4. Generation-to-lead funnel.**

| Stage | Input → output | Key metric |
|---|---|---|
| ChEMBL curation | 821 activities → 725 unique ligands → 68 drug-like anchors | pChEMBL median 6.71 |
| Route A (BRICS, CPU) | 400 novel drug-like candidates | scaffold diversity 11% |
| Route B batch 1 (DiffSBDD) | 272 generated → 37 GNINA-scored (33 developable + 4 TIER-1) | diversity 89%, best GNINA 6.56 |
| Route B batch 2 (G1 scale-up) | 271 generated (salvaged mid-run) → 24 developable scored | diversity 91%, best GNINA 6.35 |
| G2 round 1 (post seed-dominance fix) | 8 GA generations → top 15 GNINA-validated | best GNINA 7.06 (1/3 ≥ 7.0) |
| G2 round 2 (reseeded; proxy rose 7.34→7.78, GNINA-validated anyway) | 8 GA generations → top 15 GNINA-validated | best GNINA **6.64 (worse)** → safety-stop triggered, round 2 discarded |
| HQ rescore (exhaustiveness 8, round-1-era shortlist only) | top 10 re-ranked with deeper sampling | best GNINA 7.21 (2/3 ≥ 7.0) |
| G3 MD filter | 3 candidates simulated | 3/3 stay bound (max RMSD ≤ 2.9 Å) |
| G4 retrosynthesis | 3 leads + 37-candidate pool screened | 0/3 leads solved; 7/37 pool solved, capped at GNINA 6.34 |
| **Final leads** | **3 (L1, L2, L3)** | GNINA 6.95–7.21, all MD-stable, all novel |

The funnel shows roughly three orders of magnitude of winnowing (hundreds of generated
molecules down to 3 leads) driven overwhelmingly by the developability and
applicability-domain gates, not by the affinity gate alone — consistent with the
affinity↔synthesizability tension formalized in §2.5.

## 3. Discussion

On a shallow hydrophobic pocket with a near-identical paralog, two workhorse methods fail
in instructive ways: docking-*score* ranking (ROC 0.455) and ligand-based selectivity
prediction (ρ ≈ 0). What remains effective is a layered strategy — CNN rescoring for
affinity, explicit-solvent MD for pose reality, and structural analysis for selectivity —
which together still yield a defensible, testable shortlist. We argue that publishing this
capability map, including the failures, is more useful to the field than another
success-only generative case study, because the failures are target-class properties, not
project-specific accidents.

**Why this generalizes beyond FKBP12.** The two failures are not idiosyncrasies of this
protein; they are predictable consequences of two structural properties that recur across
many drug targets. A shallow, near-apolar pocket (druggability 0.517, polar SASA ≈ 0) will
defeat any scoring function whose terms are dominated by hydrophobic contact area,
independent of the specific protein — the same failure mode is expected on other shallow
immunophilin-, chaperone-, or protein–protein-interaction pockets. Likewise, a
near-identical paralog (83% sequence identity, differing at only 2 of ~13 pocket residues)
will defeat any ligand-based classifier or regressor trained on whole-molecule fingerprints,
because the training signal simply does not encode a 2-residue structural difference — this
generalizes to any target with a close paralog family (kinases, other immunophilins, ion
channel subtypes). The actionable takeaway for practitioners on structurally similar
targets is therefore to benchmark docking-score ranking and ligand-ML selectivity
*before* trusting either, rather than assuming a pipeline validated on an easier target
class will transfer.

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
assay has been performed. Every candidate is a hypothesis for experimental testing.

**Experimental validation roadmap.** We propose a staged, go/no-go cascade rather than a
single confirmatory assay, because the value of the program hinges on more than potency
alone. *Tier 1* (FK506-competition fluorescence polarization and/or PPIase/rotamase
inhibition, SPR or ITC for K_d) asks the minimal question — does the molecule bind FKBP12
at all — and is the gate for everything downstream. *Tier 2*, the decisive tier, runs two
counter-screens in parallel: a calcineurin-phosphatase/NFAT-luciferase (or IL-2) assay
testing the central design hypothesis that avoiding the 40s loop avoids calcineurin
recruitment, and an FKBP12.6 binding counter-screen testing the Trp59/Ile90 selectivity
hypothesis directly. *Tier 3* (BMP/SMAD1/5 phosphorylation or a BRE-luciferase reporter in
pulmonary-artery smooth-muscle or C2C12 cells) asks whether the mechanism is functional,
not merely structural. *Tier 4* (ADMET/PK, then rodent PAH models) is downstream of a
positive Tier 1–3 result. We prioritize **L3** for first synthesis — the best MD stability
(1.05/1.60 Å mean/max) and both selectivity handles (Trp59 + Ile90) via the shortest
predicted route (4 steps) — followed by **L2** (smallest, most drug-like, both handles),
keeping **L1** (highest affinity, one handle) as a comparator; this ordering deliberately
weights the selectivity thesis over raw predicted affinity, since selectivity, not potency,
is this program's differentiator from FK506 itself. Order-of-magnitude estimates for
synthesis and the Tier 1–2 cascade on two leads are in the low tens of thousands of USD
over several months — modest relative to the decisiveness of the resulting data, which
either supports or falsifies the core non-immunosuppressive-selectivity hypothesis before
any further computational or synthetic investment.

**Risks and open uncertainties, stated plainly.** Beyond the proxy-vs-measurement caveats
above, four risks are specific to this program and are not resolved by more computation.
(i) The magnitude of the FKBP12/FKBP12.6 selectivity window is unknown — only two
conservative pocket residues differ, so the achievable window could be narrow even if the
sign is favorable; relative FEP could narrow this uncertainty computationally but is itself
resource-gated, and the assay (Tier 2) is the more decisive test. (ii) No compound in this
program has yet been synthesized or measured; every number reported here is a validated
*proxy*, and proxies can be wrong in ways benchmarking cannot fully anticipate. (iii) The
leads' novel cores are not solved by automated retrosynthesis, which is a cost-and-schedule
risk for a synthesis campaign, not evidence the cores are unmakeable — a synthetic chemist
routinely closes gaps that template-based retrosynthesis tools cannot see, but this adds
real time and expense. (iv) The pocket's low druggability score (0.517) caps how much
potency this chemotype class can plausibly deliver; if Tier 1 shows only weak binding,
the honest conclusion is a pocket-class limitation, not a failure of the design method.

**A note on intellectual property.** The three leads are novel relative to the curated
ChEMBL FKBP12 ligand set used for training and generation (Tanimoto < 0.4 to any ligand in
that set), which is the narrow, correctly-scoped sense in which "novel" is used throughout
this paper. That is a necessary but not sufficient condition for patentability, and a
preliminary keyword-level patent-literature search (Google Patents/USPTO, not a
structure-based professional freedom-to-operate search) changes our assessment for two of
the three leads specifically. L1 and L3 share a fused proline-derived diaza-bicyclic
sulfonamide core, and this general design principle — a rigid bicyclic proline/pipecolic
core bearing a sulfonamide arm as an FKBP-binding pharmacophore — closely matches a
chemotype pioneered by Guilford Pharmaceuticals in the 1990s (e.g., US 5,721,256, on
rotamase-inhibiting sulfonamide compounds) and claimed more broadly in later filings
(e.g., CA 2,446,795A1, covering sulfonamide/sulfamide/urea/heterobicycle FKBP-binding
ligands), with a structurally related diazabicyclo[4.3.1]decane-sulfonamide FKBP-ligand
family patented as recently as 2016–2017 for the paralogs FKBP51/52 (EP 3,097,103 B1).
L1/L3's specific exemplars may or may not be individually disclosed in this prior art —
that determination requires a professional structure-based search — but the genus-level
novelty claim for their core chemotype is now in real doubt. L2's oligo-pyrrolidine core
did not return a comparable match in this preliminary pass and is the comparatively
cleaner candidate, though this is not a clean bill of health, only the absence of a hit in
a limited, non-exhaustive search. No compound has been synthesized or assayed, which
remains a separate, independent prerequisite for any filing. We report this finding rather
than the earlier, more optimistic framing because it is decision-relevant and directly
testable by any reader: at present the contribution is a methodology and a hypothesis set,
explicitly not an IP-cleared candidate, and — for L1/L3 specifically — not even confidently
a novel one at the chemotype level.

## 4. Methods

**4.1 Target definition.** The original BMPR2-pocket premise was rejected after a
literature check (BMPR2 haploinsufficiency; sotatercept approval) and the target redefined
to FKBP12 (recorded as ADR-001). Structures 1FKF, 1B6C, 1TCO, 2PPN, 1F40 were retrieved
from the RCSB and identities verified from headers.

**4.2 Pocket and interfaces.** fpocket [7] on 1FKF (FK506 removed) quantified the pocket;
interfaces were mapped with BioPython at a 5 Å cutoff.

**4.3 Ligand data.** FKBP12 activities were collected from ChEMBL (CHEMBL1902) [8];
drug-like anchors were filtered at MW ≤ 500, QED ≥ 0.5, pChEMBL ≥ 6.

**4.4 Scoring.** Docking used smina [9] (AutoDock Vina scoring [10], exhaustiveness 8, seed
42); redock RMSD and a 20-active/40-decoy enrichment set assessed pose and ranking. Binder
triage (RandomForestClassifier, 400 trees, balanced class weight), a pChEMBL regressor
(RandomForestRegressor, 400 trees), and an analogous FKBP5 selectivity probe were trained
on 2048-bit ECFP4 fingerprints [27] (fixed random\_state 42) and evaluated by 5-fold
scaffold-grouped cross-validation (GroupKFold). Developability used QED [22], Lipinski [23]
and Veber [24] rule compliance, and PAINS structural alerts [26] (SA ≤ 4, QED ≥ 0.6, 1 ≤
logP ≤ 4). GNINA CNNaffinity [11] rescoring was validated on a 60-molecule enrichment set
and adopted with a Tanimoto ≥ 0.30 applicability-domain gate.

**4.5 Generation and optimization.** BRICS recombination [25] provided a CPU baseline;
DiffSBDD [12] performed pocket-conditioned diffusion generation on a GPU. Screening-stage
GNINA rescoring used exhaustiveness 4 / 1 pose for throughput; the final ranking of
optimization winners used a deliberately more thorough exhaustiveness 8 / 5 poses (both
seed 42) — a resampling change, not a threshold adjustment. The genetic algorithm ran for 8
generations per round, maximizing the potency regressor (pred\_pchembl) subject to hard
constraints (0.30 ≤ AD Tanimoto ≤ 0.65, QED ≥ 0.55, 250 ≤ MW ≤ 500 Da, RDKit-valid), using
BRICS recombination of a fragment-donor pool (developable generated candidates + drug-like
ChEMBL anchors); anchors were never placed in the seed population itself (§2.3). Round 2
reseeded from the GNINA-confirmed round-1 winners plus the round-1 GA population; its top 15
were GNINA-validated exactly as round 1's were, and rejected by a pre-registered
safety-stop rule when the validated score fell rather than rose (§2.3).

**4.6 Molecular dynamics.** Complexes were prepared with PDBFixer (missing hydrogens added
at pH 7.0); ligands were docked into the pocket (autobox on FK506). Systems (Amber ff14SB
[13] + GAFF-2.11 [14] + TIP3P [15], 1.0 nm solvent padding, 0.15 M NaCl) were built with
OpenMM/openmmforcefields [16] SystemGenerator, energy-minimized, and run for 1 ns
(LangevinMiddleIntegrator, 300 K, 1 ps⁻¹ friction, 4 fs timestep enabled by hydrogen-mass
repartitioning to 3 amu, HBonds + rigid-water constraints, Monte Carlo barostat at 1 bar/300
K) using OpenCL (a CUDA/driver PTX incompatibility on the cloud GPU forced the OpenCL
platform; a CUDA→OpenCL→CPU fallback is implemented for reproducibility on heterogeneous
hardware). Ligand heavy-atom RMSD versus the initial pose was computed with MDTraj
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

## Supplementary Methods

**S1. Key resources table.** Software and data resources, with the pinned versions used.

| Resource | Version / identifier | Role |
|---|---|---|
| Python (local) | 3.11 (conda-forge) | main `lumen` environment |
| RDKit | ≥ 2024.03 | cheminformatics (BRICS, ECFP4, QED, PAINS, SA score) |
| Biopython | ≥ 1.83 | interface mapping, sequence alignment |
| fpocket | conda-forge build | pocket detection/druggability |
| scikit-learn | conda-forge build | RandomForest classifier/regressor |
| smina / OpenBabel | conda-forge build (`lumen-dock` env, isolated to avoid an ICU dependency conflict) | docking, format conversion |
| GNINA | v1.0 binary (self-installed per Colab session) | CNN affinity rescoring |
| DiffSBDD | GitHub, full-atom conditional checkpoint `crossdocked_fullatom_cond` (Zenodo record 8183747) | pocket-conditioned generation |
| AiZynthFinder | 4.4.1, public USPTO expansion/ring-breaker/filter models, ZINC stock (17.4 M) | retrosynthesis |
| OpenMM + openmmforcefields | Colab-session install (Python 3.12, CUDA 12 build) | MD system building and integration |
| PDBFixer, MDTraj | Colab-session install | structure repair; trajectory analysis |
| Structures | RCSB PDB: 1FKF, 1B6C, 1TCO, 2PPN, 1F40 | see §4.1 |
| Ligand data | ChEMBL, target CHEMBL1902 | see §4.3 |

**S2. GPU execution model.** All GPU-dependent steps (DiffSBDD generation, GNINA rescoring,
OpenMM MD) were written as self-contained scripts (`colab_*.py` in `src/generation/`,
`src/scoring/`, `src/md/`) executed on transient Google Colab T4 GPU instances via a
CLI wrapper (`colab run --gpu T4/L4 --timeout <s> script.py`); each script installs its own
dependencies, downloads its own structures/checkpoints, and prints delimited marker lines
(e.g., `===LUMEN_SMILES_BEGIN/END===`, `STEP ...`) that are parsed back into `results/`
logs after the transient VM is torn down. This was a deliberate frugal-compute design (no
persistent GPU infrastructure), not an artifact of the local macOS development environment
lacking a GPU.

**S3. Notable reproducibility fixes.** For transparency and to save future implementers the
same debugging cycles, the nontrivial cross-platform issues encountered and fixed were: (i)
macOS-side SSL certificate verification failures against the Colab kernel WebSocket,
resolved by pointing `SSL_CERT_FILE` at the `certifi` bundle; (ii) silent loss of
Colab-side subprocess stdout/stderr under `ipykernel`, resolved by explicit capture-and-print
wrapping; (iii) a `torch.load(weights_only=True)` default (PyTorch ≥ 2.6) rejecting the
DiffSBDD checkpoint's `argparse.Namespace` payload, resolved via a `sitecustomize.py`
override (justified because the checkpoint originates from a trusted, pinned Zenodo record);
(iv) `Bio.SeqUtils.three_to_one` removed in Biopython ≥ 1.80, replaced with
`protein_letters_3to1`; (v) an Amber/OpenMM CUDA-driver PTX incompatibility on the cloud
GPU for the MD step, resolved with an explicit CUDA→OpenCL→CPU platform fallback (§4.6); and
(vi) a scaffold-grouping bug in the potency regressor's cross-validation (Python `set()`
iteration order is not guaranteed stable), which silently made the reported Spearman ρ
non-deterministic across runs until fixed to `sorted(set(...))` — after the fix, ρ = 0.748
(R² = 0.562) reproduces exactly across independent runs.

**S4. Docking-score reproducibility check.** Because docking score is known to have weak
ranking power (§2.2), an early TIER-1 candidate set (4 molecules meeting dock ≤ −7 kcal/mol,
developable, and novel) was re-docked at three additional random seeds beyond the default
(seed 42): score standard deviation across seeds ranged 0.02–0.48 kcal/mol, and one
candidate's default-seed score (−8.54) was found to be seed-optimistic relative to its
seed-averaged value (≈ −7.5). This motivated reporting seed-averaged rather than
best-single-seed docking scores throughout the rest of the pipeline (§2.3).

## Data and code availability
All code (`src/`), gate logs (`docs/research-log/`), the master record
(`docs/PROJECT_LOG.md`), result tables (`results/`), and figures (`docs/paper/figures/`)
are available in the project repository. Structures are from the RCSB PDB; ligand data from
ChEMBL; retrosynthesis models/stock from the AiZynthFinder public data release.
**Supplementary Information** (Tables S1–S7: the full 37-candidate ranked pool, docking
and MD raw data, selectivity-contact geometry, structural/synthesizability triage, the
developability panel, and extended reproducibility notes) is provided in
`docs/paper/manuscript_SI_EN.md`.

## Author contributions
Y.C.: conceptualization, methodology, software, formal analysis, investigation, data
curation, writing — original draft, writing — review & editing. AI-assisted tooling
(Claude Code, Anthropic) was used throughout for pipeline scripting, literature/reference
verification, and manuscript drafting under the author's direction; the author verified
all reported numbers against the underlying result files and is responsible for all
scientific claims and conclusions.

## Competing interests
The authors declare no competing interests.

## Acknowledgements
No external funding was received for this work. Compute costs were limited to free/low-cost
cloud GPU sessions (Google Colab), consistent with the frugal-compute design documented in
§4.9 and the Supplementary Methods.

## References
*(Verification pass, this revision: every reference below was independently checked via
live literature search — author list, journal, volume/pages, and DOI — against the primary
source or its indexing record (PubMed/DOI registry). All 28 references now carry a
DOI and matched on author/journal/pages. References 1, 3, and 6, previously flagged as
reconstructed from memory, are now confirmed to match the primary source exactly and the
flag is removed. Reference 28 (RDKit) is software without a fixed paper DOI; a specific
release/Zenodo DOI should still be pinned at copy-editing to match the exact version used.)*

1. Morrell NW, Adnot S, Archer SL, et al. Cellular and molecular basis of pulmonary
   arterial hypertension. *J Am Coll Cardiol* 2009;54(1 Suppl):S20–S31.
   doi:10.1016/j.jacc.2009.04.018.
2. Lane KB, Machado RD, Pauciulo MW, et al. Heterozygous germline mutations in BMPR2,
   encoding a TGF-β receptor, cause familial primary pulmonary hypertension. *Nat Genet*
   2000;26:81–84. doi:10.1038/79226.
3. Machado RD, Eickelberg O, Elliott CG, et al. Genetics and genomics of pulmonary arterial
   hypertension. *J Am Coll Cardiol* 2009;54(1 Suppl):S32–S42. doi:10.1016/j.jacc.2009.04.015
   (NMD/haploinsufficiency mechanism for truncating BMPR2 variants).
4. Hoeper MM, Badesch DB, Ghofrani HA, et al. Phase 3 trial of sotatercept for treatment
   of pulmonary arterial hypertension (STELLAR). *N Engl J Med* 2023;388(16):1478–1490.
   doi:10.1056/NEJMoa2213558.
5. Spiekerkoetter E, et al. FK506 activates BMPR2, rescues endothelial dysfunction, and
   reverses pulmonary hypertension. *J Clin Invest* 2013;123(8):3600–3613. doi:10.1172/JCI65592.
6. Gold BG, Zeleny-Pooley M, Wang M-S, Chaturvedi P, Armistead DM. A nonimmunosuppressant
   FKBP-12 ligand increases nerve regeneration. *Exp Neurol* 1997;147(2):269–278.
   doi:10.1006/exnr.1997.6630 (the GPI-1046 neuroimmunophilin precedent for
   pocket-occupancy without calcineurin recruitment).
7. Le Guilloux V, Schmidtke P, Tuffery P. Fpocket: an open source platform for ligand
   pocket detection. *BMC Bioinformatics* 2009;10:168. doi:10.1186/1471-2105-10-168.
8. Mendez D, Gaulton A, Bento AP, et al. ChEMBL: towards direct deposition of bioassay
   data. *Nucleic Acids Res* 2019;47(D1):D930–D940. doi:10.1093/nar/gky1075.
9. Koes DR, Baumgartner MP, Camacho CJ. Lessons learned in empirical scoring with smina
   from the CSAR 2011 benchmarking exercise. *J Chem Inf Model* 2013;53(8):1893–1904.
   doi:10.1021/ci300604z.
10. Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking with a
    new scoring function, efficient optimization, and multithreading. *J Comput Chem*
    2010;31(2):455–461. doi:10.1002/jcc.21334.
11. McNutt AT, Francoeur P, Aggarwal R, et al. GNINA 1.0: molecular docking with deep
    learning. *J Cheminform* 2021;13:43. doi:10.1186/s13321-021-00522-2.
12. Schneuing A, Harris C, Du Y, et al. Structure-based drug design with equivariant
    diffusion models (DiffSBDD). *Nat Comput Sci* 2024;4(12):899–909.
    doi:10.1038/s43588-024-00737-x (preprint arXiv:2210.13695).
13. Maier JA, Martinez C, Kasavajhala K, Wickstrom L, Hauser KE, Simmerling C. ff14SB:
    improving the accuracy of protein side chain and backbone parameters from ff99SB.
    *J Chem Theory Comput* 2015;11(8):3696–3713. doi:10.1021/acs.jctc.5b00255.
14. Wang J, Wolf RM, Caldwell JW, Kollman PA, Case DA. Development and testing of a
    general Amber force field (GAFF). *J Comput Chem* 2004;25(9):1157–1174.
    doi:10.1002/jcc.20035.
15. Jorgensen WL, Chandrasekhar J, Madura JD, Impey RW, Klein ML. Comparison of simple
    potential functions for simulating liquid water (TIP3P). *J Chem Phys*
    1983;79(2):926–935. doi:10.1063/1.445869.
16. Eastman P, Galvelis R, Peláez RP, et al. OpenMM 8: molecular dynamics simulation with
    machine-learning potentials. *J Phys Chem B* 2024;128(1):109–116.
    doi:10.1021/acs.jpcb.3c06662.
17. Genheden S, Thakkar A, Chadimová V, et al. AiZynthFinder: a fast, robust and flexible
    open-source software for retrosynthetic planning. *J Cheminform* 2020;12:70.
    doi:10.1186/s13321-020-00472-1.
18. Ertl P, Schuffenhauer A. Estimation of synthetic accessibility score of drug-like
    molecules based on molecular complexity and fragment contributions. *J Cheminform*
    2009;1:8. doi:10.1186/1758-2946-1-8.
19. Van Duyne GD, Standaert RF, Karplus PA, Schreiber SL, Clardy J. Atomic structure of
    FKBP–FK506, an immunophilin–immunosuppressant complex. *Science* 1991;252:839–842.
    doi:10.1126/science.1709302.
20. Olivecrona M, Blaschke T, Engkvist O, Chen H. Molecular de-novo design through deep
    reinforcement learning. *J Cheminform* 2017;9:48. doi:10.1186/s13321-017-0235-x.
21. Jensen JH. A graph-based genetic algorithm and generative model/Monte Carlo tree search
    for the exploration of chemical space. *Chem Sci* 2019;10(12):3567–3572.
    doi:10.1039/C8SC05372C.
22. Bickerton GR, Paolini GV, Besnard J, Muresan S, Hopkins AL. Quantifying the chemical
    beauty of drugs. *Nat Chem* 2012;4(2):90–98. doi:10.1038/nchem.1243.
23. Lipinski CA, Lombardo F, Dominy BW, Feeney PJ. Experimental and computational
    approaches to estimate solubility and permeability in drug discovery and development
    settings. *Adv Drug Deliv Rev* 1997;23(1–3):3–25. doi:10.1016/S0169-409X(96)00423-1.
24. Veber DF, Johnson SR, Cheng HY, Smith BR, Ward KW, Kopple KD. Molecular properties that
    influence the oral bioavailability of drug candidates. *J Med Chem* 2002;45(12):2615–
    2623. doi:10.1021/jm020017n.
25. Degen J, Wegscheid-Gerlach C, Zaliani A, Rarey M. On the art of compiling and using
    'drug-like' chemical fragment spaces. *ChemMedChem* 2008;3(10):1503–1507.
    doi:10.1002/cmdc.200800178.
26. Baell JB, Holloway GA. New substructure filters for removal of pan assay interference
    compounds (PAINS) from screening libraries and for their exclusion in bioassays.
    *J Med Chem* 2010;53(7):2719–2740. doi:10.1021/jm901137j.
27. Rogers D, Hahn M. Extended-connectivity fingerprints. *J Chem Inf Model*
    2010;50(5):742–754. doi:10.1021/ci100050t.
28. Landrum G, et al. RDKit: open-source cheminformatics software (version used throughout
    this work; exact release and Zenodo DOI to be pinned at copy-editing).

*Citation numbers reflect order of addition, not order of first appearance in text; a
copy-editing pass should resequence [20]–[28] into the main numbering by first citation.*

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
