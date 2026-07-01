# Final Candidate Ranking (GNINA-validated) — Track A milestone

Date: 2026-07-02 · Related: gnina-gate-PASS, phase2b
Result: `results/phase2b/FINAL_candidates.csv` (37 candidates ranked by GNINA CNNaffinity)

## What this is
The 33 developable + 4 TIER-1 de novo candidates, ranked by the ONLY affinity
signal we validated (GNINA CNNaffinity; gate ROC-AUC 0.799 vs docking 0.455).
Reference from the gate: known actives 6.34 ± 0.36, decoys 5.80 ± 0.53. So
CNNaffinity ≥ 6.3 ≈ "known-active level".

## Result (honest)
- CNNaffinity range across candidates: 3.09 – 6.56.
- **Only 3/37 reach ≥ 6.3 (active-like).** Most novel de novo molecules score at
  decoy level (~5.8) or below — expected: pocket-conditioned generation yields
  diverse/novel matter but most are not strong binders.
- **Top candidate converges across ALL independent signals:**
  `CCC[C@@H]1CC(=O)N[C@@]2(C)[C@H](C1)[C@@H](Cc1ccc(F)cc1N)CC[C@@H]2O`
  - GNINA CNNaffinity 6.56 (above the active mean 6.34)
  - docking −8.5 (seed-avg ~−7.5), developable (QED 0.72, SA 4.39, logP 3.4, MW 362)
  - novel (Tanimoto ~0.20 to known FKBP12 ligands)
  This is a TIER-1 molecule: independently flagged by docking + developability +
  now GNINA. The strongest single hypothesis from the whole pipeline.
- Next two (CNNaffinity 6.34, 6.34): a small acid/amino-phenyl developable
  molecule (SA 2.28!, very synthesizable) and another TIER-1.

## Honest limits (unchanged)
- CNNaffinity is a validated-but-imperfect proxy (gate 0.799, not 1.0); NOT
  measured binding. These are the best *hypotheses to test*, not confirmed binders.
- Only 3 candidates reach active-level score — a small, honest shortlist.
- Selectivity (FKBP12.6/calcineurin) still not addressed here (FEP step pending).
- Real validation = MD stability (next) + ultimately wet-lab.

## Deliverable
A defensible, GNINA-validated shortlist of de novo FKBP12 candidates — the core
computational output for the methodology paper. Top ~3 are the priority molecules
for any downstream MD / synthesis / assay.

## Track A status
step1 MM-GBSA/GNINA rescoring: DONE (gate passed, candidates ranked).
Next gates: MD pose-stability of top candidates; selectivity FEP; retrosynthesis.
