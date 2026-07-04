# G7 (precursor) — FKBP12 vs FKBP12.6 selectivity feasibility

Date: 2026-07-04 · Roadmap G7 · Before committing to expensive relative FEP, answer
the prerequisite: **do the FKBP12 and FKBP12.6 binding pockets even differ?** If the
pocket is identical, no method (FEP included) can find ligand selectivity, and that is
itself the decisive, honest answer. Frugal, CPU-only.

## Method
UniProt sequences P62942 (FKBP12 / FKBP1A) and P68106 (FKBP12.6 / FKBP1B), both 108 aa,
positional comparison (they align 1:1, no indels). Differences intersected with the
FK506-binding pocket residues (Tyr26, Phe36, Asp37, Arg42, Phe46, Val55, Ile56, Trp59,
Tyr82, His87, Ile90/91, Phe99 — mature numbering; target dossier §4) and the calcineurin
40s loop (32–45). Full table: `results/selectivity/fkbp12_vs_126_pocket_diff.csv`.

## Result
- **83.3 % identical — 18 differences over 108 aa.**
- **Only 2 of the ~13 binding-pocket residues differ, and both are conservative:**
  - **Trp59 → Phe** (aromatic→aromatic; FKBP12.6 loses the indole's second ring →
    smaller, less aromatic surface)
  - **Ile90 → Val** (aliphatic→aliphatic; FKBP12.6 loses one methyl → slightly smaller)
- 40s loop: Asp32→Asn (relevant to calcineurin/immunosuppression selectivity, a
  separate axis). The rest of the 18 differences are surface/distal, outside the pocket.

## Interpretation — decisive and honest
The FKBP12 vs FKBP12.6 pockets are **nearly identical**, differing only by two
conservative, slightly-smaller substitutions. This directly explains the two earlier
negatives: ligand-based ML could not predict FKBP12/12.6 selectivity (ρ ≈ 0) and known
ligands barely discriminate (1/40) — there is almost no chemical difference to exploit.

But it is **not hopeless, and it is now specific**: the only structural handles are
- **Trp59 (FKBP12) vs Phe (FKBP12.6):** a ligand that packs against / π-stacks with the
  extra indole ring present only in FKBP12,
- **Ile90 (FKBP12) vs Val (FKBP12.6):** a ligand substituent that fills the extra methyl
  pocket present only in FKBP12 (would clash slightly less favorably in 12.6).
A selective inhibitor must gain contact energy specifically at these two positions.

## Implication for the FEP step (G7)
This is the classic small-ΔΔG regime (cf. FEP resolving FKBP51/52, also ~85% identical),
so relative FEP is the *right* tool — but the selectivity window will be narrow (2
conservative differences), so any predicted selectivity will be small and demands the
retrospective control before trust. Feasibility verdict: selectivity is **physically
possible but hard**; it is an explicit design objective (exploit Trp59/Ile90), not an
assumption. Full FEP remains **resource-gated** (free-tier Colab cannot complete it; the
1 ns MD already lost one VM to reclaim). This structural analysis is the honest,
zero-cost first pass that the roadmap's "retrospective control before expensive FEP"
principle calls for.

## Next
- Check whether L1–L3 poses actually contact Trp59 / Ile90 (per-residue contact analysis
  of the docked/MD poses) — turns "possible in principle" into "possible for our leads".
- If pursued: relative FEP FKBP12↔FKBP12.6 on reference ligands (retrospective control)
  then leads — on paid/borrowed GPU, or flag as an open experimental risk (FKBP12.6
  counter-screen) for the wet-lab stage.
