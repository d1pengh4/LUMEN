# G4 — Synthetic accessibility + retrosynthesis (AiZynthFinder)

Date: 2026-07-04 · Roadmap G4 · Runs after G3 MD PASS. Gate: each lead has
SA ≤ ~4 AND a plausible retrosynthetic route to purchasable building blocks.
Leads = the 3 MD-stable candidates L1/L2/L3 (see 2026-07-02-G3-MD).

## Method
- SA + structural triage: `src/scoring/synth_feasibility.py` →
  `results/synthesis/leads_synth_feasibility.csv`.
- Retrosynthesis: AiZynthFinder 4.4.1 (local, CPU; conda env `aizynth`, onnxruntime
  — no GPU/tensorflow). Public USPTO expansion + ringbreaker + filter models and the
  ZINC purchasable stock (17.4 M compounds), `download_public_data` →
  `data/external/aizynth/` (gitignored). Two searches per lead: default (~100 iter)
  and deep (`config_deep.yml`: 500 iter cap, time 600 s, max_transforms 8; each
  reached ~3.6–5.2 k nodes in ~220–260 s).

## SA / structural triage
| Lead | SA | stereocenters | max ring | arom rings | synthesizable_easy |
|------|----|---------------|----------|------------|--------------------|
| L1 | 4.34 | 3 | 7 | 2 | no |
| L2 | 4.05 | 3 | 6 | 0 | no |
| L3 | 4.39 | 3 | 7 | 2 | no |
SA is moderate (~4.0–4.4, borderline vs the ≤4 target); each has 3 stereocenters and
a fused/bridged ring (max 7). None clear the strict "easy" triage (SA ≤ 3.5).

## Retrosynthesis result (deep search) — `results/synthesis/aizynth/leads_retro_summary.csv`
| Lead | solved | best route steps | precursors in stock | irreducible not-in-stock precursor |
|------|--------|------------------|---------------------|------------------------------------|
| L1 | **no** | 5 | 5/6 | `C=C[C@H]1CNC(=O)[C@@H]2CCC[C@H]1N2` (bare diaza-bicyclic core) |
| L2 | **no** | 6 | 4/5 | `CCOC(=O)C[C@@H]1CCCN1C1CCCN1C1CCCCN1` (oligo-pyrrolidine chain) |
| L3 | **no** | 4 | 4/5 | `C=C[C@H]1CNC(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1cccc(C)c1` |

## Verdict — G4 CONDITIONAL (not cleanly passed; honest limitation)
AiZynthFinder finds **no fully purchasable route** for any lead, even with deep search.
In every case the peripheral chemistry is standard and templated (sulfonylation,
amide coupling, reductive amination, Suzuki) — 4/5–5/6 precursors are catalog ZINC
compounds — but the **irreducible precursor is the novel de novo core itself** (the
fused proline-derived diaza-bicyclic sulfonamide for L1/L3; the oligo-pyrrolidine for
L2). USPTO templates do not disconnect that core to stock.

Interpretation (honest): "not solved by USPTO+ZINC" ≠ "unsynthesizable." The bicyclic
proline-derived scaffold resembles known conformationally-constrained peptidomimetic
cores and a chemist could likely build it via dedicated asymmetric routes — but it is
**not one-click from catalog**, so custom core synthesis is a real cost. This is a
**synthesizability risk flagged for the wet-lab/CRO stage**, not a silent pass and not
a hard stop (cf. the resource-gated selectivity FEP). It also points the generator at a
concrete future objective: bias toward catalog-reachable scaffolds.

## Reproduce
```
conda create -n aizynth python=3.11 && conda install -n aizynth -c conda-forge \
    numba llvmlite rdkit onnxruntime pandas scipy && pip install aizynthfinder   # arm64: llvmlite must come from conda, not pip
download_public_data data/external/aizynth
aizynthcli --config results/synthesis/aizynth/config_deep.yml \
    --smiles leads.smi --output leads_routes_deep.json.gz
```

## Addendum — pool synthesizability screen (2026-07-04)
Motivated by the leads' un-synthesizable cores, screened all 37 GNINA-ranked
candidates (`FINAL_candidates.csv`) through AiZynthFinder (default search).
Cross-table: `results/synthesis/aizynth/pool_retro_x_gnina.csv`.

- **7/37 are fully solved** (complete purchasable route). Clear
  **affinity ↔ synthesizability trade-off**: the top-GNINA candidates (≥6.5,
  incl. the L1–L3 leads at 6.95–7.21) are the un-synthesizable novel-core ones;
  the solvable set tops out at GNINA ~6.34 (≈ known-active mean 6.34 ± 0.36).
- **Best synthesizable candidate — a pragmatic backup lead:**
  `O=C(O)CCCCN(CCO)c1cccc(F)c1F` — GNINA **6.34**, SA **2.28**, QED 0.71,
  **1-step** route from 2 catalog precursors. (This is the "practical lead" first
  flagged in phase2b, now with a confirmed complete route.) Next: MD-test it so we
  have a lead that is simultaneously affinity-plausible, MD-stable AND synthesizable.
- Second: `CC1([C@H](Cc2cccc(Cl)c2)OCCNC(=O)O)CCC1` — GNINA 6.34, SA 3.10, 4-step, 4/4 in stock.

Implication for the endgame: two distinct lead classes — (A) high-affinity novel
scaffolds needing custom synthesis (L1–L3), (B) readily-synthesizable ~known-active-
level binders. A wet-lab program could de-risk with (B) while pursuing (A)'s potency.

### MD test of the synthesizable backup lead (2026-07-04) — FAIL
Ran 1 ns MD (same protocol as G3) on `O=C(O)CCCCN(CCO)c1cccc(F)c1F`.
Ligand heavy-atom RMSD **mean 2.93 / max 6.04 / final 5.45 Å** (vs L1–L3 all ≤ 2.9 max)
→ **does not hold the docked pose**; the small flexible acyclic ligand drifts out.
Log: `results/md/S1_synthesizable_1ns.log`; row in `results/md/g3_md_rmsd.csv`.

Sharpened picture (honest): the trade-off is not just affinity↔synthesizability but
also **rigidity↔synthesizability**. The MD-stable leads (L1–L3) are rigid novel-core
scaffolds (hard to make); the readily-synthesizable candidate is a floppy small acid
that does not stay bound. So there is **no clean "synthesizable + stable + high-affinity"
lead in the current pool** — the de-risk-via-B route is weaker than hoped. L1–L3
(high-affinity, MD-stable, custom synthesis) remain the real leads.
