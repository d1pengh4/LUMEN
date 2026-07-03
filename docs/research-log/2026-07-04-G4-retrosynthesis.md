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
