"""Phase 2 PoC (CPU): fragment-based de novo generation + scoring + ranking.

This is the honest CPU stand-in for the GPU diffusion model (DiffSBDD):
  1. decompose validated drug-like FKBP12 binders into BRICS fragments
  2. recombine fragments into novel molecules (RDKit BRICSBuild)
  3. filter: valid, drug-like (MW 250-500, QED>=0.5), novel (not a known ligand)
  4. score with the validated FKBP12Scorer (binder_prob + potency + applicability)
  5. rank and write a candidate table

Honest limits: candidates are recombinations of known-binder fragments (modest
novelty); the scorer is an FKBP-family triage + potency model and does NOT score
selectivity. Outputs are HYPOTHESES TO TEST, not drug candidates.

Output: results/phase2/ranked_candidates.csv
"""
from __future__ import annotations

import itertools
import random
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import BRICS, Descriptors, QED

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "scoring"))
from score_molecules import FKBP12Scorer  # noqa: E402

RDLogger.DisableLog("rdApp.*")
random.seed(42)
PROC = ROOT / "data" / "processed"
OUT = ROOT / "results" / "phase2"
MAX_ATTEMPTS = 8000
TARGET_N = 400


def fragment_pool() -> list:
    anchors = pd.read_csv(PROC / "fkbp12_druglike_anchors.csv")
    frags = set()
    for smi in anchors.canonical_smiles:
        m = Chem.MolFromSmiles(smi)
        if m:
            frags |= set(BRICS.BRICSDecompose(m))
    mols = [Chem.MolFromSmiles(f) for f in frags]
    return [m for m in mols if m is not None]


def generate(known: set[str]) -> list[str]:
    frags = fragment_pool()
    print(f"fragment pool: {len(frags)} BRICS fragments")
    random.seed(42)  # BRICSBuild uses random when scrambleReagents=True
    builder = BRICS.BRICSBuild(frags, onlyCompleteMols=True, scrambleReagents=True)
    cands, seen = [], set()
    for m in itertools.islice(builder, MAX_ATTEMPTS):
        if m is None:
            continue
        try:
            Chem.SanitizeMol(m)
        except Exception:
            continue
        smi = Chem.MolToSmiles(m)
        if smi in seen or smi in known:
            continue
        seen.add(smi)
        mw = Descriptors.MolWt(m)
        if 250 <= mw <= 500 and QED.qed(m) >= 0.5:
            cands.append(smi)
        if len(cands) >= TARGET_N:
            break
    return cands


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    known = set(pd.read_csv(PROC / "fkbp12_ligands.csv").canonical_smiles)
    cands = generate(known)
    print(f"generated {len(cands)} novel drug-like candidates")
    if not cands:
        return 1

    scored = FKBP12Scorer().score(cands)
    scored = scored[scored.valid].copy()
    # rank: prefer in-domain, then binder probability, then potency, then QED
    scored = scored.sort_values(
        ["in_domain", "binder_prob", "pred_pchembl", "qed"],
        ascending=[False, False, False, False]).reset_index(drop=True)
    scored.to_csv(OUT / "ranked_candidates.csv", index=False)

    print(f"\nWrote {OUT / 'ranked_candidates.csv'} ({len(scored)} scored)")
    nd = scored[scored.in_domain]
    print(f"in-domain candidates: {len(nd)}; "
          f"binder_prob>=0.5 & in-domain: {len(nd[nd.binder_prob>=0.5])}")
    cols = ["binder_prob", "pred_pchembl", "ad_tanimoto", "qed", "mw"]
    print("\nTop 10 (in-domain):")
    print(nd.head(10)[cols + ["smiles"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
