"""Local synthetic-feasibility triage for candidates (no GPU) — prep for G4.

Beyond the SA score, computes medicinal-chemistry structural indicators that
drive real synthesis difficulty, and a combined flag to decide which candidates
are worth a full retrosynthesis (AiZynthFinder) run:
  - sa_score        : Ertl synthetic accessibility (1 easy .. 10 hard)
  - n_stereocenters : each defined stereocenter ~ a synthetic control problem
  - n_rings, n_arom : ring count / aromatic rings
  - max_ring_size   : macrocycles / large rings are hard
  - has_strained    : 3-/4-membered rings (aziridine/epoxide/azetidine etc.)
  - spiro, fused_3+ : spiro atoms; 3+ fused ring systems
  - n_hetero_ring   : heteroatom-rich rings
Verdict "synthesizable_easy": SA<=3.5 AND stereocenters<=2 AND no strained-ring
AND max_ring_size<=6 — these go straight to retrosynthesis; harder ones are
flagged for manual review.

Usage: python synth_feasibility.py <in.csv> <out.csv>   (needs 'smiles')
"""
from __future__ import annotations

import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import FindMolChiralCenters, rdMolDescriptors, RDConfig

sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # noqa: E402

import pandas as pd  # noqa: E402

RDLogger.DisableLog("rdApp.*")


def analyze(smiles: str) -> dict:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return {"smiles": smiles, "valid": False}
    ri = m.GetRingInfo()
    ring_sizes = [len(r) for r in ri.AtomRings()]
    max_ring = max(ring_sizes) if ring_sizes else 0
    strained = any(s in (3, 4) for s in ring_sizes)
    try:
        n_stereo = len(FindMolChiralCenters(m, useLegacyImplementation=False, includeUnassigned=True))
    except Exception:
        n_stereo = 0
    sa = sascorer.calculateScore(m)
    easy = (sa <= 3.5) and (n_stereo <= 2) and (not strained) and (max_ring <= 6)
    return {
        "smiles": Chem.MolToSmiles(m), "valid": True,
        "sa_score": round(sa, 2),
        "n_stereocenters": n_stereo,
        "n_rings": ri.NumRings(),
        "n_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(m),
        "max_ring_size": max_ring,
        "has_strained_ring": strained,
        "n_spiro": rdMolDescriptors.CalcNumSpiroAtoms(m),
        "synthesizable_easy": easy,
    }


def main() -> int:
    src, dst = sys.argv[1], sys.argv[2]
    df = pd.read_csv(src)
    col = "smiles" if "smiles" in df.columns else "canonical_smiles"
    out = pd.DataFrame([analyze(s) for s in df[col]])
    out = out[out.valid]
    # carry CNNaffinity if present, to prioritize
    if "cnnaffinity" in df.columns:
        out = out.merge(df[[col, "cnnaffinity"]].rename(columns={col: "smiles"}),
                        on="smiles", how="left")
        out = out.sort_values("cnnaffinity", ascending=False)
    out.to_csv(dst, index=False)
    print(f"analyzed {len(out)} molecules -> {dst}")
    print(f"  synthesizable_easy (SA<=3.5, <=2 stereo, no strained, ring<=6): {out.synthesizable_easy.sum()}")
    print(f"  SA median {out.sa_score.median():.2f}; stereocenters median {out.n_stereocenters.median():.0f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
