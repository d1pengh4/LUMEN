"""Score a generated pool (from a colab_diffsbdd_run log) and apply the
scaled-generation STAGE GATE from docs/ROADMAP.md.

Pipeline: parse SMILES -> FKBP12Scorer + developability -> gate metrics.
The gate decides whether this step is solid enough to advance to MM-GBSA.

GATE (scaled generation pool) — STRENGTHENED (2026-07-01, user: raise the bar):
  G1 scaffold diversity  > 70%          (genuinely diverse, not decorated cores)
  G2 novelty (median max Tanimoto to known) < 0.35
  G3 developable count   >= 15          (SA<=4.0, QED>=0.6, 1<=logP<=4, Lipinski+Veber, no PAINS)
  G4 PAINS fraction      < 5%
PASS requires all four. Potency/affinity is deliberately NOT a gate here: the
scorer's affinity estimate is unreliable for out-of-domain novel molecules, so
efficacy is verified later by physics (MM-GBSA/docking), not by a weak filter.
On PASS, writes a high-quality developable shortlist for downstream MM-GBSA.
On FAIL, the step is fixed/redesigned before advancing.

Usage: python src/generation/score_and_gate.py <colab.log>
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Crippen, Descriptors, QED, RDConfig, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.Scaffolds import MurckoScaffold

sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "scoring"))
from score_molecules import FKBP12Scorer  # noqa: E402

RDLogger.DisableLog("rdApp.*")
OUT = ROOT / "results" / "phase2b"
_FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
_PAINS = FilterCatalog(
    (lambda p: (p.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS), p)[1])(FilterCatalogParams()))


def parse_smiles(log: Path) -> list[str]:
    on, out = False, []
    for l in log.read_text().splitlines():
        if "===LUMEN_SMILES_BEGIN===" in l:
            on = True; continue
        if "===LUMEN_SMILES_END===" in l:
            break
        if on and l.strip():
            out.append(l.strip())
    return out


def main() -> int:
    log = Path(sys.argv[1])
    raw = parse_smiles(log)
    smis = list(dict.fromkeys(Chem.MolToSmiles(m) for s in raw if (m := Chem.MolFromSmiles(s))))
    print(f"parsed {len(raw)} lines -> {len(smis)} unique valid molecules")
    if not smis:
        print("no molecules; check the log"); return 1

    scored = FKBP12Scorer().score(smis)
    scored = scored[scored.valid].copy()
    mols = [Chem.MolFromSmiles(s) for s in scored.smiles]
    scored["sa_score"] = [round(sascorer.calculateScore(m), 2) for m in mols]
    scored["logp"] = [round(Crippen.MolLogP(m), 2) for m in mols]
    scored["pains"] = [_PAINS.HasMatch(m) for m in mols]
    scored["lipinski"] = [(Descriptors.MolWt(m) <= 500) + (Crippen.MolLogP(m) <= 5) +
                          (rdMolDescriptors.CalcNumHBD(m) <= 5) + (rdMolDescriptors.CalcNumHBA(m) <= 10) == 4
                          for m in mols]
    scored["veber"] = [(rdMolDescriptors.CalcNumRotatableBonds(m) <= 10) and
                       (rdMolDescriptors.CalcTPSA(m) <= 140) for m in mols]
    # STRENGTHENED developable: SA<=4.0, QED>=0.6, 1<=logP<=4, Lipinski+Veber, no PAINS
    scored["developable"] = (scored.sa_score <= 4.0) & (scored.qed >= 0.6) & \
                            (scored.logp >= 1) & (scored.logp <= 4) & \
                            scored.lipinski & scored.veber & (~scored.pains)
    scored.to_csv(OUT / "gen_pool_scored.csv", index=False)

    # gate metrics
    scafs = {MurckoScaffold.MurckoScaffoldSmiles(mol=m) for m in mols}
    known = [_FP.GetFingerprint(Chem.MolFromSmiles(s))
             for s in pd.read_csv(ROOT / "data/processed/fkbp12_ligands.csv").canonical_smiles
             if Chem.MolFromSmiles(s)]
    nov = float(np.median([max(DataStructs.BulkTanimotoSimilarity(_FP.GetFingerprint(m), known)) for m in mols]))
    g1 = len(scafs) / len(mols)
    g3 = int(scored.developable.sum())
    g4 = float(scored.pains.mean())

    dev = scored[scored.developable].sort_values(["binder_prob", "pred_pchembl"], ascending=False)
    dev.to_csv(OUT / "gen_pool_developable.csv", index=False)

    checks = [
        ("G1 scaffold diversity > 70%", g1 > 0.70, f"{g1*100:.0f}%"),
        ("G2 novelty (median) < 0.35", nov < 0.35, f"{nov:.2f}"),
        ("G3 developable >= 15 (SA<=4,QED>=0.6,logP1-4)", g3 >= 15, f"{g3}"),
        ("G4 PAINS fraction < 5%", g4 < 0.05, f"{g4*100:.0f}%"),
    ]
    print("\n================ SCALED-GENERATION GATE ================")
    for name, ok, val in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}  = {val}")
    verdict = all(ok for _, ok, _ in checks)
    print(f"\n  GATE: {'PASS -> advance to MM-GBSA' if verdict else 'FAIL -> fix before advancing'}")
    print(f"  developable shortlist: {g3} -> results/phase2b/gen_pool_developable.csv")
    return 0 if verdict else 2


if __name__ == "__main__":
    raise SystemExit(main())
