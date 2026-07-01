"""Local developability panel for candidate molecules (no GPU, no web service).

Computes a transparent, reproducible ADMET-proxy + synthetic-accessibility +
drug-likeness profile from RDKit only. These are RULE-BASED proxies (not a
trained ADMET model), stated honestly:
  - SA_score        : synthetic accessibility (1 easy .. 10 hard; Ertl 2009)
  - QED             : quantitative drug-likeness
  - Lipinski / Veber: oral-druglike rule flags
  - PAINS           : pan-assay-interference structural alerts
  - logP, TPSA, HBD/HBA, rotatable bonds, aromatic rings, fraction Csp3
  - BBB/absorption heuristics (TPSA/logP based)

Usage: python developability_panel.py <in.csv> <out.csv>  (needs 'smiles' col)
"""
from __future__ import annotations

import os
import sys

from rdkit import Chem, RDConfig, RDLogger
from rdkit.Chem import Crippen, Descriptors, QED, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # noqa: E402

import pandas as pd  # noqa: E402

RDLogger.DisableLog("rdApp.*")


def _pains_catalog():
    p = FilterCatalogParams()
    p.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    return FilterCatalog(p)


PAINS = _pains_catalog()


def profile(smiles: str) -> dict:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return {"smiles": smiles, "valid": False}
    mw = Descriptors.MolWt(m)
    logp = Crippen.MolLogP(m)
    hbd = rdMolDescriptors.CalcNumHBD(m)
    hba = rdMolDescriptors.CalcNumHBA(m)
    tpsa = rdMolDescriptors.CalcTPSA(m)
    rot = rdMolDescriptors.CalcNumRotatableBonds(m)
    lipinski = (mw <= 500) + (logp <= 5) + (hbd <= 5) + (hba <= 10)
    veber = (rot <= 10) and (tpsa <= 140)
    return {
        "smiles": Chem.MolToSmiles(m), "valid": True,
        "mw": round(mw, 1), "logp": round(logp, 2), "tpsa": round(tpsa, 1),
        "hbd": hbd, "hba": hba, "rot_bonds": rot,
        "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(m),
        "fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(m), 2),
        "qed": round(QED.qed(m), 3),
        "sa_score": round(sascorer.calculateScore(m), 2),
        "lipinski_pass": lipinski == 4,
        "veber_pass": bool(veber),
        "pains_alert": PAINS.HasMatch(m),
        # crude oral-CNS heuristics (honest proxies, not predictions)
        "gi_absorption_ok": tpsa <= 140 and logp <= 5,
        "bbb_likely": tpsa < 90 and 1 <= logp <= 4,
    }


def main() -> int:
    src, dst = sys.argv[1], sys.argv[2]
    df = pd.read_csv(src)
    col = "smiles" if "smiles" in df.columns else "canonical_smiles"
    out = pd.DataFrame([profile(s) for s in df[col]])
    out = out[out.valid]
    # developability flag: synthesizable + drug-like + clean
    out["developable"] = (out.sa_score <= 4.5) & (out.qed >= 0.5) & \
                         out.lipinski_pass & out.veber_pass & (~out.pains_alert)
    out.to_csv(dst, index=False)
    print(f"profiled {len(out)} molecules -> {dst}")
    print(f"  developable (SA<=4.5, QED>=0.5, Lipinski+Veber, no PAINS): {out.developable.sum()}")
    print(f"  SA_score: median {out.sa_score.median():.2f} (lower=easier)")
    print(f"  QED median {out.qed.median():.2f}; PAINS alerts {out.pains_alert.sum()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
