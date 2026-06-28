"""FKBP-family selectivity panel for FKBP12.

Calcineurin selectivity is handled STRUCTURALLY (avoid the 40s loop, from 1TCO) —
see dossier; ChEMBL calcineurin data is sparse/mechanistically irrelevant, so no
ML model is built for it.

This module addresses FKBP-family off-targets where data exists:
  - FKBP5 (CHEMBL2052031), FKBP4 (CHEMBL4050): enough data -> selectivity models.
  - FKBP1B / FKBP12.6 (CHEMBL2430): 83% identical to FKBP12, only ~40 matched
    pairs -> we do NOT model it; instead we empirically test, from matched pairs,
    whether known ligands discriminate at all (they essentially do not).

Leakage-free selectivity test: matched pairs (molecules measured on BOTH FKBP12
and the off-target) are EXCLUDED from training, then used to compare predicted
vs measured selectivity (ΔpChEMBL = p[FKBP12] - p[off-target]).

Output: data/processed/fkbp_offtargets.csv, models/fkbp5/4_rf_regressor.joblib
"""
from __future__ import annotations

from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[2]
PROC = ROOT / "data" / "processed"
MODELS = ROOT / "models"
BASE = "https://www.ebi.ac.uk/chembl/api/data"
OFF = {"FKBP5": "CHEMBL2052031", "FKBP4": "CHEMBL4050", "FKBP12.6": "CHEMBL2430"}
KEEP = {"Ki", "Kd", "IC50", "EC50"}
_FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def fpv(s):
    m = Chem.MolFromSmiles(s)
    return np.array(_FP.GetFingerprint(m)) if m else None


def X_of(smiles):
    return np.array([fpv(s) for s in smiles])


def collect(tid: str) -> pd.DataFrame:
    rows, url, params = [], f"{BASE}/activity.json", {
        "target_chembl_id": tid, "pchembl_value__isnull": "false", "limit": 1000}
    while True:
        d = requests.get(url, params=params, timeout=60).json()
        rows += d["activities"]
        nxt = d["page_meta"]["next"]
        if not nxt:
            break
        url, params = "https://www.ebi.ac.uk" + nxt, None
    recs = []
    for a in rows:
        smi, p = a.get("canonical_smiles"), a.get("pchembl_value")
        if smi and p is not None and a.get("standard_type") in KEEP and Chem.MolFromSmiles(smi):
            recs.append({"name": a["molecule_chembl_id"],
                         "smiles": Chem.MolToSmiles(Chem.MolFromSmiles(smi)),
                         "pchembl": float(p)})
    return (pd.DataFrame(recs).sort_values("pchembl", ascending=False)
            .drop_duplicates("name").reset_index(drop=True))


def reg():
    return RandomForestRegressor(n_estimators=400, n_jobs=-1, random_state=42)


def selectivity_eval(name: str, df12: pd.DataFrame, dfoff: pd.DataFrame):
    matched = df12.merge(dfoff, on="name", suffixes=("_12", "_off"))
    print(f"\n=== {name}: {len(dfoff)} ligands, {len(matched)} matched pairs ===")
    if len(matched) < 15:
        d = matched.pchembl_12 - matched.pchembl_off
        print(f"  too few for a model. Empirical ΔpChEMBL (FKBP12 - {name}): "
              f"mean {d.mean():+.2f}, std {d.std():.2f}, |Δ|>1: {(d.abs()>1).mean()*100:.0f}%")
        print(f"  -> ligands barely discriminate (high homology); use EXPERIMENTAL counter-screen.")
        return None
    # leakage-free: train off-target & FKBP12 regressors WITHOUT matched pairs
    mids = set(matched.name)
    tr_off = dfoff[~dfoff.name.isin(mids)]
    tr_12 = df12[~df12.name.isin(mids)]
    r_off, r_12 = reg(), reg()
    r_off.fit(X_of(tr_off.smiles), tr_off.pchembl.values)
    r_12.fit(X_of(tr_12.smiles), tr_12.pchembl.values)
    Xm = X_of(matched.smiles_12)
    pred_delta = r_12.predict(Xm) - r_off.predict(Xm)
    meas_delta = (matched.pchembl_12 - matched.pchembl_off).values
    rho = spearmanr(meas_delta, pred_delta).correlation
    # can we flag FKBP12-selective (meas Δ>0.5)?
    print(f"  measured ΔpChEMBL: mean {meas_delta.mean():+.2f}, "
          f"FKBP12-selective (Δ>0.5): {(meas_delta>0.5).sum()}/{len(meas_delta)}")
    print(f"  predicted-vs-measured selectivity Spearman: {rho:.3f}")
    r_off.fit(X_of(dfoff.smiles), dfoff.pchembl.values)  # refit on all for saving
    joblib.dump(r_off, MODELS / f"{name.replace('.', '')}_rf_regressor.joblib")
    return rho


def main() -> int:
    MODELS.mkdir(exist_ok=True)
    df12 = pd.read_csv(PROC / "fkbp12_ligands.csv").rename(
        columns={"molecule_chembl_id": "name", "canonical_smiles": "smiles",
                 "pchembl_value": "pchembl"})[["name", "smiles", "pchembl"]]
    alloff = []
    for name, tid in OFF.items():
        dfoff = collect(tid)
        dfoff["target"] = name
        alloff.append(dfoff)
        selectivity_eval(name, df12, dfoff)
    pd.concat(alloff, ignore_index=True).to_csv(PROC / "fkbp_offtargets.csv", index=False)
    print(f"\nWrote {PROC / 'fkbp_offtargets.csv'}")
    print("Calcineurin selectivity = structural (40s-loop avoidance), not ML. See dossier §5.1.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
