"""Train + HONESTLY validate ML scorers for FKBP12 to replace the uninformative
docking score (enrichment ROC-AUC 0.455).

Two models:
  (1) Binder triage classifier  : known FKBP12 ligands vs random drug-like ChEMBL.
  (2) Potency regressor          : pChEMBL among the 725 binders.

Validation is deliberately strict (medical-grade rigor):
  - random-split AND scaffold-split CV (scaffold split = generalization to new chemotypes)
  - external test on the held-out enrichment set (same protocol as docking)
  - similarity audit (max Tanimoto of enrichment actives/decoys to train positives)
    to expose whether the enrichment test is "easy".

Leakage control: enrichment actives/decoys are excluded from training.
Known limitation (documented, not hidden): selectivity (calcineurin, FKBP12.6) is
NOT modeled here; predictions outside the applicability domain are unreliable.

Output: models/fkbp12_rf_classifier.joblib, models/fkbp12_rf_regressor.joblib
"""
from __future__ import annotations

import random
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import requests
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import mean_absolute_error, r2_score, roc_auc_score
from sklearn.model_selection import GroupKFold, cross_val_predict, cross_val_score

RDLogger.DisableLog("rdApp.*")
random.seed(42)

ROOT = Path(__file__).resolve().parents[2]
PROC = ROOT / "data" / "processed"
MODELS = ROOT / "models"
BASE = "https://www.ebi.ac.uk/chembl/api/data"
N_NEG = 1200
_FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def fp(s):
    m = Chem.MolFromSmiles(s)
    return _FP.GetFingerprint(m) if m else None


def arr(s):
    f = fp(s)
    return np.array(f) if f is not None else None


def scaffold(s):
    m = Chem.MolFromSmiles(s)
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=m) or s
    except Exception:
        return s


def featurize(smiles):
    X, keep = [], []
    for i, s in enumerate(smiles):
        v = arr(s)
        if v is not None:
            X.append(v); keep.append(i)
    return np.array(X), keep


def fetch_negatives(exclude: set[str]) -> pd.DataFrame:
    cache = PROC / "chembl_negatives.csv"
    if cache.exists():
        return pd.read_csv(cache)[lambda d: ~d.name.isin(exclude)]
    rows, seen = [], set()
    for off in random.sample(range(0, 400000), 60):
        r = requests.get(f"{BASE}/molecule.json", timeout=60, params={
            "molecule_properties__full_mwt__gte": 250, "molecule_properties__full_mwt__lte": 550,
            "molecule_properties__qed_weighted__gte": 0.4, "limit": 50, "offset": off})
        if r.status_code != 200:
            continue
        for m in r.json()["molecules"]:
            cid = m["molecule_chembl_id"]
            smi = (m.get("molecule_structures") or {}).get("canonical_smiles")
            if smi and cid not in seen and cid not in exclude and Chem.MolFromSmiles(smi):
                seen.add(cid); rows.append({"name": cid, "smiles": smi})
        if len(rows) >= N_NEG:
            break
    df = pd.DataFrame(rows).drop_duplicates("smiles").head(N_NEG)
    df.to_csv(cache, index=False)
    return df


def grp_ids(smiles):
    g = [scaffold(s) for s in smiles]
    u = {s: i for i, s in enumerate(sorted(set(g)))}  # sorted -> deterministic folds
    return np.array([u[s] for s in g])


def main() -> int:
    MODELS.mkdir(exist_ok=True)
    enrich = pd.read_csv(PROC / "enrichment_set.csv")
    held = set(enrich.name)
    allpos = pd.read_csv(PROC / "fkbp12_ligands.csv")
    pos = allpos[~allpos.molecule_chembl_id.isin(held)]
    neg = fetch_negatives(exclude=set(allpos.molecule_chembl_id) | held)

    # ---------- (1) Binder triage classifier ----------
    Xp, kp = featurize(pos.canonical_smiles.tolist())
    Xn, kn = featurize(neg.smiles.tolist())
    Xc = np.vstack([Xp, Xn]); yc = np.array([1] * len(Xp) + [0] * len(Xn))
    gc = grp_ids(list(pos.canonical_smiles.iloc[kp]) + list(neg.smiles.iloc[kn]))
    clf = RandomForestClassifier(n_estimators=400, class_weight="balanced", n_jobs=-1, random_state=42)
    print(f"[classifier] positives={len(Xp)} negatives={len(Xn)}")
    print(f"  random-split   CV ROC-AUC: {cross_val_score(clf, Xc, yc, cv=5, scoring='roc_auc', n_jobs=-1).mean():.3f}")
    print(f"  scaffold-split CV ROC-AUC: {cross_val_score(clf, Xc, yc, cv=GroupKFold(5), groups=gc, scoring='roc_auc', n_jobs=-1).mean():.3f}")
    clf.fit(Xc, yc)
    joblib.dump(clf, MODELS / "fkbp12_rf_classifier.joblib")

    # external test on held-out enrichment set
    Xe, ke = featurize(enrich.smiles.tolist()); ye = enrich.label.values[ke]
    pe = clf.predict_proba(Xe)[:, 1]
    print(f"  external enrichment ROC-AUC: {roc_auc_score(ye, pe):.3f}  (docking was 0.455)")

    # similarity audit: is the enrichment test easy?
    pos_fps = [fp(s) for s in pos.canonical_smiles if fp(s) is not None]
    sim = lambda s: max(DataStructs.BulkTanimotoSimilarity(fp(s), pos_fps)) if fp(s) else np.nan
    ma = enrich[enrich.label == 1].smiles.map(sim).median()
    md = enrich[enrich.label == 0].smiles.map(sim).median()
    print(f"  similarity audit: actives median Tanimoto {ma:.2f} vs decoys {md:.2f} "
          f"-> {'EASY test (actives are analogues)' if ma - md > 0.3 else 'balanced'}")

    # ---------- (2) Potency regressor (the hard, useful task) ----------
    Xr, kr = featurize(allpos.canonical_smiles.tolist())
    yr = allpos.pchembl_value.values[kr]
    gr = grp_ids(list(allpos.canonical_smiles.iloc[kr]))
    reg = RandomForestRegressor(n_estimators=400, n_jobs=-1, random_state=42)
    pred = cross_val_predict(reg, Xr, yr, cv=GroupKFold(5), groups=gr, n_jobs=-1)
    print(f"\n[potency regressor] n={len(yr)} (scaffold-split CV)")
    print(f"  Spearman rho: {spearmanr(yr, pred).correlation:.3f}")
    print(f"  R2: {r2_score(yr, pred):.3f}   MAE: {mean_absolute_error(yr, pred):.2f} log units")
    reg.fit(Xr, yr)
    joblib.dump(reg, MODELS / "fkbp12_rf_regressor.joblib")

    print("\nSaved models to models/. NOTE: no selectivity model yet (calcineurin/FKBP12.6);"
          " apply an applicability-domain (Tanimoto) gate before trusting de novo predictions.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
