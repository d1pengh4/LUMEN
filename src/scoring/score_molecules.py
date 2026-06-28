"""Unified FKBP12 scoring interface for candidate molecules (the function the
generation phase will call).

For each molecule returns:
  binder_prob   : P(FKBP-pocket binder) [classifier]. NOTE: this detects the
                  FKBP-family pocket-binder chemotype, NOT FKBP12 specificity
                  (it scores FKBP5 ligands ~0.79). Use for triage vs random
                  molecules, not as a selectivity signal. (verification 2026-06-28)
  pred_pchembl  : predicted FKBP12 potency (pChEMBL)          [regressor]
  ad_tanimoto   : max ECFP4 Tanimoto to FKBP12 training actives (applicability domain)
  in_domain     : ad_tanimoto >= AD_THRESHOLD (else predictions are unreliable)
  qed, mw       : RDKit drug-likeness

Deliberately does NOT emit a selectivity score: the selectivity panel showed
ligand-based selectivity is not predictable here (see selectivity-panel log).
Selectivity is enforced structurally (calcineurin 40s loop) + experimentally
(FKBP12.6 counter-screen).

CLI: python score_molecules.py <smiles.csv|->  (column 'smiles')  > scored.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, QED, rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[2]
MODELS = ROOT / "models"
PROC = ROOT / "data" / "processed"
AD_THRESHOLD = 0.30  # below this, model is extrapolating -> untrusted
_FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def _bitvect(mol):
    return _FP.GetFingerprint(mol)


class FKBP12Scorer:
    def __init__(self):
        self.clf = joblib.load(MODELS / "fkbp12_rf_classifier.joblib")
        self.reg = joblib.load(MODELS / "fkbp12_rf_regressor.joblib")
        acts = pd.read_csv(PROC / "fkbp12_ligands.csv").canonical_smiles
        self.ref_fps = [_bitvect(m) for s in acts if (m := Chem.MolFromSmiles(s))]

    def score(self, smiles: list[str]) -> pd.DataFrame:
        out = []
        for s in smiles:
            mol = Chem.MolFromSmiles(s)
            if mol is None:
                out.append({"smiles": s, "valid": False})
                continue
            bv = _bitvect(mol)
            x = np.array(bv).reshape(1, -1)
            ad = max(DataStructs.BulkTanimotoSimilarity(bv, self.ref_fps))
            out.append({
                "smiles": Chem.MolToSmiles(mol), "valid": True,
                "binder_prob": round(float(self.clf.predict_proba(x)[0, 1]), 3),
                "pred_pchembl": round(float(self.reg.predict(x)[0]), 2),
                "ad_tanimoto": round(float(ad), 3),
                "in_domain": bool(ad >= AD_THRESHOLD),
                "qed": round(QED.qed(mol), 3),
                "mw": round(Descriptors.MolWt(mol), 1),
            })
        return pd.DataFrame(out)


def main() -> int:
    src = sys.argv[1] if len(sys.argv) > 1 else "-"
    df = pd.read_csv(sys.stdin if src == "-" else src)
    col = "smiles" if "smiles" in df.columns else "canonical_smiles"
    scored = FKBP12Scorer().score(df[col].tolist())
    scored.to_csv(sys.stdout, index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
