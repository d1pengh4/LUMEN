"""Build a docking enrichment set: known FKBP12 actives vs property-matched decoys.

Actives  = top drug-like FKBP12 binders (pChEMBL >= 6).
Decoys   = random drug-like ChEMBL molecules (assumed non-binders to FKBP12),
           filtered to a similar MW/QED window and excluding any known actives.

A docking method with real discriminative power should rank actives above decoys
(ROC-AUC > 0.5, positive early enrichment). Output: data/processed/enrichment_set.csv
"""
from __future__ import annotations

import random
from pathlib import Path

import pandas as pd
import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, QED

RDLogger.DisableLog("rdApp.*")
random.seed(42)

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "processed"
N_ACTIVES = 20
N_DECOYS = 40
BASE = "https://www.ebi.ac.uk/chembl/api/data"


def load_actives() -> pd.DataFrame:
    df = pd.read_csv(OUT / "fkbp12_druglike_anchors.csv")
    df = df.sort_values("pchembl_value", ascending=False).head(N_ACTIVES)
    return pd.DataFrame({"name": df.molecule_chembl_id, "smiles": df.canonical_smiles,
                         "label": 1, "pchembl": df.pchembl_value})


def fetch_decoys(exclude: set[str]) -> pd.DataFrame:
    """Random drug-like ChEMBL molecules, MW 300-500, as presumed non-binders."""
    rows, seen = [], set()
    # sample several random offsets to avoid one contiguous block
    for off in random.sample(range(0, 200000), 12):
        r = requests.get(f"{BASE}/molecule.json", timeout=60, params={
            "molecule_properties__full_mwt__gte": 300,
            "molecule_properties__full_mwt__lte": 500,
            "molecule_properties__qed_weighted__gte": 0.5,
            "limit": 50, "offset": off})
        if r.status_code != 200:
            continue
        for m in r.json()["molecules"]:
            cid = m["molecule_chembl_id"]
            struct = m.get("molecule_structures") or {}
            smi = struct.get("canonical_smiles")
            if not smi or cid in exclude or cid in seen:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            seen.add(cid)
            rows.append({"name": cid, "smiles": Chem.MolToSmiles(mol), "label": 0,
                         "pchembl": None, "mw": round(Descriptors.MolWt(mol), 1),
                         "qed": round(QED.qed(mol), 3)})
        if len(rows) >= N_DECOYS * 2:
            break
    df = pd.DataFrame(rows).drop_duplicates("smiles").head(N_DECOYS)
    return df[["name", "smiles", "label", "pchembl"]]


def main() -> int:
    actives = load_actives()
    all_known = set(pd.read_csv(OUT / "fkbp12_ligands.csv").molecule_chembl_id)
    decoys = fetch_decoys(exclude=all_known)
    out = pd.concat([actives, decoys], ignore_index=True)
    out.to_csv(OUT / "enrichment_set.csv", index=False)
    print(f"actives={len(actives)}  decoys={len(decoys)}  total={len(out)}")
    print(f"Wrote {OUT / 'enrichment_set.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
