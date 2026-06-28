"""Collect FKBP12 (human FKBP1A, CHEMBL1902) binding ligands from ChEMBL.

Produces a cleaned, deduplicated dataset of bioactive molecules with RDKit
descriptors -- the training/benchmark set for the affinity model and a source
of real chemical matter for validating the de novo pipeline.

Output: data/processed/fkbp12_ligands.csv (+ .summary.json)
See docs/target-dossier-FKBP12.md.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, QED

RDLogger.DisableLog("rdApp.*")

TARGET = "CHEMBL1902"  # human FKBP1A, single protein
BASE = "https://www.ebi.ac.uk/chembl/api/data"
OUT = Path(__file__).resolve().parents[2] / "data" / "processed"
KEEP_TYPES = {"Ki", "Kd", "IC50", "EC50"}


def fetch_activities(target: str) -> list[dict]:
    """Page through all activities with a pChEMBL value for the target."""
    rows, url = [], f"{BASE}/activity.json"
    params = {"target_chembl_id": target, "pchembl_value__isnull": "false", "limit": 1000}
    while True:
        r = requests.get(url, params=params, timeout=60)
        r.raise_for_status()
        data = r.json()
        rows.extend(data["activities"])
        nxt = data["page_meta"]["next"]
        if not nxt:
            break
        url, params = "https://www.ebi.ac.uk" + nxt, None  # next already has query string
    return rows


def clean(rows: list[dict]) -> pd.DataFrame:
    recs = []
    for a in rows:
        smi = a.get("canonical_smiles")
        pchembl = a.get("pchembl_value")
        if not smi or pchembl is None or a.get("standard_type") not in KEEP_TYPES:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        recs.append({
            "molecule_chembl_id": a["molecule_chembl_id"],
            "canonical_smiles": Chem.MolToSmiles(mol),
            "standard_type": a["standard_type"],
            "pchembl_value": float(pchembl),
            "mw": round(Descriptors.MolWt(mol), 1),
            "logp": round(Descriptors.MolLogP(mol), 2),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "tpsa": round(Descriptors.TPSA(mol), 1),
            "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "qed": round(QED.qed(mol), 3),
        })
    df = pd.DataFrame(recs)
    if df.empty:
        return df
    # Dedupe by canonical SMILES, keep the strongest measured affinity
    df = (df.sort_values("pchembl_value", ascending=False)
            .drop_duplicates("canonical_smiles")
            .reset_index(drop=True))
    return df


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"Fetching activities for {TARGET} ...")
    rows = fetch_activities(TARGET)
    print(f"  {len(rows)} raw activity records")
    df = clean(rows)
    if df.empty:
        print("No usable ligands.")
        return 1
    csv = OUT / "fkbp12_ligands.csv"
    df.to_csv(csv, index=False)

    summary = {
        "target": TARGET,
        "n_unique_ligands": len(df),
        "pchembl": {"min": df.pchembl_value.min(), "max": df.pchembl_value.max(),
                    "median": round(df.pchembl_value.median(), 2)},
        "by_type": df.standard_type.value_counts().to_dict(),
        "lipinski_pass": int(((df.mw <= 500) & (df.logp <= 5) &
                              (df.hbd <= 5) & (df.hba <= 10)).sum()),
        "qed_median": round(df.qed.median(), 3),
    }
    (OUT / "fkbp12_ligands.summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nWrote {len(df)} unique ligands -> {csv}")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
