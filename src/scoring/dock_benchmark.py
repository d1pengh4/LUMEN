"""Docking benchmark to validate the scoring pipeline against known binders.

Prepares the FKBP12 receptor from 1FKF (FK506 removed), defines the box from the
native FK506 position, then redocks FK506 (positive control) and docks the top
drug-like ChEMBL anchors. A working pipeline should score these strongly negative.

Engine: smina (env `lumen-dock`). Run with the `lumen` python (needs pandas).
Output: results/phase1/docking/docking_scores.csv
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DOCK_BIN = Path("/usr/local/Caskroom/miniforge/base/envs/lumen-dock/bin")
OBABEL, SMINA = str(DOCK_BIN / "obabel"), str(DOCK_BIN / "smina")
SRC_PDB = ROOT / "data" / "raw" / "structures" / "1FKF.pdb"
WORK = ROOT / "results" / "phase1" / "docking"
N_ANCHORS = 6
BOX = 24.0  # Angstrom (FK506 is large)


def sh(cmd: list[str]) -> str:
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def prep_receptor() -> Path:
    prot = WORK / "receptor_protein.pdb"
    prot.write_text("\n".join(l for l in SRC_PDB.read_text().splitlines()
                              if l.startswith("ATOM")) + "\n")
    pdbqt = WORK / "receptor.pdbqt"
    sh([OBABEL, str(prot), "-O", str(pdbqt), "-xr", "-p", "7.4"])
    return pdbqt


def box_center() -> tuple[float, float, float]:
    xs = ys = zs = n = 0
    for l in SRC_PDB.read_text().splitlines():
        if l.startswith("HETATM") and l[17:20].strip() == "FK5":
            xs += float(l[30:38]); ys += float(l[38:46]); zs += float(l[46:54]); n += 1
    return xs / n, ys / n, zs / n


def native_fk506_pdbqt() -> Path:
    nat = WORK / "fk506_native.pdb"
    nat.write_text("\n".join(l for l in SRC_PDB.read_text().splitlines()
                             if l.startswith("HETATM") and l[17:20].strip() == "FK5") + "\n")
    out = WORK / "fk506_native.pdbqt"
    sh([OBABEL, str(nat), "-O", str(out), "-p", "7.4"])
    return out


def smiles_to_pdbqt(smiles: str, name: str) -> Path | None:
    out = WORK / f"{name}.pdbqt"
    try:
        sh([OBABEL, f"-:{smiles}", "-O", str(out), "--gen3d", "-p", "7.4"])
    except subprocess.CalledProcessError:
        return None
    return out if out.exists() and out.stat().st_size > 0 else None


def dock(lig: Path, cx, cy, cz, name: str) -> float | None:
    out = WORK / f"{name}_docked.pdbqt"
    try:
        sh([SMINA, "-r", str(WORK / "receptor.pdbqt"), "-l", str(lig),
            "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
            "--size_x", str(BOX), "--size_y", str(BOX), "--size_z", str(BOX),
            "--exhaustiveness", "8", "--seed", "42", "--cpu", "4", "-o", str(out)])
    except subprocess.CalledProcessError:
        return None
    m = re.search(r"minimizedAffinity\s+(-?\d+\.\d+)", out.read_text())
    return float(m.group(1)) if m else None


def main() -> int:
    WORK.mkdir(parents=True, exist_ok=True)
    prep_receptor()
    cx, cy, cz = box_center()
    print(f"Box center (native FK506): ({cx:.1f}, {cy:.1f}, {cz:.1f})  size {BOX} A")

    results = []
    # Positive control: redock native FK506
    fk = native_fk506_pdbqt()
    score = dock(fk, cx, cy, cz, "fk506")
    print(f"  FK506 (redock control): {score} kcal/mol")
    results.append({"name": "FK506", "molecule_chembl_id": "native", "smiles": "",
                    "pchembl_value": None, "dock_score": score})

    # Drug-like anchors
    anchors = pd.read_csv(ROOT / "data" / "processed" / "fkbp12_druglike_anchors.csv").head(N_ANCHORS)
    for _, row in anchors.iterrows():
        cid = row["molecule_chembl_id"]
        lig = smiles_to_pdbqt(row["canonical_smiles"], cid)
        score = dock(lig, cx, cy, cz, cid) if lig else None
        print(f"  {cid} (pKi {row['pchembl_value']}): {score} kcal/mol")
        results.append({"name": cid, "molecule_chembl_id": cid,
                        "smiles": row["canonical_smiles"],
                        "pchembl_value": row["pchembl_value"], "dock_score": score})

    df = pd.DataFrame(results)
    out_csv = WORK / "docking_scores.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nWrote {out_csv}")
    ok = df.dock_score.dropna()
    print(f"Docked {len(ok)}/{len(df)}; mean score {ok.mean():.2f} kcal/mol")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
