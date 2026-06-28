"""Structure-based validation of DiffSBDD (B) molecules.

The ligand-ML scorer cannot judge B's molecules (all out of applicability domain,
Tanimoto < 0.3 to known binders). So we validate them the right way: dock them
into the FKBP12 pocket and compare the score distribution to known binders.
Docking SCORE ranks weakly (enrichment 0.455) but a pose + a score in the range
of real binders is genuine structural evidence the novel scaffold can occupy the
pocket. Reuses the receptor/box from src/scoring/dock_benchmark.py.

Output: results/phase2b/diffsbdd_docked.csv
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DOCK = Path("/usr/local/Caskroom/miniforge/base/envs/lumen-dock/bin")
OBABEL, SMINA = str(DOCK / "obabel"), str(DOCK / "smina")
DOCKDIR = ROOT / "results" / "phase1" / "docking"
RECEPTOR = DOCKDIR / "receptor.pdbqt"
SRC_PDB = ROOT / "data" / "raw" / "structures" / "1FKF.pdb"
WORK = ROOT / "results" / "phase2b"
LIGDIR = WORK / "dock_ligs"
BOX = 24.0


def sh(cmd):
    return subprocess.run(cmd, capture_output=True, text=True).returncode == 0


def box_center():
    xs = ys = zs = n = 0
    for l in SRC_PDB.read_text().splitlines():
        if l.startswith("HETATM") and l[17:20].strip() == "FK5":
            xs += float(l[30:38]); ys += float(l[38:46]); zs += float(l[46:54]); n += 1
    return xs / n, ys / n, zs / n


def dock(smiles, name, cx, cy, cz):
    lig = LIGDIR / f"{name}.pdbqt"
    if not sh([OBABEL, f"-:{smiles}", "-O", str(lig), "--gen3d", "-p", "7.4"]):
        return None
    out = LIGDIR / f"{name}_d.pdbqt"
    if not sh([SMINA, "-r", str(RECEPTOR), "-l", str(lig),
               "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
               "--size_x", str(BOX), "--size_y", str(BOX), "--size_z", str(BOX),
               "--exhaustiveness", "8", "--seed", "42", "--cpu", "4", "-o", str(out)]):
        return None
    m = re.search(r"minimizedAffinity\s+(-?\d+\.\d+)", out.read_text())
    return float(m.group(1)) if m else None


def main() -> int:
    LIGDIR.mkdir(parents=True, exist_ok=True)
    cx, cy, cz = box_center()
    df = pd.read_csv(WORK / "diffsbdd_scored.csv")
    scores = []
    for i, row in df.iterrows():
        s = dock(row["smiles"], f"b{i}", cx, cy, cz)
        scores.append(s)
        print(f"[{i+1}/{len(df)}] b{i}: {s}", flush=True)
    df["dock_score"] = scores
    df.to_csv(WORK / "diffsbdd_docked.csv", index=False)
    ok = df.dock_score.dropna()
    print(f"\ndocked {len(ok)}/{len(df)}")
    print(f"B dock score: min {ok.min():.2f}, median {ok.median():.2f}, max {ok.max():.2f}")
    print("reference: FK506 -10.3 ; known drug-like binders -6.7..-7.6")
    print(f"B molecules scoring <= -7.0 (binder-like): {(ok<=-7.0).sum()}/{len(ok)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
