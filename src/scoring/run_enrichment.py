"""Dock the enrichment set and quantify discriminative power of the docking score.

Metrics:
  ROC-AUC  = P(a random active scores better than a random decoy). 0.5 = random.
  EF@k%    = enrichment factor in the top k% by score.
Output: results/phase1/docking/enrichment_results.csv (+ printed metrics)
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DOCK = Path("/usr/local/Caskroom/miniforge/base/envs/lumen-dock/bin")
OBABEL, SMINA = str(DOCK / "obabel"), str(DOCK / "smina")
WORK = ROOT / "results" / "phase1" / "docking"
LIGDIR = WORK / "enrichment_ligs"
SRC_PDB = ROOT / "data" / "raw" / "structures" / "1FKF.pdb"
BOX = 24.0


def sh(cmd):
    return subprocess.run(cmd, capture_output=True, text=True).returncode == 0


def box_center():
    xs = ys = zs = n = 0
    for l in SRC_PDB.read_text().splitlines():
        if l.startswith("HETATM") and l[17:20].strip() == "FK5":
            xs += float(l[30:38]); ys += float(l[38:46]); zs += float(l[46:54]); n += 1
    return xs / n, ys / n, zs / n


def dock_smiles(smiles, name, cx, cy, cz):
    lig = LIGDIR / f"{name}.pdbqt"
    if not sh([OBABEL, f"-:{smiles}", "-O", str(lig), "--gen3d", "-p", "7.4"]):
        return None
    out = LIGDIR / f"{name}_docked.pdbqt"
    ok = sh([SMINA, "-r", str(WORK / "receptor.pdbqt"), "-l", str(lig),
             "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
             "--size_x", str(BOX), "--size_y", str(BOX), "--size_z", str(BOX),
             "--exhaustiveness", "8", "--seed", "42", "--cpu", "4", "-o", str(out)])
    if not ok or not out.exists():
        return None
    m = re.search(r"minimizedAffinity\s+(-?\d+\.\d+)", out.read_text())
    return float(m.group(1)) if m else None


def roc_auc(labels, scores):
    # higher 'goodness' = more negative score -> use -score
    g = -np.asarray(scores)
    y = np.asarray(labels)
    pos, neg = g[y == 1], g[y == 0]
    # Mann-Whitney U / (n_pos*n_neg)
    wins = sum((p > neg).sum() + 0.5 * (p == neg).sum() for p in pos)
    return wins / (len(pos) * len(neg))


def ef(labels, scores, k):
    df = pd.DataFrame({"y": labels, "s": scores}).sort_values("s")  # most negative first
    n_top = max(1, int(len(df) * k))
    top = df.head(n_top)
    return (top.y.mean()) / (df.y.mean())


def main():
    LIGDIR.mkdir(parents=True, exist_ok=True)
    cx, cy, cz = box_center()
    df = pd.read_csv(ROOT / "data" / "processed" / "enrichment_set.csv")
    scores = []
    for i, row in df.iterrows():
        s = dock_smiles(row["smiles"], row["name"], cx, cy, cz)
        scores.append(s)
        print(f"[{i+1}/{len(df)}] {row['name']} label={row['label']} score={s}")
    df["dock_score"] = scores
    df = df.dropna(subset=["dock_score"])
    df.to_csv(WORK / "enrichment_results.csv", index=False)

    auc = roc_auc(df.label.values, df.dock_score.values)
    print("\n================ ENRICHMENT ================")
    print(f"docked: {len(df)} ({(df.label==1).sum()} actives, {(df.label==0).sum()} decoys)")
    print(f"ROC-AUC : {auc:.3f}   (0.5 = random, >0.7 useful)")
    print(f"EF@10%  : {ef(df.label.values, df.dock_score.values, 0.10):.2f}")
    print(f"EF@20%  : {ef(df.label.values, df.dock_score.values, 0.20):.2f}")
    print(f"active mean score: {df[df.label==1].dock_score.mean():.2f}  "
          f"decoy mean score: {df[df.label==0].dock_score.mean():.2f}")


if __name__ == "__main__":
    main()
