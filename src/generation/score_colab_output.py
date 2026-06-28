"""Parse SMILES from a saved `colab run` log and score them locally with the
validated FKBP12Scorer + the same physchem/diversity audit applied to Phase 2-A,
so DiffSBDD (B) can be compared head-to-head against the BRICS baseline (A).

Usage:
    colab run ... colab_diffsbdd_run.py 200 | tee results/phase2b/colab.log
    python src/generation/score_colab_output.py results/phase2b/colab.log
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "scoring"))
from score_molecules import FKBP12Scorer  # noqa: E402

RDLogger.DisableLog("rdApp.*")
OUT = ROOT / "results" / "phase2b"
_FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def extract(log: Path) -> list[str]:
    lines, on = [], False
    for l in log.read_text().splitlines():
        if "===LUMEN_SMILES_BEGIN===" in l:
            on = True; continue
        if "===LUMEN_SMILES_END===" in l:
            break
        if on and l.strip():
            lines.append(l.strip())
    return lines


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    log = Path(sys.argv[1])
    smis = extract(log)
    smis = [Chem.MolToSmiles(m) for s in smis if (m := Chem.MolFromSmiles(s))]
    smis = list(dict.fromkeys(smis))  # dedupe, keep order
    if not smis:
        print("No SMILES found between markers. Check the colab log.")
        return 1
    print(f"parsed {len(smis)} unique valid molecules from {log}")

    scored = FKBP12Scorer().score(smis)
    scored = scored[scored.valid].sort_values(
        ["in_domain", "binder_prob", "pred_pchembl"], ascending=False).reset_index(drop=True)
    scored.to_csv(OUT / "diffsbdd_scored.csv", index=False)

    # same audit as A
    mols = [Chem.MolFromSmiles(s) for s in scored.smiles]
    scafs = {MurckoScaffold.MurckoScaffoldSmiles(mol=m) for m in mols}
    known = [_FP.GetFingerprint(Chem.MolFromSmiles(s))
             for s in pd.read_csv(ROOT / "data/processed/fkbp12_ligands.csv").canonical_smiles
             if Chem.MolFromSmiles(s)]
    nov = np.median([max(DataStructs.BulkTanimotoSimilarity(_FP.GetFingerprint(m), known)) for m in mols])
    logp = np.median([Descriptors.MolLogP(m) for m in mols])

    print(f"\n=== DiffSBDD (B) vs BRICS (A) audit ===")
    print(f"  scaffold diversity: {len(scafs)}/{len(mols)} = {len(scafs)/len(mols)*100:.0f}%  (A was 11%)")
    print(f"  novelty (median max Tanimoto to known): {nov:.2f}  (A was 0.46)")
    print(f"  median logP: {logp:.1f}  (A was 5.3)")
    print(f"  in-domain & binder_prob>=0.5: {len(scored[(scored.in_domain)&(scored.binder_prob>=0.5)])}")
    print(f"\nWrote {OUT / 'diffsbdd_scored.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
