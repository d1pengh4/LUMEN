"""G2: goal-directed evolutionary optimization (local, CPU, frugal).

Rationale: brute-force generation plateaued (~6.5 CNNaffinity, see scaleup log).
RL (REINVENT4) with a GNINA reward is infeasible on free GPU (thousands of dockings).
Instead: a BRICS fragment-recombination genetic algorithm that OPTIMIZES the fast
validated ML scorer, keeping molecules developable and IN the applicability domain
(hard constraints, to resist reward-hacking). The optimized top molecules are then
GNINA-validated separately (the expensive but validated check).

Objective (maximize): pred_pchembl (potency regressor), subject to hard gates
  in_domain (ad_tanimoto >= 0.30) AND QED >= 0.55 AND 250 <= MW <= 500 AND valid.
Seeds: our developable generated candidates + drug-like ChEMBL binders (fragment donors).

Output: results/phase2b/ga_optimized.csv (top optimized molecules, ML-scored)
"""
from __future__ import annotations

import random
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import BRICS

RDLogger.DisableLog("rdApp.*")
random.seed(42)
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "scoring"))
from score_molecules import FKBP12Scorer  # noqa: E402

PROC = ROOT / "data" / "processed"
N_GEN = 8
POP = 60
OFFSPRING = 400
# in-domain but NOT a near-duplicate of a known binder (prevents seed-dominance /
# regurgitating training molecules — see attempt 1 in the G2 log).
AD_MIN, AD_MAX, QED_MIN = 0.30, 0.65, 0.55


def seed_smiles() -> list[str]:
    # Seed the POPULATION only with NOVEL generated candidates (NOT anchors).
    s = []
    for p in [ROOT / "results" / "phase2b" / "gen600_developable.csv",
              ROOT / "results" / "phase2b" / "gen_pool_developable.csv"]:
        if p.exists():
            s += pd.read_csv(p)["smiles"].tolist()
    out = []
    seen = set()
    for smi in s:
        m = Chem.MolFromSmiles(smi)
        if m:
            c = Chem.MolToSmiles(m)
            if c not in seen:
                seen.add(c); out.append(c)
    return out


def fragments(pop: list[str]) -> list:
    frags = set()
    for smi in pop:
        m = Chem.MolFromSmiles(smi)
        if m:
            try:
                frags |= set(BRICS.BRICSDecompose(m))
            except Exception:
                pass
    mols = [Chem.MolFromSmiles(f) for f in frags]
    return [m for m in mols if m is not None]


def anchor_donor_fragments() -> list:
    p = PROC / "fkbp12_druglike_anchors.csv"
    if not p.exists():
        return []
    return fragments(pd.read_csv(p)["canonical_smiles"].tolist())


_ANCHOR_FRAGS = None


def breed(pop: list[str], n: int) -> list[str]:
    global _ANCHOR_FRAGS
    if _ANCHOR_FRAGS is None:
        _ANCHOR_FRAGS = anchor_donor_fragments()  # anchors contribute FRAGMENTS only
    frags = fragments(pop) + _ANCHOR_FRAGS
    if len(frags) < 2:
        return []
    random.seed()
    builder = BRICS.BRICSBuild(frags, onlyCompleteMols=True, scrambleReagents=True)
    kids, seen = [], set()
    import itertools
    for m in itertools.islice(builder, n * 8):
        if m is None:
            continue
        try:
            Chem.SanitizeMol(m)
        except Exception:
            continue
        smi = Chem.MolToSmiles(m)
        if smi not in seen:
            seen.add(smi); kids.append(smi)
        if len(kids) >= n:
            break
    return kids


def main() -> int:
    scorer = FKBP12Scorer()
    pop = seed_smiles()
    print(f"seeds: {len(pop)}")

    def evaluate(smis):
        df = scorer.score(smis)
        df = df[df.valid]
        # hard constraints, then objective = pred_pchembl
        ok = df[(df.ad_tanimoto >= AD_MIN) & (df.ad_tanimoto <= AD_MAX) &
                (df.qed >= QED_MIN) & (df.mw >= 250) & (df.mw <= 500)].copy()
        return ok.sort_values("pred_pchembl", ascending=False)

    cur = evaluate(pop)
    best_hist = []
    all_seen = {r.smiles: r for r in cur.itertuples()}
    for gen in range(N_GEN):
        kids = breed(list(cur.smiles.head(POP)) if len(cur) else pop, OFFSPRING)
        scored_kids = evaluate(kids)
        merged = pd.concat([cur, scored_kids], ignore_index=True).drop_duplicates("smiles")
        cur = merged.sort_values("pred_pchembl", ascending=False).head(POP).reset_index(drop=True)
        best = cur.pred_pchembl.max()
        best_hist.append(best)
        print(f"gen {gen+1}: pop {len(cur)}, best pred_pchembl {best:.3f}, "
              f"kids kept {len(scored_kids)}")

    cur.to_csv(ROOT / "results" / "phase2b" / "ga_optimized.csv", index=False)
    print(f"\nbest pred_pchembl trajectory: {[round(x,2) for x in best_hist]}")
    print(f"wrote {len(cur)} optimized molecules -> results/phase2b/ga_optimized.csv")
    print("\nTop 8 (ML-optimized, in-domain, developable):")
    print(cur.head(8)[["pred_pchembl", "binder_prob", "qed", "ad_tanimoto", "mw", "smiles"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
