"""G4-follow-up: synthesizability-constrained GA optimization (local, CPU, frugal).

Central finding so far: an affinity<->synthesizability trade-off — the high-GNINA
leads (L1-L3) have novel cores AiZynthFinder cannot route to catalog stock, while the
synthesizable pool caps at ~known-active affinity (and its best member fails MD).

Question this script answers: can directed optimization escape that trade-off — i.e.
how high can predicted potency go while the molecule stays SYNTHESIZABLE-EASY?

Same BRICS-GA machinery and hard developability/domain gates as ga_optimize_r2, plus a
hard SYNTHESIZABILITY gate (from synth_feasibility.analyze): SA <= 3.5 AND
stereocenters <= 2 AND max ring <= 6 AND no strained ring. Objective stays pred_pchembl.
The winners are then AiZynthFinder-checked separately (does a real catalog route exist?).

Output: results/phase2b/ga_synth_optimized.csv
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
from synth_feasibility import analyze  # noqa: E402

PROC = ROOT / "data" / "processed"
N_GEN = 8
POP = 60
OFFSPRING = 400
AD_MIN, AD_MAX, QED_MIN = 0.30, 0.65, 0.55
# synthesizability gate (the new constraint)
SA_MAX, STEREO_MAX, RING_MAX = 3.5, 2, 6


def seed_smiles() -> list[str]:
    s = []
    for p in [ROOT / "results" / "phase2b" / "gaopt_full.csv",
              ROOT / "results" / "phase2b" / "ga_optimized_r2.csv",
              ROOT / "results" / "phase2b" / "gen600_developable.csv",
              ROOT / "results" / "phase2b" / "gen_pool_developable.csv"]:
        if p.exists():
            s += pd.read_csv(p)["smiles"].tolist()
    out, seen = [], set()
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
    return [m for m in (Chem.MolFromSmiles(f) for f in frags) if m is not None]


def anchor_donor_fragments() -> list:
    p = PROC / "fkbp12_druglike_anchors.csv"
    return fragments(pd.read_csv(p)["canonical_smiles"].tolist()) if p.exists() else []


_ANCHOR_FRAGS = None


def breed(pop: list[str], n: int) -> list[str]:
    global _ANCHOR_FRAGS
    if _ANCHOR_FRAGS is None:
        _ANCHOR_FRAGS = anchor_donor_fragments()
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


def synth_ok(smi: str) -> bool:
    a = analyze(smi)
    if not a.get("valid"):
        return False
    return (a["sa_score"] <= SA_MAX and a["n_stereocenters"] <= STEREO_MAX
            and a["max_ring_size"] <= RING_MAX and not a["has_strained_ring"])


def main() -> int:
    scorer = FKBP12Scorer()
    pop = seed_smiles()
    print(f"seeds: {len(pop)}")

    def evaluate(smis):
        df = scorer.score(smis)
        df = df[df.valid]
        ok = df[(df.ad_tanimoto >= AD_MIN) & (df.ad_tanimoto <= AD_MAX) &
                (df.qed >= QED_MIN) & (df.mw >= 250) & (df.mw <= 500)].copy()
        if len(ok):
            ok = ok[ok.smiles.map(synth_ok)]   # the added synthesizability gate
        return ok.sort_values("pred_pchembl", ascending=False)

    cur = evaluate(pop)
    best_hist = []
    for gen in range(N_GEN):
        kids = breed(list(cur.smiles.head(POP)) if len(cur) else pop, OFFSPRING)
        scored_kids = evaluate(kids)
        merged = pd.concat([cur, scored_kids], ignore_index=True).drop_duplicates("smiles")
        cur = merged.sort_values("pred_pchembl", ascending=False).head(POP).reset_index(drop=True)
        best = cur.pred_pchembl.max() if len(cur) else float("nan")
        best_hist.append(best)
        print(f"gen {gen+1}: pop {len(cur)}, best pred_pchembl {best:.3f}, kids kept {len(scored_kids)}")

    out = ROOT / "results" / "phase2b" / "ga_synth_optimized.csv"
    cur.to_csv(out, index=False)
    print(f"\nbest pred_pchembl trajectory: {[round(x,2) for x in best_hist]}")
    print(f"wrote {len(cur)} synthesizable-constrained molecules -> {out}")
    print("\nTop 8 (synthesizable-easy, in-domain, developable):")
    if len(cur):
        print(cur.head(8)[["pred_pchembl", "binder_prob", "qed", "ad_tanimoto", "mw", "smiles"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
