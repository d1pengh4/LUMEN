#!/usr/bin/env python
"""GNINA candidate scoring: rank the generated candidates (33 developable + TIER-1)
by CNNaffinity. Same robust batched method as the gate (RDKit 3D + one gnina call).
Assumes gnina binary + rec.pdb/ref.pdb already present at /content/gnina/ from the
gate run (reuse the live session). Rebuilds them if missing.
"""
import subprocess, urllib.request, json, os
from pathlib import Path
WORK=Path("/content/gnina"); WORK.mkdir(parents=True, exist_ok=True)
DATA=json.loads('[["C=C(S)[C@@H]1[C@@H](CO)CCCN1C(=O)O", "gen600_dev"], ["NCC=C=CCN1CCCC1c1ccccc1", "gen600_dev"], ["CCn1c(C)c(C(=O)O)c(O)c1Sc1ccccc1", "gen600_dev"], ["O=C(O)CCO[C@H]1CCC[C@@H]2CCCC(=O)N2C1", "gen600_dev"], ["OCCCN1CCC[C@@H]2CCC[C@H]21", "gen600_dev"], ["CCCCC1=C[C@@H](C)CC[C@@H](O)C1", "gen600_dev"], ["NC(=O)NC(=O)N1CCC(F)(F)[C@@H](Cc2ccccc2)C1", "gen600_dev"], ["CC(=O)CC[C@@H]1CCCCS[C@@]1(C)O", "gen600_dev"], ["CCCCN1C[C@H](C)[C@@](C[C@H]2CCNC2)(c2ccccc2)C1", "gen600_dev"], ["O=C(O)[C@@H]1CC[C@@]2(CCCC[C@H]2Oc2ccccc2)O1", "gen600_dev"], ["Cc1nc2c(s1)[C@@H](C(=O)O)CC=C2", "gen600_dev"], ["CC[C@@H]1NC([C@H]2CN2)=N[C@H]1c1ccccc1", "gen600_dev"], ["O[C@@H](CCC1(c2cc(Cl)ccc2F)CC1)C(O)(F)F", "gen600_dev"], ["OC1=Cc2ccccc2CC=C1", "gen600_dev"], ["CCNC(=O)CCOC1=CC=CCC(C(=O)NCC)=C1", "gen600_dev"], ["CCCc1nc2c(s1)N1C[C@]1(F)c1ccc(F)cc1OC2", "gen600_dev"], ["CC[C@H](C)NCC[C@@H]1[C@H](O)c2c(OCCOC)cccc2[C@@H]1N", "gen600_dev"], ["O[C@H]1CCCCCc2ncccc21", "gen600_dev"], ["N[C@H](CSc1cccs1)C(=O)O", "gen600_dev"], ["CC(=O)[C@H]1CCN(C(=O)[C@@]2(C)CCc3ccccc32)C1", "gen600_dev"], ["OS(O)(CCF)OCc1ccccc1", "gen600_dev"], ["NCCC(F)(F)Cc1cnccn1", "gen600_dev"], ["C[C@@H](CP(=O)(O)O)Nc1ncccc1C(F)(F)F", "gen600_dev"], ["N[C@]1(CCO)SCc2cc(F)c(F)cc21", "gen600_dev"]]')
def sh(cmd, check=False, timeout=None):
    r=subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
    if r.stderr and check: print("[err]",r.stderr[-500:],flush=True)
    return r
def main():
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")
    if not (WORK/"rec.pdb").exists():
        raw=urllib.request.urlopen("https://files.rcsb.org/download/1FKF.pdb",timeout=60).read().decode()
        (WORK/"rec.pdb").write_text("\n".join(l for l in raw.splitlines() if l.startswith("ATOM"))+"\n")
        (WORK/"ref.pdb").write_text("\n".join(l for l in raw.splitlines() if l.startswith("HETATM") and l[17:20].strip()=="FK5")+"\n")
    if not (WORK/"gnina").exists():
        sh("wget -q https://github.com/gnina/gnina/releases/download/v1.3/gnina -O /content/gnina/gnina && chmod +x /content/gnina/gnina")
    print("STEP embedding candidates", flush=True)
    w=Chem.SDWriter(str(WORK/"cand.sdf")); order=[]
    for i,(smi,src) in enumerate(DATA):
        m=Chem.MolFromSmiles(smi)
        if m is None: continue
        m=Chem.AddHs(m)
        if AllChem.EmbedMolecule(m, randomSeed=42, maxAttempts=50, useRandomCoords=True)!=0: continue
        try: AllChem.MMFFOptimizeMolecule(m, maxIters=200)
        except Exception: pass
        m.SetProp("_Name", f"c{i}"); w.write(m); order.append((i,smi,src))
    w.close()
    print(f"embedded {len(order)}/{len(DATA)}", flush=True)
    out=WORK/"cand_docked.sdf"
    print("STEP gnina docking candidates", flush=True)
    sh(f"timeout 3000 /content/gnina/gnina -r {WORK}/rec.pdb -l {WORK}/cand.sdf --autobox_ligand {WORK}/ref.pdb "
       f"-o {out} --cnn_scoring rescore --exhaustiveness 4 --num_modes 1 --seed 42 2>&1 | tail -3")
    print("STEP done, parsing", flush=True)
    best={}
    if out.exists():
        for mol in Chem.SDMolSupplier(str(out)):
            if mol is None: continue
            nm=mol.GetProp("_Name"); a=float(mol.GetProp("CNNaffinity")) if mol.HasProp("CNNaffinity") else None
            if a is not None and (nm not in best or a>best[nm]): best[nm]=a
    print("\n===CAND_BEGIN===", flush=True)
    for i,smi,src in order:
        print(f"{best.get(f'c{i}')}\t{src}\t{smi}", flush=True)
    print("===CAND_END===", flush=True)
if __name__=="__main__": raise SystemExit(main())
