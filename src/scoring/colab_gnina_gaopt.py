#!/usr/bin/env python
"""GNINA candidate scoring: rank the generated candidates (33 developable + TIER-1)
by CNNaffinity. Same robust batched method as the gate (RDKit 3D + one gnina call).
Assumes gnina binary + rec.pdb/ref.pdb already present at /content/gnina/ from the
gate run (reuse the live session). Rebuilds them if missing.
"""
import subprocess, urllib.request, json, os
from pathlib import Path
WORK=Path("/content/gnina"); WORK.mkdir(parents=True, exist_ok=True)
DATA=json.loads('[["C=C[C@H]1CN(CC(=O)O)C(=O)[C@@H]2CCC[C@H]1N2Cc1cccc(CC(N)=O)c1", "ga_opt"], ["C=C[C@@H]1CN(CC(N)=O)C(=O)[C@@H]2CCC[C@H]1N2Cc1cccc(CC(N)=O)c1", "ga_opt"], ["C=C[C@@H]1CN(CC)C(=O)[C@@H]2CCC[C@H]1N2Cc1cccc(CC(N)=O)c1", "ga_opt"], ["C=C[C@H]1CN(CCF)C(=O)[C@@H]2CCC[C@H]1N2Cc1cccc(CC(N)=O)c1", "ga_opt"], ["C=C[C@@H]1CN(CC#N)C(=O)[C@@H]2CCC[C@H]1N2Cc1cccc(CC(N)=O)c1", "ga_opt"], ["C=C[C@H]1CN(c2cccc(NCC)c2)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1ccc(C)cc1", "ga_opt"], ["C=C[C@H]1CN(c2cccc(NCC)c2)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1cccc(C)c1", "ga_opt"], ["C=C[C@H]1CN(c2cccc(NCC)c2)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1ccc(C#N)cc1", "ga_opt"], ["C#CCN1C[C@H](C=C)[C@H]2CCC[C@@H](C1=O)N2Cc1cccc(CC(N)=O)c1", "ga_opt"], ["C=C[C@H]1CN(c2cccc(NCC)c2)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1ccccc1", "ga_opt"], ["CC(C)(C)C(=O)C(=O)N1CCC[C@H]1c1cccc(-c2ccc(Cl)nc2)c1", "ga_opt"], ["CC(C)(C)C(=O)C(=O)N1CCCCC1N1CCCC1N1CCC[C@H]1CC(N)=O", "ga_opt"], ["C=C[C@H]1CN(C)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)Cc1cccc(NCC)c1", "ga_opt"], ["C=C[C@H]1CN(c2cccc(NCC)c2)C(=O)[C@@H]2CCC[C@H]1N2CC(C)CC", "ga_opt"], ["CC(C)(C)C(=O)C(=O)N1CCC[C@H]1N1CCCC1N1CCC[C@H]1CC(N)=O", "ga_opt"]]')
def sh(cmd, check=False, timeout=None):
    r=subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
    if r.stderr and check: print("[err]",r.stderr[-500:],flush=True)
    return r
def main():
    sh("pip -q install rdkit 2>/dev/null || true", timeout=300)  # fresh session needs rdkit
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
