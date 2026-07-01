#!/usr/bin/env python
"""GNINA candidate scoring: rank the generated candidates (33 developable + TIER-1)
by CNNaffinity. Same robust batched method as the gate (RDKit 3D + one gnina call).
Assumes gnina binary + rec.pdb/ref.pdb already present at /content/gnina/ from the
gate run (reuse the live session). Rebuilds them if missing.
"""
import subprocess, urllib.request, json, os
from pathlib import Path
WORK=Path("/content/gnina"); WORK.mkdir(parents=True, exist_ok=True)
DATA=json.loads('[["COC(=O)N1CCC[C@@H]1C1=CCC=CC1=O", "developable"], ["O=C1[C@@H](Cc2ccccc2)C[C@H]2C[C@H]2N1CC1CCOCC1", "developable"], ["Cc1ccc(CN2CCC[C@@H](CCCC(=O)O)C2)cc1", "developable"], ["CCSC1=NCC[C@@H](CN(CF)C(=O)NCCc2ccccc2)N1", "developable"], ["NCCCCS[C@H]1[C@H](N2CCCCC2)CCCS1(=O)=O", "developable"], ["CCCOC[C@@H]1CC[C@H](O)C[C@@H]1NCC(=O)C1CCCC1", "developable"], ["CNC(=O)[C@H]1CCO[C@H](C[C@H]2CCCO2)[C@@H]1c1ccc(C#N)cc1", "developable"], ["NC(=O)[C@@H]1CCCCCCc2c(OOF)cccc2O1", "developable"], ["NC(=O)CCO[C@@H]1CC[C@H]2C[C@@H]2CC(=O)C1", "developable"], ["C[C@H]1CC[C@@H]([C@@]2(O)CCC[C@H](N)C2)CC1", "developable"], ["C[C@]1(c2cccc(O)c2)CCOC[C@H](CC#N)N1", "developable"], ["Cc1ccc(Cn2ccnn2)cc1NC=C1CC1", "developable"], ["CCCCCN[C@@H]1[C@@H](c2ccccc2)CCCC[C@@H]1O", "developable"], ["Nc1cc(F)ccc1[C@@H]1CCCC[C@@H](C=CC(=O)O)CO1", "developable"], ["CC[C@@H](CCC(=O)O)N(C(=O)O)[C@@H]1C(=O)[C@H]1C", "developable"], ["C[C@H]1CCC[C@@](Oc2ccccc2)(C(=O)NCC(=O)O)[C@@H]1O", "developable"], ["O=C(O)CCCN1Cc2cc(F)ccc2[C@H](F)C1", "developable"], ["CCCC[C@H]1CC=C[C@H]2CC[C@@H](C(=O)O)C[C@@H]12", "developable"], ["CN1COc2ccc(O)cc21", "developable"], ["C[C@@H]1C(=O)COC(F)(F)[C@H]1c1ccccc1", "developable"], ["C1CCC(C[C@H]2CCCS2)CC1", "developable"], ["Nc1ccc(Cl)c(CC2NOCS2)c1", "developable"], ["CCCCCN1OC(=O)CCC12CC2", "developable"], ["C[C@@H]1SCC[C@@H]1C(=O)O", "developable"], ["CC[C@@H]1C[C@@]1(C)c1ccc(F)cc1N", "developable"], ["C=C(C)N(CCC)C(=O)O", "developable"], ["FC(F)(F)O[C@@H](c1ccccc1)C1CC1", "developable"], ["NCCc1cccc2c1C(=O)SC2=O", "developable"], ["CCC(=O)CC[C@H](C)C(=O)O", "developable"], ["C=C(C)NNS(=O)(=O)c1ccccc1Cc1ncccn1", "developable"], ["O=C(O)Cc1ccccc1", "developable"], ["O=C(O)CCCCN(CCO)c1cccc(F)c1F", "developable"], ["O=C(O)CNc1ccccc1", "developable"], ["CCC[C@@H]1CC(=O)N[C@@]2(C)[C@H](C1)[C@@H](Cc1ccc(F)cc1N)CC[C@@H]2O", "tier1"], ["CC1([C@H](Cc2cccc(Cl)c2)OCCNC(=O)O)CCC1", "tier1"], ["OCCCC[C@H]1CCCC[C@H](CC2c3ccccc32)N1", "tier1"], ["C[C@H]1CCCN(c2ccc(F)cc2)[C@@H]1CCCCC(=O)O", "tier1"]]')
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
