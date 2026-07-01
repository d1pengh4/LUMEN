#!/usr/bin/env python
"""GNINA enrichment GATE: score the 60-mol enrichment set (actives vs decoys)
on a Colab GPU, so we can compute ROC-AUC and decide if GNINA beats docking
(0.455). Data embedded for a self-contained ephemeral VM run.
"""
import subprocess, sys, urllib.request, re, json
from pathlib import Path
WORK=Path("/content/gnina"); WORK.mkdir(parents=True, exist_ok=True)
DATA=json.loads('[["COC[C@H]1CN(Cc2ccccn2)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1cc(Cl)cc(Cl)c1", 1], ["CCC(C)(C)C(=O)C(=O)N1CCC[C@H]1C(=O)NCCCc1ccccc1", 1], ["C=C[C@H]1CN(CC(=O)O)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1cc(Cl)cc(Cl)c1", 1], ["C=C[C@H]1CN(CCNC)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1cc(Cl)cc(Cl)c1", 1], ["C=C[C@@H]1CN(CCN(C)C)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1cc(Cl)cc(Cl)c1", 1], ["C=C[C@H]1CN(CC(N)=O)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1cc(Cl)cc(Cl)c1", 1], ["CCC(C)(C)C(=O)C(=O)N1CCC[C@H]1C(=O)SCCCc1cccnc1", 1], ["CCC(C)(C)C(=O)C(=O)N1CCC[C@H]1C(=O)SCCCc1ccccc1", 1], ["C=C[C@H]1CN(Cc2ncccc2F)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1cc(Cl)cc(Cl)c1", 1], ["CCC(C)(C)C(=O)C(=O)N1CCC[C@H]1C(=O)SCCc1ccccc1", 1], ["C=C[C@H]1CN(CCOC)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1cc(Cl)cc(Cl)c1", 1], ["COc1cc2c(c(OC)c1)C1C3CCCC(C(=O)N1CC2)N3C(=O)C(=O)c1cccc(Cl)c1", 1], ["C=C[C@@H]1CN(CCF)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1cc(Cl)cc(Cl)c1", 1], ["C=C[C@H]1CN(CC#N)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1cc(Cl)cc(Cl)c1", 1], ["C#CCN1C[C@H](C=C)[C@H]2CCC[C@@H](C1=O)N2S(=O)(=O)c1cc(Cl)cc(Cl)c1", 1], ["C=C[C@H]1CN(C)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1cc(Cl)cc(Cl)c1", 1], ["C=C[C@H]1CNC(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1cc(Cl)cc(Cl)c1", 1], ["CCC(C)CNC(=O)N1CCCCC1C(=O)OCCCc1cccnc1", 1], ["O=C(OCCCc1cccnc1)C1CCCN1S(=O)(=O)Cc1ccccc1", 1], ["COc1cc2c(c(OC)c1)C1C3CCCC(C(=O)N1CC2)N3C(=O)C(=O)c1ccc(Cl)cc1", 1], ["CCCCCN1CCc2c(nc(CC(C)(C)C)n2CC2CC2)C1", 0], ["CS(=O)(=O)c1ccc(-c2cc(C(F)(F)F)nc(NC3CCCCC3)n2)cc1", 0], ["COc1ccccc1/C=N/c1nnc(Cn2c(C)nc3ccccc32)s1", 0], ["Cc1nc2ccccc2n1Cc1nnc(/N=C/c2ccc(O)cc2)s1", 0], ["COc1ccc(/C=N/c2nnc(Cn3c(C)nc4ccccc43)s2)cc1", 0], ["COc1ccc(-c2nnc(CCC(=O)c3ccc(C)cc3)o2)cc1", 0], ["Cc1nc(S)n(Nc2cccc(Br)c2)c1-c1ccccc1", 0], ["Cc1ccc(C2CC(=NOCc3ccccc3)CC(c3ccc(C)cc3)N2)cc1", 0], ["COc1ccc(C2CC(=NOCc3ccccc3)CC(c3ccc(OC)cc3)N2)cc1", 0], ["O=S(=O)(c1ccccc1)N1N=C2CCCCC2C1c1ccc(Cl)cc1", 0], ["CCc1cc2c(cc1OC1CCOC1)Cc1c(-c3cccnc3)n[nH]c1-2", 0], ["Nc1nc2c(ncn2/C=C2\\\\CC23CCP(=O)(O)OC3)c(=O)[nH]1", 0], ["COc1cc2nc(N3CCN(C(=O)[C@@H]4Cc5ccccc5CN4)CC3)nc(N)c2cc1OC", 0], ["COc1cc(/C=C/c2n[nH]c3c2Cc2ccccc2-3)ccc1O", 0], ["Cn1c2c(c(=O)n(C)c1=O)C(c1cccc(C(F)(F)F)c1F)C1=C(N2)c2ccccc2C1=O", 0], ["O=C(O)C(=O)N(c1cc(Cl)ccc1C(=O)O)c1cccc2ccc(O)cc12", 0], ["Oc1ccc(C/N=N/c2ccc3nnc(C(F)(F)F)n3n2)c(O)c1", 0], ["Cc1cccc2c1c(=O)c1ccccc1n2CC(O)CN1CCCCC1", 0], ["O=C(O)C(=O)N(c1ccc2c(c1)OCCO2)c1ccccc1C(=O)O", 0], ["O=C(O)c1cccc2c(=O)c3ccccc3n(CC(O)CN3CCOCC3)c12", 0], ["COC(=O)c1c(O)cccc1OC[C@H]1C[C@H]1c1cccc(-c2cc(C(=O)O)no2)c1", 0], ["CCCC(=O)Nc1cc(-c2cc(C(=O)O)no2)ccc1NC(C)=O", 0], ["Cn1cnc(C(=O)N(Cc2cccc(OC(F)(F)F)c2)C[C@H]2[C@@H]3CNC[C@@H]32)c1", 0], ["Cn1cnc(C(=O)N(Cc2ccc(F)c(Cl)c2)C[C@H]2[C@@H]3CNC[C@@H]32)c1", 0], ["CC(=O)N1CCC(C(=O)Nc2cccc(-c3cc(C(=O)O)no3)c2)CC1", 0], ["CC1=C(C(=O)O)N2C(=O)[C@@H](NC(=O)[C@H](N)c3ccccc3)[C@H]2SC1", 0], ["C=CCn1c(=O)ccc2cc(F)c(-n3c(=O)cc(C(F)(F)F)n(C)c3=O)cc21", 0], ["COc1ccc(Nc2cc(C)c3c(C)c(C)sc3c2C)cc1", 0], ["C#CCOc1ccc(CCNC(=O)C(OCCF)c2ccc(Cl)cc2)cc1OC", 0], ["COc1ccc(Nc2cc(C)c3c(C)c(C)sc3c2C)c(Br)c1", 0], ["C#CCOc1ccc(CCNC(=O)C(OC(C)(C)C#C)c2ccc(Cl)cc2)cc1OC", 0], ["N=C(N)c1ccc2cc(CCc3ccc4c(c3)CCNC4)ccc2c1", 0], ["Cc1cccc(OC(C)(C)C(=O)N/N=C/c2ccc([N+](=O)[O-])cc2)c1", 0], ["C=CCN(Cc1cccc(Cl)c1)CC(O)(Cn1cncn1)c1ccc(F)cc1F", 0], ["COC(=O)c1sc2ccccc2c1C#Cc1ccc(F)cc1", 0], ["COc1cc(-c2cc(N3CCCCC3)nc(N)n2)cc(OC)c1OC", 0], ["COC(=O)c1sc2cc(OC)ccc2c1C#Cc1ccc(F)cc1", 0], ["CC(C)=C1CCC2=C(CC[C@H]3[C@](C)(C(=O)O)CCC[C@]23C)C1", 0], ["O=C(N/N=C/c1ccco1)c1cc([N+](=O)[O-])cc([N+](=O)[O-])c1", 0], ["O=C(N[C@@H](Cc1ccc(O)cc1)C(=O)O)OCC1c2ccccc2-c2ccccc21", 0]]')

def sh(cmd, check=True):
    r=subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.stderr and check: print("[err]", r.stderr[-800:], flush=True)
    if check and r.returncode: raise RuntimeError(f"rc={r.returncode}: {cmd}")
    return r

def main():
    sh("nvidia-smi -L", check=False)
    raw=urllib.request.urlopen("https://files.rcsb.org/download/1FKF.pdb",timeout=60).read().decode()
    (WORK/"rec.pdb").write_text("\n".join(l for l in raw.splitlines() if l.startswith("ATOM"))+"\n")
    (WORK/"ref.pdb").write_text("\n".join(l for l in raw.splitlines() if l.startswith("HETATM") and l[17:20].strip()=="FK5")+"\n")
    sh("apt-get -qq install -y openbabel >/dev/null 2>&1 || true", check=False)
    sh("wget -q https://github.com/gnina/gnina/releases/download/v1.3/gnina -O /content/gnina/gnina && chmod +x /content/gnina/gnina")
    results=[]
    for i,(smi,lab) in enumerate(DATA):
        f=WORK/f"m{i}.smi"; f.write_text(f"{smi} m{i}\n")
        sdf=WORK/f"m{i}.sdf"
        if sh(f"obabel {f} -O {sdf} --gen3d 2>/dev/null", check=False).returncode or not sdf.exists():
            results.append((i,lab,None)); continue
        out=WORK/f"m{i}_d.sdf"
        sh(f"/content/gnina/gnina -r {WORK}/rec.pdb -l {sdf} --autobox_ligand {WORK}/ref.pdb "
           f"-o {out} --cnn_scoring rescore --num_modes 3 --seed 42 >/dev/null 2>&1", check=False)
        aff=None
        if out.exists():
            m=re.findall(r"CNNaffinity[>\s]*\n?\s*(-?\d+\.\d+)", out.read_text())
            if m: aff=max(float(x) for x in m)
        results.append((i,lab,aff))
        print(f"[{i+1}/{len(DATA)}] label={lab} CNNaffinity={aff}", flush=True)
    print("\n===GNINA_GATE_BEGIN===", flush=True)
    for i,lab,aff in results: print(f"{i}\t{lab}\t{aff}", flush=True)
    print("===GNINA_GATE_END===", flush=True)

if __name__=="__main__":
    raise SystemExit(main())
