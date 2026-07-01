"""Smoke test: does GNINA run on a Colab GPU and separate known FKBP12 binders
from an obvious non-binder? Prerequisite before the full enrichment GATE.

Run: src/generation/colabx.sh run --gpu T4 --timeout 1800 src/scoring/colab_gnina_smoke.py
Prints CNNaffinity per molecule between markers. Two known drug-like FKBP12
actives should score higher (more favorable) than aspirin.
"""
import subprocess
import sys
import urllib.request
from pathlib import Path

WORK = Path("/content/gnina")
WORK.mkdir(parents=True, exist_ok=True)

# (name, smiles, expected)
MOLS = [
    ("active1", "COC[C@H]1CN(Cc2ccccn2)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1cc(Cl)cc(Cl)c1", "binder"),
    ("active2", "CCC(C)(C)C(=O)C(=O)N1CCC[C@H]1C(=O)NCCCc1ccccc1", "binder"),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O", "non-binder"),
]


def sh(cmd, check=True):
    print(f"$ {cmd}", flush=True)
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.stdout:
        print(r.stdout[-2000:], flush=True)
    if r.stderr:
        print("[err]", r.stderr[-1500:], flush=True)
    if check and r.returncode:
        raise RuntimeError(f"failed rc={r.returncode}: {cmd}")
    return r


def main() -> int:
    sh("nvidia-smi -L")
    # receptor + reference ligand from 1FKF
    raw = urllib.request.urlopen("https://files.rcsb.org/download/1FKF.pdb", timeout=60).read().decode()
    (WORK / "rec.pdb").write_text(
        "\n".join(l for l in raw.splitlines() if l.startswith("ATOM")) + "\n")
    (WORK / "ref.pdb").write_text(
        "\n".join(l for l in raw.splitlines()
                  if l.startswith("HETATM") and l[17:20].strip() == "FK5") + "\n")

    # tools: openbabel CLI for 3D, gnina static binary
    sh("apt-get -qq update >/dev/null 2>&1 || true", check=False)
    sh("apt-get -qq install -y openbabel >/dev/null 2>&1 || true", check=False)
    sh("wget -q https://github.com/gnina/gnina/releases/download/v1.3/gnina "
       "-O /content/gnina/gnina && chmod +x /content/gnina/gnina")
    sh("/content/gnina/gnina --version 2>&1 | head -3", check=False)

    results = []
    for name, smi, exp in MOLS:
        smi_f = WORK / f"{name}.smi"
        smi_f.write_text(f"{smi} {name}\n")
        sdf = WORK / f"{name}.sdf"
        sh(f"obabel {smi_f} -O {sdf} --gen3d 2>/dev/null", check=False)
        out = WORK / f"{name}_docked.sdf"
        r = sh(f"/content/gnina/gnina -r {WORK}/rec.pdb -l {sdf} "
               f"--autobox_ligand {WORK}/ref.pdb -o {out} "
               f"--cnn_scoring rescore --num_modes 3 --seed 42 2>&1 | tail -20", check=False)
        # parse CNNaffinity from output sdf tags or stdout
        aff = None
        if out.exists():
            txt = out.read_text()
            import re
            m = re.findall(r"CNNaffinity[>\s]*\n?\s*(-?\d+\.\d+)", txt)
            if m:
                aff = max(float(x) for x in m)  # higher CNNaffinity = better
        results.append((name, exp, aff))

    print("\n===GNINA_SMOKE_BEGIN===", flush=True)
    for name, exp, aff in results:
        print(f"{name}\t{exp}\tCNNaffinity={aff}", flush=True)
    print("===GNINA_SMOKE_END===", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
