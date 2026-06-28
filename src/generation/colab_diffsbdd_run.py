#!/usr/bin/env -S colab run --gpu L4 --timeout 2400
"""Self-contained DiffSBDD generation job for a Colab GPU VM (LUMEN Phase 2-B).

Run from the LUMEN repo on this Mac:
    colab auth                       # ONE TIME, interactive Google login (user does this)
    colab run --gpu L4 --timeout 2400 src/generation/colab_diffsbdd_run.py 200

The script is fully self-contained so the ephemeral VM needs no uploads:
  1. fetch 1FKF from RCSB, build the FKBP12 pocket (chain A + FK506 A:108)
  2. clone DiffSBDD, best-effort install deps, download the full-atom conditional
     checkpoint from Zenodo (record 8183747)
  3. generate N pocket-conditioned molecules
  4. print the SMILES to stdout (captured by `colab run`) and save an SDF

Honest risks (smoke-test first):
  - DiffSBDD pins old torch/CUDA; we install torch matching the VM's CUDA and
    compatible torch-scatter, which MAY need adjustment. If imports fail, the
    error is printed and the job exits non-zero — that's expected on run #1.
SMILES between the markers below are the deliverable to bring back and score
locally with src/scoring/score_molecules.py.
"""
import subprocess
import sys
import urllib.request
from pathlib import Path

N = int(sys.argv[1]) if len(sys.argv) > 1 else 100
WORK = Path("/content/lumen") if Path("/content").exists() else Path.cwd() / "lumen_run"
WORK.mkdir(parents=True, exist_ok=True)


def sh(cmd, **kw):
    print(f"$ {cmd}", flush=True)
    return subprocess.run(cmd, shell=True, check=True, **kw)


def build_pocket() -> Path:
    raw = urllib.request.urlopen("https://files.rcsb.org/download/1FKF.pdb", timeout=60).read().decode()
    keep = [l for l in raw.splitlines()
            if l.startswith("ATOM") or (l.startswith("HETATM") and l[17:20].strip() == "FK5")]
    pdb = WORK / "fkbp12_fk506.pdb"
    pdb.write_text("\n".join(keep) + "\nEND\n")
    print(f"pocket pdb: {sum(l.startswith('ATOM') for l in keep)} protein + "
          f"{sum(l.startswith('HETATM') for l in keep)} FK506 atoms", flush=True)
    return pdb


def setup():
    sh("nvidia-smi -L")
    repo = WORK / "DiffSBDD"
    if not repo.exists():
        sh(f"git clone --depth 1 https://github.com/arneschneuing/DiffSBDD {repo}")
    # best-effort deps matched to the VM's torch/CUDA (not DiffSBDD's pinned old ones)
    sh(f'{sys.executable} -m pip -q install "rdkit>=2022.03" biopython scipy networkx')
    sh(f'{sys.executable} -m pip -q install torch-scatter '
       f'-f https://data.pyg.org/whl/torch-$(python -c "import torch;print(torch.__version__.split(\'+\')[0])")'
       f'+cu$(python -c "import torch;print(torch.version.cuda.replace(\'.\',\'\'))").html || true')
    ckpt = repo / "crossdocked_fullatom_cond.ckpt"
    if not ckpt.exists():
        sh(f'wget -q -O {ckpt} '
           f'"https://zenodo.org/record/8183747/files/crossdocked_fullatom_cond.ckpt?download=1"')
    return repo, ckpt


def generate(repo: Path, ckpt: Path, pocket: Path) -> Path:
    out = WORK / "fkbp12_generated.sdf"
    sh(f'cd {repo} && {sys.executable} generate_ligands.py {ckpt} '
       f'--pdbfile {pocket} --ref_ligand A:108 --outfile {out} '
       f'--n_samples {N} --sanitize')
    return out


def emit_smiles(sdf: Path):
    from rdkit import Chem
    smis = [Chem.MolToSmiles(m) for m in Chem.SDMolSupplier(str(sdf)) if m is not None]
    print(f"\n===LUMEN_SMILES_BEGIN=== ({len(smis)} molecules)", flush=True)
    for s in smis:
        print(s, flush=True)
    print("===LUMEN_SMILES_END===", flush=True)


def main() -> int:
    pocket = build_pocket()
    repo, ckpt = setup()
    sdf = generate(repo, ckpt, pocket)
    emit_smiles(sdf)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
