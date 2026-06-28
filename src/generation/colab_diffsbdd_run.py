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


def sh(cmd, check=True):
    # Capture child output explicitly: ipykernel does NOT forward an OS-level
    # subprocess's stdout/stderr, so we must print it ourselves to see errors.
    print(f"$ {cmd}", flush=True)
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.stdout:
        print(r.stdout, flush=True)
    if r.stderr:
        print("[stderr]", r.stderr, flush=True)
    if check and r.returncode:
        raise RuntimeError(f"command failed (rc={r.returncode}): {cmd}")
    return r


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
    import torch
    print(f"TORCH {torch.__version__} CUDA {torch.version.cuda} avail {torch.cuda.is_available()}", flush=True)
    repo = WORK / "DiffSBDD"
    if not repo.exists():
        sh(f"git clone --depth 1 https://github.com/arneschneuing/DiffSBDD {repo}")
    # Patch DiffSBDD for modern Biopython (three_to_one removed in 1.80+; 1.79 has
    # no py3.12 wheel). Keep current biopython, shim three_to_one via protein_letters_3to1.
    shim = ('from Bio.PDB.Polypeptide import protein_letters_3to1 as _P3; '
            'three_to_one = lambda r: _P3.get(r.upper(), "X")')
    sh(f"cd {repo} && grep -rl 'from Bio.PDB.Polypeptide import three_to_one' --include='*.py' . | "
       f'xargs -r sed -i "s|from Bio.PDB.Polypeptide import three_to_one|{shim}|g"')
    # deps matched to the VM's torch/CUDA (not DiffSBDD's pinned old ones)
    sh(f'{sys.executable} -m pip -q install "rdkit>=2022.03" biopython scipy networkx pytorch-lightning openbabel-wheel')
    tv = torch.__version__.split("+")[0]
    cu = torch.version.cuda.replace(".", "")
    sh(f'{sys.executable} -m pip -q install torch-scatter '
       f'-f https://data.pyg.org/whl/torch-{tv}+cu{cu}.html')
    # torch>=2.6 defaults torch.load(weights_only=True), which rejects the
    # checkpoint's argparse.Namespace. We trust the Zenodo checkpoint -> force
    # weights_only=False via an auto-imported sitecustomize (no DiffSBDD edits).
    (WORK / "sitecustomize.py").write_text(
        "import torch as _t\n"
        "_orig = _t.load\n"
        "_t.load = lambda *a, **k: _orig(*a, **{**k, 'weights_only': False})\n")
    ckpt = repo / "crossdocked_fullatom_cond.ckpt"
    if not ckpt.exists():
        sh(f'wget -q -O {ckpt} '
           f'"https://zenodo.org/record/8183747/files/crossdocked_fullatom_cond.ckpt?download=1"')
        sh(f'ls -lh {ckpt}')
    return repo, ckpt


def generate(repo: Path, ckpt: Path, pocket: Path) -> Path:
    out = WORK / "fkbp12_generated.sdf"
    sh(f'cd {repo} && PYTHONPATH={WORK}:$PYTHONPATH {sys.executable} generate_ligands.py {ckpt} '
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
