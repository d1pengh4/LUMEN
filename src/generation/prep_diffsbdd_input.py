"""Prepare the pocket input for DiffSBDD (pocket-conditioned 3D diffusion, GPU).

DiffSBDD defines the binding pocket from a protein PDB + a reference ligand.
We use 1FKF: FKBP12 (chain A, res 1-107) + FK506 (FK5, A:108) as the reference
ligand. Waters are stripped. The reference ligand id for DiffSBDD is "A:108".

Output: data/processed/diffsbdd_input/fkbp12_fk506.pdb  (+ README with the cmd)
"""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "data" / "raw" / "structures" / "1FKF.pdb"
OUT = ROOT / "data" / "processed" / "diffsbdd_input"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    keep = []
    for l in SRC.read_text().splitlines():
        if l.startswith("ATOM"):
            keep.append(l)
        elif l.startswith("HETATM") and l[17:20].strip() == "FK5":
            keep.append(l)
    keep.append("END")
    pdb = OUT / "fkbp12_fk506.pdb"
    pdb.write_text("\n".join(keep) + "\n")

    n_prot = sum(1 for l in keep if l.startswith("ATOM"))
    n_lig = sum(1 for l in keep if l.startswith("HETATM"))
    print(f"Wrote {pdb}: {n_prot} protein atoms + {n_lig} FK506 atoms (waters removed)")
    print("DiffSBDD reference ligand id: A:108")
    (OUT / "README.txt").write_text(
        "DiffSBDD input for FKBP12 (LUMEN project)\n"
        "Protein+ligand: fkbp12_fk506.pdb (FKBP12 chain A 1-107 + FK506 A:108)\n"
        "Reference ligand for --ref_ligand: A:108\n\n"
        "Generation command (run in Colab GPU, in the DiffSBDD repo):\n"
        "  python generate_ligands.py <checkpoint.ckpt> \\\n"
        "    --pdbfile fkbp12_fk506.pdb --ref_ligand A:108 \\\n"
        "    --outfile fkbp12_generated.sdf --n_samples 500 --sanitize\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
