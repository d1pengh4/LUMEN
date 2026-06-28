"""Fetch FKBP12-related experimental structures from the RCSB PDB.

Target rationale: see docs/decisions/ADR-001-target-redefinition.md
Each structure anchors a specific design objective:
  - drug pocket (where the de novo molecule binds)
  - ON-target interface (FKBP12 -- type I receptor; signal to activate)
  - OFF-target interface (FKBP12 -- calcineurin; selectivity to avoid)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import requests

RAW = Path(__file__).resolve().parents[2] / "data" / "raw" / "structures"

# id -> (role, description)
STRUCTURES = {
    "1FKF": ("drug_pocket", "FKBP12 + FK506 complex; defines the small-molecule binding pocket"),
    "1B6C": ("on_target_interface", "FKBP12 + TGF-beta type I receptor (GS domain); ON-target interface analogous to BMP type I receptors"),
    "1TCO": ("off_target_interface", "FKBP12 + FK506 + calcineurin; OFF-target interface to AVOID (immunosuppression)"),
    "2PPN": ("apo", "Apo FKBP12, high resolution; ligand-free reference"),
}

PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"


def fetch(pdb_id: str) -> bytes:
    r = requests.get(PDB_URL.format(pdb_id=pdb_id), timeout=30)
    r.raise_for_status()
    return r.content


def main() -> int:
    RAW.mkdir(parents=True, exist_ok=True)
    manifest = []
    for pdb_id, (role, desc) in STRUCTURES.items():
        out = RAW / f"{pdb_id}.pdb"
        try:
            data = fetch(pdb_id)
        except requests.HTTPError as e:
            print(f"[FAIL] {pdb_id}: {e}", file=sys.stderr)
            return 1
        out.write_bytes(data)
        n_atoms = sum(1 for line in data.splitlines() if line.startswith((b"ATOM", b"HETATM")))
        print(f"[OK]   {pdb_id}  {role:22s} {len(data):>8d} B  {n_atoms:>5d} atom-records")
        manifest.append({"pdb_id": pdb_id, "role": role, "description": desc,
                         "file": str(out.relative_to(RAW.parents[1])), "atom_records": n_atoms})

    (RAW / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"\nWrote {len(manifest)} structures + manifest to {RAW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
