"""Map FKBP12 interface residues for the ON-target (receptor) and OFF-target
(calcineurin) surfaces. Defines the selectivity basis for de novo design.

ON-target  (1B6C): FKBP12 chain A vs TGF-beta type I receptor chain B  -> want to DISRUPT
OFF-target (1TCO): FKBP12 chain C vs calcineurin chains A/B             -> want to AVOID recreating
See docs/target-dossier-FKBP12.md / ADR-001.
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

from Bio.PDB import NeighborSearch, PDBParser

warnings.simplefilter("ignore")

ROOT = Path(__file__).resolve().parents[2]
STRUCT = ROOT / "data" / "raw" / "structures"
OUT = ROOT / "results" / "phase1" / "interfaces"
CUTOFF = 5.0  # Angstrom, heavy-atom contact

# (pdb, fkbp_chain, partner_chains, label)
JOBS = [
    ("1B6C", "A", ["B"], "on_target_receptor"),
    ("1TCO", "C", ["A", "B"], "off_target_calcineurin"),
]


def interface_residues(pdb_id, fkbp_chain, partner_chains):
    s = PDBParser(QUIET=True).get_structure(pdb_id, STRUCT / f"{pdb_id}.pdb")[0]
    partner_atoms = [a for cid in partner_chains for a in s[cid].get_atoms()
                     if a.get_parent().id[0] == " "]
    ns = NeighborSearch(partner_atoms)
    hits = {}
    for res in s[fkbp_chain]:
        if res.id[0] != " ":
            continue
        for atom in res:
            if ns.search(atom.coord, CUTOFF):
                hits[res.id[1]] = res.resname
                break
    return dict(sorted(hits.items()))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    summary = {}
    for pdb_id, fk, partners, label in JOBS:
        res = interface_residues(pdb_id, fk, partners)
        summary[label] = {"pdb": pdb_id, "fkbp_chain": fk, "partners": partners,
                          "n_residues": len(res),
                          "residues": [f"{n}:{r}" for n, r in res.items()]}
        print(f"{label} ({pdb_id}, FKBP12 chain {fk} vs {partners}): "
              f"{len(res)} interface residues")
        print("  " + ", ".join(f"{r}{n}" for n, r in res.items()))

    on = set(summary["on_target_receptor"]["residues"])
    off = set(summary["off_target_calcineurin"]["residues"])
    summary["analysis"] = {
        "shared": sorted(on & off),
        "on_only": sorted(on - off),
        "off_only": sorted(off - on),
    }
    print(f"\nShared (both interfaces): {sorted(on & off)}")
    print(f"OFF-only (calcineurin-specific, AVOID engaging): {sorted(off - on)}")

    (OUT / "interface_residues.json").write_text(json.dumps(summary, indent=2))
    print(f"\nWrote {OUT / 'interface_residues.json'}")


if __name__ == "__main__":
    main()
