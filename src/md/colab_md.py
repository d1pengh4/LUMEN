#!/usr/bin/env python
"""G3: MD pose-stability test on Colab GPU (detached). Decisive check that a
GNINA-top lead actually STAYS bound in the FKBP12 pocket (vs a lucky docking pose).

Pipeline (self-contained on the VM):
  1. install openmm + openmmforcefields + openff-toolkit + pdbfixer + mdtraj + rdkit
  2. fetch 1FKF, protein = chain A (ATOM); clean with PDBFixer (add H, pH 7)
  3. download gnina, dock the lead SMILES into the pocket (autobox on FK506) -> top pose
  4. SystemGenerator (ff14SB + GAFF-2.11 + TIP3P), build protein+ligand complex, solvate
  5. minimize; short NPT MD (default ~1 ns); record ligand heavy-atom RMSD vs start
  6. print RMSD summary (mean/max/final) between markers -> PASS if ligand stays bound

Honest: heavy dependency stack; expect iteration. Every step prints a STEP marker.
Arg: SMILES of the lead. Env-frugal: single lead, short run.
"""
import subprocess, sys, urllib.request, json
from pathlib import Path

LIG_SMILES = sys.argv[1] if len(sys.argv) > 1 else "C=C[C@H]1CN(c2cccc(NCC)c2)C(=O)[C@@H]2CCC[C@H]1N2S(=O)(=O)c1ccc(C)cc1"
NS = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0  # ns
WORK = Path("/content/md"); WORK.mkdir(parents=True, exist_ok=True)


def sh(cmd, check=True, timeout=None):
    print(f"$ {cmd}", flush=True)
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
    if r.stdout: print(r.stdout[-1500:], flush=True)
    if r.stderr: print("[stderr]", r.stderr[-1500:], flush=True)
    if check and r.returncode:
        raise RuntimeError(f"failed rc={r.returncode}: {cmd}")
    return r


def main():
    print("STEP nvidia-smi", flush=True); sh("nvidia-smi -L", check=False)
    # deps come from the conda env this script is launched with (see bootstrap);
    # no pip install here (openff-toolkit / pdbfixer / ambertools are conda-only).

    # ---- protein prep ----
    print("STEP fetch+prep protein (1FKF chain A)", flush=True)
    raw = urllib.request.urlopen("https://files.rcsb.org/download/1FKF.pdb", timeout=60).read().decode()
    (WORK/"prot_raw.pdb").write_text("\n".join(l for l in raw.splitlines() if l.startswith("ATOM")) + "\nEND\n")
    (WORK/"ref.pdb").write_text("\n".join(l for l in raw.splitlines() if l.startswith("HETATM") and l[17:20].strip()=="FK5") + "\nEND\n")
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    fixer = PDBFixer(filename=str(WORK/"prot_raw.pdb"))
    fixer.findMissingResidues(); fixer.findMissingAtoms(); fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(WORK/"protein.pdb","w"))
    print("  protein.pdb written", flush=True)

    # ---- dock lead to get a pocket pose ----
    print("STEP dock lead (gnina)", flush=True)
    if not (WORK/"gnina").exists():
        sh("wget -q https://github.com/gnina/gnina/releases/download/v1.3/gnina -O /content/md/gnina && chmod +x /content/md/gnina")
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")
    m = Chem.AddHs(Chem.MolFromSmiles(LIG_SMILES))
    AllChem.EmbedMolecule(m, randomSeed=42, maxAttempts=50, useRandomCoords=True)
    AllChem.MMFFOptimizeMolecule(m, maxIters=200)
    Chem.MolToMolFile(m, str(WORK/"lig_init.sdf"))
    sh(f"cd {WORK} && ./gnina -r protein.pdb -l lig_init.sdf --autobox_ligand ref.pdb "
       f"-o lig_docked.sdf --cnn_scoring rescore --exhaustiveness 8 --num_modes 1 --seed 42 2>&1 | tail -3", check=False)

    # ---- parameterize + build complex ----
    print("STEP parameterize + solvate (SystemGenerator)", flush=True)
    from openff.toolkit import Molecule
    from openmmforcefields.generators import SystemGenerator
    from openmm import app, unit, LangevinMiddleIntegrator, Platform, MonteCarloBarostat
    import openmm as mm
    import mdtraj as md
    import numpy as np

    lig = Molecule.from_file(str(WORK/"lig_docked.sdf"))
    protein = app.PDBFile(str(WORK/"protein.pdb"))
    ff_kwargs = dict(constraints=app.HBonds, rigidWater=True, removeCMMotion=True,
                     hydrogenMass=3*unit.amu)
    sysgen = SystemGenerator(forcefields=["amber/ff14SB.xml", "amber/tip3p_standard.xml"],
                             small_molecule_forcefield="gaff-2.11",
                             molecules=[lig], forcefield_kwargs=ff_kwargs)
    modeller = app.Modeller(protein.topology, protein.positions)
    lig_top = lig.to_topology().to_openmm(); lig_pos = lig.conformers[0].to_openmm()
    modeller.add(lig_top, lig_pos)
    modeller.addSolvent(sysgen.forcefield, model="tip3p", padding=1.0*unit.nanometer,
                        ionicStrength=0.15*unit.molar)
    system = sysgen.create_system(modeller.topology, molecules=[lig])
    system.addForce(MonteCarloBarostat(1*unit.bar, 300*unit.kelvin))
    integ = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.004*unit.picoseconds)
    sim = app.Simulation(modeller.topology, system, integ,
                         Platform.getPlatformByName("CUDA"))
    sim.context.setPositions(modeller.positions)

    print("STEP minimize", flush=True)
    sim.minimizeEnergy()
    sim.context.setVelocitiesToTemperature(300*unit.kelvin)

    print(f"STEP MD {NS} ns", flush=True)
    steps = int(NS * 1e6 / 4)  # ns -> number of 4 fs steps
    dcd = str(WORK/"traj.dcd")
    sim.reporters.append(app.DCDReporter(dcd, 2500))  # every 10 ps
    top_pdb = str(WORK/"complex.pdb")
    app.PDBFile.writeFile(modeller.topology, modeller.positions, open(top_pdb,"w"))
    sim.step(steps)
    print("  MD done", flush=True)

    # ---- ligand RMSD ----
    print("STEP RMSD analysis", flush=True)
    t = md.load(dcd, top=top_pdb)
    lig_atoms = t.topology.select("resname UNK or resname LIG or resname MOL")
    if len(lig_atoms) == 0:
        # ligand is the non-protein, non-water/ion residue; pick heavy atoms not protein/water
        lig_atoms = t.topology.select("not protein and not water and element != Na and element != Cl and element H")
    prot_ca = t.topology.select("protein and name CA")
    t.superpose(t, 0, atom_indices=prot_ca)
    import numpy as np
    lig_heavy = [a for a in lig_atoms if t.topology.atom(a).element.symbol != "H"]
    ref = t.xyz[0, lig_heavy]
    rmsd = np.sqrt(((t.xyz[:, lig_heavy] - ref)**2).sum(axis=2).mean(axis=1)) * 10  # nm->A
    print("\n===MD_RESULT_BEGIN===", flush=True)
    print(f"ns={NS} frames={t.n_frames} lig_heavy_atoms={len(lig_heavy)}", flush=True)
    print(f"ligand_RMSD_A mean={rmsd.mean():.2f} max={rmsd.max():.2f} final={rmsd[-1]:.2f}", flush=True)
    print(f"smiles={LIG_SMILES}", flush=True)
    print("===MD_RESULT_END===", flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
