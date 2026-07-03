#!/usr/bin/env python
"""G3 MD bootstrap (the "see bootstrap" referenced by colab_md.py).

Colab VMs are ephemeral and ship no conda; the MD stack (openmm, openff-toolkit,
pdbfixer, ambertools for GAFF) is conda-only, not pip. This installs a static
micromamba and builds a self-contained env at /content/mm/envs/md so the run is
reproducible instead of set up by hand in a since-dead session.

Run once per fresh session (colab exec -s <sess> -f colab_md_setup.py --timeout 1800),
then run the science with that env's python:
    micromamba run -p /content/mm/envs/md python colab_md.py "<SMILES>" <ns>
"""
import subprocess

ENV = "/content/mm/envs/md"
PKGS = ("python=3.11 openmm openmmforcefields openff-toolkit "
        "pdbfixer ambertools mdtraj rdkit")


def sh(cmd, timeout=None):
    print(f"$ {cmd}", flush=True)
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
    if r.stdout:
        print(r.stdout[-2000:], flush=True)
    if r.stderr:
        print("[stderr]", r.stderr[-2000:], flush=True)
    return r.returncode


def main():
    print("STEP install micromamba", flush=True)
    sh("cd /content && curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest "
       "| tar -xj bin/micromamba && /content/bin/micromamba --version", timeout=300)

    print("STEP create env (conda solve; slow)", flush=True)
    rc = sh(f"export MAMBA_ROOT_PREFIX=/content/mm && "
            f"/content/bin/micromamba create -y -p {ENV} -c conda-forge {PKGS}",
            timeout=1700)
    if rc != 0:
        raise SystemExit("env create failed")

    print("STEP verify imports", flush=True)
    rc = sh(f"/content/bin/micromamba run -p {ENV} python -c "
            f"'import openmm, openff.toolkit, pdbfixer, openmmforcefields, mdtraj, rdkit;"
            f"print(\"OK openmm\", openmm.__version__)'", timeout=300)
    print("===SETUP_OK===" if rc == 0 else "===SETUP_FAIL===", flush=True)
    raise SystemExit(rc)


if __name__ == "__main__":
    main()
