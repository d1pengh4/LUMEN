DiffSBDD input for FKBP12 (LUMEN project)
Protein+ligand: fkbp12_fk506.pdb (FKBP12 chain A 1-107 + FK506 A:108)
Reference ligand for --ref_ligand: A:108

Generation command (run in Colab GPU, in the DiffSBDD repo):
  python generate_ligands.py <checkpoint.ckpt> \
    --pdbfile fkbp12_fk506.pdb --ref_ligand A:108 \
    --outfile fkbp12_generated.sdf --n_samples 500 --sanitize
