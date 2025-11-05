"""
3D molecular structure generation using RDKit.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from .config import RDKIT_RANDOM_SEED


def generate_3d_structure(smiles, output_xyz):
    """
    Generate 3D structure using RDKit ETKDG + MMFF optimization.
    
    Args:
        smiles: SMILES string representation of the molecule
        output_xyz: Path to output XYZ file
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates using ETKDG
        params = AllChem.ETKDGv3()
        params.randomSeed = RDKIT_RANDOM_SEED
        if AllChem.EmbedMolecule(mol, params) == -1:
            print("Failed to embed molecule")
            return False
        
        # Optimize geometry
        if AllChem.MMFFOptimizeMolecule(mol) == 1:
            print("MMFF optimization failed, trying UFF")
            AllChem.UFFOptimizeMolecule(mol)
        
        # Write XYZ file
        _write_xyz_file(mol, smiles, output_xyz)
        return True
    except Exception as e:
        print(f"Error generating 3D structure: {e}")
        return False


def _write_xyz_file(mol, smiles, output_xyz):
    """
    Write molecule coordinates to XYZ file format.
    
    Args:
        mol: RDKit molecule object with 3D coordinates
        smiles: SMILES string for comment line
        output_xyz: Path to output XYZ file
    """
    conf = mol.GetConformer()
    with open(output_xyz, 'w') as f:
        f.write(f"{mol.GetNumAtoms()}\n")
        f.write(f"Generated from SMILES: {smiles}\n")
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            f.write(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")

