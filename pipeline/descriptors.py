"""
2D molecular descriptor calculations using RDKit.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, MolSurf, Crippen, rdMolDescriptors


def calculate_2d_descriptors(smiles):
    """
    Calculate RDKit 2D descriptors from SMILES string.
    
    Args:
        smiles: SMILES string representation of the molecule
        
    Returns:
        dict: Dictionary containing calculated descriptors, or empty dict on error
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}
        
        descriptors = {
            'MolWt': Descriptors.MolWt(mol),
            'MolLogP': Crippen.MolLogP(mol),
            'RingCount': rdMolDescriptors.CalcNumRings(mol),
            'HBD': rdMolDescriptors.CalcNumHBD(mol),
            'HBA': rdMolDescriptors.CalcNumHBA(mol),
            'TPSA': MolSurf.TPSA(mol),
            'Fsp3': rdMolDescriptors.CalcFractionCSP3(mol),
            'SP2_Carbons': _count_sp2_carbons(mol),
            'SP3_Carbons': _count_sp3_carbons(mol),
            'ChiralCenters': len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        }
        return descriptors
    except Exception as e:
        print(f"Error calculating 2D descriptors: {e}")
        return {}


def _count_sp2_carbons(mol):
    """Count SP2 hybridized carbon atoms."""
    return sum(1 for atom in mol.GetAtoms() 
               if atom.GetHybridization() == Chem.HybridizationType.SP2 
               and atom.GetAtomicNum() == 6)


def _count_sp3_carbons(mol):
    """Count SP3 hybridized carbon atoms."""
    return sum(1 for atom in mol.GetAtoms() 
               if atom.GetHybridization() == Chem.HybridizationType.SP3 
               and atom.GetAtomicNum() == 6)

