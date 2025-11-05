"""
Configuration constants for the molecular property pipeline.
"""

# Default file paths
DEFAULT_INPUT_CSV = 'data/test.csv'
DEFAULT_SMILES_COLUMN = 'smiles'
DEFAULT_OUTPUT_CSV = 'data/results/molecular_properties_complete.csv'

# CREST settings
CREST_TIMEOUT = 1200  # seconds
CREST_METHOD = 'gfn2'
CREST_SOLVENT = 'water'
CREST_EWIN = '6'
CREST_CORES = '1'

# xTB settings
XTB_TIMEOUT = 600  # seconds
XTB_METHOD = '2'  # GFN2
XTB_SOLVENT = 'water'
XTB_OPTIMIZE = True

# RDKit settings
RDKIT_RANDOM_SEED = 42

