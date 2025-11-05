# Molecular Property Calculation Pipeline

A modular Python pipeline for calculating molecular properties using RDKit, CREST, and xTB.

## Project Structure

```
pipeline/
├── pipeline/                # Main package directory
│   ├── __init__.py          # Package initialization
│   ├── config.py            # Configuration constants
│   ├── descriptors.py       # 2D molecular descriptor calculations
│   ├── structure_generation.py  # 3D structure generation using RDKit
│   ├── external_tools.py    # Interface for CREST and xTB
│   ├── data_extraction.py   # Extract properties from xTB logs
│   ├── pipeline.py          # Main pipeline orchestration
│   └── utils.py             # Utility functions
├── data/                    # Data directory
│   ├── test.csv             # Input CSV file with SMILES
│   └── results/             # Output directory for results
├── main.py                  # Main entry point
├── requirements.txt         # Python dependencies
├── .gitignore              # Git ignore rules
└── README.md               # This file
```

## Modules

### `config.py`
Contains all configuration constants including:
- File paths and column names
- CREST settings (timeout, method, solvent)
- xTB settings (timeout, method, solvent)
- RDKit settings

### `descriptors.py`
Calculates 2D molecular descriptors from SMILES using RDKit:
- Molecular weight
- LogP
- Ring count
- Hydrogen bond donors/acceptors
- TPSA
- Fsp3
- SP2/SP3 carbon counts
- Chiral centers

### `structure_generation.py`
Generates 3D molecular structures:
- Uses RDKit ETKDG for initial 3D embedding
- MMFF optimization (falls back to UFF if needed)
- Exports to XYZ format

### `external_tools.py`
Interfaces with external computational chemistry tools:
- `run_crest()`: Runs CREST conformer search
- `run_xtb()`: Runs xTB quantum chemistry calculations

### `data_extraction.py`
Extracts properties from xTB log files:
- HOMO/LUMO energies
- HOMO-LUMO gap
- Total energy
- Dipole moment
- Mulliken charges

### `pipeline.py`
Main orchestration class (`MolecularPropertyPipeline`):
- Coordinates all steps of the calculation pipeline
- Manages molecule processing workflow
- Handles result collection and saving

### `utils.py`
Utility functions:
- Temporary directory management
- File system operations

## Setup

### 1. Create Virtual Environment

```bash
python -m venv venv
```

### 2. Activate Virtual Environment

**Windows:**
```bash
venv\Scripts\activate
```

**Linux/Mac:**
```bash
source venv/bin/activate
```

### 3. Install Dependencies

```bash
pip install -r requirements.txt
```

### 4. Install External Tools

The pipeline requires CREST and xTB to be installed and available in your PATH:

- **CREST**: Download from [https://github.com/grimme-lab/crest](https://github.com/grimme-lab/crest)
- **xTB**: Download from [https://github.com/grimme-lab/xtb](https://github.com/grimme-lab/xtb)

Make sure both `crest` and `xtb` commands are available in your terminal.

### 5. Prepare Input Data

Place your input CSV file in the `data/` directory. The file should contain a column named `smiles` with SMILES strings.

Example `data/test.csv`:
```csv
smiles,type
O=S(C=C)C1=CC=CC=C1,SUL
NC1=CC=C(OC)C=C1,An
```

## Usage

### Command Line

Simply run:
```bash
python main.py
```

This will:
1. Read SMILES from `data/test.csv`
2. Process each molecule through the pipeline
3. Save results to `data/results/molecular_properties_complete.csv`

### Python API

```python
from pipeline import MolecularPropertyPipeline

pipeline = MolecularPropertyPipeline(
    input_csv='data/test.csv',
    smiles_column='smiles',
    output_csv='data/results/my_results.csv'
)

pipeline.run_pipeline()
```

## Dependencies

### Python Packages
- pandas >= 2.3.1
- numpy >= 2.3.2
- rdkit-pypi >= 2025.3.5
- pillow >= 11.3.0

### External Binaries
- **CREST**: Conformer search tool (must be in PATH)
- **xTB**: Quantum chemistry calculator (must be in PATH)

## Configuration

Edit `pipeline/config.py` to customize:
- Input/output file paths
- CREST and xTB parameters
- Timeout values
- Solvent models
- RDKit random seed

## Output

The pipeline generates a CSV file with the following columns:
- **Index**: Molecule index
- **SMILES**: Input SMILES string
- **Status**: Processing status (Success/Failed)
- **2D Descriptors**: MolWt, MolLogP, RingCount, HBD, HBA, TPSA, Fsp3, etc.
- **Quantum Properties**: HOMO/LUMO energies, gap, total energy, dipole moment, Mulliken charges

