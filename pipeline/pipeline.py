"""
Main pipeline orchestration for molecular property calculations.
"""

import os
import pandas as pd
from .descriptors import calculate_2d_descriptors
from .structure_generation import generate_3d_structure
from .external_tools import run_crest, run_xtb
from .data_extraction import extract_xtb_data
from .utils import setup_temp_directory, cleanup_temp_directory, ensure_directory_exists


class MolecularPropertyPipeline:
    """
    Complete pipeline for molecular property calculation:
    1. Read SMILES from CSV
    2. Generate 3D structures (RDKit ETKDG + MMFF)
    3. Run CREST for conformer search
    4. Run xTB with solvent correction
    5. Extract all properties and save to CSV
    """
    
    def __init__(self, input_csv, smiles_column='smiles', output_csv='molecular_properties_results.csv'):
        """
        Initialize the pipeline.
        
        Args:
            input_csv: Path to input CSV file containing SMILES
            smiles_column: Name of the column containing SMILES strings
            output_csv: Path to output CSV file for results
        """
        self.input_csv = input_csv
        self.smiles_column = smiles_column
        self.output_csv = output_csv
        self.results = []
        self.temp_dir = None
    
    def process_molecule(self, idx, smiles):
        """
        Process a single molecule through the entire pipeline.
        
        Args:
            idx: Index of the molecule
            smiles: SMILES string representation
            
        Returns:
            dict: Dictionary containing all calculated properties
        """
        print(f"\n{'='*60}")
        print(f"Processing molecule {idx + 1}: {smiles}")
        print(f"{'='*60}")
        
        result = {
            'Index': idx,
            'SMILES': smiles,
            'Status': 'Failed'
        }
        
        try:
            # Step 1: Calculate 2D descriptors
            print("Step 1: Calculating 2D descriptors...")
            result.update(calculate_2d_descriptors(smiles))
            
            # Step 2: Generate 3D structure
            mol_dir = os.path.join(self.temp_dir, f'molecule_{idx}')
            ensure_directory_exists(mol_dir)
            
            print("Step 2: Generating 3D structure...")
            initial_xyz = os.path.join(mol_dir, 'initial.xyz')
            if not generate_3d_structure(smiles, initial_xyz):
                return result
            
            # Step 3: Run CREST
            print("Step 3: Running CREST...")
            crest_xyz = os.path.join(mol_dir, 'crest_best.xyz')
            if not run_crest(initial_xyz, crest_xyz):
                crest_xyz = initial_xyz
            
            # Step 4: Run xTB
            print("Step 4: Running xTB...")
            xtb_log = os.path.join(mol_dir, 'xtb.log')
            if not run_xtb(crest_xyz, xtb_log):
                return result
            
            # Step 5: Extract xTB data
            print("Step 5: Extracting xTB data...")
            result.update(extract_xtb_data(xtb_log))
            
            result['Status'] = 'Success'
        
        except Exception as e:
            result['Status'] = f"Error: {e}"
        
        return result
    
    def run_pipeline(self):
        """Execute the complete pipeline for all molecules in the input CSV."""
        print("\n" + "="*70)
        print("MOLECULAR PROPERTY CALCULATION PIPELINE")
        print("="*70)
        
        try:
            df = pd.read_csv(self.input_csv)
        except Exception as e:
            print(f"Error reading CSV: {e}")
            return
        
        if self.smiles_column not in df.columns:
            print(f"Column '{self.smiles_column}' not found.")
            return
        
        print(f"Found {len(df)} molecules.")
        
        # Setup temporary directory
        self.temp_dir = setup_temp_directory()
        
        # Process each molecule
        for idx, smiles in enumerate(df[self.smiles_column]):
            result = self.process_molecule(idx, smiles)
            self.results.append(result)
            
            # Save intermediate results every 10 molecules
            if (idx + 1) % 10 == 0:
                self.save_results(self.output_csv + ".tmp")
        
        # Save final results
        self.save_results(self.output_csv)
        cleanup_temp_directory(self.temp_dir)
        self.print_summary()
    
    def save_results(self, output_file):
        """
        Save results to CSV file.
        
        Args:
            output_file: Path to output CSV file
        """
        df_results = pd.DataFrame(self.results)
        df_results.to_csv(output_file, index=False)
        print(f"Results saved to: {output_file}")
    
    def print_summary(self):
        """Print summary statistics of the pipeline run."""
        print("\n" + "="*70)
        print("PIPELINE SUMMARY")
        print("="*70)
        
        total = len(self.results)
        success = sum(1 for r in self.results if r['Status'] == 'Success')
        
        print(f"Total: {total}")
        print(f"Successful: {success}")
        print(f"Failed: {total - success}")

