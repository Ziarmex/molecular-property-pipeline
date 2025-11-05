"""
Main entry point for the molecular property calculation pipeline.
"""

from pipeline import MolecularPropertyPipeline
from pipeline.config import DEFAULT_INPUT_CSV, DEFAULT_SMILES_COLUMN, DEFAULT_OUTPUT_CSV


def main():
    """Run the molecular property calculation pipeline."""
    pipeline = MolecularPropertyPipeline(
        input_csv=DEFAULT_INPUT_CSV,
        smiles_column=DEFAULT_SMILES_COLUMN,
        output_csv=DEFAULT_OUTPUT_CSV
    )
    
    pipeline.run_pipeline()


if __name__ == "__main__":
    main()

