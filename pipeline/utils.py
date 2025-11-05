"""
Utility functions for file management and temporary directories.
"""

import os
import tempfile
import shutil
from pathlib import Path


def setup_temp_directory(prefix='mol_calc_'):
    """
    Create a temporary directory for calculations.
    
    Args:
        prefix: Prefix for the temporary directory name
        
    Returns:
        str: Path to the created temporary directory
    """
    temp_dir = tempfile.mkdtemp(prefix=prefix)
    print(f"Working directory: {temp_dir}")
    return temp_dir


def cleanup_temp_directory(temp_dir):
    """
    Clean up temporary files and directory.
    
    Args:
        temp_dir: Path to the temporary directory to clean up
    """
    if temp_dir and os.path.exists(temp_dir):
        # Uncomment the line below to actually delete temp files
        # shutil.rmtree(temp_dir)
        print("Cleaned up temporary files")


def ensure_directory_exists(directory_path):
    """
    Ensure a directory exists, creating it if necessary.
    
    Args:
        directory_path: Path to the directory
    """
    os.makedirs(directory_path, exist_ok=True)

