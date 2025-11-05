"""
Interface for external computational chemistry tools (CREST and xTB).
"""

import os
import shutil
import subprocess
from .config import (
    CREST_TIMEOUT, CREST_METHOD, CREST_SOLVENT, CREST_EWIN, CREST_CORES,
    XTB_TIMEOUT, XTB_METHOD, XTB_SOLVENT, XTB_OPTIMIZE
)


def run_crest(input_xyz, output_xyz):
    """
    Run CREST conformer search (corrected for CREST 3.0.x).
    
    Args:
        input_xyz: Path to input XYZ file
        output_xyz: Path to output XYZ file (crest_best.xyz will be copied here)
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        workdir = os.path.dirname(input_xyz)

        cmd = [
            "crest",
            input_xyz,
            f"-{CREST_METHOD}",
            "-alpb", CREST_SOLVENT,
            "-ewin", CREST_EWIN,
            "-c", CREST_CORES
        ]

        result = subprocess.run(
            cmd,
            cwd=workdir,
            capture_output=True,
            text=True,
            timeout=CREST_TIMEOUT
        )

        if result.returncode != 0:
            print(f"CREST failed:\n{result.stderr}")
            shutil.copy(input_xyz, output_xyz)
            return False

        crest_best = os.path.join(workdir, "crest_best.xyz")

        if os.path.exists(crest_best):
            if os.path.abspath(crest_best) != os.path.abspath(output_xyz):
                shutil.copy(crest_best, output_xyz)
            else:
                shutil.copy(crest_best, output_xyz + "_copy.xyz")
            return True

        print("CREST did NOT produce crest_best.xyz — using initial structure")
        shutil.copy(input_xyz, output_xyz)
        return False

    except subprocess.TimeoutExpired:
        print("CREST timed out — using initial structure")
        shutil.copy(input_xyz, output_xyz)
        return False

    except Exception as e:
        print(f"Unexpected CREST error: {e}")
        shutil.copy(input_xyz, output_xyz)
        return False


def run_xtb(input_xyz, output_log):
    """
    Run xTB calculation with solvent correction.
    
    Args:
        input_xyz: Path to input XYZ file
        output_log: Path to output log file
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        cmd = [
            'xtb',
            input_xyz,
            '--gfn', XTB_METHOD,
            '--alpb', XTB_SOLVENT
        ]
        
        if XTB_OPTIMIZE:
            cmd.append('--opt')
        
        result = subprocess.run(
            cmd,
            cwd=os.path.dirname(input_xyz),
            capture_output=True,
            text=True,
            timeout=XTB_TIMEOUT
        )
        
        with open(output_log, 'w') as f:
            f.write(result.stdout)
        
        if result.returncode != 0:
            print(f"xTB failed: {result.stderr}")
            return False
        
        return True
    except Exception as e:
        print(f"Error running xTB: {e}")
        return False

