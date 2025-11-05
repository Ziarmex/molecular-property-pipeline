"""
Extract molecular properties from xTB log files.
"""

import re
import numpy as np


def extract_xtb_data(log_file):
    """
    Extract all properties from xTB log file.
    
    Args:
        log_file: Path to xTB log file
        
    Returns:
        dict: Dictionary containing extracted properties
    """
    try:
        with open(log_file, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        
        data = {}
        
        # Extract HOMO/LUMO energies
        _extract_homo_lumo(content, data)
        
        # Extract HOMO-LUMO gap
        _extract_hl_gap(content, data)
        
        # Extract total energy
        _extract_total_energy(content, data)
        
        # Extract dipole moment
        _extract_dipole_moment(content, data)
        
        # Extract Mulliken charges
        _extract_mulliken_charges(content, data)
        
        return data
    
    except Exception as e:
        print(f"Error extracting xTB data: {e}")
        return {}


def _extract_homo_lumo(content, data):
    """Extract HOMO and LUMO energies from xTB output."""
    homo_lumo_pattern = (
        r'(\d+)\s+2\.0000\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
        r'\s+\(HOMO\)\s*\n\s*(\d+)\s+(-?\d+\.\d+)'
        r'\s+(-?\d+\.\d+)\s+\(LUMO\)'
    )
    matches = re.findall(homo_lumo_pattern, content)
    if matches:
        m = matches[-1]  # Use last match
        data['HOMO_Eh'] = float(m[1])
        data['HOMO_eV'] = float(m[2])
        data['LUMO_Eh'] = float(m[4])
        data['LUMO_eV'] = float(m[5])


def _extract_hl_gap(content, data):
    """Extract HOMO-LUMO gap from xTB output."""
    gap_pattern = r'HL-Gap\s+(-?\d+\.\d+)\s+Eh\s+(-?\d+\.\d+)\s+eV'
    matches = re.findall(gap_pattern, content)
    if matches:
        m = matches[-1]  # Use last match
        data['HLGap_Eh'] = float(m[0])
        data['HLGap_eV'] = float(m[1])


def _extract_total_energy(content, data):
    """Extract total energy from xTB output."""
    energy_pattern = r'::\s+total energy\s+(-?\d+\.\d+)\s+Eh\s+::'
    matches = re.findall(energy_pattern, content)
    if matches:
        data['TotalEnergy_Eh'] = float(matches[-1])


def _extract_dipole_moment(content, data):
    """Extract dipole moment from xTB output."""
    dipole_pattern = (
        r'molecular dipole:\s*\n\s*x\s+y\s+z\s+tot\s+\(Debye\)'
        r'\s*\n\s*q only:.*?\n\s*full:\s+(-?\d+\.\d+)'
        r'\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)'
    )
    match = re.search(dipole_pattern, content, re.DOTALL)
    if match:
        data['Dipole_X'] = float(match.group(1))
        data['Dipole_Y'] = float(match.group(2))
        data['Dipole_Z'] = float(match.group(3))
        data['Dipole_Total_Debye'] = float(match.group(4))


def _extract_mulliken_charges(content, data):
    """Extract Mulliken charges from xTB output (α(0) version)."""
    charges_pattern = (
        r'#\s+Z\s+covCN\s+q\s+C6AA\s+α\(0\)\s*\n'
        r'((?:\s*\d+\s+\d+\s+[\d.]+\s+([-\d.]+)\s+[\d.]+\s+[\d.]+\s*\n)+)'
    )
    match = re.search(charges_pattern, content)
    if match:
        block = match.group(1)
        charges = []
        for line in block.strip().split('\n'):
            parts = line.split()
            if len(parts) >= 4:
                try:
                    charges.append(float(parts[3]))  # q-column
                except:
                    pass
        
        if charges:
            data['Mulliken_Avg'] = np.mean(charges)
            data['Mulliken_Max'] = max(charges)
            data['Mulliken_Min'] = min(charges)
            data['Mulliken_StdDev'] = np.std(charges)

