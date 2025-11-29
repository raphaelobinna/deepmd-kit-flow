#!/usr/bin/env python3
"""
Prepare ABACUS input files from conformation PDB files.

This script processes each conformation PDB file and creates ABACUS input files
(INPUT, STRU, KPT) and copies pseudopotentials/orbitals.

Usage:
    python 3_prepare_abacus_inputs.py --conformations_dir data/conformations --output_dir data/abacus_inputs --pp_dir resources
"""

import os
import sys
import shutil
from pathlib import Path
from collections import defaultdict

# Element mapping
ELEMENT_MAP = {
    'H': 'H', 'C': 'C', 'N': 'N', 'O': 'O', 'S': 'S',
    'P': 'P', 'F': 'F', 'CL': 'Cl', 'BR': 'Br', 'I': 'I'
}

# Atomic masses
ATOMIC_MASSES = {
    'H': 1.008, 'C': 12.01, 'N': 14.01, 'O': 16.00, 'S': 32.07,
    'P': 30.97, 'F': 19.00, 'Cl': 35.45, 'Br': 79.90, 'I': 126.90
}

# Pseudopotential and orbital file names
PP_FILES = {
    'H': 'H_ONCV_PBE-1.0.upf',
    'C': 'C_ONCV_PBE-1.0.upf',
    'N': 'N_ONCV_PBE-1.0.upf',
    'O': 'O_ONCV_PBE-1.0.upf',
    'S': 'S_ONCV_PBE-1.0.upf',
}

ORB_FILES = {
    'H': 'H_gga_6au_100Ry_2s1p.orb',
    'C': 'C_gga_7au_100Ry_2s2p1d.orb',
    'N': 'N_gga_7au_100Ry_2s2p1d.orb',
    'O': 'O_gga_7au_100Ry_2s2p1d.orb',
    'S': None,
}


def parse_pdb(pdb_file):
    """Parse PDB file and extract elements and coordinates."""
    elements = []
    coords = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                parts = line.split()
                atom_name = parts[2] if len(parts) > 2 else ''
                element = parts[-1] if len(parts) > 10 else atom_name[0]
                
                try:
                    x = float(parts[5] if len(parts) > 5 else parts[6])
                    y = float(parts[6] if len(parts) > 6 else parts[7])
                    z = float(parts[7] if len(parts) > 7 else parts[8])
                except (ValueError, IndexError):
                    continue
                
                elements.append(element)
                coords.append((x, y, z))
    
    return elements, coords


def get_unique_elements(elements):
    """Get unique elements in order."""
    seen = set()
    unique = []
    for elem in elements:
        if elem not in seen:
            seen.add(elem)
            unique.append(elem)
    return unique


def generate_stru_file(elements, coords, output_file, box_size=20.0):
    """Generate ABACUS STRU file from elements and coordinates."""
    unique_elements = get_unique_elements(elements)
    element_coords = defaultdict(list)
    for elem, coord in zip(elements, coords):
        element_coords[elem].append(coord)
    
    with open(output_file, 'w') as f:
        f.write("ATOMIC_SPECIES\n")
        for elem in unique_elements:
            mass = ATOMIC_MASSES.get(elem, 1.0)
            pp_file = PP_FILES.get(elem, f"{elem}_ONCV_PBE-1.2.upf")
            f.write(f"{elem:2s} {mass:6.2f} {pp_file}\n")
        f.write("\n")
        
        f.write("NUMERICAL_ORBITAL\n")
        for elem in unique_elements:
            orb_file = ORB_FILES.get(elem)
            if orb_file:
                f.write(f"{orb_file}\n")
        f.write("\n")
        
        f.write("LATTICE_CONSTANT\n")
        f.write("1.8897261254578281\n")
        f.write("\n")
        
        f.write("LATTICE_VECTORS\n")
        f.write(f"{box_size:6.1f}  0  0\n")
        f.write(f"0  {box_size:6.1f}  0\n")
        f.write(f"0  0  {box_size:6.1f}\n")
        f.write("\n")
        
        f.write("ATOMIC_POSITIONS\n")
        f.write("Cartesian    # Cartesian(Unit is LATTICE_CONSTANT)\n")
        
        for elem in unique_elements:
            f.write(f"{elem}\n")
            f.write("0.0\n")
            f.write(f"{len(element_coords[elem])}\n")
            for x, y, z in element_coords[elem]:
                f.write(f"{x:12.6f} {y:12.6f} {z:12.6f} 1 1 1\n")


def generate_input_file(output_file, ntype, calculation='md', md_nstep=10, md_restart=None):
    """
    Generate ABACUS INPUT file.
    
    If md_restart is None, automatically detects if restart files exist.
    """
    # Auto-detect restart if not specified
    if md_restart is None:
        output_dir = Path(output_file).parent
        restart_file = output_dir / "OUT.ABACUS" / "Restart_md.dat"
        md_restart = 1 if restart_file.exists() else 0
    with open(output_file, 'w') as f:
        f.write("INPUT_PARAMETERS\n\n")
        f.write("#Parameters (1.General)\n")
        f.write(f"ntype                   {ntype}\n")
        f.write("symmetry                 0\n")
        f.write("nspin                    1\n\n")
        f.write("#Parameters (2.Iteration)\n")
        f.write("ecutwfc                 100\n")
        f.write("scf_thr                 1e-7\n")
        f.write("scf_nmax                120\n\n")
        f.write("#Parameters (3.Basis)\n")
        f.write("basis_type              lcao\n\n")
        f.write("#Parameters (4.Smearing)\n")
        f.write("smearing_method         gaussian\n")
        f.write("smearing_sigma          0.002\n\n")
        f.write("#Parameters (5.Mixing)\n")
        f.write("mixing_type             pulay\n")
        f.write("mixing_beta             0.4\n\n")
        
        if calculation == 'md':
            f.write("#Parameters (6.md)\n")
            f.write("calculation             md\n")
            f.write("cal_force               1\n")
            f.write("cal_stress              1\n")
            f.write(f"md_nstep                {md_nstep}\n")
            f.write("md_type                 nvt\n")
            f.write("md_tfirst               300\n")
            f.write(f"md_restart              {md_restart}\n")
            f.write("md_dumpfreq             1\n")
            f.write("out_stru                1\n")
            f.write("dump_force              1\n")
            f.write("dump_vel                1\n")
            f.write("dump_virial             1\n")
        else:
            f.write("#Parameters (6.scf)\n")
            f.write("calculation             scf\n")
            f.write("cal_force               1\n")
            f.write("cal_stress              0\n")


def generate_kpt_file(output_file):
    """Generate ABACUS KPT file."""
    with open(output_file, 'w') as f:
        f.write("K_POINTS\n")
        f.write("0\n")
        f.write("Gamma\n")
        f.write("1 1 1 0 0 0\n")


def prepare_abacus_inputs(conformations_dir, output_dir, pp_dir=None, calculation='md', md_nstep=10):
    """
    Prepare ABACUS inputs for all conformation PDB files.
    
    Parameters:
    -----------
    conformations_dir : str
        Directory containing conformation PDB files
    output_dir : str
        Output directory for ABACUS calculations
    pp_dir : str
        Directory containing pseudopotentials and orbitals
    calculation : str
        'md' for MD or 'scf' for single-point
    md_nstep : int
        Number of MD steps (default: 10 for testing)
    """
    conf_path = Path(conformations_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Find all PDB files (recursively, in case conformations are in fragment subdirectories)
    pdb_files = sorted(conf_path.rglob("*.pdb"))
    if not pdb_files:
        print(f"Error: No PDB files found in {conf_path}")
        return
    
    print(f"Preparing ABACUS inputs for {len(pdb_files)} conformations...")
    print(f"Output directory: {output_dir}\n")
    
    # Get unique elements from first conformation
    first_conf = pdb_files[0]
    elements, coords = parse_pdb(first_conf)
    if not elements:
        print(f"Error: Could not parse {first_conf}")
        return
    
    unique_elements = get_unique_elements(elements)
    ntype = len(unique_elements)
    print(f"Elements found: {unique_elements}")
    print(f"Number of types (ntype): {ntype}\n")
    
    # Copy pseudopotentials and orbitals if provided
    if pp_dir:
        pp_path = Path(pp_dir)
        print(f"Copying pseudopotentials and orbitals from {pp_dir}...")
        for element in unique_elements:
            if element in PP_FILES:
                src_pp = pp_path / PP_FILES[element]
                if src_pp.exists():
                    shutil.copy2(src_pp, output_path / PP_FILES[element])
                    print(f"  ✓ {PP_FILES[element]}")
            
            if element in ORB_FILES and ORB_FILES[element]:
                src_orb = pp_path / ORB_FILES[element]
                if src_orb.exists():
                    shutil.copy2(src_orb, output_path / ORB_FILES[element])
                    print(f"  ✓ {ORB_FILES[element]}")
        print()
    
    # Prepare inputs for each conformation
    success_count = 0
    for i, pdb_file in enumerate(pdb_files, 1):
        conf_name = pdb_file.stem
        conf_dir = output_path / conf_name
        conf_dir.mkdir(exist_ok=True)
        
        if i % 10 == 0 or i == len(pdb_files):
            print(f"[{i}/{len(pdb_files)}] Processing {conf_name}...")
        
        # Parse PDB
        elements, coords = parse_pdb(pdb_file)
        if not elements:
            continue
        
        # Generate STRU file
        generate_stru_file(elements, coords, conf_dir / "STRU")
        
        # Generate INPUT file (auto-detect restart if previous run exists)
        generate_input_file(conf_dir / "INPUT", ntype, calculation, md_nstep)
        
        # Generate KPT file
        generate_kpt_file(conf_dir / "KPT")
        
        # Copy pseudopotentials and orbitals
        if pp_dir:
            for element in unique_elements:
                if element in PP_FILES:
                    src_pp = output_path / PP_FILES[element]
                    if src_pp.exists():
                        shutil.copy2(src_pp, conf_dir / PP_FILES[element])
                
                if element in ORB_FILES and ORB_FILES[element]:
                    src_orb = output_path / ORB_FILES[element]
                    if src_orb.exists():
                        shutil.copy2(src_orb, conf_dir / ORB_FILES[element])
        
        success_count += 1
    
    print(f"\n✓ Complete! Prepared {success_count}/{len(pdb_files)} conformations")
    print(f"  Each conformation has INPUT, STRU, KPT, and PP/ORB files")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Prepare ABACUS input files from conformation PDB files"
    )
    parser.add_argument(
        "--conformations_dir",
        type=str,
        default="data/conformations",
        help="Directory containing conformation PDB files (default: data/conformations)"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="data/abacus_inputs",
        help="Output directory for ABACUS calculations (default: data/abacus_inputs)"
    )
    parser.add_argument(
        "--pp_dir",
        type=str,
        default="resources",
        help="Directory containing pseudopotentials and orbitals (default: resources)"
    )
    parser.add_argument(
        "--calculation",
        type=str,
        choices=['md', 'scf'],
        default='md',
        help="Calculation type: 'md' for MD or 'scf' for single-point (default: md)"
    )
    parser.add_argument(
        "--md_nstep",
        type=int,
        default=10,
        help="Number of MD steps (default: 10 for testing, use 50-200 for production)"
    )
    
    args = parser.parse_args()
    
    prepare_abacus_inputs(
        args.conformations_dir,
        args.output_dir,
        args.pp_dir,
        args.calculation,
        args.md_nstep
    )


if __name__ == "__main__":
    main()

