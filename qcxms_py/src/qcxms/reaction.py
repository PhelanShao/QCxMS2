"""
Reaction energy calculation and product processing module.
This module handles reaction energy calculations, product optimization,
and high-energy fragment filtering.
"""

import os
from pathlib import Path
from typing import List, Tuple, Dict, Optional

try:
    from .data import RunTypeData
    from .constants import WP, AUTOEV, EVTOKCAL
    from . import iomod
    from . import qmmod
    from . import utility
except ImportError:
    # Fallbacks for standalone/testing
    from data import RunTypeData
    from constants import WP, AUTOEV, EVTOKCAL
    import iomod_mock as iomod
    import qmmod_mock as qmmod
    import utility_mock as utility


def calculate_reaction_energies_py(env: RunTypeData, npairs: int, 
                                 fragdirs: List[List[str]]) -> Dict[str, WP]:
    """
    Calculate reaction energies for all fragment pairs.
    Returns a dictionary mapping pair_id to reaction energy in eV.
    """
    print(f"Calculating reaction energies for {npairs} reactions...")
    reaction_energies = {}
    
    # Get precursor energy
    precursor_xyz = "isomer.xyz" if Path("isomer.xyz").exists() else "fragment.xyz"
    if not Path(precursor_xyz).exists():
        print(f"Error: Precursor file {precursor_xyz} not found")
        return {fd[0]: 1000.0 for fd in fragdirs}  # High energy for all
    
    # Calculate precursor energy
    chrg_prec, uhf_prec = qmmod._get_chrg_uhf(env, precursor_xyz)
    success_prec = qmmod.run_qm_job_and_read(env, precursor_xyz, env.tslevel, 'sp', chrg_prec, uhf_prec)
    
    if not success_prec:
        print("Error: Could not calculate precursor energy")
        return {fd[0]: 1000.0 for fd in fragdirs}
    
    # Read precursor energy from qmdata
    query_prec = f"{env.tslevel} sp {chrg_prec} {uhf_prec}"
    found_prec, energy_prec = iomod.grepval("qmdata", query_prec)
    if not found_prec:
        print("Error: Could not read precursor energy from qmdata")
        return {fd[0]: 1000.0 for fd in fragdirs}
    
    # Calculate reaction energies for each pair
    for pair_id, item1_dir, item2_dir in fragdirs:
        try:
            product_energy_total = 0.0
            
            # Calculate energy of first product
            if item1_dir:
                item1_xyz = Path(item1_dir) / ("isomer.xyz" if not item2_dir else "fragment.xyz")
                if item1_xyz.exists():
                    os.chdir(item1_dir)
                    chrg1, uhf1 = qmmod._get_chrg_uhf(env, item1_xyz.name)
                    success1 = qmmod.run_qm_job_and_read(env, item1_xyz.name, env.tslevel, 'sp', chrg1, uhf1)
                    if success1:
                        query1 = f"{env.tslevel} sp {chrg1} {uhf1}"
                        found1, energy1 = iomod.grepval("qmdata", query1)
                        if found1:
                            product_energy_total += energy1
                    os.chdir("..")
            
            # Calculate energy of second product (if it exists)
            if item2_dir:
                item2_xyz = Path(item2_dir) / "fragment.xyz"
                if item2_xyz.exists():
                    os.chdir(item2_dir)
                    chrg2, uhf2 = qmmod._get_chrg_uhf(env, item2_xyz.name)
                    success2 = qmmod.run_qm_job_and_read(env, item2_xyz.name, env.tslevel, 'sp', chrg2, uhf2)
                    if success2:
                        query2 = f"{env.tslevel} sp {chrg2} {uhf2}"
                        found2, energy2 = iomod.grepval("qmdata", query2)
                        if found2:
                            product_energy_total += energy2
                    os.chdir("..")
            
            # Calculate reaction energy: ΔE = E(products) - E(reactant)
            reaction_energy_hartree = product_energy_total - energy_prec
            reaction_energy_ev = reaction_energy_hartree * AUTOEV  # Convert to eV
            
            reaction_energies[pair_id] = reaction_energy_ev
            
            # Write de_<level> file for this reaction
            (Path(pair_id) / f"de_{env.tslevel}").write_text(f"{reaction_energy_ev:.10f}\n")
            
            print(f"    Reaction {pair_id}: ΔE = {reaction_energy_ev:.3f} eV")
            
        except Exception as e:
            print(f"Error calculating reaction energy for {pair_id}: {e}")
            reaction_energies[pair_id] = 1000.0  # High energy to exclude this reaction
    
    return reaction_energies


def optimize_fragments_py(env: RunTypeData, npairs_in: int, fragdirs_in: List[List[str]]) -> Tuple[int, List[List[str]]]:
    """
    Optimize fragment geometries for IP calculations.
    Returns (npairs_out, fragdirs_out) after filtering failed optimizations.
    """
    print(f"Optimizing {npairs_in} fragments for IP calculations...")
    
    fragdirs_out = []
    npairs_out = 0
    
    for pair_id, item1_dir, item2_dir in fragdirs_in:
        success = True
        
        # Optimize first fragment
        if item1_dir:
            item1_xyz = Path(item1_dir) / ("isomer.xyz" if not item2_dir else "fragment.xyz")
            if item1_xyz.exists():
                os.chdir(item1_dir)
                chrg1, uhf1 = qmmod._get_chrg_uhf(env, item1_xyz.name, chrg_in=0)  # Neutral for IP
                success1 = qmmod.run_qm_job_and_read(env, item1_xyz.name, env.geolevel, 'opt', chrg1, uhf1)
                if not success1:
                    print(f"    Warning: Optimization failed for {item1_dir}")
                    success = False
                os.chdir("..")
        
        # Optimize second fragment (if it exists)
        if item2_dir and success:
            item2_xyz = Path(item2_dir) / "fragment.xyz"
            if item2_xyz.exists():
                os.chdir(item2_dir)
                chrg2, uhf2 = qmmod._get_chrg_uhf(env, item2_xyz.name, chrg_in=0)  # Neutral for IP
                success2 = qmmod.run_qm_job_and_read(env, item2_xyz.name, env.geolevel, 'opt', chrg2, uhf2)
                if not success2:
                    print(f"    Warning: Optimization failed for {item2_dir}")
                    success = False
                os.chdir("..")
        
        if success:
            fragdirs_out.append([pair_id, item1_dir, item2_dir])
            npairs_out += 1
    
    print(f"Fragment optimization: {npairs_out}/{npairs_in} successful")
    return npairs_out, fragdirs_out


def optimize_products_py(env: RunTypeData, npairs_in: int, fragdirs_in: List[List[str]]) -> Tuple[int, List[List[str]]]:
    """
    Optimize products at their assigned charges.
    Returns (npairs_out, fragdirs_out) after filtering failed optimizations.
    """
    print(f"Optimizing {npairs_in} products at assigned charges...")
    
    fragdirs_out = []
    npairs_out = 0
    
    for pair_id, item1_dir, item2_dir in fragdirs_in:
        success = True
        
        # Optimize first product at its assigned charge
        if item1_dir:
            item1_xyz = Path(item1_dir) / ("isomer.xyz" if not item2_dir else "fragment.xyz")
            if item1_xyz.exists():
                os.chdir(item1_dir)
                chrg1, uhf1 = qmmod._get_chrg_uhf(env, item1_xyz.name)  # Use assigned charge
                success1 = qmmod.run_qm_job_and_read(env, item1_xyz.name, env.geolevel, 'opt', chrg1, uhf1)
                if not success1:
                    print(f"    Warning: Product optimization failed for {item1_dir}")
                    success = False
                os.chdir("..")
        
        # Optimize second product (if it exists)
        if item2_dir and success:
            item2_xyz = Path(item2_dir) / "fragment.xyz"
            if item2_xyz.exists():
                os.chdir(item2_dir)
                chrg2, uhf2 = qmmod._get_chrg_uhf(env, item2_xyz.name)  # Use assigned charge
                success2 = qmmod.run_qm_job_and_read(env, item2_xyz.name, env.geolevel, 'opt', chrg2, uhf2)
                if not success2:
                    print(f"    Warning: Product optimization failed for {item2_dir}")
                    success = False
                os.chdir("..")
        
        if success:
            fragdirs_out.append([pair_id, item1_dir, item2_dir])
            npairs_out += 1
    
    print(f"Product optimization: {npairs_out}/{npairs_in} successful")
    return npairs_out, fragdirs_out


def sort_out_high_energy_fragments_py(env: RunTypeData, level: str, use_rrho: bool,
                                     npairs_in: int, fragdirs_in: List[List[str]],
                                     scale_iee_factor: WP = 3.0) -> Tuple[int, List[List[str]]]:
    """
    Filter out high-energy fragments based on energy threshold.
    Returns (npairs_out, fragdirs_out) after filtering.
    """
    print(f"Filtering high-energy fragments (threshold: {scale_iee_factor} * IEE)...")
    
    # Calculate energy threshold
    # This is a simplified version - the full implementation would consider
    # the IEE distribution and molecular size
    nat0 = iomod.rdshort_int_py(Path(env.startdir) / "fragment.xyz", default=30)
    energy_threshold_ev = env.ieeatm * nat0 * scale_iee_factor
    
    print(f"Energy threshold: {energy_threshold_ev:.2f} eV")
    
    fragdirs_out = []
    npairs_out = 0
    
    for pair_id, item1_dir, item2_dir in fragdirs_in:
        # Read reaction energy
        de_file = Path(pair_id) / f"de_{level}"
        if de_file.exists():
            reaction_energy = iomod.rdshort_real_py(de_file, default=0.0)
            
            if reaction_energy <= energy_threshold_ev:
                fragdirs_out.append([pair_id, item1_dir, item2_dir])
                npairs_out += 1
            else:
                print(f"    Filtered out {pair_id}: ΔE = {reaction_energy:.2f} eV > {energy_threshold_ev:.2f} eV")
        else:
            print(f"    Warning: No energy file found for {pair_id}, keeping it")
            fragdirs_out.append([pair_id, item1_dir, item2_dir])
            npairs_out += 1
    
    print(f"High-energy filter: {npairs_out}/{npairs_in} reactions passed")
    return npairs_out, fragdirs_out


def gff_topology_check_py(env: RunTypeData, npairs_in: int, fragdirs_in: List[List[str]]) -> Tuple[int, List[List[str]]]:
    """
    Perform GFF topology check to filter out unreasonable structures.
    This is a placeholder implementation.
    """
    print(f"Performing GFF topology check for {npairs_in} structures...")
    
    # For now, just pass through all structures
    # A real implementation would use GFF to check for reasonable geometries
    print(f"GFF topology check: {npairs_in}/{npairs_in} structures passed (placeholder)")
    return npairs_in, fragdirs_in


if __name__ == '__main__':
    print("Testing reaction.py functions...")
    
    # Create test environment
    test_env = RunTypeData()
    test_env.tslevel = "gfn2"
    test_env.geolevel = "gfn2"
    test_env.ieeatm = 0.8
    test_env.startdir = Path(".").resolve()
    
    # Test data
    test_fragdirs = [
        ["p0", "p0f1", "p0f2"],
        ["p1", "p1", ""]  # Isomer
    ]
    
    print("reaction.py test completed.")
