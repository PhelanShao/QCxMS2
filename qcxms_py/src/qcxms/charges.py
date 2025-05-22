import math
from pathlib import Path
from typing import Tuple, Optional, List, Dict, Union

try:
    from .data import RunTypeData, WP
    from . import iomod
    from . import qmmod
    from . import utility
except ImportError:
    # Fallbacks for standalone/testing
    print("Attempting to import dummy/mock modules for charges.py standalone run.")
    from data import RunTypeData, WP # type: ignore
    import iomod_mock as iomod # type: ignore
    import qmmod_mock as qmmod # type: ignore
    import utility_mock as utility # type: ignore

# Constants (ensure these are consistent with the project's definitions)
HARTREE_TO_EV = 27.211386245988  # Conversion factor for Hartrees to eV
KB_EV_K = 8.6173332621415e-5   # Boltzmann constant in eV/K


def _calculate_ip_for_fragment(
    env: RunTypeData,
    fragment_xyz_path: Path, # Absolute path to fragment.xyz
    level: str,
    original_charge: int # The charge of the system these fragments came from (env.chrg)
) -> Tuple[Optional[WP], bool]:
    """
    Calculates Ionization Potential for a single fragment.
    Manages QM calculations for neutral and charged species.
    Returns (IP_in_eV, ip_ok_bool).
    This function operates in the CWD, which should be the fragment's specific directory.
    """
    ip_value_ev: Optional[WP] = None
    ip_ok = False
    
    # --- Neutral Species Calculation (Charge 0) ---
    # Path for neutral calculation (e.g., ./charge0/)
    neutral_calc_dir = Path("charge0")
    neutral_calc_dir.mkdir(exist_ok=True)
    # Copy fragment XYZ to neutral calc dir
    neutral_xyz_path = neutral_calc_dir / fragment_xyz_path.name
    iomod.copy_file(fragment_xyz_path, neutral_xyz_path)
    
    cwd_orig = Path.cwd()
    os.chdir(neutral_calc_dir) # Change to neutral directory

    # Prep QM for neutral
    cmd_neut, out_neut, patt_neut, clean_neut, cached_neut, val_neut = \
        qmmod.prepqm_py(env, neutral_xyz_path.name, level, 'sp', chrg_in=0, uhf_in=None, restart=False) # UHF will be auto-determined for neutral

    e_neutral_au: Optional[WP] = None
    failed_neutral = False

    if cached_neut and val_neut is not None:
        e_neutral_au = val_neut
    elif cmd_neut:
        result_neut = iomod.execute_command(cmd_neut, shell=True if isinstance(cmd_neut, str) else False) # shell=True for string command
        if result_neut.returncode == 0:
            e_neutral_au, failed_neutral = qmmod.readoutqm_py(env, neutral_xyz_path.name, level, 'sp', out_neut, patt_neut)
        else:
            failed_neutral = True
        if clean_neut: iomod.execute_command(clean_neut.split(), shell=True) # Cleanup uses shell
    else: # prepqm returned no command and not cached (error)
        failed_neutral = True
    
    os.chdir(cwd_orig) # Return to fragment's main directory
    if failed_neutral or e_neutral_au is None:
        print(f"Error: Neutral SP calculation failed for {fragment_xyz_path.parent.name}")
        return None, False

    # --- Charged Species Calculation (Charge original_charge) ---
    # This happens in the fragment's main directory (current CWD)
    cmd_chrg, out_chrg, patt_chrg, clean_chrg, cached_chrg, val_chrg = \
        qmmod.prepqm_py(env, fragment_xyz_path.name, level, 'sp', chrg_in=original_charge, uhf_in=None, restart=False)

    e_charged_au: Optional[WP] = None
    failed_charged = False

    if cached_chrg and val_chrg is not None:
        e_charged_au = val_chrg
    elif cmd_chrg:
        result_chrg = iomod.execute_command(cmd_chrg, shell=True if isinstance(cmd_chrg, str) else False)
        if result_chrg.returncode == 0:
            e_charged_au, failed_charged = qmmod.readoutqm_py(env, fragment_xyz_path.name, level, 'sp', out_chrg, patt_chrg)
        else:
            failed_charged = True
        if clean_chrg: iomod.execute_command(clean_chrg.split(), shell=True)
    else:
        failed_charged = True

    if failed_charged or e_charged_au is None:
        print(f"Error: Charged SP calculation failed for {fragment_xyz_path.parent.name}")
        return None, False

    # Calculate IP
    ip_value_au = e_charged_au - e_neutral_au
    ip_value_ev = ip_value_au * HARTREE_TO_EV

    # Validity check (e.g., positive IP for positive ion mode)
    if original_charge >= 0 and ip_value_ev < 0:
        print(f"Warning: Negative IP ({ip_value_ev:.2f} eV) calculated for {fragment_xyz_path.parent.name} with charge {original_charge}. Treating as problematic.")
        ip_ok = False 
        # In Fortran, this was an error that stopped processing for the pair.
        # Here, returning False allows caller to handle.
    else:
        ip_ok = True
    
    # Write IP to file (ip_<level>)
    iomod.wrshort_real(f"ip_{level}", ip_value_ev)
    
    return ip_value_ev, ip_ok


def _boltzmann_distribution_py(
    env: RunTypeData, 
    ip1_ev: WP, 
    ip2_ev: WP, 
    pair_base_dir_str: str # e.g. "pair_1" -> writes to "pair_1f1/statchrg", "pair_1f2/statchrg"
):
    """Calculates Boltzmann distribution for charge between two fragments."""
    temp_k = float(env.temp)
    # kB * T in eV units
    factor_kT_ev = KB_EV_K * temp_k
    if factor_kT_ev == 0: # Avoid division by zero if T=0
        statchrg1 = 1.0 if ip1_ev <= ip2_ev else 0.0
        statchrg2 = 1.0 - statchrg1
    else:
        try:
            # Shift IPs to prevent overflow with large positive IPs (exp becomes tiny)
            # And to handle large negative IPs (exp becomes huge)
            # Smallest IP will be reference (0)
            min_ip = min(ip1_ev, ip2_ev)
            shifted_ip1 = ip1_ev - min_ip
            shifted_ip2 = ip2_ev - min_ip

            q1 = math.exp(-shifted_ip1 / factor_kT_ev)
            q2 = math.exp(-shifted_ip2 / factor_kT_ev)
            
            sum_q = q1 + q2
            if sum_q == 0: # Both probabilities are zero (e.g. very high IPs relative to kT)
                # Default to lower IP getting the charge, or 50/50 if equal
                statchrg1 = 0.5 if ip1_ev == ip2_ev else (1.0 if ip1_ev < ip2_ev else 0.0)
            else:
                statchrg1 = q1 / sum_q
        except OverflowError: # Should be rare with IP shifting
            statchrg1 = 0.5 if ip1_ev == ip2_ev else (1.0 if ip1_ev < ip2_ev else 0.0)
            print(f"Warning: Overflow during Boltzmann calculation for IPs {ip1_ev}, {ip2_ev}. Defaulting charge assignment.")

        statchrg2 = 1.0 - statchrg1

    # Write to statchrg files in respective fragment subdirectories (f1, f2)
    # Assuming pair_base_dir_str is like "p1" and fragdirs are "p1/f1", "p1/f2"
    # This needs careful path handling. Fortran was "dir//f1"
    # If pair_base_dir_str is "p1", then it wrote to "p1f1/statchrg"
    # This implies fragdirs_in[i,1] was the base for the pair.
    
    # The Fortran code wrote to fragdirs_in(i,1) + "f1" / statchrg
    # fragdirs_in(i,1) is the path to the pair directory, e.g. "calc/run1/p01"
    # Fragment 1 data is in fragdirs_in(i,2), e.g. "calc/run1/p01/f1"
    # Fragment 2 data is in fragdirs_in(i,3), e.g. "calc/run1/p01/f2"
    # The statchrg is written *into* the fragment data directories.

    frag1_data_dir = Path(pair_base_dir_str) / "f1" # This assumes pair_base_dir_str is like "pair_results/p1"
                                                # And f1, f2 are subdirs of that.
                                                # This is inconsistent with how fragdirs_in was used for IP calcs (where fragdirs_in(i,2) was the full path to f1).
                                                # Let's assume pair_base_dir_str is the *common parent path* for f1 and f2.
                                                # And fragdirs_in(i,2) and fragdirs_in(i,3) are full paths to f1 and f2.
                                                # The statchrg should be written *inside* fragdirs_in(i,2) and fragdirs_in(i,3).
                                                # The Fortran `call chdir(trim(tmpdir))` where `tmpdir = trim(dir)//'f1'` suggests `dir` is the pair identifier (e.g. "p0")
                                                # and the CWD is the main run directory.
                                                # So, if fragdirs_in(i,2) is "p0/f1", statchrg goes into "p0/f1/statchrg".

    # Corrected logic: pair_base_dir_str argument for _boltzmann_py should actually be fragdirs_in[i,1] (the pair identifier string)
    # and the actual writing should use the full paths fragdirs_in(i,2) and fragdirs_in(i,3)
    # For now, this function will just calculate. Writing will be done by caller.
    # No, Fortran `boltzmann` takes `dir` which is `fragdirs_in(i,1)` (e.g., "p0") and writes to "p0f1/statchrg" and "p0f2/statchrg"
    # This means the `f1` and `f2` are subdirectories of the path indicated by `dir`.
    # This implies `fragdirs_in(i,2)` and `fragdirs_in(i,3)` are actually `dir + '/f1'` and `dir + '/f2'`.

    # Let's adjust to take the actual fragment data directories
    # This function will be called with frag1_dir_path and frag2_dir_path
    # No, the Fortran code writes to subdirs of `dir` argument.
    # `dir` is `fragdirs_in(i,1)`. So it writes to `fragdirs_in(i,1)/f1/statchrg` and `fragdirs_in(i,1)/f2/statchrg`.
    # This means `fragdirs_in(i,2)` and `fragdirs_in(i,3)` are NOT the paths used for `statchrg` writing.
    # This is confusing. Let's assume `pair_base_dir_str` is the directory that CONTAINS f1 and f2.
    
    # Re-evaluating Fortran: `write (tmpdir, '(a)') trim(dir)//'f1'` then `chdir(tmpdir)`.
    # So, `dir` is the parent of `f1` and `f2`.
    # If `dir` = "pair01", then it writes to "pair01/f1/statchrg".
    # This makes `pair_base_dir_str` the correct argument name.

    frag1_statchrg_path = Path(pair_base_dir_str) / "f1" / "statchrg"
    frag2_statchrg_path = Path(pair_base_dir_str) / "f2" / "statchrg"
    
    frag1_statchrg_path.parent.mkdir(parents=True, exist_ok=True)
    frag2_statchrg_path.parent.mkdir(parents=True, exist_ok=True)

    iomod.wrshort_real(frag1_statchrg_path, statchrg1)
    iomod.wrshort_real(frag2_statchrg_path, statchrg2)
    

def assign_charges_py(env: RunTypeData, npairs_in: int, fragdirs_in: List[List[str]]) -> Tuple[int, List[List[str]]]:
    """
    Assigns statistical charges to fragment pairs based on Ionization Potentials.
    fragdirs_in: list of [pair_ID_str, fullpath_to_frag1_dir, fullpath_to_frag2_dir]
    """
    print("Assigning statistical charges via Delta SCF procedure...")
    original_path = Path.cwd()
    ip_threshold_ev = 2.0  # eV, for refining with ip2level

    fragpairs_data = [] # Store (pair_id, frag1_dir, frag2_dir, ip1, ip2, ip1_ok, ip2_ok)

    # --- Stage 1: Calculate IPs at env.iplevel ---
    print(f"Calculating IPs at level: {env.iplevel}")
    for i in range(npairs_in):
        pair_id_str, frag1_dir_str, frag2_dir_str = fragdirs_in[i]
        if not (frag1_dir_str and frag2_dir_str): # Skip if not a full pair
            fragpairs_data.append((pair_id_str, frag1_dir_str, frag2_dir_str, None, None, False, False))
            continue

        ips_for_pair = [None, None] # [IP_frag1, IP_frag2]
        ips_ok_for_pair = [False, False]

        for frag_idx, frag_dir_abs_str in enumerate([frag1_dir_str, frag2_dir_str]):
            frag_dir_path = Path(frag_dir_abs_str)
            frag_xyz_file = frag_dir_path / "fragment.xyz" # Standard name

            # Check if IP already calculated and stored for this level
            ip_file = frag_dir_path / f"ip_{env.iplevel}"
            if env.restartrun and ip_file.exists():
                # print(f"IP for {frag_dir_path.name} at {env.iplevel} already exists (restart mode). Reading.")
                ips_for_pair[frag_idx] = iomod.rdshort_real(ip_file)
                # Assume it was OK if it exists in restart mode, or re-evaluate based on value
                ips_ok_for_pair[frag_idx] = True if ips_for_pair[frag_idx] is not None and (env.chrg < 0 or ips_for_pair[frag_idx] >= 0) else False # Basic check
                if not ips_ok_for_pair[frag_idx]: print(f"Warning: Pre-existing IP for {frag_dir_path.name} at {env.iplevel} is invalid ({ips_for_pair[frag_idx]}). Will attempt recompute.")
                if ips_ok_for_pair[frag_idx]: continue # Skip recompute

            print(f"  Processing IP for fragment: {frag_dir_path.name} (Pair {pair_id_str})")
            os.chdir(frag_dir_path) # Change to fragment directory for QM calculation
            ip_val, ip_ok = _calculate_ip_for_fragment(env, frag_xyz_file, env.iplevel, env.chrg)
            os.chdir(original_path) # Change back
            
            ips_for_pair[frag_idx] = ip_val
            ips_ok_for_pair[frag_idx] = ip_ok
        
        fragpairs_data.append((pair_id_str, frag1_dir_str, frag2_dir_str, 
                               ips_for_pair[0], ips_for_pair[1], 
                               ips_ok_for_pair[0], ips_ok_for_pair[1]))

    # --- Stage 2: Refine IPs with env.ip2level if specified and IPs are close ---
    if env.ip2level and env.ip2level.strip() and env.ip2level.strip() != env.iplevel.strip():
        print(f"\nRefining close IPs at level: {env.ip2level}")
        for i in range(len(fragpairs_data)):
            pair_id, f1d, f2d, ip1, ip2, ip1_ok, ip2_ok = fragpairs_data[i]
            if not (ip1_ok and ip2_ok and ip1 is not None and ip2 is not None):
                continue # Cannot refine if initial IPs are not OK

            if abs(ip1 - ip2) < ip_threshold_ev:
                print(f"  IPs for pair {pair_id} are close ({ip1:.2f} eV, {ip2:.2f} eV). Refining with {env.ip2level}.")
                refined_ips = [ip1, ip2] # Start with previous IPs as fallback
                refined_ips_ok = [True, True]

                for frag_idx, frag_dir_abs_str in enumerate([f1d, f2d]):
                    frag_dir_path = Path(frag_dir_abs_str)
                    frag_xyz_file = frag_dir_path / "fragment.xyz"
                    
                    # Check cache for ip2level
                    ip_file_l2 = frag_dir_path / f"ip_{env.ip2level}"
                    if env.restartrun and ip_file_l2.exists():
                        # print(f"IP for {frag_dir_path.name} at {env.ip2level} already exists (restart mode). Reading.")
                        refined_ips[frag_idx] = iomod.rdshort_real(ip_file_l2)
                        refined_ips_ok[frag_idx] = True if refined_ips[frag_idx] is not None and (env.chrg < 0 or refined_ips[frag_idx] >= 0) else False
                        if not refined_ips_ok[frag_idx]: print(f"Warning: Pre-existing IP for {frag_dir_path.name} at {env.ip2level} is invalid ({refined_ips[frag_idx]}). Will attempt recompute.")
                        if refined_ips_ok[frag_idx]: continue

                    print(f"    Refining IP for fragment: {frag_dir_path.name} (Pair {pair_id})")
                    os.chdir(frag_dir_path)
                    ip_val_l2, ip_ok_l2 = _calculate_ip_for_fragment(env, frag_xyz_file, env.ip2level, env.chrg)
                    os.chdir(original_path)

                    if ip_ok_l2 and ip_val_l2 is not None:
                        refined_ips[frag_idx] = ip_val_l2
                    else: # Refinement failed, keep original IP for this fragment
                        print(f"    Refinement failed for {frag_dir_path.name}, using IP from {env.iplevel}.")
                        refined_ips_ok[frag_idx] = False # Mark refinement as failed for this one
                
                # Update with refined IPs, falling back to original if refinement failed for one
                fragpairs_data[i] = (pair_id, f1d, f2d, 
                                     refined_ips[0], refined_ips[1], 
                                     refined_ips_ok[0], refined_ips_ok[1])


    # --- Stage 3: Calculate Boltzmann distribution and write final charges ---
    print("\nCalculating Boltzmann charge distributions and assigning final charges...")
    valid_fragdirs_out_list: List[List[str]] = []

    for pair_id, f1d_str, f2d_str, ip1_ev, ip2_ev, ip1_ok, ip2_ok in fragpairs_data:
        if not (f1d_str and f2d_str): # Not a pair
            if f1d_str: # Single fragment (e.g. from initial molecule if no dissociation)
                 valid_fragdirs_out_list.append([pair_id, f1d_str, ""]) # Keep it
            continue

        if ip1_ok and ip2_ok and ip1_ev is not None and ip2_ev is not None:
            # The `pair_id_str` (fragdirs_in[i,1]) is the base directory for f1, f2 subdirs in Fortran's boltzmann.
            # Example: if pair_id_str is "p0", it writes to "p0/f1/statchrg", "p0/f2/statchrg".
            # This implies f1d_str is "p0/f1" and f2d_str is "p0/f2".
            # The _boltzmann_py function needs the parent directory of f1 and f2.
            common_parent_for_pair = Path(f1d_str).parent 
            _boltzmann_py(env, ip1_ev, ip2_ev, str(common_parent_for_pair))
            
            # Read statchrg and assign final .CHRG
            statchrg1_path = Path(f1d_str) / "statchrg"
            statchrg1 = iomod.rdshort_real(statchrg1_path, default=0.5) # Default to 0.5 if file missing

            final_chrg1, final_chrg2 = (env.chrg, 0) if statchrg1 >= 0.5 else (0, env.chrg)
            
            iomod.wrshort_int(Path(f1d_str) / ".CHRG", final_chrg1)
            iomod.remove_file(Path(f1d_str) / ".UHF")
            iomod.wrshort_int(Path(f2d_str) / ".CHRG", final_chrg2)
            iomod.remove_file(Path(f2d_str) / ".UHF")
            
            valid_fragdirs_out_list.append([pair_id, f1d_str, f2d_str])
        else:
            print(f"  Skipping charge assignment for pair {pair_id} due to failed/invalid IP calculation.")
            if env.removedirs:
                print(f"    Removing directories for failed pair {pair_id}: {f1d_str}, {f2d_str}")
                if f1d_str: iomod.rmrf(f1d_str)
                if f2d_str: iomod.rmrf(f2d_str)
                # Also remove the parent pair directory if it's now empty or only contains f1/f2
                # This logic is complex and depends on directory structure; for now, only removing f1/f2 dirs.
                # common_parent = Path(f1d_str).parent
                # if common_parent.exists() and not any(common_parent.iterdir()):
                #     iomod.rmrf(common_parent)


    npairs_out = len(valid_fragdirs_out_list)
    print(f"Charge assignment complete. {npairs_out} pairs processed successfully.")
    return npairs_out, valid_fragdirs_out_list


if __name__ == '__main__':
    # This block is for testing the charges.py module functionality.
    # It would require setting up a dummy RunTypeData, mock qmmod/iomod/utility functions,
    # and dummy file structures.
    print("Testing charges.py functions...")

    # Dummy RunTypeData
    env_test = RunTypeData()
    env_test.temp = 298 # Kelvin
    env_test.chrg = 1
    env_test.iplevel = "gfn2"
    env_test.ip2level = "r2scan3c" # Test refinement
    # env_test.ip2level = "" # Test no refinement
    env_test.restartrun = False
    env_test.removedirs = True

    # Mock file system
    Path("test_run").mkdir(exist_ok=True)
    os.chdir("test_run")
    
    pair1_id = "p0"
    pair1_f1_dir = Path(pair1_id) / "f1"
    pair1_f2_dir = Path(pair1_id) / "f2"
    pair1_f1_dir.mkdir(parents=True, exist_ok=True)
    pair1_f2_dir.mkdir(parents=True, exist_ok=True)
    (pair1_f1_dir / "fragment.xyz").write_text("1\nH\nH 0 0 0")
    (pair1_f2_dir / "fragment.xyz").write_text("1\nHe\nHe 0 0 1")

    pair2_id = "p1" # This one will simulate failed IP
    pair2_f1_dir = Path(pair2_id) / "f1"
    pair2_f2_dir = Path(pair2_id) / "f2"
    pair2_f1_dir.mkdir(parents=True, exist_ok=True)
    pair2_f2_dir.mkdir(parents=True, exist_ok=True)
    (pair2_f1_dir / "fragment.xyz").write_text("1\nLi\nLi 0 0 0")
    (pair2_f2_dir / "fragment.xyz").write_text("1\nBe\nBe 0 0 1")


    fragdirs_in_test = [
        [pair1_id, str(pair1_f1_dir.resolve()), str(pair1_f2_dir.resolve())],
        [pair2_id, str(pair2_f1_dir.resolve()), str(pair2_f2_dir.resolve())]
    ]
    npairs_in_test = len(fragdirs_in_test)

    # Mocking _calculate_ip_for_fragment behavior for the test
    def mock_calc_ip(env_mock, frag_xyz_path_mock, level_mock, chrg_mock) -> Tuple[Optional[WP], bool]:
        # Simulate IP calculation based on fragment name or path
        # print(f"Mock IP calc for: {frag_xyz_path_mock.parent.name} level {level_mock}")
        if "p0" in str(frag_xyz_path_mock) and "f1" in str(frag_xyz_path_mock): # H
            val = 13.6 if level_mock == "gfn2" else 13.8 
            ok = True
        elif "p0" in str(frag_xyz_path_mock) and "f2" in str(frag_xyz_path_mock): # He
            val = 24.5 if level_mock == "gfn2" else 24.2 # Make it close for p0 to test refinement
            ok = True
        elif "p1" in str(frag_xyz_path_mock) and "f1" in str(frag_xyz_path_mock): # Li
            val = 5.0
            ok = True
        elif "p1" in str(frag_xyz_path_mock) and "f2" in str(frag_xyz_path_mock): # Be - simulate failure
            val = None
            ok = False
        else:
            val = 10.0; ok = True
        
        if ok and val is not None:
            iomod.wrshort_real(Path(f"ip_{level_mock}"), val)
        return val, ok

    # Replace the actual function with the mock for the duration of the test
    original_calc_ip = _calculate_ip_for_fragment
    _calculate_ip_for_fragment = mock_calc_ip # type: ignore

    npairs_out_test, fragdirs_out_test = assign_charges_py(env_test, npairs_in_test, fragdirs_in_test)

    print(f"\nTest Results: Npairs_out = {npairs_out_test}")
    for i, p_data in enumerate(fragdirs_out_test):
        print(f"Pair {i}: ID={p_data[0]}")
        if p_data[1]: print(f"  Frag1 CHRG: {iomod.rdshort_int(Path(p_data[1]) / '.CHRG', -99)}")
        if p_data[2]: print(f"  Frag2 CHRG: {iomod.rdshort_int(Path(p_data[2]) / '.CHRG', -99)}")
        # Check statchrg files if they exist
        statchrg1_path = Path(p_data[0]) / "f1" / "statchrg" # Path construction needs care
        statchrg2_path = Path(p_data[0]) / "f2" / "statchrg"
        # This path construction for statchrg is based on _boltzmann_py, assuming p_data[0] is the PAIR_ID string.
        # However, _boltzmann_py was called with Path(f1d_str).parent.
        # So, if f1d_str = "test_run/p0/f1", parent is "test_run/p0".
        # Then _boltzmann_py writes to "test_run/p0/f1/statchrg".
        # This means statchrg files are inside the fragment directories themselves.
        statchrg1_path_corrected = Path(p_data[1]) / "statchrg" if p_data[1] else None
        statchrg2_path_corrected = Path(p_data[2]) / "statchrg" if p_data[2] else None

        if statchrg1_path_corrected and statchrg1_path_corrected.exists():
            print(f"  Frag1 statchrg: {iomod.rdshort_real(statchrg1_path_corrected):.3f}")
        if statchrg2_path_corrected and statchrg2_path_corrected.exists():
            print(f"  Frag2 statchrg: {iomod.rdshort_real(statchrg2_path_corrected):.3f}")


    # Restore original function
    _calculate_ip_for_fragment = original_calc_ip

    os.chdir("..")
    if Path("test_run").exists():
        shutil.rmtree("test_run")
    print("Test run finished and cleaned up.")

```
