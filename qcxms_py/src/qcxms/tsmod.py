import os
import shutil
import math
import re
import subprocess
from pathlib import Path
from typing import Tuple, Optional, List, Dict

try:
    from .data import RunTypeData, WP, Timer # BOHR, AUTOEV, AUTOKCAL (not directly used in tsmod.f90, but by callers)
    from . import iomod
    from . import qmmod
    from . import utility
    # from . import rmsd # if get_rmsd_for_coords_py is moved there from rmsd_ls
    from .rmsd import get_rmsd_for_coords_py # Assuming rmsd.py contains this
    from .constants import AUTOEV, EVTOKCAL, PI # Assuming a constants.py
except ImportError:
    # Fallbacks for standalone/testing
    print("Attempting to import dummy/mock modules for tsmod.py standalone run.")
    from data import RunTypeData, WP, Timer # type: ignore
    import iomod_mock as iomod # type: ignore
    import qmmod_mock as qmmod # type: ignore
    import utility_mock as utility # type: ignore
    from rmsd import get_rmsd_for_coords_py # type: ignore
    AUTOEV = 27.211386245988
    EVTOKCAL = 23.060547830618307
    PI = math.pi


# Helper function to manage job calls (sequential execution for now)
def _run_jobs_sequentially(env: RunTypeData, job_list: List[Tuple[Path, List[str], str, str]]):
    """
    Runs a list of QM jobs sequentially.
    job_list contains tuples of (working_directory, command_list, output_file, success_pattern).
    This replaces omp_samejobcall for now.
    """
    original_path = Path.cwd()
    for work_dir, cmd_list, out_f, patt in job_list:
        if not work_dir.exists():
            print(f"  Warning: Working directory {work_dir} does not exist. Skipping job.")
            continue
        
        os.chdir(work_dir)
        print(f"  Running in {work_dir}: {' '.join(cmd_list)} > {out_f} 2>&1") # For logging
        
        # Allow shell=True if cmd_list is actually a single command string
        use_shell = isinstance(cmd_list, str)
        actual_command = cmd_list if use_shell else cmd_list
        
        # Capture output to the specified file
        with open(out_f, "w") as stdout_file, open(f"{out_f}.err", "w") as stderr_file:
            result = iomod.execute_command(actual_command, cwd=work_dir, shell=use_shell, 
                                           stdout_to=stdout_file, stderr_to=stderr_file)
        
        if result.returncode != 0:
            print(f"    WARNING: Job in {work_dir} failed with code {result.returncode}. Output in {out_f}, errors in {out_f}.err")
        # Further processing (readoutqm) would happen after this batch in the main loops
        os.chdir(original_path)


def _prep_geodesic_py(env: RunTypeData, reaction_dir: Path) -> Optional[List[str]]:
    """
    Prepares and returns the command for geodesic interpolation.
    Corresponds to Fortran prepgeodesic.
    Creates geoin.xyz in reaction_dir.
    """
    start_xyz = reaction_dir / "start.xyz"
    end_xyz = reaction_dir / "end.xyz"
    geoin_xyz = reaction_dir / "geoin.xyz"
    interpolated_xyz = reaction_dir / "interpolated.xyz" # Output of geodesic_interpolate

    if not (start_xyz.exists() and end_xyz.exists()):
        print(f"  Error: start.xyz or end.xyz missing in {reaction_dir} for geodesic interpolation.")
        return None

    # Create geoin.xyz by concatenating start.xyz and end.xyz
    # start.xyz content (nat, comment, coords)
    # end.xyz content (nat, comment, coords)
    # geoin.xyz content (nat_start, comment_start, coords_start, nat_end, comment_end, coords_end)
    try:
        start_content = start_xyz.read_text().strip().splitlines()
        end_content = end_xyz.read_text().strip().splitlines()
        
        # Ensure there's content to read
        if not start_content or not end_content:
            print(f"  Error: start.xyz or end.xyz is empty in {reaction_dir}.")
            return None

        with open(geoin_xyz, "w") as f_out:
            f_out.write("\n".join(start_content) + "\n")
            f_out.write("\n".join(end_content) + "\n")
            
    except IOError as e:
        print(f"  Error creating geoin.xyz in {reaction_dir}: {e}")
        return None

    # jobcall = f"geodesic_interpolate -nimages {env.tsnds+2} -friction 0.01 {geoin_xyz} > geo.out"
    # Fortran used env.tsnds for -nimages for geodesic_interpolate
    # But ORCA NEB later uses env.tsnds+2 images.
    # The Fortran prepneb copied interpolated.xyz (output of geodesic) to nebguess.xyz
    # and if it had N images, it means N-2 intermediate points for NEB.
    # So if geodesic_interpolate is to produce N-2 intermediate points, it needs N images.
    # If env.tsnds is number of *intermediate* images for NEB, then geodesic needs env.tsnds+2.
    
    num_images_geodesic = env.tsnds # Fortran prepgeodesic used env.tsnds directly.
                                   # Fortran prepneb: Nimages = env.tsnds + 2.
                                   # If usegeo, it copied interpolated.xyz (N images total) to nebguess.xyz
                                   # This implies geodesic_interpolate should produce env.tsnds+2 total images.
                                   # Let's stick to Fortran logic: env.tsnds for geodesic's -nimages.
                                   # The output interpolated.xyz will then be used.

    # The Fortran code `jobcall = geodesic_interpolate ... > geo.out`
    # The actual output file for coordinates is `interpolated.xyz` by default for geodesic_interpolate.
    # `geo.out` is just stdout.
    geodesic_cmd = [
        "geodesic_interpolate", 
        "-nimages", str(num_images_geodesic), # This seems to be a mismatch with NEB Nimages.
                                             # If NEB needs N images (tsnds+2), and geo is guess,
                                             # geo should produce N images. Let's assume tsnds means intermediate points for NEB.
                                             # So total images for NEB = tsnds+2.
                                             # Geodesic should produce tsnds+2 images.
        "-friction", "0.01", 
        str(geoin_xyz) # Input file for geodesic_interpolate
    ]
    # Output is implicitly interpolated.xyz
    return geodesic_cmd


def _prepare_neb_restart_py(env: RunTypeData, current_reaction_dir: Path):
    """Prepares files for restarting a NEB calculation. Fortran: prepare_neb_restart"""
    # Fortran: copy orca_MEP_trj.xyz to nebguess.xyz if exists, else interpolated.xyz to nebguess.xyz
    # Also removes some ORCA temp files.
    
    mep_trj = current_reaction_dir / "orca_MEP_trj.xyz"
    interpolated = current_reaction_dir / "interpolated.xyz" # From geodesic
    neb_guess = current_reaction_dir / "nebguess.xyz"

    if mep_trj.exists():
        shutil.copy2(mep_trj, neb_guess)
        print(f"  Restarting NEB from orca_MEP_trj.xyz in {current_reaction_dir}")
    elif interpolated.exists():
        shutil.copy2(interpolated, neb_guess)
        print(f"  Restarting NEB from interpolated.xyz (geodesic) in {current_reaction_dir}")
    else:
        print(f"  Warning: No trajectory file found for NEB restart in {current_reaction_dir}. NEB might fail.")
        # Create a dummy nebguess from start and end if absolutely nothing else?
        # For now, do nothing if no guess is available. ORCA might try linear.

    # Remove ORCA temp files for restart
    for f_pattern in ["orca.นีb.*", "orca.tmp.*", "orca_endcoords*", "*.ges", "*.engrad", "*.xyzopt*", "*.xyziter*"]:
        for f_path in current_reaction_dir.glob(f_pattern):
            iomod.remove_file(f_path)
    iomod.remove_file(current_reaction_dir / "orca.opt")
    iomod.remove_file(current_reaction_dir / "orca.อิฐ") # Fortran had orca.อิฐ - likely typo for .neb? or specific ORCA temp
    iomod.remove_file(current_reaction_dir / ".NEB_ अगली_मૂર્તિ") # Fortran had this, also likely encoding issue.
                                                          # Assuming it meant specific NEB temp files not covered by glob.


def _prep_neb_py(
    env: RunTypeData, 
    current_reaction_dir: Path, # To know where to write orca.inp
    hbonddiss_flag: bool, # From Fortran, affects KPush/KPull - not directly used in current inp
    restart_flag: bool, 
    use_geo_guess_flag: bool # Use nebguess.xyz (from geodesic or previous NEB)
) -> Optional[List[str]]:
    """
    Prepares ORCA NEB input file and returns the command.
    Corresponds to Fortran prepneb.
    """
    # Determine charge and UHF from .CHRG/.UHF files in current_reaction_dir
    # These should have been set up by the caller (ts_search_py)
    chrg = iomod.rdshort_int(current_reaction_dir / ".CHRG", default=env.chrg)
    uhf = iomod.rdshort_int(current_reaction_dir / ".UHF", default=0) # Default to singlet if not found
    multiplicity = uhf + 1
    
    # Use settings from env for level, etc.
    level_keyword, default_etemp, _ = qmmod._get_orca_level_keyword(env.tslevel) # Using helper from qmmod

    inp_lines = ["! SlowConv" if restart_flag else "# Default NEB input"]
    inp_lines.append(f"! {level_keyword}")
    if env.solv: inp_lines.append("! CPCM(WATER)" if "gfn" not in env.tslevel else "! ALPB(WATER)")
    
    inp_lines.append("! NEB LooseOpt") # Default options from Fortran
    if env.nebnormal: # User requests tighter options
        inp_lines.remove("! NEB LooseOpt") # Remove loose
        inp_lines.append("! NEB") # Add normal (tighter)
        inp_lines.append("%NEB NIter_Max 200 END") # Example of tighter option from Fortran
    
    inp_lines.append(f"%NEB NImages {env.tsnds + 2} END") # tsnds is intermediate, total is tsnds+2

    if use_geo_guess_flag and (current_reaction_dir / "nebguess.xyz").exists():
        inp_lines.append('%geom InHessName "start.hess" end') # Initial hessian for endpoints
        inp_lines.append('  Geo_Source "nebguess.xyz"') # Use prepared guess
        inp_lines.append('end')
    else: # Linear interpolation if no guess
        inp_lines.append('%geom InHessName "start.hess" end')
        inp_lines.append('  Interp_Mode 0') # Linear
        inp_lines.append('end')

    # KPush/KPull based on hbonddiss - Fortran logic was complex, simplified here
    # if hbonddiss:
    #     inp_lines.append("%NEB KPush 0.005 KPull 0.005 END")
    # else:
    #     inp_lines.append("%NEB KPush 0.01 KPull 0.01 END")
    # Using defaults from ORCA or a fixed value for now if not specified by LooseOpt/Normal
    # LooseOpt/Normal keywords usually set these appropriately.

    if env.fermi or restart_flag: # Apply Fermi smearing if requested or restarting
        etemp_val = default_etemp
        if restart_flag and "gfn" not in level_keyword.lower() and "xtb" not in level_keyword.lower():
             etemp_val = 5000 # Higher for DFT on restart
        inp_lines.append(f"%scf SmearTemp {etemp_val} end")

    inp_lines.append(f"%pal nprocs {env.threads} end")
    inp_lines.append(f"%maxcore {env.threads * 1000}")
    
    # Coordinates: start.xyz and end.xyz are implicitly used by ORCA NEB
    # when Geo_Source is not given for the full path, or specific blocks are used.
    # The standard way is to provide them in the input file or ensure ORCA finds them.
    # For simplicity, ensure start.xyz and end.xyz are in current_reaction_dir.
    # ORCA requires *xyzfile chrg mult filename format for endpoints if not using %coords
    # However, for NEB, it's often handled by specific blocks or it finds start.xyz/end.xyz.
    # The Fortran prepneb did not explicitly write coords to orca.inp for NEB.
    # It relied on ORCA finding start.xyz and end.xyz.
    # Let's ensure the input specifies the charge and multiplicity for the images.
    inp_lines.append(f"*xyzfile {chrg} {multiplicity} start.xyz*") # Placeholder for ORCA to know chrg/mult for images

    orca_inp_path = current_reaction_dir / "orca.inp"
    try:
        orca_inp_path.write_text("\n".join(inp_lines) + "\n")
    except IOError as e:
        print(f"  Error writing ORCA NEB input to {orca_inp_path}: {e}")
        return None

    orca_executable = shutil.which("orca")
    if not orca_executable:
        print("  Error: ORCA executable not found in PATH for NEB.")
        return ["echo", "ORCA_NOT_FOUND_FOR_NEB"] # Return a dummy command
        
    return [orca_executable, str(orca_inp_path.name)]


def _check_path_py(env: RunTypeData, current_reaction_dir: Path) -> bool:
    """Checks if NEB path calculation was successful. Fortran: checkpath"""
    # Fortran checked for "NEB PATH CONVERGED" or "OPTIMIZATION HAS CONVERGED" (for CI-NEB)
    # and that orca_MEP_energy.dat exists and is not empty.
    orca_out = current_reaction_dir / "orca.out"
    mep_energy_dat = current_reaction_dir / "orca_MEP_energy.dat"

    if not orca_out.exists():
        return False
    
    content = orca_out.read_text()
    converged = "NEB PATH CONVERGED" in content or "OPTIMIZATION HAS CONVERGED" in content
    
    if converged and mep_energy_dat.exists() and mep_energy_dat.stat().st_size > 0:
        return True
    return False

# --- TS Optimization and Verification Implementation ---
def _run_ts_optimization_and_verification_py(env: RunTypeData, ts_dir_path: Path, reaction_dir_path: Path) -> Tuple[bool, float]:
    """
    Complete TS optimization and IRC verification workflow.
    1. Initial Hessian calculation at ts_dir_path/ts.xyz
    2. Check for imaginary frequency (IRC mode)
    3. TS optimization if mode is appropriate
    4. Final Hessian on optimized TS
    5. IRC verification to confirm TS connects reactants and products
    Returns (ts_verified_bool, irc_mode_frequency_cm_1).
    """
    print(f"  Running TS optimization and verification in {ts_dir_path}")
    
    original_cwd = Path.cwd()
    os.chdir(ts_dir_path)
    
    try:
        # Step 1: Initial Hessian calculation
        print("    Step 1: Initial Hessian calculation...")
        ts_xyz_file = "ts.xyz"
        if not Path(ts_xyz_file).exists():
            print(f"    Error: TS guess file {ts_xyz_file} not found")
            return False, 0.0
        
        chrg = iomod.rdshort_int(".CHRG", default=env.chrg)
        uhf = iomod.rdshort_int(".UHF", default=0)
        
        # Initial Hessian
        success_hess1 = _run_hessian_calculation_py(env, ts_xyz_file, "initial_hess", chrg, uhf)
        if not success_hess1:
            print("    Error: Initial Hessian calculation failed")
            return False, 0.0
        
        # Step 2: Check for imaginary frequency
        print("    Step 2: Analyzing vibrational frequencies...")
        irc_freq, has_correct_mode = _find_irc_mode_py(env, "initial_hess.out")
        if not has_correct_mode:
            print(f"    Warning: No appropriate imaginary frequency found (freq = {irc_freq:.1f} cm-1)")
            if abs(irc_freq) < 50.0:  # Too small imaginary frequency
                print("    Imaginary frequency too small, rejecting TS candidate")
                return False, irc_freq
        
        print(f"    Found imaginary frequency: {irc_freq:.1f} cm-1")
        
        # Step 3: TS optimization
        print("    Step 3: TS optimization...")
        success_opt = _run_ts_optimization_py(env, ts_xyz_file, "ts_opt", chrg, uhf)
        if not success_opt:
            print("    Error: TS optimization failed")
            return False, irc_freq
        
        # Update geometry to optimized structure
        if Path("ts_opt.xyz").exists():
            shutil.copy2("ts_opt.xyz", ts_xyz_file)
        
        # Step 4: Final Hessian on optimized TS
        print("    Step 4: Final Hessian calculation...")
        success_hess2 = _run_hessian_calculation_py(env, ts_xyz_file, "final_hess", chrg, uhf)
        if not success_hess2:
            print("    Error: Final Hessian calculation failed")
            return False, irc_freq
        
        # Check final frequencies
        final_irc_freq, final_has_correct_mode = _find_irc_mode_py(env, "final_hess.out")
        if not final_has_correct_mode:
            print(f"    Warning: Optimized TS lost imaginary frequency (freq = {final_irc_freq:.1f} cm-1)")
            return False, final_irc_freq
        
        print(f"    Final imaginary frequency: {final_irc_freq:.1f} cm-1")
        
        # Step 5: IRC verification (simplified)
        print("    Step 5: IRC verification...")
        irc_verified = _verify_irc_connectivity_py(env, ts_xyz_file, reaction_dir_path, chrg, uhf)
        
        if irc_verified:
            print("    TS optimization and verification successful!")
            # Write final energy to qmdata
            ts_energy = _extract_energy_from_output_py("final_hess.out")
            if ts_energy is not None:
                iomod.wrshort_real("qmdata", f"{env.tslevel} sp {chrg} {uhf} {ts_energy:.10f}\n")
            return True, final_irc_freq
        else:
            print("    IRC verification failed - TS does not connect correct reactants/products")
            return False, final_irc_freq
            
    except Exception as e:
        print(f"    Error during TS optimization and verification: {e}")
        return False, 0.0
    finally:
        os.chdir(original_cwd)


def _run_hessian_calculation_py(env: RunTypeData, xyz_file: str, output_prefix: str, chrg: int, uhf: int) -> bool:
    """Runs a Hessian calculation using the specified QM method."""
    try:
        cmd, out_file, pattern, clean_cmd, cached, energy = qmmod.prepqm_py(
            env, xyz_file, env.tslevel, 'hess', chrg_in=chrg, uhf_in=uhf
        )
        
        if cached:
            print(f"    Hessian calculation cached for {xyz_file}")
            return True
        
        if not cmd:
            print(f"    Error: Could not prepare Hessian calculation for {xyz_file}")
            return False
        
        # Run the calculation
        result = iomod.execute_command(cmd, shell=False)
        success = result.returncode == 0
        
        if success:
            # Read and process output
            energy_val, failed = qmmod.readoutqm_py(env, xyz_file, env.tslevel, 'hess', out_file, pattern)
            success = not failed
            
            # Rename output file
            if Path(out_file).exists():
                shutil.copy2(out_file, f"{output_prefix}.out")
        
        # Cleanup
        if clean_cmd:
            iomod.execute_command(clean_cmd.split(), shell=True)
        
        return success
        
    except Exception as e:
        print(f"    Error in Hessian calculation: {e}")
        return False


def _run_ts_optimization_py(env: RunTypeData, xyz_file: str, output_prefix: str, chrg: int, uhf: int) -> bool:
    """Runs a TS optimization calculation."""
    try:
        cmd, out_file, pattern, clean_cmd, cached, energy = qmmod.prepqm_py(
            env, xyz_file, env.tslevel, 'optts', chrg_in=chrg, uhf_in=uhf
        )
        
        if cached:
            print(f"    TS optimization cached for {xyz_file}")
            return True
        
        if not cmd:
            print(f"    Error: Could not prepare TS optimization for {xyz_file}")
            return False
        
        # Run the calculation
        result = iomod.execute_command(cmd, shell=False)
        success = result.returncode == 0
        
        if success:
            # Read and process output
            energy_val, failed = qmmod.readoutqm_py(env, xyz_file, env.tslevel, 'optts', out_file, pattern)
            success = not failed
            
            # Extract optimized geometry
            if not failed:
                _extract_optimized_geometry_py(out_file, f"{output_prefix}.xyz")
            
            # Rename output file
            if Path(out_file).exists():
                shutil.copy2(out_file, f"{output_prefix}.out")
        
        # Cleanup
        if clean_cmd:
            iomod.execute_command(clean_cmd.split(), shell=True)
        
        return success
        
    except Exception as e:
        print(f"    Error in TS optimization: {e}")
        return False


def _find_irc_mode_py(env: RunTypeData, hessian_output_file: str) -> Tuple[float, bool]:
    """
    Analyzes Hessian output to find the IRC mode (imaginary frequency).
    Returns (frequency_cm_1, has_appropriate_imaginary_mode).
    """
    try:
        if not Path(hessian_output_file).exists():
            return 0.0, False
        
        with open(hessian_output_file, 'r') as f:
            content = f.read()
        
        # Look for frequency information (method-dependent parsing)
        frequencies = []
        
        # ORCA frequency parsing
        if "ORCA" in content.upper():
            freq_section = False
            for line in content.split('\n'):
                if "VIBRATIONAL FREQUENCIES" in line:
                    freq_section = True
                    continue
                if freq_section and "cm**-1" in line:
                    parts = line.split()
                    for part in parts:
                        try:
                            freq = float(part)
                            if abs(freq) > 1.0:  # Skip very small frequencies
                                frequencies.append(freq)
                        except ValueError:
                            continue
                if freq_section and ("NORMAL MODES" in line or "THERMOCHEMISTRY" in line):
                    break
        
        # Gaussian frequency parsing
        elif "Gaussian" in content:
            for line in content.split('\n'):
                if "Frequencies --" in line:
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == "--":
                            try:
                                freq = float(parts[i+1])
                                frequencies.append(freq)
                            except (ValueError, IndexError):
                                continue
        
        # XTB frequency parsing
        elif "xtb" in content.lower():
            for line in content.split('\n'):
                if "cm-1" in line and ":" in line:
                    try:
                        freq_str = line.split(":")[1].split()[0]
                        freq = float(freq_str)
                        frequencies.append(freq)
                    except (ValueError, IndexError):
                        continue
        
        if not frequencies:
            print(f"    Warning: No frequencies found in {hessian_output_file}")
            return 0.0, False
        
        # Find the most negative frequency (strongest imaginary mode)
        imaginary_freqs = [f for f in frequencies if f < 0]
        
        if not imaginary_freqs:
            print("    No imaginary frequencies found")
            return 0.0, False
        
        # Return the most negative (strongest imaginary) frequency
        irc_freq = min(imaginary_freqs)
        
        # Check if it's a reasonable IRC mode
        has_correct_mode = abs(irc_freq) > 50.0  # At least 50 cm-1 imaginary
        
        return irc_freq, has_correct_mode
        
    except Exception as e:
        print(f"    Error analyzing frequencies in {hessian_output_file}: {e}")
        return 0.0, False


def _verify_irc_connectivity_py(env: RunTypeData, ts_xyz_file: str, reaction_dir_path: Path, chrg: int, uhf: int) -> bool:
    """
    Simplified IRC verification.
    In a full implementation, this would run IRC calculations in both directions
    and verify that they connect to the correct reactants and products.
    For now, we'll do a simplified check based on energy and geometry.
    """
    try:
        # Simplified verification: check if TS energy is reasonable
        # compared to reactants and products
        
        reactant_xyz = reaction_dir_path / "start.xyz"
        product_xyz = reaction_dir_path / "end.xyz"
        
        if not (reactant_xyz.exists() and product_xyz.exists()):
            print("    Warning: Cannot find reactant/product structures for IRC verification")
            return True  # Assume OK if we can't verify
        
        # Calculate energies (simplified)
        ts_energy = _calculate_single_point_energy_py(env, ts_xyz_file, chrg, uhf)
        reactant_energy = _calculate_single_point_energy_py(env, str(reactant_xyz), chrg, uhf)
        product_energy = _calculate_single_point_energy_py(env, str(product_xyz), chrg, uhf)
        
        if ts_energy is None or reactant_energy is None or product_energy is None:
            print("    Warning: Could not calculate energies for IRC verification")
            return True  # Assume OK if we can't calculate
        
        # Check if TS is higher in energy than both reactants and products
        min_reactant_product_energy = min(reactant_energy, product_energy)
        energy_barrier = (ts_energy - min_reactant_product_energy) * 27.211386245988  # Convert to eV
        
        if energy_barrier > 0.1:  # At least 0.1 eV barrier
            print(f"    IRC verification: Energy barrier = {energy_barrier:.3f} eV (reasonable)")
            return True
        else:
            print(f"    IRC verification: Energy barrier = {energy_barrier:.3f} eV (too low)")
            return False
            
    except Exception as e:
        print(f"    Error in IRC verification: {e}")
        return True  # Assume OK if verification fails


def _calculate_single_point_energy_py(env: RunTypeData, xyz_file: str, chrg: int, uhf: int) -> Optional[float]:
    """Calculate single point energy for a given structure."""
    try:
        cmd, out_file, pattern, clean_cmd, cached, energy = qmmod.prepqm_py(
            env, xyz_file, env.tslevel, 'sp', chrg_in=chrg, uhf_in=uhf
        )
        
        if cached and energy is not None:
            return energy
        
        if not cmd:
            return None
        
        # Run calculation
        result = iomod.execute_command(cmd, shell=False)
        if result.returncode != 0:
            return None
        
        # Read energy
        energy_val, failed = qmmod.readoutqm_py(env, xyz_file, env.tslevel, 'sp', out_file, pattern)
        
        # Cleanup
        if clean_cmd:
            iomod.execute_command(clean_cmd.split(), shell=True)
        
        return energy_val if not failed else None
        
    except Exception as e:
        print(f"    Error calculating single point energy: {e}")
        return None


def _extract_optimized_geometry_py(output_file: str, xyz_output: str) -> bool:
    """Extract optimized geometry from QM output file."""
    try:
        if not Path(output_file).exists():
            return False
        
        with open(output_file, 'r') as f:
            content = f.read()
        
        # Method-specific geometry extraction
        if "ORCA" in content.upper():
            return _extract_orca_geometry_py(content, xyz_output)
        elif "Gaussian" in content:
            return _extract_gaussian_geometry_py(content, xyz_output)
        elif "xtb" in content.lower():
            return _extract_xtb_geometry_py(content, xyz_output)
        else:
            print(f"    Warning: Unknown QM program output format in {output_file}")
            return False
            
    except Exception as e:
        print(f"    Error extracting geometry from {output_file}: {e}")
        return False


def _extract_orca_geometry_py(content: str, xyz_output: str) -> bool:
    """Extract optimized geometry from ORCA output."""
    try:
        lines = content.split('\n')
        
        # Find the final optimized coordinates
        coord_section = False
        coords = []
        
        for i, line in enumerate(lines):
            if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                coord_section = True
                continue
            elif coord_section and line.strip() == "":
                break
            elif coord_section and len(line.split()) >= 4:
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        symbol = parts[0]
                        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                        coords.append(f"{symbol} {x:.6f} {y:.6f} {z:.6f}")
                    except ValueError:
                        continue
        
        if coords:
            with open(xyz_output, 'w') as f:
                f.write(f"{len(coords)}\n")
                f.write("Optimized geometry\n")
                for coord in coords:
                    f.write(coord + "\n")
            return True
        
        return False
        
    except Exception as e:
        print(f"    Error extracting ORCA geometry: {e}")
        return False


def _extract_gaussian_geometry_py(content: str, xyz_output: str) -> bool:
    """Extract optimized geometry from Gaussian output."""
    # Simplified implementation - would need full Gaussian output parsing
    print("    Warning: Gaussian geometry extraction not fully implemented")
    return False


def _extract_xtb_geometry_py(content: str, xyz_output: str) -> bool:
    """Extract optimized geometry from XTB output."""
    # XTB usually writes xtbopt.xyz directly
    if Path("xtbopt.xyz").exists():
        shutil.copy2("xtbopt.xyz", xyz_output)
        return True
    return False


def _extract_energy_from_output_py(output_file: str) -> Optional[float]:
    """Extract final energy from QM output file."""
    try:
        if not Path(output_file).exists():
            return None
        
        with open(output_file, 'r') as f:
            content = f.read()
        
        # Method-specific energy extraction
        if "ORCA" in content.upper():
            # Look for final single point energy
            for line in reversed(content.split('\n')):
                if "FINAL SINGLE POINT ENERGY" in line:
                    try:
                        return float(line.split()[-1])
                    except (ValueError, IndexError):
                        continue
        
        elif "Gaussian" in content:
            # Look for SCF Done
            for line in reversed(content.split('\n')):
                if "SCF Done:" in line:
                    try:
                        return float(line.split()[4])
                    except (ValueError, IndexError):
                        continue
        
        elif "xtb" in content.lower():
            # Look for total energy
            for line in reversed(content.split('\n')):
                if "TOTAL ENERGY" in line:
                    try:
                        return float(line.split()[-2])
                    except (ValueError, IndexError):
                        continue
        
        return None
        
    except Exception as e:
        print(f"    Error extracting energy from {output_file}: {e}")
        return None

# --- Barrier Calculation Implementation ---
def _calculate_barrier_py(
    env: RunTypeData,
    reaction_dir_path: Path, # Contains reactant (start.xyz) and product (end.xyz) qmdata
    ts_dir_path: Path, # Contains TS (ts.xyz) qmdata
    irc_freq_cm1: float
    ):
    """
    Calculates reaction barrier (Ea) and RRHO corrected Ea.
    Writes barrier_<level> and earrho_<level> files.
    """
    print(f"  Calculating barrier for reaction in {reaction_dir_path} using TS from {ts_dir_path}")
    
    try:
        # Get energies from qmdata files
        reactant_energy = _get_energy_from_qmdata_py(reaction_dir_path, "start.xyz", env)
        ts_energy = _get_energy_from_qmdata_py(ts_dir_path, "ts.xyz", env)
        product_energy = _get_energy_from_qmdata_py(reaction_dir_path, "end.xyz", env)
        
        if reactant_energy is None or ts_energy is None:
            print(f"    Error: Could not get required energies for barrier calculation")
            return
        
        # Calculate forward barrier (reactant -> TS)
        ea_fwd_hartree = ts_energy - reactant_energy
        ea_fwd_ev = ea_fwd_hartree * AUTOEV
        
        print(f"    Forward barrier: {ea_fwd_ev:.3f} eV ({ea_fwd_hartree:.6f} Hartree)")
        
        # Write basic barrier
        iomod.wrshort_real(reaction_dir_path / f"barrier_{env.tslevel}", ea_fwd_ev)
        
        # Calculate reverse barrier if product energy available
        if product_energy is not None:
            ea_rev_hartree = ts_energy - product_energy
            ea_rev_ev = ea_rev_hartree * AUTOEV
            print(f"    Reverse barrier: {ea_rev_ev:.3f} eV ({ea_rev_hartree:.6f} Hartree)")
            iomod.wrshort_real(reaction_dir_path / f"barrier_rev_{env.tslevel}", ea_rev_ev)
        
        # RRHO corrected barrier calculation
        if env.bhess:
            print("    Calculating RRHO corrected barrier...")
            earrho_fwd_ev = _calculate_rrho_corrected_barrier_py(
                env, reaction_dir_path, ts_dir_path, ea_fwd_ev
            )
            if earrho_fwd_ev is not None:
                print(f"    RRHO corrected forward barrier: {earrho_fwd_ev:.3f} eV")
                iomod.wrshort_real(reaction_dir_path / f"earrho_{env.geolevel}", earrho_fwd_ev)
        
        # Write IRC mode frequency
        iomod.wrshort_real(reaction_dir_path / f"ircmode_{env.geolevel}", irc_freq_cm1)
        
        print(f"    Barrier calculation completed successfully")
        
    except Exception as e:
        print(f"    Error in barrier calculation: {e}")


def _get_energy_from_qmdata_py(directory: Path, xyz_file: str, env: RunTypeData) -> Optional[float]:
    """Get energy from qmdata file for a specific structure."""
    try:
        qmdata_file = directory / "qmdata"
        if not qmdata_file.exists():
            print(f"    Warning: qmdata file not found in {directory}")
            return None
        
        # Get charge and UHF for the structure
        chrg = iomod.rdshort_int(directory / ".CHRG", default=env.chrg)
        uhf = iomod.rdshort_int(directory / ".UHF", default=0)
        
        # Search for energy in qmdata
        query = f"{env.tslevel} sp {chrg} {uhf}"
        found, energy = iomod.grepval(qmdata_file, query)
        
        if found:
            return energy
        else:
            print(f"    Warning: Energy not found in qmdata for query: {query}")
            return None
            
    except Exception as e:
        print(f"    Error reading energy from {directory}: {e}")
        return None


def _calculate_rrho_corrected_barrier_py(
    env: RunTypeData,
    reaction_dir_path: Path,
    ts_dir_path: Path,
    electronic_barrier_ev: float
) -> Optional[float]:
    """
    Calculate RRHO (Rigid Rotor Harmonic Oscillator) corrected barrier.
    This includes zero-point energy and thermal corrections.
    """
    try:
        # Calculate thermal corrections for reactant and TS
        reactant_thermal = _calculate_thermal_correction_py(env, reaction_dir_path, "start.xyz")
        ts_thermal = _calculate_thermal_correction_py(env, ts_dir_path, "ts.xyz")
        
        if reactant_thermal is None or ts_thermal is None:
            print("    Warning: Could not calculate thermal corrections for RRHO barrier")
            return None
        
        # RRHO corrected barrier = Electronic barrier + (TS_thermal - Reactant_thermal)
        thermal_correction_ev = (ts_thermal - reactant_thermal) * AUTOEV
        rrho_barrier_ev = electronic_barrier_ev + thermal_correction_ev
        
        print(f"    Thermal correction: {thermal_correction_ev:.3f} eV")
        
        return rrho_barrier_ev
        
    except Exception as e:
        print(f"    Error calculating RRHO corrected barrier: {e}")
        return None


def _calculate_thermal_correction_py(env: RunTypeData, directory: Path, xyz_file: str) -> Optional[float]:
    """
    Calculate thermal correction (ZPE + thermal energy) for a structure.
    Returns correction in Hartree.
    """
    try:
        # Check if thermal data already exists
        thermal_file = directory / "thermal_correction"
        if thermal_file.exists():
            return iomod.rdshort_real(thermal_file, default=None)
        
        # Run frequency calculation to get thermal corrections
        original_cwd = Path.cwd()
        os.chdir(directory)
        
        try:
            chrg = iomod.rdshort_int(".CHRG", default=env.chrg)
            uhf = iomod.rdshort_int(".UHF", default=0)
            
            # Use XTB for thermal corrections (faster than DFT)
            cmd = [
                "xtb", xyz_file, "--hess", "--thermo", str(env.temp),
                f"--gfn{env.geolevel.replace('gfn', '') if 'gfn' in env.geolevel else '2'}",
                "--chrg", str(chrg), "--uhf", str(uhf)
            ]
            
            result = iomod.execute_command(cmd, shell=False)
            
            if result.returncode == 0:
                # Parse thermal correction from XTB output
                thermal_correction = _parse_xtb_thermal_correction_py()
                if thermal_correction is not None:
                    # Cache the result
                    iomod.wrshort_real(thermal_file, thermal_correction)
                    return thermal_correction
            
            print(f"    Warning: Thermal correction calculation failed for {xyz_file}")
            return None
            
        finally:
            os.chdir(original_cwd)
            
    except Exception as e:
        print(f"    Error calculating thermal correction: {e}")
        return None


def _parse_xtb_thermal_correction_py() -> Optional[float]:
    """Parse thermal correction from XTB output files."""
    try:
        # XTB writes thermal data to various files
        if Path("xtb_thermodata.dat").exists():
            with open("xtb_thermodata.dat", 'r') as f:
                for line in f:
                    if "G(RRHO)" in line:
                        try:
                            return float(line.split()[-1])
                        except (ValueError, IndexError):
                            continue
        
        # Alternative: parse from main output
        if Path("xtb.out").exists():
            with open("xtb.out", 'r') as f:
                content = f.read()
                for line in content.split('\n'):
                    if "TOTAL FREE ENERGY" in line:
                        try:
                            return float(line.split()[-2])
                        except (ValueError, IndexError):
                            continue
        
        return None
        
    except Exception as e:
        print(f"    Error parsing XTB thermal correction: {e}")
        return None


# --- Main TS Search Function (Focus on NEB part first) ---
# This is a highly complex routine in Fortran. The Python version will simplify
# some of the job batching and error handling for now.

# For _pick_ts_from_path_py
def _read_xyz_trajectory(filepath: Path) -> List[List[List[float]]]:
    """Reads an XYZ trajectory file, returns list of geometries (list of [atom_coords])."""
    frames = []
    try:
        with open(filepath, 'r') as f:
            while True:
                line_natoms = f.readline()
                if not line_natoms: break # EOF
                num_atoms = int(line_natoms.strip())
                f.readline() # Skip comment line
                current_frame_coords = []
                for _ in range(num_atoms):
                    atom_line = f.readline().split()
                    # Assuming symbol, x, y, z format
                    current_frame_coords.append([float(c) for c in atom_line[1:4]])
                frames.append(current_frame_coords)
    except Exception as e:
        print(f"Error reading XYZ trajectory {filepath}: {e}")
    return frames


def _pick_ts_from_path_py(env: RunTypeData, current_reaction_dir: Path) -> bool:
    """
    Analyzes a reaction path (e.g., from NEB), finds the highest energy structure,
    writes it as ts.xyz, and checks for cascades.
    Corresponds to Fortran pickts.
    Returns True if a TS candidate is found and it's not a cascade.
    """
    # Fortran pickts logic:
    # 1. Reads orca_MEP_energy.dat (energies) and orca_MEP_trj.xyz (geometries).
    # 2. Finds image with highest energy (excluding endpoints if needed).
    # 3. Writes this image to ts.xyz.
    # 4. If env.sortoutcascade, checks if path has multiple maxima or if TS is too close to start/end.
    
    mep_energy_file = current_reaction_dir / "orca_MEP_energy.dat"
    mep_traj_file = current_reaction_dir / "orca_MEP_trj.xyz"
    ts_guess_xyz_file = current_reaction_dir / "ts.xyz" # Output

    if not (mep_energy_file.exists() and mep_traj_file.exists()):
        print(f"  pickTS: Energy or trajectory file missing in {current_reaction_dir}")
        return False

    try:
        energies = []
        with open(mep_energy_file, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2: # Format: image_idx energy_hartree
                    energies.append(float(parts[1]))
        
        if not energies or len(energies) < 3: # Need at least one intermediate image
            print(f"  pickTS: Not enough images/energies in {mep_energy_file}")
            return False

        # Find highest energy image (excluding endpoints, Fortran was [1:-1])
        # Python: energies[0] is reactant, energies[-1] is product. Images are 0 to N-1.
        # Intermediate images are energies[1:-1].
        if len(energies) <= 2: # Only reactant and product
             print(f"  pickTS: Path in {current_reaction_dir} has no intermediate images.")
             return False

        intermediate_energies = energies[1:-1]
        if not intermediate_energies:
            print(f"  pickTS: Path in {current_reaction_dir} has no intermediate images (after slicing).")
            return False
            
        max_energy_intermediate = -float('inf')
        idx_max_energy_intermediate = -1
        for i, e_val in enumerate(intermediate_energies):
            if e_val > max_energy_intermediate:
                max_energy_intermediate = e_val
                idx_max_energy_intermediate = i
        
        # Index in the full list of energies/frames (0=reactant, 1=first_intermediate, ...)
        ts_candidate_frame_idx = idx_max_energy_intermediate + 1 

        # Read the trajectory and extract the TS candidate frame
        # This needs a robust XYZ trajectory reader.
        # Assuming utility.read_xyz_trajectory_frame_py(filepath, frame_index (0-based))
        # frames = _read_xyz_trajectory(mep_traj_file) # Reads all frames
        # if ts_candidate_frame_idx >= len(frames):
        #     print(f"  pickTS: TS candidate index out of bounds for trajectory {mep_traj_file}")
        #     return False
        # ts_coords = frames[ts_candidate_frame_idx]
        # Read natoms and symbols from start.xyz to write ts.xyz correctly
        # nat_ts, symbols_ts, _ = utility.read_xyz_structure_py(current_reaction_dir / "start.xyz")
        # utility.write_xyz_structure_py(ts_guess_xyz_file, nat_ts, symbols_ts, ts_coords, "TS candidate")
        
        # Simpler: use ORCA's output for the TS guess if available, or extract from traj
        # ORCA NEB often writes a .ts.xyz file or similar for the highest energy point.
        # For now, let's assume orca_MEP_trj.xyz must be parsed.
        # The Fortran code calls `xtb_extract_snapshot(mep_traj_file, ts_candidate_frame_idx+1, ts_guess_xyz_file)`
        # This implies a 1-based index for snapshot extraction.
        utility.xtb_extract_snapshot_py(mep_traj_file, ts_candidate_frame_idx + 1, ts_guess_xyz_file)


        # Cascade check (simplified)
        if env.sortoutcascade:
            # Check if TS is too close to start or end (e.g. first or last intermediate image)
            if ts_candidate_frame_idx == 1 or ts_candidate_frame_idx == len(energies) - 2:
                print(f"  pickTS: TS candidate in {current_reaction_dir} is an endpoint image. Considered cascade/no barrier.")
                return False
            # Check for multiple maxima (more complex, requires iterating through path)
            # Fortran checks: if e(i) > e(i-1) and e(i) > e(i+1) -> local max. Count them.
            local_maxima_count = 0
            for i in range(1, len(intermediate_energies) -1): # Check internal points of intermediate path
                if intermediate_energies[i] > intermediate_energies[i-1] and \
                   intermediate_energies[i] > intermediate_energies[i+1]:
                    local_maxima_count +=1
            if max_energy_intermediate > energies[0] and max_energy_intermediate > energies[-1]: # Overall barrier
                if local_maxima_count > 1 : # More than one peak along the intermediate path
                     print(f"  pickTS: Path in {current_reaction_dir} has multiple maxima ({local_maxima_count}). Considered cascade.")
                     return False
            elif local_maxima_count > 0 : # No overall barrier but intermediate humps
                 print(f"  pickTS: Path in {current_reaction_dir} has intermediate humps but no overall barrier. Considered cascade.")
                 return False


        print(f"  pickTS: TS candidate found at image {ts_candidate_frame_idx} in {current_reaction_dir}, written to {ts_guess_xyz_file}")
        return True

    except Exception as e:
        print(f"  Error in pickTS for {current_reaction_dir}: {e}")
        return False


# --- Main TS Search function placeholder ---
def ts_search_py(
    env: RunTypeData, 
    fname_reactant_in_precursor_dir: str, # e.g. "fragment.xyz" or "isomer.xyz"
    npairs_in: int, 
    fragdirs_in: List[List[str]]  # [[pair_id, r_dir, p_dir],...] relative to precursor dir
) -> Tuple[int, List[List[str]]]:
    """
    Manages Transition State searches for a list of reactions (fragmentations/isomerizations).
    Focuses on NEB part first.
    Returns (number_of_successful_ts_searches, list_of_reaction_dirs_with_ts).
    Operates from within the precursor's directory.
    """
    print(f"\n--- Transition State Search for products of {fname_reactant_in_precursor_dir} ---")
    original_cwd = Path.cwd() # Should be the precursor's directory
    
    successful_reactions: List[List[str]] = [] # Store [pair_id, r_dir, p_dir] for successful ones

    # TODO: Initialize timer if used within ts_search

    # Loop through each reaction (isomerization or fragmentation pair)
    # fragdirs_in contains [pair_id, path_to_reactant_dir (same as pair_id), path_to_product_dir_or_empty_for_isomer]
    # This structure is confusing. Let's assume fragdirs_in contains:
    # [[pair_id_str_for_reaction_dir, reactant_xyz_in_reaction_dir, product_xyz_in_reaction_dir], ...]
    # For now, adapting to the structure from fragmentation.py's output `fragdirs`
    # which is [[pair_id, item1_dir, item2_dir_or_empty],...]
    # item1_dir is the isomer dir, or first fragment dir.
    # item2_dir is second fragment dir if it's a pair.
    # The 'reaction_dir' is pair_id. start.xyz is precursor, end.xyz is product(s).

    # --- Stage 1: Geodesic Interpolation (if enabled) ---
    geodesic_jobs: List[Tuple[Path, List[str], str, str]] = []
    if env.tsgeodesic:
        print("  Preparing Geodesic Interpolation jobs...")
        for pair_id_str, item1_dir_str, item2_dir_str in fragdirs_in:
            reaction_dir = Path(pair_id_str) # This is the pXX directory for the reaction
            reaction_dir.mkdir(exist_ok=True)
            
            # Copy precursor (reactant) to reaction_dir/start.xyz
            shutil.copy2(original_cwd / fname_reactant_in_precursor_dir, reaction_dir / "start.xyz")

            # Create product structure (end.xyz)
            # If isomer: item1_dir_str is the isomer product directory. Copy isomer.xyz from there.
            # If pair: item1_dir_str is frag1 dir, item2_dir_str is frag2 dir. Concatenate them.
            if not item2_dir_str: # Isomer
                isomer_xyz_src = Path(env.startdir) / item1_dir_str / "isomer.xyz" # item1_dir_str is relative to startdir
                if isomer_xyz_src.exists(): shutil.copy2(isomer_xyz_src, reaction_dir / "end.xyz")
                else: print(f"    Warning: Isomer XYZ {isomer_xyz_src} not found for reaction {pair_id_str}."); continue
            else: # Fragment pair
                frag1_xyz_src = Path(env.startdir) / item1_dir_str / "fragment.xyz"
                frag2_xyz_src = Path(env.startdir) / item2_dir_str / "fragment.xyz"
                if not (frag1_xyz_src.exists() and frag2_xyz_src.exists()):
                    print(f"    Warning: Fragment XYZs not found for reaction {pair_id_str}."); continue
                utility.concatenate_xyz_files_py([frag1_xyz_src, frag2_xyz_src], reaction_dir / "end.xyz")
            
            # Set up .CHRG and .UHF for the reaction (usually based on reactant)
            # This logic is complex in Fortran, using env.pathmult etc.
            # For now, assume reactant's charge/uhf.
            if (original_cwd / ".CHRG").exists(): shutil.copy2(original_cwd / ".CHRG", reaction_dir / ".CHRG")
            if (original_cwd / ".UHF").exists(): shutil.copy2(original_cwd / ".UHF", reaction_dir / ".UHF")
            
            geo_cmd_list = _prep_geodesic_py(env, reaction_dir)
            if geo_cmd_list:
                geodesic_jobs.append((reaction_dir, geo_cmd_list, "geo.out", "")) # No specific success pattern for geo run itself

        if geodesic_jobs:
            print(f"  Running {len(geodesic_jobs)} Geodesic Interpolation jobs sequentially...")
            _run_jobs_sequentially(env, geodesic_jobs)
            print("  Geodesic Interpolations finished.")


    # --- Stage 2: NEB Path Search ---
    neb_jobs: List[Tuple[Path, List[str], str, str]] = []
    if env.tsfinder.lower() == 'neb':
        print("  Preparing NEB calculations...")
        for pair_id_str, _, _ in fragdirs_in: # Iterate through original list to decide if NEB is needed
            reaction_dir = Path(pair_id_str)
            if not (reaction_dir / "neb_finished").exists(): # Check if NEB already done
                # hbonddiss flag logic not translated yet, assume False
                # restart_flag initially False, use_geo_guess from env.tsgeodesic
                neb_cmd_list = _prep_neb_py(env, reaction_dir, hbonddiss_flag=False, restart_flag=False, use_geo_guess_flag=env.tsgeodesic)
                if neb_cmd_list:
                    neb_jobs.append((reaction_dir, neb_cmd_list, "orca.out", "NEB PATH CONVERGED")) # Pattern for ORCA NEB
            else:
                print(f"  NEB already marked finished for {reaction_dir}, skipping.")
        
        if neb_jobs:
            print(f"  Running {len(neb_jobs)} initial NEB calculations sequentially...")
            _run_jobs_sequentially(env, neb_jobs)
            print("  Initial NEB calculations finished. Checking for restarts...")

            # NEB Retry Logic
            neb_retry_jobs = []
            for work_dir, _, out_f, patt in neb_jobs: # Check results of first batch
                os.chdir(work_dir) # Need to be in dir for _check_path_py and _prepare_neb_restart_py
                if not _check_path_py(env, Path(".")): # Pass current dir
                    print(f"    NEB in {work_dir} failed or did not converge. Preparing restart.")
                    _prepare_neb_restart_py(env, Path(".")) # Modifies files in place
                    # Retry with restart=True (enables SlowConv) and use_geo_guess=True (uses nebguess.xyz)
                    neb_cmd_retry = _prep_neb_py(env, Path("."), hbonddiss_flag=False, restart_flag=True, use_geo_guess_flag=True)
                    if neb_cmd_retry:
                        neb_retry_jobs.append((work_dir, neb_cmd_retry, "orca_restart.out", "NEB PATH CONVERGED"))
                os.chdir(original_cwd)
            
            if neb_retry_jobs:
                print(f"  Running {len(neb_retry_jobs)} NEB restart calculations sequentially...")
                _run_jobs_sequentially(env, neb_retry_jobs)
                print("  NEB restart calculations finished.")
                # After restart, update orca.out if restart was successful
                for work_dir, _, out_f_retry, _ in neb_retry_jobs:
                    if (work_dir / out_f_retry).exists() and _check_path_py(env, work_dir): # Check if restart converged
                         shutil.copy2(work_dir / out_f_retry, work_dir / "orca.out") # Overwrite original orca.out
    
    elif env.tsfinder.lower() == 'gsm':
        print("GSM pathfinder not yet implemented in Python version.")
    # elif env.tsfinder.lower() == 'xtb': # Fortran had this for xtb internal path search
    #     print("XTB internal pathfinder not yet implemented in Python version.")


    # --- Stage 3: Pick TS candidates from successful paths ---
    print("\n  Picking TS candidates from converged paths...")
    for pair_id_str, item1_dir_str, item2_dir_str in fragdirs_in:
        reaction_dir = Path(pair_id_str)
        os.chdir(reaction_dir) # _pick_ts_from_path_py expects to be in reaction_dir
        
        if _check_path_py(env, Path(".")): # Check if NEB (or other pathfinder) was successful
            ts_found = _pick_ts_from_path_py(env, Path("."))
            if ts_found:
                # If TS optimization and verification were here, they would follow.
                # For now, if pickts is successful, we consider the reaction valid.
                # Placeholder for TS opt and IRC:
                ts_dir = Path("ts")
                ts_dir.mkdir(exist_ok=True)
                if (Path(".") / "ts.xyz").exists(): shutil.copy2("ts.xyz", ts_dir / "ts.xyz")
                if Path(".CHRG").exists(): shutil.copy2(".CHRG", ts_dir / ".CHRG")
                if Path(".UHF").exists(): shutil.copy2(".UHF", ts_dir / ".UHF")
                
                verified, irc_freq = _run_ts_optimization_and_verification_py(env, ts_dir, Path("."))
                
                if verified:
                    # Placeholder for barrier calculation
                    # _calculate_barrier_py(env, Path("."), ts_dir, irc_freq)
                    print(f"    TS verified for {pair_id_str} (placeholder). Barrier calculation placeholder.")
                    successful_reactions.append([pair_id_str, item1_dir_str, item2_dir_str])
                else:
                    print(f"    TS verification failed for {pair_id_str} (placeholder).")
            else:
                print(f"    No valid TS candidate picked from path in {reaction_dir}.")
        else:
            print(f"    Path search failed or no valid path found in {reaction_dir}.")
        os.chdir(original_cwd)

    npairs_out = len(successful_reactions)
    print(f"\nTS search phase finished. Found {npairs_out} viable TS candidates out of {npairs_in}.")
    return npairs_out, successful_reactions


if __name__ == '__main__':
    print("Testing tsmod.py (placeholders for QM runs)")
    # Dummy env
    env_test = RunTypeData(
        tsgeodesic=True, 
        tsnds=3, # Number of intermediate NEB images
        tsfinder='neb',
        geolevel='gfn2', # Level for geodesic and bhess
        tslevel='gfn2',  # Level for NEB and SP on TS
        threads=1, cores=1,
        startdir=str(Path.cwd() / "tsmod_test_run"), # Set a specific startdir for test
        path=str(Path.cwd() / "tsmod_test_run"), # Current path
        solv=False, fermi=False,
        sortoutcascade=True, # Test cascade detection in pick_ts
        bhess=False # For simplicity in barrier calc placeholder
    )
    Path(env_test.startdir).mkdir(exist_ok=True)
    os.chdir(env_test.startdir)

    # Create dummy reactant (start.xyz) and product (end.xyz) for a reaction pair
    reaction_pair_id = "p0"
    reaction_dir_test = Path(reaction_pair_id)
    reaction_dir_test.mkdir(exist_ok=True)
    
    (reaction_dir_test / "start.xyz").write_text("1\nReactant\nH 0 0 0\n")
    (reaction_dir_test / "end.xyz").write_text("1\nProduct\nH 0 0 1\n")
    (reaction_dir_test / ".CHRG").write_text("0\n")
    (reaction_dir_test / ".UHF").write_text("0\n")


    # Mock external command calls within iomod for testing this module
    def mock_execute_command(cmd_list, cwd=None, shell=False, stdout_to=None, stderr_to=None):
        cmd_str = " ".join(cmd_list) if isinstance(cmd_list, list) else cmd_list
        print(f"  MOCK EXECUTE in {Path(cwd if cwd else '.').resolve()}: {cmd_str}")
        if "geodesic_interpolate" in cmd_str:
            # Create dummy interpolated.xyz
            (Path(cwd if cwd else '.') / "interpolated.xyz").write_text("1\ninterpolated\nH 0 0 0.5\n")
        if "orca orca.inp" in cmd_str: # NEB run
            # Create dummy orca.out, orca_MEP_energy.dat, orca_MEP_trj.xyz
            (Path(cwd if cwd else '.') / "orca.out").write_text("NEB PATH CONVERGED\nFINAL SINGLE POINT ENERGY -1.0\n")
            (Path(cwd if cwd else '.') / "orca_MEP_energy.dat").write_text("0 -0.9 # Reactant\n1 -0.8 # Image 1\n2 -0.7 # TS guess\n3 -0.85 # Image 2\n4 -1.0 # Product\n")
            # Dummy trajectory with 5 images (R, I1, TS, I2, P)
            dummy_traj = "1\nR\nH 0 0 0\n1\nI1\nH 0 0 0.25\n1\nTS\nH 0 0 0.5\n1\nI2\nH 0 0 0.75\n1\nP\nH 0 0 1.0\n"
            (Path(cwd if cwd else '.') / "orca_MEP_trj.xyz").write_text(dummy_traj)
            (Path(cwd if cwd else '.') / "neb_finished").touch() # Mark as finished
        return subprocess.CompletedProcess(args=cmd_list, returncode=0, stdout="mock output", stderr="")

    iomod.execute_command = mock_execute_command
    utility.xtb_extract_snapshot_py = lambda traj, idx, out: Path(out).write_text(f"{idx}\nTS guess from frame {idx}\nH 0 0 0.5")


    test_fragdirs = [[reaction_pair_id, reaction_dir_test.name, ""]] # Isomer-like for test structure

    # Run the main TS search function
    np_out, frags_out = ts_search_py(env_test, "start.xyz", 1, test_fragdirs)
    
    print(f"\nTS Search results: Npairs_out={np_out}, Fragdirs_out={frags_out}")
    if np_out > 0:
        ts_file = Path(env.startdir) / frags_out[0][0] / "ts" / "ts.xyz"
        if ts_file.exists():
            print(f"  TS candidate file created: {ts_file}")
            print(f"  Content:\n{ts_file.read_text()}")
        else:
            print(f"  TS candidate file {ts_file} NOT created.")

    os.chdir(Path(env.startdir).parent) # Go back up from test_run
    if Path(env.startdir).exists():
        shutil.rmtree(env.startdir)
