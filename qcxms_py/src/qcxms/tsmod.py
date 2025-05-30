import os
import shutil
import math
import re
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

# --- Placeholder for TS Optimization and Verification ---
def _run_ts_optimization_and_verification_py(env: RunTypeData, ts_dir_path: Path, reaction_dir_path: Path) -> Tuple[bool, float]:
    """
    Placeholder for the complex TS optimization and IRC verification workflow.
    This would involve:
    1. Initial Hessian at ts_dir_path/ts.xyz.
    2. find_irc_py to check for correct mode.
    3. OptTS calculation if mode is good.
    4. Final Hessian on optimized TS.
    5. find_irc_py again to confirm TS connects reactants and products.
    Returns (ts_verified_bool, irc_mode_frequency_cm_1).
    """
    print(f"  PLACEHOLDER: Running TS optimization and verification in {ts_dir_path}")
    # Simulate finding a relevant imaginary frequency
    # Write dummy files that subsequent steps might expect
    (ts_dir_path / "g98.out").write_text("Imaginary Freq: -150.0 cm-1 connecting reactants and products.\n")
    (ts_dir_path / "qmdata").write_text(f"{env.tslevel} sp {iomod.rdshort_int(ts_dir_path / '.CHRG')} {iomod.rdshort_int(ts_dir_path / '.UHF')} -100.0\n") # Dummy TS energy
    return True, -150.0 

# --- Placeholder for Barrier Calculation ---
def _calculate_barrier_py(
    env: RunTypeData, 
    reaction_dir_path: Path, # Contains reactant (start.xyz) and product (end.xyz) qmdata
    ts_dir_path: Path, # Contains TS (ts.xyz) qmdata
    irc_freq_cm1: float
    ):
    """
    Placeholder for calculating reaction barrier (Ea) and RRHO corrected Ea.
    Writes barrier_<level> and earrho_<level> files.
    """
    print(f"  PLACEHOLDER: Calculating barrier for reaction in {reaction_dir_path} using TS from {ts_dir_path}")
    # Dummy values
    ea_fwd_ev = 1.0 # eV
    earrho_fwd_ev = 0.8 # eV
    
    iomod.wrshort_real(reaction_dir_path / f"barrier_{env.tslevel}", ea_fwd_ev)
    if env.bhess:
        iomod.wrshort_real(reaction_dir_path / f"earrho_{env.geolevel}", earrho_fwd_ev)
    iomod.wrshort_real(reaction_dir_path / f"ircmode_{env.geolevel}", irc_freq_cm1)


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
```
