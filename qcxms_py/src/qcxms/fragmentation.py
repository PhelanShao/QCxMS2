import os
import shutil
from pathlib import Path
from typing import Tuple, Optional, List, Dict

try:
    from .data import RunTypeData, WP
    from . import iomod
    from . import qmmod
    from . import charges
    from . import reaction
    from . import utility
    from . import tsmod # Import the new tsmod
    # from . import mcsimu # To be imported when available
    # from . import structools # To be imported when available
    from .constants import AUTOEV, EVTOKCAL # Assuming these are in a constants module
except ImportError:
    # Fallbacks for standalone/testing
    print("Attempting to import dummy/mock modules for fragmentation.py standalone run.")
    from data import RunTypeData, WP # type: ignore
    import iomod_mock as iomod # type: ignore
    import qmmod_mock as qmmod # type: ignore
    import charges_mock as charges # type: ignore
    import reaction_mock as reaction # type: ignore
    import utility_mock as utility # type: ignore
    import tsmod_mock as tsmod # type: ignore # Add mock for tsmod if testing standalone
    AUTOEV = 27.211386245988
    EVTOKCAL = 23.060547830618307


def _generate_fragments_with_crest_py(
    env: RunTypeData, 
    fname_xyz: str, 
    precursor_charge: int, # Actual charge of the precursor for CREST
    precursor_uhf: int     # Actual UHF of the precursor for CREST
) -> Tuple[int, bool]:
    """
    Generates fragments using CREST --msreact.
    Returns (number_of_pairs_generated, success_bool).
    Operates in the current working directory.
    """
    print(f"Generating fragments for {fname_xyz} with CREST --msreact...")
    original_path = Path.cwd()

    # Determine energy window for CREST
    # Fortran: nat0 = atoms in initial molecule of the whole cascade (from env.startdir/in.xyz)
    #          sumreac = energy to form current precursor from initial molecule (read from sumreac_tslevel)
    #          ewin_kcal = (env.ieeatm * nat0 * 3.0 * EVTOKCAL) - (sumreac * EVTOKCAL)
    
    # This requires knowing the path to the absolute start of the cascade.
    # For now, let's use a simplified ewin or make it a large default if nat0/sumreac are hard to get here.
    # Placeholder for nat0 and sumreac_to_current_precursor_kcal
    nat0_initial_molecule = iomod.rdshort_int(Path(env.startdir) / "in.xyz", default=30) # Default if not found
    sumreac_to_current_precursor_kcal = iomod.rdshort_real(original_path / f"sumreac_{env.tslevel}", default=0.0) * EVTOKCAL
    
    ewin_kcal = (env.ieeatm * nat0_initial_molecule * 3.0 * EVTOKCAL) - sumreac_to_current_precursor_kcal
    if ewin_kcal <= 0: # Ensure positive energy window
        ewin_kcal = 100.0 # Default to a reasonable window like 100 kcal/mol if calculated is too small/negative
        print(f"  Warning: Calculated ewin for CREST was <=0. Using default: {ewin_kcal:.2f} kcal/mol")
    else:
        print(f"  Calculated CREST ewin: {ewin_kcal:.2f} kcal/mol (based on nat0={nat0_initial_molecule}, sumreac={sumreac_to_current_precursor_kcal/EVTOKCAL:.2f} eV)")

    iomod.copy_file(fname_xyz, "infrag.xyz") # CREST input name

    # Backup .CHRG and .UHF as CREST might delete them
    chrg_temp_path = original_path / "chrgtemp"
    uhf_temp_path = original_path / "uhftemp"
    if (original_path / ".CHRG").exists(): shutil.move(original_path / ".CHRG", chrg_temp_path)
    if (original_path / ".UHF").exists(): shutil.move(original_path / ".UHF", uhf_temp_path)

    utility.set_omp_to_one_py() # Ensure XTB runs with 1 thread if CREST calls it that way

    crest_cmd_parts = [
        "crest", "infrag.xyz", "--msreact", "--mslargeprint",
        "--msnbonds", str(env.msnbonds),
        "--chrg", str(precursor_charge),
        "--uhf", str(precursor_uhf),
        "--T", str(env.cores), # Total cores for CREST internal parallelization
        f"--{env.geolevel if env.geolevel in ['gfn2', 'gfn1'] else 'gfn2'}", # CREST msreact usually with GFNn
        "--ewin", f"{ewin_kcal:.5f}" # kcal/mol for CREST
    ]

    if env.msnoiso or (env.nfrag_lvl > 1 and not env.msfulliso): # nfrag_lvl needs to be added to env from main loop
        crest_cmd_parts.append("--msnoiso")
    if env.msiso:
        crest_cmd_parts.append("--msiso")
    if env.msnshifts > 0:
        crest_cmd_parts.extend(["--msnshifts", str(env.msnshifts)])
    if env.msnshifts2 > 0:
        crest_cmd_parts.extend(["--msnshifts2", str(env.msnshifts2)])
    
    crest_cmd_parts.append("--msattrh") # Default in Fortran

    if env.msmolbar: crest_cmd_parts.append("--msmolbar")
    if env.msinchi: crest_cmd_parts.append("--msinchi")
    
    if env.msfragdist > 0.0:
        Path("crestms.inp").write_text(f"fragdist {env.msfragdist:.6f}\n")
        crest_cmd_parts.extend(["--msinput", "crestms.inp"])

    # Execute CREST
    print(f"  Running CREST: {' '.join(crest_cmd_parts)}")
    # Redirect output to msreact.out and errors to cresterror.out
    with open("msreact.out", "w") as f_out, open("cresterror.out", "w") as f_err:
        result = iomod.execute_command(crest_cmd_parts, stdout_to=f_out, stderr_to=f_err) # Adapted execute_command

    # Restore .CHRG, .UHF
    if chrg_temp_path.exists(): shutil.move(chrg_temp_path, original_path / ".CHRG")
    if uhf_temp_path.exists(): shutil.move(uhf_temp_path, original_path / ".UHF")

    if result.returncode != 0:
        print(f"  CREST --msreact failed with code {result.returncode}. Check msreact.out and cresterror.out.")
        return 0, False
    
    npairs = iomod.rdshort_int("npairs", default=0)
    print(f"  CREST generated {npairs} initial fragment pairs/isomers.")

    # Cleanup some CREST files
    files_to_remove = ["crest_msreact_products.xyz", "infrag.xyz", "isomers.xyz", "crestms.inp", "crest.restart"]
    if env.printlevel < 2: # Keep if debug
        files_to_remove.extend(["crest_best.xyz", "crest_input_copy.xyz", "crest_unique_products.xyz", 
                                "crest_unique_products.sorted", "fragmentpairs.xyz", "coord", 
                                "scoord.1", "coordprot.0", "lmocent.coord", "cresterror.out"])
    for f_name in files_to_remove:
        iomod.remove_file(Path(f_name))
    if env.mskeepdir:
        print("  Keeping MSDIR directory as per --mskeepdir.")
    else:
        iomod.rmrf(Path("MSDIR"))
        
    return npairs, True


def _read_crest_fragment_list_py(
    env: RunTypeData, 
    npairs_from_crest: int, # Value read from "npairs" file
    current_fragmentation_level: int
) -> Tuple[int, List[List[str]], List[WP]]:
    """
    Parses msreact.out to get fragment directory info and relative energies.
    Returns (actual_npairs, fragdirs_list, rel_energies_kcal_list).
    fragdirs_list: [[pair_id, frag1_dir, frag2_dir_or_empty],...]
    """
    fragdirs_list: List[List[str]] = []
    rel_energies_kcal_list: List[WP] = []
    
    if not Path("msreact.out").exists():
        print("Error: msreact.out not found for parsing fragment list.")
        return 0, [], []

    nisomers = 0
    n_actual_products = 0

    with open("msreact.out", "r") as f:
        found_header = False
        for line in f:
            if "directory | fragment type | rel. energy [kcal/mol]" in line:
                found_header = True
                if env.msiso: # Fortran logic: read one more line if msiso
                    next(f, None) 
                continue
            
            if found_header:
                parts = line.split()
                if not parts or len(parts) < 4: # Expecting at least pXX isomer rel_E unit
                    if "===================================================" in line: # End of table
                        break
                    continue # Skip empty or malformed lines

                pair_id_str = parts[0]
                frag_kind = parts[1].lower()
                try:
                    rel_e_kcal = float(parts[2])
                except ValueError:
                    print(f"Warning: Could not parse energy from msreact.out line: {line.strip()}")
                    continue # Skip this entry

                if frag_kind == "isomer":
                    nisomers += 1
                    if env.msnoiso or (current_fragmentation_level > 1 and not env.msfulliso):
                        print(f"  Skipping isomer {pair_id_str} due to msnoiso/msfulliso logic.")
                        continue # Skip isomer based on settings
                    
                    fragdirs_list.append([pair_id_str, pair_id_str, ""]) # Isomer: item1_dir is pair_id, no item2
                    rel_energies_kcal_list.append(rel_e_kcal)
                    n_actual_products += 1
                
                elif frag_kind == "fragmentpair":
                    # Next two lines should be the fragment details
                    try:
                        line_f1 = next(f).split()
                        line_f2 = next(f).split()
                        frag1_dir = line_f1[0]
                        frag2_dir = line_f2[0]
                        fragdirs_list.append([pair_id_str, frag1_dir, frag2_dir])
                        rel_energies_kcal_list.append(rel_e_kcal)
                        n_actual_products += 1
                    except StopIteration:
                        print("Error: Unexpected end of file while reading fragment pair details in msreact.out.")
                        break
                    except IndexError:
                        print(f"Error: Malformed fragment pair entry in msreact.out for pair {pair_id_str}.")
                        continue
                
                if n_actual_products >= npairs_from_crest : # Stop if we've read all expected pairs
                    break
    
    print(f"Read {n_actual_products} products ({nisomers} isomers before filtering) from msreact.out.")
    return n_actual_products, fragdirs_list, rel_energies_kcal_list


def _collect_fragments_py(
    env: RunTypeData,
    npairs: int,
    fragdirs: List[List[str]], # [[pair_id, frag1_dir, frag2_dir_or_empty],...]
    output_list_fname: str, # e.g., "fragments"
    precursor_dir_in_cascade: str # e.g., "" for initial, "p0/f1" for sub-fragment
) -> Tuple[List[str], int]:
    """
    Collects important fragment information into a list and writes to a file.
    This is a simplified version; Fortran `collectfrags` might be more complex.
    Returns (list_of_significant_fragment_paths, count).
    """
    # Fortran `collectfrags` logic is not fully available.
    # This placeholder assumes we list all valid, processed fragment directories.
    # "Significance" would typically be determined by mcsimu results (intensities).
    
    significant_fragments: List[str] = []
    for pair_id, item1_dir, item2_dir in fragdirs:
        if item1_dir: # Isomer or first fragment of a pair
            # Construct full path relative to the start of the cascade for uniqueness
            # full_path1 = Path(env.startdir) / precursor_dir_in_cascade / item1_dir
            # For now, just use the relative path from current precursor
            rel_path1 = Path(precursor_dir_in_cascade) / item1_dir
            significant_fragments.append(str(rel_path1))
        if item2_dir: # Second fragment of a pair
            # full_path2 = Path(env.startdir) / precursor_dir_in_cascade / item2_dir
            rel_path2 = Path(precursor_dir_in_cascade) / item2_dir
            significant_fragments.append(str(rel_path2))
            
    try:
        with open(output_list_fname, "w") as f:
            for frag_path in significant_fragments:
                f.write(frag_path + "\n")
    except IOError:
        print(f"Error writing fragment list to {output_list_fname}")

    return significant_fragments, len(significant_fragments)


# --- Main Fragmentation Workflow ---
def calculate_fragments_py(
    env: RunTypeData, 
    precursor_xyz_fname: str, # Filename of the precursor (e.g., "fragment.xyz" or "isomer.xyz")
    # eiee_dist: List[WP], # IEE distribution energies
    # piee_dist: List[WP], # IEE distribution probabilities/weights
    # These are now assumed to be accessible via env or re-read if needed by mcsimu
    fragmentation_level: int, # Depth of fragmentation (nfragl in Fortran)
    precursor_dir_in_cascade: str, # Path string like "p0/f1", or "" for initial molecule
) -> Tuple[List[str], int]: # Returns (list_of_next_fragments_to_process, count)
    """
    Main fragmentation routine for a given precursor molecule.
    Orchestrates fragment generation, QM calculations, charge assignment, filtering, and simulation.
    """
    
    # This function operates within the specific directory of the precursor molecule.
    # precursor_dir_in_cascade is its path relative to env.startdir (the overall run root).
    # current_abs_path = Path(env.startdir) / precursor_dir_in_cascade
    # os.chdir(current_abs_path) # Ensure we are in the right directory

    print(f"\n--- Fragmentation Level {fragmentation_level} for Precursor: {Path.cwd() / precursor_xyz_fname} ---")

    # Check if precursor is a single atom
    is_atom = utility.check_atom_py(precursor_xyz_fname)
    if is_atom:
        if env.printlevel >= 2: print(f"  Fragment {precursor_xyz_fname} is an atom, cannot be further fragmented.")
        return [], 0

    # If restarting, check if fragments were already generated
    # This is simplified; a more robust restart would check progress at multiple stages.
    npairs_generated: int
    fragdirs_from_crest: List[List[str]]
    # rel_energies_kcal: List[WP] # Relative energies from CREST output

    if Path("msreact.out").exists() and env.restartrun:
        print("  Found existing msreact.out, attempting to read fragment list (restart mode).")
        npairs_generated = iomod.rdshort_int("npairs", default=0)
        if npairs_generated > 0:
            npairs_generated, fragdirs_from_crest, _ = _read_crest_fragment_list_py(env, npairs_generated, fragmentation_level)
            # TODO: Logic to move fragment directories back if they were globally stored (Fortran code)
            # This part is complex and depends on how fragdirs are stored globally vs locally.
            # For now, assume if msreact.out exists, subdirs p0, p1 etc. also exist locally.
        else: # msreact.out exists but npairs is 0 or unreadable
            print("  msreact.out found but no pairs listed or readable. Regenerating fragments.")
            npairs_generated, success = _generate_fragments_with_crest_py(env, precursor_xyz_fname, env.chrg, utility.get_uhf_from_file_py())
            if not success: return [], 0
            npairs_generated, fragdirs_from_crest, _ = _read_crest_fragment_list_py(env, npairs_generated, fragmentation_level)
    else:
        precursor_chrg_val, precursor_uhf_val = qmmod._get_chrg_uhf(env, precursor_xyz_fname)
        npairs_generated, success = _generate_fragments_with_crest_py(env, precursor_xyz_fname, precursor_chrg_val, precursor_uhf_val)
        if not success: return [], 0
        npairs_generated, fragdirs_from_crest, _ = _read_crest_fragment_list_py(env, npairs_generated, fragmentation_level)

    if npairs_generated == 0:
        print("  No fragment pairs or isomers generated by CREST.")
        return [], 0
    
    current_fragdirs = fragdirs_from_crest
    current_npairs = npairs_generated

    # GFF Topology Check (only if new fragments, i.e., not restart with existing msreact.out)
    # The Fortran code did this if not restart.
    # if not (Path("msreact.out").exists() and env.restartrun):
    #    current_npairs, current_fragdirs = reaction.gff_topology_check_py(env, current_npairs, current_fragdirs)
    #    if current_npairs == 0: print("  No products left after GFF topology check."); return [],0
    # This function is not fully translated in reaction.py yet. Skipping for now.
    print("  Skipping GFF topology check (placeholder).")


    # Optimize fragments for IP calculation (if not hot IPs)
    if not env.hotip:
        # current_npairs, current_fragdirs = reaction.optimize_fragments_py(env, current_npairs, current_fragdirs)
        # This function is not fully translated in reaction.py yet. Skipping for now.
        print("  Skipping fragment optimization for IPs (assuming hotip or placeholder).")
        if current_npairs == 0: print("  No products left after fragment optimization for IP."); return [],0

    # Assign charges
    current_npairs, current_fragdirs = charges.assign_charges_py(env, current_npairs, current_fragdirs)
    if current_npairs == 0: print("  No products left after charge assignment."); return [],0

    # Optimize products at their assigned charge
    # current_npairs, current_fragdirs = reaction.optimize_products_py(env, current_npairs, current_fragdirs)
    # This function is not fully translated in reaction.py yet. Skipping for now.
    print("  Skipping product optimization at assigned charge (placeholder).")
    if current_npairs == 0: print("  No products left after product optimization."); return [],0

    # Calculate reaction energies
    # This returns a dict, but subsequent steps need e_pairs list matching fragdirs.
    # The Fortran calcreactionenergy directly filled an e_pairs array.
    # reaction_energy_dict = reaction.calculate_reaction_energies_py(env, current_npairs, current_fragdirs)
    # e_pairs_for_sorting = [reaction_energy_dict.get(fd[0], 10000.0) for fd in current_fragdirs]
    # This function is not fully translated in reaction.py yet. Using dummy energies.
    print("  Skipping reaction energy calculation (placeholder). Assigning dummy DE=0 for all.")
    e_pairs_for_sorting = [0.0] * current_npairs 
    # Also, need to write dummy de_<level> and sumreac_<level> files for sortouthighe to read.
    for p_id, _, _ in current_fragdirs:
        (Path(p_id) / f"de_{env.tslevel}").write_text("0.0\n")
    (Path(".") / f"sumreac_{env.tslevel}").write_text("0.0\n") # Precursor sumreac


    # Sort out high energy fragments
    current_npairs, current_fragdirs = reaction.sort_out_high_energy_fragments_py(
        env, env.tslevel, use_rrho=env.bhess, 
        npairs_in=current_npairs, fragdirs_in=current_fragdirs, 
        scale_iee_factor=3.0 # Fortran used 3.0 here
    )
    if current_npairs == 0: print("  No products left after high energy filter."); return [],0

    # Transition State Search (if not env.nots)
    if not env.nots:
        print(f"  Starting Transition State search for {current_npairs} potential reactions...")
        # precursor_xyz_fname is the name of the reactant file within the current working directory (precursor_dir)
        current_npairs, current_fragdirs = tsmod.ts_search_py(
            env,
            precursor_xyz_fname, 
            current_npairs,
            current_fragdirs 
        )
        if current_npairs == 0: 
            print("  No products left after TS search filter.")
            return [], 0
    else:
        print("  Skipping Transition State search due to -nots flag.")
    
    # Monte Carlo Simulation
    # This requires eiee/piee to be passed or read by mcsimu
    # mcsimu_results = mcsimu.run_monte_carlo_py(env, current_npairs, current_fragdirs, eiee_dist, piee_dist, fragmentation_level)
    # This function from mcsimu.py is not translated yet.
    print("  Skipping Monte Carlo simulation (placeholder).")

    # Collect important fragments (based on mcsimu_results, which we don't have)
    # For now, assume all remaining fragments are "important" enough to be listed.
    # The path stored in fraglist needs to be relative to env.startdir for uniqueness across levels.
    final_fragment_list, num_final_frags = _collect_fragments_py(
        env, current_npairs, current_fragdirs, 
        "fragments", precursor_dir_in_cascade
    )
    
    # File management for multi-level fragmentation
    if fragmentation_level > 1: # Fortran used nfragl (fragmentation level)
        if Path("fragments").exists():
            with open(original_path / ".." / "allfragments", "a") as outfile, open("fragments", "r") as infile:
                outfile.write(infile.read())
        
        # Move processed fragment/isomer directories to a global location relative to env.startdir
        # Fortran: mv fragdirs(i,j) env%startdir/startdir_fragdirs(i,j)
        # Example: current precursor is "p0/f1". A product is "p0/f1/p0".
        # It should be moved to "p0f1p0" at the top level (env.startdir).
        for pair_id, item1_dir_str, item2_dir_str_or_empty in current_fragdirs:
            items_to_move = [Path(item1_dir_str)]
            if item2_dir_str_or_empty: items_to_move.append(Path(item2_dir_str_or_empty))

            for local_item_path in items_to_move: # e.g. local_item_path is "p0" relative to "p0/f1"
                # Global name should incorporate the precursor path.
                # precursor_dir_in_cascade might be "p0/f1"
                # local_item_path.name is "p0"
                # global_name = precursor_dir_in_cascade.replace("/", "") + local_item_path.name 
                # This logic needs to exactly match Fortran's naming for recursive calls.
                # Fortran: dir2 = trim(startdir)//trim(fragdirs(i, j))
                #          jobcall = 'mv '//trim(fragdirs(i, j))//' '//trim(env%startdir)//'/'//trim(dir2)
                # If startdir = "p0f1", fragdirs(i,j) = "p0", then target is env.startdir / "p0f1p0"
                
                # Construct the global unique name
                global_unique_name_parts = []
                if precursor_dir_in_cascade:
                    global_unique_name_parts.append(precursor_dir_in_cascade.replace(os.sep, ""))
                global_unique_name_parts.append(local_item_path.name)
                global_target_name = "".join(global_unique_name_parts)
                global_target_path = Path(env.startdir) / global_target_name

                if local_item_path.exists(): # local_item_path is relative to current precursor dir
                    print(f"  Moving {local_item_path} to {global_target_path}")
                    shutil.move(str(local_item_path), str(global_target_path))
    else: # fragmentation_level == 1
        if Path("fragments").exists():
            iomod.copy_file("fragments", "allfragments") # For the first level

    return final_fragment_list, num_final_frags


def calc_start_fragment_py(env: RunTypeData, fname_xyz: str, fragmentation_level: int):
    """
    Performs initial QM calculations (opt, sp, bhess) for a precursor molecule.
    """
    print(f"\n--- Initial calculations for precursor: {fname_xyz} (Level {fragmentation_level}) ---")
    # utility.set_temperature_py(env, fname_xyz) # Needs translation

    # Get current charge/UHF or set defaults for the precursor
    # This function operates in the precursor's own directory.
    chrg, uhf = qmmod._get_chrg_uhf(env, fname_xyz, chrg_in=env.chrg if fragmentation_level==1 else None)


    # 1. Optimization
    print(f"  Optimizing at {env.geolevel}...")
    cmd_opt, out_opt, patt_opt, clean_opt, cached_opt, _ = \
        qmmod.prepqm_py(env, fname_xyz, env.geolevel, 'opt', chrg_in=chrg, uhf_in=uhf)
    
    failed_opt = False
    if not cached_opt and cmd_opt:
        res_opt = iomod.execute_command(cmd_opt, shell=False) # prepqm_py returns list
        if res_opt.returncode == 0:
            _, failed_opt = qmmod.readoutqm_py(env, fname_xyz, env.geolevel, 'opt', out_opt, patt_opt)
        else: failed_opt = True
        if clean_opt: iomod.execute_command(clean_opt.split(), shell=True)
    
    if failed_opt: # Retry logic from Fortran
        print(f"  Retrying optimization for {fname_xyz}...")
        cmd_opt, out_opt, patt_opt, clean_opt, _, _ = \
            qmmod.prepqm_py(env, fname_xyz, env.geolevel, 'opt', chrg_in=chrg, uhf_in=uhf, restart=True)
        if cmd_opt: # Should always get a command unless error in prepqm
            res_opt = iomod.execute_command(cmd_opt, shell=False)
            if res_opt.returncode == 0:
                _, failed_opt = qmmod.readoutqm_py(env, fname_xyz, env.geolevel, 'opt', out_opt, patt_opt)
            else: failed_opt = True # Still failed
            if clean_opt: iomod.execute_command(clean_opt.split(), shell=True)
    
    if failed_opt:
        msg = f"CRITICAL: Initial geometry optimization failed for {fname_xyz}."
        if fragmentation_level == 1: sys.exit(msg)
        else: print(msg + " Skipping this fragment."); return False # False indicates failure to proceed

    # 2. Single Point / TDDFT
    qm_job_type_sp = "tddft" if env.exstates > 0 else "sp"
    print(f"  Calculating SP/TDDFT at {env.tslevel}...")
    cmd_sp, out_sp, patt_sp, clean_sp, cached_sp, _ = \
        qmmod.prepqm_py(env, fname_xyz, env.tslevel, qm_job_type_sp, chrg_in=chrg, uhf_in=uhf)
        
    failed_sp = False
    if not cached_sp and cmd_sp:
        res_sp = iomod.execute_command(cmd_sp, shell=False)
        if res_sp.returncode == 0:
            _, failed_sp = qmmod.readoutqm_py(env, fname_xyz, env.tslevel, qm_job_type_sp, out_sp, patt_sp)
        else: failed_sp = True
        if clean_sp: iomod.execute_command(clean_sp.split(), shell=True)

    if failed_sp: # Retry
        print(f"  Retrying SP/TDDFT for {fname_xyz}...")
        cmd_sp, out_sp, patt_sp, clean_sp, _, _ = \
            qmmod.prepqm_py(env, fname_xyz, env.tslevel, qm_job_type_sp, chrg_in=chrg, uhf_in=uhf, restart=True)
        if cmd_sp:
            res_sp = iomod.execute_command(cmd_sp, shell=False)
            if res_sp.returncode == 0:
                _, failed_sp = qmmod.readoutqm_py(env, fname_xyz, env.tslevel, qm_job_type_sp, out_sp, patt_sp)
            else: failed_sp = True
            if clean_sp: iomod.execute_command(clean_sp.split(), shell=True)

    if failed_sp:
        msg = f"CRITICAL: Initial SP/TDDFT calculation failed for {fname_xyz}."
        if fragmentation_level == 1: sys.exit(msg)
        else: print(msg + " Skipping this fragment."); return False

    # 3. Hessian (bhess)
    if env.bhess:
        print(f"  Calculating Hessian (bhess) at {env.geolevel}...")
        cmd_hess, out_hess, patt_hess, clean_hess, cached_hess, _ = \
            qmmod.prepqm_py(env, fname_xyz, env.geolevel, 'bhess', chrg_in=chrg, uhf_in=uhf)
        
        failed_hess = False
        if not cached_hess and cmd_hess:
            res_hess = iomod.execute_command(cmd_hess, shell=False)
            if res_hess.returncode == 0:
                _, failed_hess = qmmod.readoutqm_py(env, fname_xyz, env.geolevel, 'bhess', out_hess, patt_hess)
            else: failed_hess = True
            if clean_hess: iomod.execute_command(clean_hess.split(), shell=True)

        if failed_hess: # Retry
            print(f"  Retrying Hessian (bhess) for {fname_xyz}...")
            cmd_hess, out_hess, patt_hess, clean_hess, _, _ = \
                qmmod.prepqm_py(env, fname_xyz, env.geolevel, 'bhess', chrg_in=chrg, uhf_in=uhf, restart=True)
            if cmd_hess:
                res_hess = iomod.execute_command(cmd_hess, shell=False)
                if res_hess.returncode == 0:
                    _, failed_hess = qmmod.readoutqm_py(env, fname_xyz, env.geolevel, 'bhess', out_hess, patt_hess)
                else: failed_hess = True
                if clean_hess: iomod.execute_command(clean_hess.split(), shell=True)
        
        if failed_hess:
            # Fortran code stops if first SPH fails.
            sys.exit(f"CRITICAL: Initial Hessian (bhess) calculation failed for {fname_xyz}.")
            # return False # Would be for non-critical failure
            
    print(f"  Initial calculations for {fname_xyz} complete.")
    return True # Success


if __name__ == '__main__':
    print("Testing fragmentation.py functions...")
    # Requires extensive setup of dummy files, RunTypeData, and mock modules.
    # Example:
    # test_env_frag = RunTypeData(startdir=str(Path(".").resolve()))
    # test_env_frag.geolevel = "gfn2"
    # test_env_frag.tslevel = "gfn2"
    # test_env_frag.iplevel = "gfn2"
    # test_env_frag.chrg = 1
    # test_env_frag.cores = 2
    # test_env_frag.printlevel = 3 # Keep files for debugging
    # test_env_frag.nfrag_lvl = 1 # Add this field to RunTypeData if it's used
    
    # # Create a dummy input XYZ for the initial molecule
    # initial_xyz_content = "2\nCO\nC 0 0 0\nO 0 0 1.12\n"
    # Path(test_env_frag.startdir / "in.xyz").write_text(initial_xyz_content)
    # Path("fragment.xyz").write_text(initial_xyz_content) # Assume it's copied here

    # # Mock utility.get_uhf_from_file_py if it's directly called
    # utility.get_uhf_from_file_py = lambda: 0

    # success_init = calc_start_fragment_py(test_env_frag, "fragment.xyz", 1)
    # if success_init:
    #     print("calc_start_fragment_py seemed to run.")
    #     # Further test calculate_fragments_py if desired, needs more mocks (CREST output etc.)
    #     # next_frags, num_next = calculate_fragments_py(test_env_frag, "fragment.xyz", 1, "")
    #     # print(f"Next fragments to process: {next_frags}, Count: {num_next}")

    # # Cleanup
    # for f in [".CHRG", ".UHF", "qmdata", "infrag.xyz", "msreact.out", "npairs", 
    #           "cresterror.out", "crestms.inp", "fragment.xyz", Path(test_env_frag.startdir) / "in.xyz"]:
    #     Path(f).unlink(missing_ok=True)

    pass
```
