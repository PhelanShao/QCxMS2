import sys
import os
import random
from pathlib import Path
from typing import List, Set, Optional # Added Optional

try:
    from .data import RunTypeData, Timer, WP
    from . import argparser as ap
    from . import iomod
    from . import utility
    from . import qmmod
    from . import charges
    from . import isotopes
    from . import iee
    from . import fragmentation
    from . import reaction
    from . import tsmod 
    from . import mcsimu # Actual import for mcsimu
    # from . import plot # Placeholder
    from .constants import AUTOEV, EVTOKCAL, KB_EV_K, PI_CONST, WP
except ImportError:
    print("Attempting to import local/mock modules for main.py standalone run.")
    from data import RunTypeData, Timer, WP # type: ignore
    import argparser_mock as ap # type: ignore
    import iomod_mock as iomod # type: ignore
    import utility_mock as utility # type: ignore
    import qmmod_mock as qmmod # type: ignore
    import charges_mock as charges # type: ignore
    import isotopes_mock as isotopes # type: ignore
    import iee_mock as iee # type: ignore
    import fragmentation_mock as fragmentation # type: ignore
    import reaction_mock as reaction # type: ignore
    import tsmod_mock as tsmod # type: ignore
    import mcsimu_mock as mcsimu # type: ignore # Mock for mcsimu
    import math 
    AUTOEV = 27.211386245988
    EVTOKCAL = 23.060547830618307
    KB_EV_K = 8.6173332621415e-5
    PI_CONST = math.pi


def qcxms2_header_py():
    print("="*70)
    print("QCxMS2 - Quantum Chemistry guided Mass Spectrometry")
    print("Version: Python Port (In Progress)")
    print("="*70)

def disclaimer_py():
    print("\nDisclaimer: This is a research code. Use at your own risk.\n")


def main(argv: Optional[List[str]] = None):
    qcxms2_header_py()
    disclaimer_py()
    print("Command line input:")
    print(f"> python {' '.join(sys.argv)}") 

    main_timer = Timer() 
    main_timer.init(20) 
    overall_timer_idx = 0 
    main_timer.start(overall_timer_idx, "Total QCxMS2 Run")

    env_startdir = Path.cwd()
    env = ap.parse_arguments(argv) 
    env.startdir = str(env_startdir) 
    env.path = str(env_startdir)     

    if env.logo: ap.qcxms2_logo_py()
    ap.citation_py(env) 

    if env.iseed == 0: random.seed()
    else: random.seed(env.iseed)
    print(f"Initialized random seed (seed: {'system' if env.iseed == 0 else env.iseed}).")

    ap.check_settings_py(env)
    ap.check_programs_py(env) 

    if env.ircrun:
        print("IRC run mode (-ircrun) selected. Not fully implemented.")
        sys.exit(0)
    if env.picktsrun:
        print("Pick TS run mode (-picktsrun) selected. Not fully implemented.")
        sys.exit(0)
    
    if env.oplot:
        print("Plotting only mode (-oplot).")
        if not (Path(env.startdir) / "allfrags.dat").exists():
            sys.exit("Error: allfrags.dat not found for -oplot mode.")
        ntot_frags = iomod.rdshort_int_py(Path(env.startdir) / "allfrags.dat", default=0)
        all_fragment_paths: List[str] = []
        if ntot_frags > 0:
            with open(Path(env.startdir) / "allfrags.dat", "r") as f:
                f.readline(); all_fragment_paths = [line.strip() for line in f if line.strip()]
        print(f"Plotting for {ntot_frags} fragments. Plotting not fully implemented.")
        main_timer.stop(overall_timer_idx); ap.eval_timer_py(main_timer); sys.exit(0)

    if env.esim:
        print("Energy simulation only mode (-esim).")
        env.cores = 1; utility.set_omp_threads_py(env, 1) 
        os.chdir(env.startdir) 
        ip_file_path = Path(f"ip_{env.iplevel}")
        if not ip_file_path.exists():
             sys.exit(f"Error: IP file {ip_file_path} not found for -esim mode.")
        env.ehomo = iomod.rdshort_real_py(ip_file_path)
        print(f"Read IP from {ip_file_path}: {env.ehomo:.3f} eV")
        eiee_dist, piee_dist = iee.get_energy_distribution_py(env, env.nsamples, plotlog=True)
        if not eiee_dist: sys.exit("Failed to generate IEE distribution for -esim mode.")
        
        print("Calling mcesim_py for energy simulation...")
        mcsimu.mcesim_py(env, eiee_dist, piee_dist) # Actual call
        print("Finished mcesim_py call.")
        
        if not (Path(env.startdir) / "allfrags.dat").exists():
             print("Warning: allfrags.dat not found for -esim mode plotting.")
        ntot_frags = iomod.rdshort_int_py(Path(env.startdir) / "allfrags.dat", default=0)
        all_fragment_paths = []
        if ntot_frags > 0:
            with open(Path(env.startdir) / "allfrags.dat", "r") as f:
                f.readline(); all_fragment_paths = [line.strip() for line in f if line.strip()]
        print(f"Plotting for {ntot_frags} fragments after esim. Plotting not fully implemented.")
        main_timer.stop(overall_timer_idx); ap.eval_timer_py(main_timer); sys.exit(0)

    if (Path(env.startdir) / "QCxMS2_finished").exists():
        print("Prior QCxMS2 calculation found. Restarting...")
        iomod.cleanup_qcxms2_for_restart_py(env) 
        env.restartrun = True
    else: env.restartrun = False

    current_precursor_xyz_fname = Path(env.infile).name 
    os.chdir(env.startdir)
    
    if env.mode.lower() == "ei":
        ip_file = Path(f"ip_{env.iplevel}")
        if env.restartrun and ip_file.exists():
            env.ehomo = iomod.rdshort_real_py(ip_file)
            print(f"Restart: Using existing IP from {ip_file}: {env.ehomo:.3f} eV")
        else:
            success_opt = qmmod.run_qm_job_and_read(env, current_precursor_xyz_fname, env.geolevel, 'opt', chrg_in=0, uhf_in=None)
            if not success_opt: sys.exit(f"Initial optimization of {current_precursor_xyz_fname} failed.")
            print(f"Computing IP at {env.iplevel} level...")
            success_ip = qmmod.calculate_ip_py(env, current_precursor_xyz_fname)
            if not success_ip:
                print("Warning: IP calculation failed, using default value of 10.0 eV")
                env.ehomo = 10.0
            else:
                env.ehomo = qmmod.read_ip_result_py(env, current_precursor_xyz_fname)
            iomod.wrshort_real_py(ip_file, env.ehomo)
        print(f"IP of input molecule is {env.ehomo:.3f} eV")
    
    success_start_frag = fragmentation.calc_start_fragment_py(env, current_precursor_xyz_fname, 1) # Level 1 for initial precursor
    if not success_start_frag: sys.exit(f"Initial calculations for {current_precursor_xyz_fname} failed.")
    iomod.copy_file(current_precursor_xyz_fname, "fragment.xyz")
    current_precursor_xyz_fname = "fragment.xyz" 

    print("\nGenerating IEE distribution...")
    eiee_dist, piee_dist = iee.get_energy_distribution_py(env, env.nsamples, plotlog=True)
    if not eiee_dist: sys.exit("Failed to generate IEE distribution.")

    print("\n" + "-"*32 + "\n| Starting first fragmentation |\n" + "-"*32 + "\n")
    main_timer.start(1, "Fragmentation Level 1")
    precursor_cascade_path = "" 
    env.path = str(Path(env.startdir) / precursor_cascade_path) 
    os.chdir(env.path)
    env.nfrag_lvl = 1 # type: ignore 
    initial_precursor_intensity = 1.0

    fragments_for_next_level, num_new_frags = fragmentation.calculate_fragments_py(
        env, current_precursor_xyz_fname, 
        eiee_dist_energies=eiee_dist, piee_dist_probs=piee_dist,
        precursor_intensity=initial_precursor_intensity, 
        fragmentation_level=1, # Precursor is level 0, its products will be level 1. Or this is "generation 1".
        precursor_dir_in_cascade=precursor_cascade_path
    )
    os.chdir(env.startdir) 
    main_timer.stop(1)

    all_fragment_paths_generated: Set[str] = set(fragments_for_next_level)
    # Add initial precursor to all_fragment_paths for collection purposes if not already there via isomer mechanism
    # This path should be relative to env.startdir
    initial_precursor_rel_path = Path(precursor_cascade_path) / current_precursor_xyz_fname if precursor_cascade_path else Path(current_precursor_xyz_fname)
    all_fragment_paths_generated.add(str(initial_precursor_rel_path))

    if not fragments_for_next_level:
        print("No fragments generated or survived from the first fragmentation level.")
    else:
        for k_level in range(2, env.nfrag + 2): # Max generations: env.nfrag. Loop k_level from 2 to env.nfrag+1.
                                                # If env.nfrag=0, no loop. If env.nfrag=1, k_level is 2.
                                                # This means products up to generation env.nfrag+1 are explored.
            print(f"\n" + "-"*50 + f"\n| Starting fragmentation level {k_level} for {len(fragments_for_next_level)} precursor(s) |\n" + "-"*50 + "\n")
            main_timer.start(k_level, f"Fragmentation Level {k_level}")
            current_precursors_this_level = list(fragments_for_next_level) 
            fragments_for_next_level = [] 
            num_total_new_frags_this_level = 0
            env.nfrag_lvl = k_level # type: ignore
            intensity_threshold = env.pthr / 100.0 # Convert % to fraction

            for precursor_rel_path_str in current_precursors_this_level:
                precursor_abs_path = Path(env.startdir) / precursor_rel_path_str
                pfrag_val = iomod.rdshort_real_py(precursor_abs_path / "pfrag", default=0.0)
                
                if pfrag_val >= intensity_threshold:
                    os.chdir(precursor_abs_path); env.path = str(precursor_abs_path)
                    current_xyz = "isomer.xyz" if Path("isomer.xyz").exists() else "fragment.xyz"
                    if not Path(current_xyz).exists():
                        print(f"  Skipping {precursor_rel_path_str}: {current_xyz} not found."); os.chdir(env.startdir); continue
                    print(f"  --- Processing precursor for level {k_level}: {precursor_rel_path_str} (Intensity: {pfrag_val:.4f}) ---")
                    success_sub_frag_init = fragmentation.calc_start_fragment_py(env, current_xyz, k_level)
                    if not success_sub_frag_init:
                        print(f"    Initial calculations failed. Skipping further fragmentation."); os.chdir(env.startdir); continue
                    current_xyz = "fragment.xyz" if Path("fragment.xyz").exists() else "isomer.xyz" # May have been renamed by calc_start_fragment
                    
                    newly_gen_paths, num_new = fragmentation.calculate_fragments_py(
                        env, current_xyz,
                        eiee_dist_energies=eiee_dist, piee_dist_probs=piee_dist,
                        precursor_intensity=pfrag_val, 
                        fragmentation_level=k_level, # This is the level of the *current precursor* being fragmented
                        precursor_dir_in_cascade=precursor_rel_path_str 
                    )
                    os.chdir(env.startdir) 
                    fragments_for_next_level.extend(newly_gen_paths)
                    all_fragment_paths_generated.update(newly_gen_paths)
                    num_total_new_frags_this_level += num_new
                else:
                    print(f"  Skipping precursor {precursor_rel_path_str}: Intensity ({pfrag_val:.4f}) < threshold ({intensity_threshold:.4f}).")
            main_timer.stop(k_level)
            if num_total_new_frags_this_level == 0:
                print(f"No new viable fragments generated at level {k_level}."); break 
            fragments_for_next_level = sorted(list(set(fragments_for_next_level)))

    print("\n--- Finalizing ---")
    allfrags_list = sorted(list(all_fragment_paths_generated))
    ntot_unique_frags = len(allfrags_list)
    try:
        with open(Path(env.startdir) / "allfrags.dat", "w") as f:
            f.write(f"{ntot_unique_frags}\n")
            for frag_path in allfrags_list: f.write(f"{frag_path}\n")
        print(f"Written {ntot_unique_frags} unique fragment paths to allfrags.dat.")
    except IOError: print("Error writing allfrags.dat.")

    # Generate final spectrum and plots
    print("\n=== Generating Final Results ===")
    
    if (Path(env.startdir) / "allfrags.dat").exists():
        try:
            # Import plotting module
            from . import plotting
            
            # Generate mass spectrum from pfrag files
            _generate_mass_spectrum_py(env)
            
            # Generate all plots
            if plotting.MATPLOTLIB_AVAILABLE:
                print("Generating plots...")
                plotting.generate_all_plots_py(env, Path(env.startdir))
            else:
                print("Matplotlib not available - skipping plot generation")
                
            # Generate simulation summary
            _generate_simulation_summary_py(env)
            
        except ImportError as e:
            print(f"Warning: Could not import plotting module: {e}")
        except Exception as e:
            print(f"Error generating final results: {e}")
    else:
        print("allfrags.dat not found, skipping plotting stage.")
    
    main_timer.stop(overall_timer_idx)
    ap.eval_timer_py(main_timer) 
    print("\nQCxMS2 Python port finished normally (or with placeholders).")
    iomod.touch(Path(env.startdir) / "QCxMS2_finished")


def _generate_mass_spectrum_py(env: RunTypeData) -> bool:
    """
    Generate mass spectrum from pfrag files.
    Collects all fragment intensities and creates spectrum data.
    """
    try:
        print("Generating mass spectrum from fragment intensities...")
        
        spectrum_data = {}  # mass -> intensity
        
        # Read allfrags.dat to get all fragment paths
        allfrags_file = Path(env.startdir) / "allfrags.dat"
        if not allfrags_file.exists():
            print("  Error: allfrags.dat not found")
            return False
        
        with open(allfrags_file, 'r') as f:
            lines = f.readlines()
            if len(lines) < 2:
                print("  Error: allfrags.dat is empty or malformed")
                return False
            
            # Process each fragment
            for line in lines[1:]:  # Skip first line (count)
                frag_path = line.strip()
                if not frag_path:
                    continue
                
                frag_dir = Path(env.startdir) / frag_path
                if not frag_dir.exists():
                    continue
                
                # Read fragment intensity
                pfrag_file = frag_dir / "pfrag"
                if pfrag_file.exists():
                    intensity = iomod.rdshort_real_py(pfrag_file, default=0.0)
                    
                    # Get fragment mass
                    xyz_file = None
                    for xyz_name in ["fragment.xyz", "isomer.xyz"]:
                        if (frag_dir / xyz_name).exists():
                            xyz_file = frag_dir / xyz_name
                            break
                    
                    if xyz_file and intensity > 0:
                        # Calculate molecular mass
                        natoms, atomic_numbers, _ = utility.get_atomic_numbers_and_coords_py(xyz_file)
                        if natoms > 0 and atomic_numbers:
                            mass = utility.get_average_mol_mass_py(atomic_numbers)
                            
                            # Add to spectrum (sum intensities for same mass)
                            if mass in spectrum_data:
                                spectrum_data[mass] += intensity
                            else:
                                spectrum_data[mass] = intensity
        
        if not spectrum_data:
            print("  Warning: No fragment data found for spectrum generation")
            return False
        
        # Normalize intensities to 100%
        max_intensity = max(spectrum_data.values())
        if max_intensity > 0:
            for mass in spectrum_data:
                spectrum_data[mass] = (spectrum_data[mass] / max_intensity) * 100.0
        
        # Write spectrum file
        spectrum_file = Path(env.startdir) / "allspec.dat"
        with open(spectrum_file, 'w') as f:
            f.write("# Mass spectrum from QCxMS2 simulation\n")
            f.write("# m/z  Intensity(%)\n")
            
            # Sort by mass
            for mass in sorted(spectrum_data.keys()):
                f.write(f"{mass:.1f}  {spectrum_data[mass]:.2f}\n")
        
        print(f"  Mass spectrum written to {spectrum_file}")
        print(f"  Spectrum contains {len(spectrum_data)} peaks")
        
        return True
        
    except Exception as e:
        print(f"  Error generating mass spectrum: {e}")
        return False


def _generate_simulation_summary_py(env: RunTypeData) -> bool:
    """
    Generate comprehensive simulation summary report.
    """
    try:
        print("Generating simulation summary report...")
        
        summary_file = Path(env.startdir) / "simulation_summary.txt"
        
        with open(summary_file, 'w') as f:
            f.write("QCxMS2 Simulation Summary Report\n")
            f.write("=" * 50 + "\n\n")
            
            # Simulation parameters
            f.write("Simulation Parameters:\n")
            f.write("-" * 20 + "\n")
            f.write(f"Input molecule: {env.infile}\n")
            f.write(f"Charge: {env.chrg}\n")
            f.write(f"Mode: {env.mode}\n")
            f.write(f"Temperature: {env.temp} K\n")
            f.write(f"Geometry level: {env.geolevel}\n")
            f.write(f"TS level: {env.tslevel}\n")
            f.write(f"IP level: {env.iplevel}\n")
            if hasattr(env, 'ip2level') and env.ip2level:
                f.write(f"IP2 level: {env.ip2level}\n")
            f.write(f"Max fragmentation levels: {env.nfrag}\n")
            f.write(f"Probability threshold: {env.pthr}%\n")
            f.write(f"Number of cores: {env.cores}\n")
            f.write(f"TS finder: {env.tsfinder}\n")
            
            if env.mode.lower() == "cid":
                f.write(f"CID collision energy: {getattr(env, 'cid_elab', 'N/A')} eV\n")
                f.write(f"CID target gas: {getattr(env, 'cid_target_gas', 'N/A')}\n")
            
            f.write("\n")
            
            # Fragment statistics
            allfrags_file = Path(env.startdir) / "allfrags.dat"
            if allfrags_file.exists():
                with open(allfrags_file, 'r') as af:
                    lines = af.readlines()
                    if lines:
                        total_frags = int(lines[0].strip())
                        f.write("Fragment Statistics:\n")
                        f.write("-" * 20 + "\n")
                        f.write(f"Total unique fragments: {total_frags}\n")
                        
                        # Count fragments by level
                        level_counts = {}
                        for line in lines[1:]:
                            frag_path = line.strip()
                            if frag_path:
                                level = frag_path.count('p')  # Simple level estimation
                                level_counts[level] = level_counts.get(level, 0) + 1
                        
                        for level in sorted(level_counts.keys()):
                            f.write(f"  Level {level}: {level_counts[level]} fragments\n")
                        
                        f.write("\n")
            
            # Spectrum statistics
            spectrum_file = Path(env.startdir) / "allspec.dat"
            if spectrum_file.exists():
                masses = []
                intensities = []
                with open(spectrum_file, 'r') as sf:
                    for line in sf:
                        if not line.startswith('#') and line.strip():
                            parts = line.split()
                            if len(parts) >= 2:
                                try:
                                    masses.append(float(parts[0]))
                                    intensities.append(float(parts[1]))
                                except ValueError:
                                    continue
                
                if masses and intensities:
                    f.write("Mass Spectrum Statistics:\n")
                    f.write("-" * 25 + "\n")
                    f.write(f"Number of peaks: {len(masses)}\n")
                    f.write(f"Mass range: {min(masses):.1f} - {max(masses):.1f} m/z\n")
                    f.write(f"Base peak: {masses[intensities.index(max(intensities))]:.1f} m/z\n")
                    f.write(f"Molecular ion: {max(masses):.1f} m/z\n")
                    
                    # Count significant peaks (>5% intensity)
                    significant_peaks = sum(1 for i in intensities if i > 5.0)
                    f.write(f"Significant peaks (>5%): {significant_peaks}\n")
                    f.write("\n")
            
            # File listing
            f.write("Generated Files:\n")
            f.write("-" * 15 + "\n")
            
            important_files = [
                "allfrags.dat", "allspec.dat", "simulation_summary.txt",
                "mass_spectrum.png", "energy_distribution.png",
                "fragmentation_tree.png", "results_summary.html"
            ]
            
            for filename in important_files:
                filepath = Path(env.startdir) / filename
                if filepath.exists():
                    size_kb = filepath.stat().st_size / 1024
                    f.write(f"  {filename} ({size_kb:.1f} KB)\n")
            
            # Count reaction directories
            reaction_dirs = list(Path(env.startdir).glob("p*"))
            reaction_dirs = [d for d in reaction_dirs if d.is_dir()]
            if reaction_dirs:
                f.write(f"  {len(reaction_dirs)} reaction directories (p*)\n")
            
            f.write("\n")
            f.write("End of Summary\n")
        
        print(f"  Simulation summary written to {summary_file}")
        return True
        
    except Exception as e:
        print(f"  Error generating simulation summary: {e}")
        return False


if __name__ == "__main__":
    main()
