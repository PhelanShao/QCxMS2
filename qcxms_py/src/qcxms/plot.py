import math
import os
from pathlib import Path
from typing import Tuple, Optional, List, Dict
from collections import defaultdict
import numpy as np # For easier array manipulations if allowed

try:
    from .data import RunTypeData, WP
    from . import iomod
    from . import utility
    from . import isotopes # For get_isotopic_pattern_py and get_average_mol_mass_py
    from .constants import EVTOKCAL # If energies from files are in eV and need kcal/mol for output
except ImportError:
    print("Attempting to import dummy/mock modules for plot.py standalone run.")
    from data import RunTypeData, WP # type: ignore
    import iomod_mock as iomod # type: ignore
    import utility_mock as utility # type: ignore
    import isotopes_mock as isotopes # type: ignore
    EVTOKCAL = 23.060547830618307


def _merge_spectrum_data(
    master_masses: List[WP], 
    master_intensities: List[WP], 
    new_masses: List[WP], 
    new_intensities: List[WP],
    mass_tolerance: WP = 1e-5 # Tolerance for considering masses identical
) -> Tuple[List[WP], List[WP]]:
    """
    Merges a new isotopic pattern into a master list of masses and intensities.
    If a new mass is very close to an existing mass, their intensities are summed.
    Otherwise, the new mass/intensity pair is added.
    This is a simplified version of Fortran's check_entries. A more robust way
    would be to collect all peaks and then bin/merge them once.
    """
    if not new_masses:
        return master_masses, master_intensities

    # For simplicity and robustness with float comparisons, this function assumes
    # that a final binning step will occur later. Here, we just append
    # and rely on the final binning to consolidate.
    # A direct translation of check_entries' loop-and-compare logic:
    
    current_master_len = len(master_masses)
    for i in range(len(new_masses)):
        nm = new_masses[i]
        ni = new_intensities[i]
        found_match = False
        for j in range(current_master_len): # Check against existing master entries
            if abs(master_masses[j] - nm) < mass_tolerance:
                master_intensities[j] += ni
                found_match = True
                break
        if not found_match:
            master_masses.append(nm)
            master_intensities.append(ni)
            # Update current_master_len for subsequent new_masses in this call,
            # though Fortran check_entries used a global list_masses and count_mass
            # which effectively grew.
            current_master_len = len(master_masses) 
            
    return master_masses, master_intensities

def _bin_spectrum(
    masses: List[WP], 
    intensities: List[WP], 
    bin_width: Optional[WP] = None, # If None, sums identical masses after potential rounding
    round_digits: Optional[int] = None # If int_masses, round_digits = 0
) -> Tuple[List[WP], List[WP]]:
    """Bins spectrum by summing intensities for masses that fall into the same bin."""
    if not masses:
        return [], []

    binned_spectrum: Dict[Union[int, float], float] = defaultdict(float)

    for m, i in zip(masses, intensities):
        bin_key: Union[int, float]
        if round_digits is not None:
            bin_key = round(m, round_digits)
            if round_digits == 0: bin_key = int(bin_key) # For integer masses
        elif bin_width is not None and bin_width > 0:
            bin_key = round(m / bin_width) * bin_width
        else: # Default: group by very close floating point numbers (e.g. after MC)
              # This requires pre-sorting and merging adjacent if close, or using a tolerance.
              # For now, if no bin_width/rounding, assume masses are distinct enough or already binned.
              # A simple way to handle this is to round to a high precision.
            bin_key = round(m, 5) # Default rounding for grouping if no other specified

        binned_spectrum[bin_key] += i
    
    sorted_masses = sorted(binned_spectrum.keys())
    final_masses = [WP(m) for m in sorted_masses]
    final_intensities = [WP(binned_spectrum[m]) for m in sorted_masses]
    
    return final_masses, final_intensities


def get_peaks_py(env: RunTypeData, num_total_frags: int, allfrag_paths_list: List[str]):
    """
    Aggregates isotopic patterns from all fragments, bins them, normalizes,
    and writes the final spectrum to peaks.dat/peaks.csv and allpeaks.dat.
    """
    original_cwd = Path.cwd()
    os.chdir(env.startdir) # Ensure we are in the main run directory

    print("\n" + "-"*37)
    print("| Reaction network generation finished |") # Fortran: Reaction network generation finished
    print("-"*37 + "\n")
    print("Normalizing Intensities for all fragments...")
    print(f"Total number of unique fragment/isomer paths considered: {num_total_frags}")

    # --- Step 1: Accumulate isotopic patterns from all fragments ---
    # Use a dictionary for accumulation, keys rounded to a certain precision
    # to handle minor float variations from MC simulations before explicit binning.
    # Precision for accumulation dictionary keys (e.g., 4-5 decimal places)
    # This replaces the Fortran check_entries logic more robustly.
    accumulated_raw_spectrum: Dict[float, float] = defaultdict(float)
    
    # Determine max_atoms for random number array (if used by isotopes.py)
    # Fortran read this from the initial input file.
    # nat_initial_mol, _, _ = utility.read_xyz_structure_py(env.infile) # Assuming this utility
    nat_initial_mol = iomod.rdshort_int(env.infile, default=0) # Placeholder if utility not ready

    # Process initial molecule (precursor of the whole cascade)
    print(f"  Processing isotopic pattern for initial molecule: {env.infile}")
    # The `charge_state` for get_isotopic_pattern_py should be the charge of the species
    # whose m/z is being calculated. For EI, this is typically +1 (or env.chrg).
    masses_init, intensities_init = isotopes.get_isotopic_pattern_py(
        env, env.infile, num_trials=env.nsamples, charge_state=env.chrg 
    )
    for m, inten_val in zip(masses_init, intensities_init):
        # Round mass to a defined precision for dictionary key to merge close peaks
        mass_key = round(m, 5) 
        accumulated_raw_spectrum[mass_key] += inten_val

    # Process each fragment/isomer from allfrags.dat
    for frag_path_str in allfrag_paths_list:
        frag_path = Path(frag_path_str)
        full_frag_path = Path(env.startdir) / frag_path
        
        # Determine if it's an isomer or fragment file
        current_xyz_fname = "isomer.xyz" if (full_frag_path / "isomer.xyz").exists() else "fragment.xyz"
        full_xyz_path = full_frag_path / current_xyz_fname

        if not full_xyz_path.exists():
            print(f"  Warning: XYZ file not found in {full_frag_path}. Skipping.")
            continue
        
        print(f"  Processing isotopic pattern for: {frag_path_str}/{current_xyz_fname}")
        # Charge for fragments should be read from their .CHRG file
        frag_chrg = iomod.rdshort_int(full_frag_path / ".CHRG", default=env.chrg)

        masses_frag, intensities_frag = isotopes.get_isotopic_pattern_py(
            env, str(full_xyz_path), num_trials=env.nsamples, charge_state=frag_chrg
        )
        for m, inten_val in zip(masses_frag, intensities_frag):
            mass_key = round(m, 5)
            accumulated_raw_spectrum[mass_key] += inten_val
            
    # Convert accumulated spectrum to lists
    temp_masses = list(accumulated_raw_spectrum.keys())
    temp_intensities = list(accumulated_raw_spectrum.values())
    
    # Sort by mass
    sorted_indices = np.argsort(temp_masses)
    binned_masses = [temp_masses[i] for i in sorted_indices]
    binned_intensities = [temp_intensities[i] for i in sorted_indices]


    # --- Step 2: Apply Integer Mass Rounding (if requested) ---
    if env.int_masses:
        print("  Rounding m/z values to integers.")
        binned_masses, binned_intensities = _bin_spectrum(binned_masses, binned_intensities, round_digits=0)

    # --- Step 3: Further Binning/Merging if a mass tolerance (env.mtol) is specified ---
    # The Fortran code's primary merging was for identical masses (after potential int rounding).
    # If env.mtol exists and is > small_epsilon, implement explicit binning here.
    # For now, assume the above rounding/dict accumulation handles most direct overlaps.
    # A more general binning:
    # if hasattr(env, 'mtol') and env.mtol > 1e-4: # Example tolerance
    #    print(f"  Binning peaks with mass tolerance: {env.mtol:.4f} Da")
    #    binned_masses, binned_intensities = _bin_spectrum(binned_masses, binned_intensities, bin_width=env.mtol)


    # --- Step 4: Intensity Normalization and Thresholding ---
    if not binned_intensities:
        print("  No peaks to process after accumulation/binning.")
        final_masses, final_intensities = [], []
    else:
        max_intensity = max(binned_intensities) if binned_intensities else 0.0
        if max_intensity > 0:
            normalized_intensities = [(i / max_intensity) * 10000.0 for i in binned_intensities]
        else:
            normalized_intensities = [0.0] * len(binned_intensities)

        final_masses: List[WP] = []
        final_intensities: List[WP] = []
        # Fortran env.intthr was % * 100. Here, if env.intthr is e.g. 0.1 for 0.1%, scale it for 10000 max.
        # Assuming env.intthr from parser is direct percentage (e.g. 1.0 for 1%)
        intensity_abs_threshold = (env.intthr / 100.0) * 10000.0 if env.intthr > 0 else 0.0
        
        for m, i_norm in zip(binned_masses, normalized_intensities):
            if m >= env.mthr and i_norm >= intensity_abs_threshold:
                final_masses.append(m)
                final_intensities.append(i_norm)
    
    # --- Step 5: Output peaks.dat and peaks.csv ---
    print("\nComputed spectrum (peaks.dat/peaks.csv):")
    print("(m/z | intensity - normalized to 10000)")
    try:
        with open("peaks.dat", "w") as f_dat, open("peaks.csv", "w") as f_csv:
            f_csv.write("m/z,Intensity\n")
            for m, i in zip(final_masses, final_intensities):
                f_dat.write(f"{m:<10.6f} {i:<10.1f}\n")
                f_csv.write(f"{m:.6f},{i:.1f}\n")
                if i / 100.0 > 10.0: # Print if > 10% of max (i.e., > 1000 if max is 10000)
                    print(f"  {m:<10.6f} {i:<10.1f}")
    except IOError as e:
        print(f"Error writing peaks files: {e}")

    # --- Step 6: Output allpeaks.dat ---
    print("\nWriting all contributing fragment/isomer peaks (monoisotopic/most_intense) to allpeaks.dat...")
    all_peaks_data: List[Tuple[str, WP, WP]] = []
    
    # Temporarily set env.noiso for get_isotopic_pattern_py
    original_noiso_setting = env.noiso
    env.noiso = True 

    # Initial molecule
    masses_init_noiso, intensities_init_noiso = isotopes.get_isotopic_pattern_py(
        env, env.infile, num_trials=env.nsamples, charge_state=env.chrg
    )
    if masses_init_noiso: # Should be a single peak
        all_peaks_data.append(("input_structure", masses_init_noiso[0], intensities_init_noiso[0]))

    # Fragments/isomers
    for frag_path_str in allfrag_paths_list:
        frag_path = Path(frag_path_str)
        full_frag_path = Path(env.startdir) / frag_path
        current_xyz_fname = "isomer.xyz" if (full_frag_path / "isomer.xyz").exists() else "fragment.xyz"
        full_xyz_path = full_frag_path / current_xyz_fname

        if full_xyz_path.exists():
            frag_chrg = iomod.rdshort_int(full_frag_path / ".CHRG", default=env.chrg)
            masses_frag_noiso, intensities_frag_noiso = isotopes.get_isotopic_pattern_py(
                env, str(full_xyz_path), num_trials=env.nsamples, charge_state=frag_chrg
            )
            if masses_frag_noiso:
                 all_peaks_data.append((frag_path_str, masses_frag_noiso[0], intensities_frag_noiso[0]))
        else:
            print(f"  Warning: XYZ for {frag_path_str} not found for allpeaks.dat generation.")

    env.noiso = original_noiso_setting # Restore setting

    # Normalize allpeaks.dat intensities
    max_allpeaks_intensity = max(item[2] for item in all_peaks_data) if all_peaks_data else 0.0
    
    try:
        with open("allpeaks.dat", "w") as f:
            f.write("Fragment_Path | Representative_Mass | Relative_Intensity_Contribution\n")
            for path_str, mass_val, inten_val in all_peaks_data:
                norm_inten = (inten_val / max_allpeaks_intensity) * 10000.0 if max_allpeaks_intensity > 0 else 0.0
                f.write(f"{path_str:<40} {mass_val:<15.6f} {norm_inten:<20.1f}\n")
                if norm_inten / 100.0 > 10.0: # Print if > 10%
                    print(f"  {path_str:<40} {mass_val:<15.6f} {norm_inten:<20.1f}")
    except IOError as e:
        print(f"Error writing allpeaks.dat: {e}")
        
    os.chdir(original_cwd) # Change back to original directory if chdir happened


def collect_fragments_info_py(
    env: RunTypeData, 
    npairs: int, 
    fragdirs: List[List[str]], # [[pair_id, frag1_dir, frag2_dir_or_empty],...]
    output_fname: str,
    precursor_cascade_path: str # Path of the current precursor, e.g. "p0/f1"
) -> Tuple[List[str], int]:
    """
    Collects information about generated fragments/isomers and writes to a summary file.
    Also prepares a list of fragment paths for the next fragmentation level.
    """
    print(f"Writing fragment/isomer summary to: {output_fname}")
    
    next_level_frag_paths: List[str] = []
    num_collected_frags = 0
    
    original_cwd = Path.cwd() # Should be precursor's directory

    try:
        with open(output_fname, "w") as f: # Overwrite if exists for this level
            f.write("Dir  Product_Type  SumReac(kcal/mol)  DeltaE(kcal/mol)  Barrier(kcal/mol)  IRC(cm-1)\n")
            f.write("  Fragment_Path  Mass  SumFormula  Rel_Intensity(%)\n")
            f.write("-" * 80 + "\n")

            for i in range(npairs):
                pair_id_str, item1_dir_str, item2_dir_str = fragdirs[i]
                
                # Paths relative to current precursor dir
                pair_dir_path_local = Path(pair_id_str) 
                item1_path_local = Path(item1_dir_str)
                
                # Paths relative to env.startdir (for global reference in fraglist)
                item1_global_path_str = str(Path(precursor_cascade_path) / item1_dir_str).replace("\\", "/")

                # Read reaction data from the pair_id directory
                # Energies are in eV, convert to kcal/mol for output
                barrier_ev = iomod.rdshort_real(pair_dir_path_local / f"barrier_{env.tslevel}", default=999.0)
                sumreac_ev = iomod.rdshort_real(pair_dir_path_local / f"sumreac_{env.tslevel}", default=0.0)
                de_ev = iomod.rdshort_real(pair_dir_path_local / f"de_{env.tslevel}", default=999.0)
                irc_cm1 = iomod.rdshort_real(pair_dir_path_local / f"ircmode_{env.geolevel}", default=0.0)

                is_fragment_pair = bool(item2_dir_str)
                product_type = "FragmentPair" if is_fragment_pair else "Isomer"
                
                f.write(f"{Path(precursor_cascade_path) / pair_id_str} {product_type} "
                        f"{sumreac_ev * EVTOKCAL:.1f} {de_ev * EVTOKCAL:.1f} "
                        f"{barrier_ev * EVTOKCAL:.1f} {irc_cm1:.1f}\n")

                # Process item1 (isomer or first fragment)
                os.chdir(item1_path_local) # Enter item1's directory
                mass1 = iomod.rdshort_real("mass", default=0.0)
                pfrag1 = iomod.rdshort_real("pfrag", default=0.0)
                # sumform1 = utility.get_sum_formula_from_xyz_py("fragment.xyz" if is_fragment_pair else "isomer.xyz")
                sumform1 = utility.get_sum_formula_from_xyz_py(item1_path_local.name + ".xyz") # Needs fix
                os.chdir(original_cwd) # Back to precursor dir
                
                f.write(f"  {item1_global_path_str} {mass1:.3f} {sumform1} {pfrag1 / env.nsamples * 100.0:.1f}\n")
                if pfrag1 / env.nsamples * 100.0 >= env.pthr: # pthr is %
                    next_level_frag_paths.append(item1_global_path_str)
                    if env.printlevel >= 1 and product_type=="Isomer": print(f"  Isomer: {item1_global_path_str}, Mass: {mass1:.3f}, Formula: {sumform1}, Intensity: {pfrag1 / env.nsamples * 100.0:.1f}%")


                if is_fragment_pair:
                    item2_path_local = Path(item2_dir_str)
                    item2_global_path_str = str(Path(precursor_cascade_path) / item2_dir_str).replace("\\", "/")
                    os.chdir(item2_path_local)
                    mass2 = iomod.rdshort_real("mass", default=0.0)
                    pfrag2 = iomod.rdshort_real("pfrag", default=0.0)
                    sumform2 = utility.get_sum_formula_from_xyz_py("fragment.xyz")
                    os.chdir(original_cwd)
                    
                    f.write(f"  {item2_global_path_str} {mass2:.3f} {sumform2} {pfrag2 / env.nsamples * 100.0:.1f}\n")
                    if pfrag2 / env.nsamples * 100.0 >= env.pthr:
                         next_level_frag_paths.append(item2_global_path_str)
                    if env.printlevel >= 1 and (pfrag1 / env.nsamples * 100.0 >= env.pthr or pfrag2 / env.nsamples * 100.0 >= env.pthr) :
                         print(f"  Pair: {Path(precursor_cascade_path) / pair_id_str}")
                         print(f"    Frag1: {item1_global_path_str}, Mass: {mass1:.3f}, Formula: {sumform1}, Intensity: {pfrag1 / env.nsamples * 100.0:.1f}%")
                         print(f"    Frag2: {item2_global_path_str}, Mass: {mass2:.3f}, Formula: {sumform2}, Intensity: {pfrag2 / env.nsamples * 100.0:.1f}%")

                f.write("-" * 80 + "\n")
                
    except IOError as e:
        print(f"Error writing summary file {output_fname}: {e}")

    num_collected_frags = len(next_level_frag_paths)
    return next_level_frag_paths, num_collected_frags


if __name__ == '__main__':
    print("Testing plot.py functions...")
    # Needs extensive mocking of env, file system, and other modules (isotopes, utility, iomod)
    # Example structure for testing get_peaks_py:
    # 1. Create dummy env
    # 2. Create dummy allfrags.dat
    # 3. Create dummy fragment directories listed in allfrags.dat
    #    - Each needs env.infile (e.g. fragment.xyz)
    #    - Each needs pfrag file
    #    - isotopes.get_isotopic_pattern_py needs to be callable (mocked or real)
    # 4. Call get_peaks_py
    # 5. Check output files peaks.dat, peaks.csv, allpeaks.dat

    # Example for _bin_spectrum
    masses = [100.01, 100.03, 100.05, 101.02, 101.04, 102.00]
    intensities = [10, 20, 5, 50, 15, 30]
    binned_m, binned_i = _bin_spectrum(masses, intensities, bin_width=0.1)
    print("\nBinned Spectrum (width 0.1):")
    for m,i_val in zip(binned_m, binned_i): print(f"  m/z: {m:.2f}, Int: {i_val:.1f}")
    
    binned_m_int, binned_i_int = _bin_spectrum(masses, intensities, round_digits=0)
    print("\nBinned Spectrum (integer masses):")
    for m,i_val in zip(binned_m_int, binned_i_int): print(f"  m/z: {m:.0f}, Int: {i_val:.1f}")
    pass

```
