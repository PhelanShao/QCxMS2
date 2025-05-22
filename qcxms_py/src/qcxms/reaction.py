import math
import os
from pathlib import Path
from typing import Tuple, Optional, List, Dict, Union

try:
    from .data import RunTypeData, WP
    from . import iomod
    from . import qmmod
    from . import utility
    from .constants import AUTOEV, KB_EV_K # Assuming these are in a constants module
except ImportError:
    # Fallbacks for standalone/testing
    print("Attempting to import dummy/mock modules for reaction.py standalone run.")
    from data import RunTypeData, WP # type: ignore
    import iomod_mock as iomod # type: ignore
    import qmmod_mock as qmmod # type: ignore
    import utility_mock as utility # type: ignore
    AUTOEV = 27.211386245988
    KB_EV_K = 8.6173332621415e-5


def calculate_reaction_energies_py(
    env: RunTypeData, 
    npairs: int, 
    fragdirs: List[List[str]], # [[pair_id, frag1_dir, frag2_dir_or_empty],...]
    # work_dir_base: Path # Base path where all these calculations are happening
) -> Dict[str, float]: # Returns a dict mapping pair_id to its reaction energy (dE_eV)
    """
    Calculates reaction energies (Delta E) for isomers or fragment pairs.
    Includes electronic energies and optional RRHO thermal corrections.
    Writes 'de_<level>' and 'sumreac_<level>' files in product directories.
    Operates relative to the current working directory, which should be the precursor's directory.
    """
    print(f"Calculating reaction energies for {npairs} products at level {env.tslevel}...")
    original_path = Path.cwd()
    reaction_energies_eV: Dict[str, float] = {}

    # 1. Get Precursor Energy
    precursor_fname = "isomer.xyz" if Path("isomer.xyz").exists() else "fragment.xyz"
    if not Path(precursor_fname).exists():
        print(f"Error: Precursor geometry file ({precursor_fname}) not found in {original_path}.")
        return reaction_energies_eV # Empty

    precursor_chrg, precursor_uhf = qmmod._get_chrg_uhf(env, precursor_fname) # Assumes _get_chrg_uhf is accessible

    qm_job_type = "tddft" if env.exstates > 0 else "sp"
    precursor_query = qmmod._construct_qmdata_query(env.tslevel, qm_job_type, precursor_chrg, precursor_uhf)
    
    cached_precursor, e_start_au = qmmod._check_qmdata_cache(precursor_query)
    if not cached_precursor or e_start_au is None:
        print(f"ERROR: Energy for precursor {precursor_fname} ({precursor_query}) not found in qmdata.")
        # Optionally, run the calculation here if missing
        # For now, assuming it must exist from a prior step.
        return reaction_energies_eV
    
    rrho_start_au: WP = 0.0
    if env.bhess:
        precursor_bhess_query = qmmod._construct_qmdata_query(env.geolevel, "bhess", precursor_chrg, precursor_uhf)
        cached_bhess, rrho_start_au_val = qmmod._check_qmdata_cache(precursor_bhess_query)
        if not cached_bhess or rrho_start_au_val is None:
            print(f"ERROR: RRHO (bhess) energy for precursor {precursor_fname} ({precursor_bhess_query}) not found.")
            return reaction_energies_eV
        rrho_start_au = rrho_start_au_val

    e_start_total_au = e_start_au + rrho_start_au

    # Read sumreac from precursor directory (sum of reaction energies up to this point)
    sumreac_previous_eV = iomod.rdshort_real(original_path / f"sumreac_{env.tslevel}", default=0.0)

    # 2. Calculate Product Energies (SP/TDDFT and optionally BHESS)
    # This involves collecting all necessary QM jobs first, then running them.
    
    jobs_to_run: List[Tuple[Path, str, str, str, int, Optional[int]]] = [] # dir, fname, level, job, chrg, uhf
    
    # Stage 1: SP/TDDFT energies for all products
    for i in range(npairs):
        pair_id, item1_dir_str, item2_dir_str = fragdirs[i]
        item_dirs_abs = [Path(item1_dir_str)]
        if item2_dir_str: item_dirs_abs.append(Path(item2_dir_str))

        for prod_dir_abs in item_dirs_abs:
            prod_fname = "isomer.xyz" if not item2_dir_str else "fragment.xyz"
            # Charge/UHF for products are determined by their own .CHRG/.UHF or default to precursor if missing
            # For products, we typically use their assigned charge after IP calculation
            # This needs careful handling if .CHRG/.UHF are not yet in product dirs.
            # Assuming they are, or _get_chrg_uhf handles it.
            os.chdir(prod_dir_abs)
            prod_chrg, prod_uhf = qmmod._get_chrg_uhf(env, prod_fname) # Reads/writes .CHRG/.UHF in prod_dir
            os.chdir(original_path)

            prod_query = qmmod._construct_qmdata_query(env.tslevel, qm_job_type, prod_chrg, prod_uhf)
            is_cached, _ = qmmod._check_qmdata_cache(prod_query) # Check in prod_dir's qmdata
            if not is_cached:
                jobs_to_run.append((prod_dir_abs, prod_fname, env.tslevel, qm_job_type, prod_chrg, prod_uhf))
    
    # Run collected SP/TDDFT jobs (sequentially for now)
    print(f"Need to run {len(jobs_to_run)} SP/TDDFT calculations for products at {env.tslevel}...")
    for work_dir, fname_in_dir, level, job_type, chrg_val, uhf_val in jobs_to_run:
        os.chdir(work_dir)
        cmd_list, out_f, patt, clean_cmd, _, _ = qmmod.prepqm_py(env, fname_in_dir, level, job_type, chrg_val, uhf_val)
        if cmd_list:
            # print(f"  Running in {work_dir}: {' '.join(cmd_list)}")
            res = iomod.execute_command(cmd_list, cwd=work_dir, shell=False) # shell=False for list
            if res.returncode == 0:
                qmmod.readoutqm_py(env, fname_in_dir, level, job_type, out_f, patt)
            else:
                print(f"    ERROR: QM calculation failed in {work_dir}. Stdout:\n{res.stdout}\nStderr:\n{res.stderr}")
            if clean_cmd: iomod.execute_command(clean_cmd.split(), cwd=work_dir, shell=True)
        os.chdir(original_path)
    
    # Stage 2: BHESS energies for all products if env.bhess
    bhess_jobs_to_run: List[Tuple[Path, str, str, str, int, Optional[int]]] = []
    if env.bhess:
        for i in range(npairs):
            pair_id, item1_dir_str, item2_dir_str = fragdirs[i]
            item_dirs_abs = [Path(item1_dir_str)]
            if item2_dir_str: item_dirs_abs.append(Path(item2_dir_str))

            for prod_dir_abs in item_dirs_abs:
                prod_fname = "isomer.xyz" if not item2_dir_str else "fragment.xyz"
                os.chdir(prod_dir_abs)
                prod_chrg, prod_uhf = qmmod._get_chrg_uhf(env, prod_fname)
                os.chdir(original_path)
                
                bhess_query = qmmod._construct_qmdata_query(env.geolevel, "bhess", prod_chrg, prod_uhf)
                is_cached, _ = qmmod._check_qmdata_cache(bhess_query) # Check in prod_dir
                if not is_cached:
                    bhess_jobs_to_run.append((prod_dir_abs, prod_fname, env.geolevel, "bhess", prod_chrg, prod_uhf))

        print(f"Need to run {len(bhess_jobs_to_run)} BHESS calculations for products at {env.geolevel}...")
        for work_dir, fname_in_dir, level, job_type, chrg_val, uhf_val in bhess_jobs_to_run:
            os.chdir(work_dir)
            cmd_list, out_f, patt, clean_cmd, _, _ = qmmod.prepqm_py(env, fname_in_dir, level, job_type, chrg_val, uhf_val)
            if cmd_list:
                # print(f"  Running BHESS in {work_dir}: {' '.join(cmd_list)}")
                res = iomod.execute_command(cmd_list, cwd=work_dir, shell=False)
                if res.returncode == 0:
                    qmmod.readoutqm_py(env, fname_in_dir, level, job_type, out_f, patt)
                else:
                    print(f"    ERROR: BHESS calculation failed in {work_dir}. Stdout:\n{res.stdout}\nStderr:\n{res.stderr}")
                if clean_cmd: iomod.execute_command(clean_cmd.split(), cwd=work_dir, shell=True)
            os.chdir(original_path)

    # Stage 3: Read all energies and calculate reaction energies
    for i in range(npairs):
        pair_id, item1_dir_str, item2_dir_str = fragdirs[i]
        item1_dir = Path(item1_dir_str)
        
        e_products_au: WP = 0.0
        rrho_products_au: WP = 0.0
        products_valid = True

        # Product 1 (isomer or first fragment of a pair)
        os.chdir(item1_dir)
        prod1_fname = "isomer.xyz" if not item2_dir_str else "fragment.xyz"
        prod1_chrg, prod1_uhf = qmmod._get_chrg_uhf(env, prod1_fname)
        
        prod1_query = qmmod._construct_qmdata_query(env.tslevel, qm_job_type, prod1_chrg, prod1_uhf)
        cached_prod1, e_f1_au = qmmod._check_qmdata_cache(prod1_query)
        if not cached_prod1 or e_f1_au is None: products_valid = False; print(f"Energy for {item1_dir.name} not found!")
        else: e_products_au += e_f1_au

        if env.bhess and products_valid:
            prod1_bhess_query = qmmod._construct_qmdata_query(env.geolevel, "bhess", prod1_chrg, prod1_uhf)
            cached_bhess1, rrho_f1_au_val = qmmod._check_qmdata_cache(prod1_bhess_query)
            if not cached_bhess1 or rrho_f1_au_val is None: products_valid = False; print(f"BHESS for {item1_dir.name} not found!")
            else: rrho_products_au += rrho_f1_au_val
        os.chdir(original_path)

        # Product 2 (if it's a pair)
        if item2_dir_str and products_valid:
            item2_dir = Path(item2_dir_str)
            os.chdir(item2_dir)
            prod2_fname = "fragment.xyz"
            prod2_chrg, prod2_uhf = qmmod._get_chrg_uhf(env, prod2_fname)

            prod2_query = qmmod._construct_qmdata_query(env.tslevel, qm_job_type, prod2_chrg, prod2_uhf)
            cached_prod2, e_f2_au = qmmod._check_qmdata_cache(prod2_query)
            if not cached_prod2 or e_f2_au is None: products_valid = False; print(f"Energy for {item2_dir.name} not found!")
            else: e_products_au += e_f2_au
            
            if env.bhess and products_valid:
                prod2_bhess_query = qmmod._construct_qmdata_query(env.geolevel, "bhess", prod2_chrg, prod2_uhf)
                cached_bhess2, rrho_f2_au_val = qmmod._check_qmdata_cache(prod2_bhess_query)
                if not cached_bhess2 or rrho_f2_au_val is None: products_valid = False; print(f"BHESS for {item2_dir.name} not found!")
                else: rrho_products_au += rrho_f2_au_val
            os.chdir(original_path)
        
        if not products_valid:
            reaction_energies_eV[pair_id] = 10000.0 # High energy for failed products
            de_eV = 10000.0
            drrho_eV = 0.0
        else:
            e_products_total_au = e_products_au + rrho_products_au
            de_eV = (e_products_total_au - e_start_total_au) * AUTOEV
            reaction_energies_eV[pair_id] = de_eV
            
            # Write de_<level> and drrho_<geolevel> (if bhess)
            # These are written in the pair_id directory (fragdirs[i][0])
            pair_dir_path = Path(pair_id) # This assumes pair_id is a path string like "p0"
            pair_dir_path.mkdir(exist_ok=True) # Ensure pair directory exists

            de_file_path = pair_dir_path / f"de_{env.tslevel}"
            iomod.wrshort_real(de_file_path, de_eV)

            if env.bhess:
                drrho_au = rrho_products_au - rrho_start_au
                drrho_eV = drrho_au * AUTOEV
                drrho_file_path = pair_dir_path / f"drrho_{env.geolevel}"
                iomod.wrshort_real(drrho_file_path, drrho_eV)
        
        # Update and write sumreac for products
        current_sumreac_eV = sumreac_previous_eV + de_eV
        
        sumreac_file_path_prod1 = item1_dir / f"sumreac_{env.tslevel}"
        iomod.wrshort_real(sumreac_file_path_prod1, current_sumreac_eV)
        if item2_dir_str:
            sumreac_file_path_prod2 = Path(item2_dir_str) / f"sumreac_{env.tslevel}"
            iomod.wrshort_real(sumreac_file_path_prod2, current_sumreac_eV)
        # Also write to pair_id directory
        iomod.wrshort_real(Path(pair_id) / f"sumreac_{env.tslevel}", current_sumreac_eV)


    print("Reaction energy calculation complete.")
    return reaction_energies_eV


def sort_out_high_energy_fragments_py(
    env: RunTypeData,
    sort_level: str, # Level of theory used for de values (usually env.tslevel)
    # etots_in: List[WP], # Not needed if reading de from files
    use_rrho: bool, # Corresponds to env.bhess
    npairs_in: int,
    fragdirs_in: List[List[str]], # [[pair_id, frag1_dir, frag2_dir_or_empty],...]
    scale_iee_factor: WP # Factor to scale average IEE for thresholding
) -> Tuple[int, List[List[str]]]:
    """
    Filters out fragment pairs that are too high in reaction energy.
    Operates in the precursor's directory.
    """
    print("Sorting out fragment pairs with excessively high reaction energies...")
    original_path = Path.cwd()
    
    valid_fragdirs_out: List[List[str]] = []
    
    prob_threshold = 0.01 # Fortran: pthr = 0.01_wp

    # Get precursor info
    precursor_fname = "isomer.xyz" if Path("isomer.xyz").exists() else "fragment.xyz"
    if not Path(precursor_fname).exists():
        print(f"Error: Precursor file {precursor_fname} not found in {original_path} for sorting.")
        return 0, []
    
    num_atoms_precursor = iomod.rdshort_int(precursor_fname) # Reads from first line of XYZ
    if num_atoms_precursor == 0:
        print(f"Error: Could not read atom count from {precursor_fname}.")
        return 0, []

    avg_internal_energy_iee = env.ieeatm * num_atoms_precursor # Average IEE of precursor
    
    num_vib_modes = 3 * num_atoms_precursor - 6
    if num_vib_modes <= 0: num_vib_modes = 1 # Avoid non-positive for diatomic/atomic

    # sumreac for the precursor (energy consumed to form it)
    sumreac_precursor_eV = iomod.rdshort_real(original_path / f"sumreac_{sort_level}", default=0.0)
    
    # Effective IEE available for further reactions
    effective_iee = avg_internal_energy_iee - sumreac_precursor_eV

    # Find H-dissociation scaling (needs findhdiss_py from utility)
    # is_h_dissociation_list = utility.findhdiss_py(env, precursor_fname, num_atoms_precursor, npairs_in, fragdirs_in)
    # For now, assume no special scaling for H-dissociation to simplify
    is_h_dissociation_list = [False] * npairs_in 
    h_diss_scale_factor = env.scaleeinthdiss

    de_min_eV = 5.0 # Default minimum reaction energy to consider (eV)
    # Find actual de_min among fragmentation reactions
    all_de_values_for_min_search: List[WP] = []
    for i in range(npairs_in):
        pair_id, item1_dir_str, item2_dir_str = fragdirs_in[i]
        if item2_dir_str: # Only consider fragmentations for de_min
            de_val = iomod.rdshort_real(Path(pair_id) / f"de_{sort_level}", default=10000.0)
            if de_val < 10000.0 : all_de_values_for_min_search.append(de_val)
    
    if all_de_values_for_min_search:
        de_min_eV = min(all_de_values_for_min_search)
        de_min_eV = max(de_min_eV, 0.0) # Ensure de_min is not negative for rate formula
        # Fortran had a hard cap: if (de_min .lt. 5.0_wp) de_min = 5.0_wp
        # This seems like a high floor if actual reactions are very exothermic.
        # Let's use a smaller floor or just max(de_min, 0.0)
        de_min_eV = max(de_min_eV, 0.1) # Ensure positive for log, if any rate formula needs it
        if env.mode == "cid" and de_min_eV < 5.0: de_min_eV = 5.0 # Specific CID floor
    
    print(f"  Effective IEE for sorting: {effective_iee:.2f} eV. Min reaction energy (de_min) considered: {de_min_eV:.2f} eV.")

    for i in range(npairs_in):
        pair_id, item1_dir_str, item2_dir_str = fragdirs_in[i]
        if not item1_dir_str : continue # Already removed

        de_eV = iomod.rdshort_real(Path(pair_id) / f"de_{sort_level}", default=10000.0)
        if de_eV >= 9999.0: # Problematic reaction, filter out
            if env.removedirs: iomod.rmrf(Path(pair_id)) # Also remove sub-fragment dirs
            continue

        iee_for_rate_calc = effective_iee * h_diss_scale_factor if is_h_dissociation_list[i] else effective_iee
        if iee_for_rate_calc <=0 : # No energy to drive reaction
            k_relative = 0.0
        elif env.eyring:
            # calc_releyring(energy, dbarrier, nvib)
            # dbarrier = de(i) - de_min
            k_relative = utility.calc_releyring(iee_for_rate_calc, de_eV - de_min_eV, num_vib_modes)
        else: # Simplified RRKM-like from Fortran
            # k_rel = EXP(-(nvib - 1)*(de(i) - de_min)/(IEE2*scaleIEE))
            # Note: scaleIEE is an input to this function.
            denom = iee_for_rate_calc * scale_iee_factor
            if denom == 0: k_relative = 0.0
            else:
                exponent = -(num_vib_modes - 1) * (de_eV - de_min_eV) / denom
                try: k_relative = math.exp(exponent)
                except OverflowError: k_relative = 0.0
        
        if k_relative < prob_threshold:
            print(f"  Filtering out pair {pair_id}: DE={de_eV:.2f} eV, k_rel={k_relative:.2e} < threshold {prob_threshold:.2e}")
            if env.removedirs:
                iomod.rmrf(Path(item1_dir_str))
                if item2_dir_str: iomod.rmrf(Path(item2_dir_str))
                iomod.rmrf(Path(pair_id)) # Remove the base pair_id dir
        else:
            valid_fragdirs_out.append(fragdirs_in[i])
            
    npairs_out = len(valid_fragdirs_out)
    print(f"  From {npairs_in} pairs, {npairs_out} remain after high energy filter.")
    return npairs_out, valid_fragdirs_out


# Other subroutines like optfragments, optproducts, checkproducttopo, etc. would follow a similar pattern
# of translating Fortran logic, using Python's file system and subprocess capabilities via iomod/qmmod,
# and managing data structures (lists of strings for directories, dictionaries for energies).
# The parallel execution aspects (omp_samejobcall) would be noted for future implementation
# using Python's multiprocessing.

if __name__ == '__main__':
    print("Testing reaction.py functions...")
    # Setup a dummy environment and file structure for testing
    # This would be more involved due to dependencies on qmdata files, .CHRG/.UHF, etc.
    # Example:
    # env_test = RunTypeData(tslevel="gfn2", geolevel="gfn2", bhess=True, ieeatm=0.8, chrg=1)
    # test_pair_dir = Path("test_reaction_dir/p0")
    # test_pair_dir.mkdir(parents=True, exist_ok=True)
    # (test_pair_dir / "de_gfn2").write_text("1.5\n") # 1.5 eV
    # ... (setup precursor qmdata, etc.) ...
    # np_out, frags_out = sort_out_high_energy_fragments_py(env_test, "gfn2", True, 1, [["p0", str(test_pair_dir/"f1"), str(test_pair_dir/"f2")]], 2.0)
    # print(f"Sort out result: Npairs={np_out}, Frags={frags_out}")
    # shutil.rmtree("test_reaction_dir", ignore_errors=True)
    pass

```
