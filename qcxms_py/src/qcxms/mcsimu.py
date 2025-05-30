import os
import math
import re
import shutil
import subprocess
from pathlib import Path
from typing import List, Tuple, Dict, Optional, Union
import random # Ensure random is imported

try:
    from .data import RunTypeData, WP
    from . import iomod
    from . import utility
    from . import qmmod 
    from . import cid 
    from .constants import AUTOEV, EVTOKCAL, HARTREE_TO_EV, KB_EV_K, PI
except ImportError:
    print("Attempting to import dummy/mock modules for mcsimu.py standalone run.")
    from data import RunTypeData, WP # type: ignore
    import iomod_mock as iomod # type: ignore
    import utility_mock as utility # type: ignore
    import qmmod_mock as qmmod # type: ignore
    import cid_mock as cid # type: ignore
    AUTOEV = 27.211386245988 
    EVTOKCAL = 23.060547830618307
    HARTREE_TO_EV = AUTOEV 
    KB_EV_K = 8.6173332621415e-5
    PI = math.pi


def _readout_thermo_py(thermo_out_path: Path, num_expected_temps: int) -> List[WP]:
    grrho_values_hartree: List[WP] = []
    if not thermo_out_path.exists():
        print(f"  Warning: Thermo output file {thermo_out_path} not found for readout.")
        return [0.0] * num_expected_temps
    try:
        with open(thermo_out_path, 'r') as f: lines = f.readlines()
        found_table = False
        temp_count = 0
        for line in lines:
            line_s = line.strip()
            if not line_s: continue
            if "T/K" in line_s and "S" in line_s and "H" in line_s and "U" in line_s and "G(RRHO)" in line_s:
                found_table = True
                continue
            if found_table:
                parts = line_s.split()
                if len(parts) >= 5 and parts[0][0].isdigit(): 
                    try:
                        grrho_values_hartree.append(float(parts[4])) 
                        temp_count +=1
                        if temp_count >= num_expected_temps: break 
                    except (ValueError, IndexError): continue 
        if temp_count < num_expected_temps:
             print(f"  Warning: Could not parse enough G(RRHO) values from {thermo_out_path}. Expected {num_expected_temps}, Found {temp_count}.")
             padding_value = grrho_values_hartree[-1] if grrho_values_hartree else 1e10 
             grrho_values_hartree.extend([padding_value] * (num_expected_temps - temp_count))
    except Exception as e:
        print(f"  Error reading or parsing thermo output file {thermo_out_path}: {e}")
        return [1e10] * num_expected_temps 
    return grrho_values_hartree

def _xtb_thermo_py(
    env: RunTypeData, num_total_increments: int, num_vib_dof: int, 
    max_internal_energy_ev: float, molecule_xyz_path: Path, 
    is_ts_or_fragment: bool, imag_freq_cutoff_override: Optional[float] = None
) -> List[WP]:
    if not molecule_xyz_path.exists(): return [1e10] * num_total_increments 
    ithr = imag_freq_cutoff_override if imag_freq_cutoff_override is not None else env.ithr
    all_grrho_values: List[WP] = []
    xtb_gfn_level = env.geolevel_for_xtb_thermo if hasattr(env, 'geolevel_for_xtb_thermo') and env.geolevel_for_xtb_thermo else ("gfn1" if env.geolevel.lower() == "gfn1" else "gfn2")
    temperatures_k: List[float] = []
    if num_total_increments == 1 and max_internal_energy_ev == 0 : temperatures_k.append(float(env.temp))
    elif num_total_increments > 0 and num_vib_dof > 0:
        energy_step_ev = max_internal_energy_ev / num_total_increments
        for i_inc in range(1, num_total_increments + 1): temperatures_k.append(max(1.0, utility.calctemp_py(num_vib_dof, i_inc * energy_step_ev)))
    else: 
        temperatures_k.append(max(1.0, float(env.temp)))
        if num_total_increments == 0 and max_internal_energy_ev > 0: temperatures_k = [max(1.0, utility.calctemp_py(num_vib_dof, max_internal_energy_ev))]
        elif num_total_increments > 1: num_total_increments = 1
    num_total_increments = len(temperatures_k)
    if num_total_increments == 0: return []
    batch_size = 100 
    num_batches = (num_total_increments + batch_size - 1) // batch_size
    current_working_dir = Path.cwd() 
    for batch_num in range(num_batches):
        start_idx, end_idx = batch_num * batch_size, min((batch_num + 1) * batch_size, num_total_increments)
        current_batch_temps = temperatures_k[start_idx:end_idx]
        if not current_batch_temps: continue
        temp_str = ",".join(f"{t:.2f}" for t in current_batch_temps)
        thermo_out_file = current_working_dir / f"thermo_batch_{molecule_xyz_path.stem}_{batch_num}.out"
        xtb_cmd = ["xtb", str(molecule_xyz_path.resolve()), "--thermo", temp_str, f"--gfn{xtb_gfn_level.replace('gfn', '')}", "--sthr", str(env.sthr), "--ithr", str(ithr)]
        sph_hessian_path = current_working_dir / "hessian_sph" 
        if is_ts_or_fragment and sph_hessian_path.exists(): xtb_cmd.extend(["--bhess", str(sph_hessian_path.resolve())])
        else: xtb_cmd.append("--hess") 
        with open(thermo_out_file, "w") as f_out, open(f"{thermo_out_file}.err", "w") as f_err: result = iomod.execute_command(xtb_cmd, stdout_to=f_out, stderr_to=f_err, cwd=current_working_dir)
        if result.returncode == 0: all_grrho_values.extend(_readout_thermo_py(thermo_out_file, len(current_batch_temps)))
        else: all_grrho_values.extend([1e10] * len(current_batch_temps)) 
        thermo_out_file.unlink(missing_ok=True); (thermo_out_file.parent / f"{thermo_out_file.name}.err").unlink(missing_ok=True)
        default_hess_file = current_working_dir / "hessian"
        if default_hess_file.exists(): 
            if is_ts_or_fragment and not sph_hessian_path.exists(): shutil.move(str(default_hess_file), str(sph_hessian_path))
            else: default_hess_file.unlink(missing_ok=True)
        for f_pattern in ["gfnff_parameters", "xtbhess.pc", "xtb_thermodata.dat", "wbo", ".xtb_thermodata.tmp*"]:
            for f_path in current_working_dir.glob(f_pattern): iomod.remove_file(f_path)
    if len(all_grrho_values) != num_total_increments:
        padding_value = all_grrho_values[-1] if all_grrho_values else 1e10 
        all_grrho_values.extend([padding_value] * (num_total_increments - len(all_grrho_values)))
    return all_grrho_values

def montecarlo_py(
    env: RunTypeData, npairs: int, fragdirs: List[List[str]], 
    eiee_dist_energies: List[WP], eiee_dist_probs: List[WP], 
    fragmentation_level: int, p0_intensity_fraction: WP
) -> bool:
    precursor_dir_path = Path.cwd()
    print(f"\n--- Monte Carlo Simulation for precursor: {precursor_dir_path.name} (Level {fragmentation_level}) ---")
    print(f"  Number of reaction channels (npairs): {npairs}")
    tav0 = iomod.rdshort_real_py(precursor_dir_path.parent / "tav", default=0.0) if fragmentation_level > 1 else 0.0
    tf_us = env.tf 
    tf_reaction_s = ((tf_us - tav0 * 1e6) / (1.0 + env.tfscale * (fragmentation_level -1) ) * 1.0e-6) if env.tfscale > 1e-6 else ((tf_us - tav0 * 1e6) * 1.0e-6)
    if tf_reaction_s <= 0: tf_reaction_s = 1e-12
    current_eiee_ev, current_piee = list(eiee_dist_energies), list(eiee_dist_probs)    
    pfrag_abs: List[List[WP]] = [[0.0, 0.0] for _ in range(npairs + 1)]
    precursor_fname = "isomer.xyz" if (precursor_dir_path / "isomer.xyz").exists() else "fragment.xyz"
    precursor_full_path = precursor_dir_path / precursor_fname
    if not precursor_full_path.exists(): return False
    num_atoms_precursor, _, _ = utility.get_atomic_numbers_and_coords_py(precursor_full_path)
    if num_atoms_precursor == 0 : return False
    is_linear_precursor = utility.is_linear_molecule_py(str(precursor_full_path))
    num_vib_dof_precursor = max(1, (3 * num_atoms_precursor - 5) if is_linear_precursor else (3 * num_atoms_precursor - 6))
    barriers_ev, delta_e_ev, irc_cm_minus_1, is_rearrangement_flags = [0.0]*npairs, [0.0]*npairs, [0.0]*npairs, [False]*npairs
    for i in range(npairs):
        pair_id_str, reaction_dir = fragdirs[i][0], precursor_dir_path / fragdirs[i][0] 
        delta_e_ev[i] = iomod.rdshort_real_py(reaction_dir / f"de_{env.tslevel}", default=1000.0)
        if env.nots: 
            barriers_ev[i], irc_cm_minus_1[i] = max(0.0, delta_e_ev[i]), 100.0 
        else:
            barriers_ev[i] = iomod.rdshort_real_py(reaction_dir / f"barrier_{env.tslevel}", default=1000.0)
            irc_cm_minus_1[i] = iomod.rdshort_real_py(reaction_dir / f"ircmode_{env.geolevel}", default=-100.0)
            if irc_cm_minus_1[i] <= 0 and not env.eyring: irc_cm_minus_1[i] = 100.0 
        is_rearrangement_flags[i] = not bool(fragdirs[i][2])
    is_h_dissociation_flags = [False] * npairs 
    scale_iee_for_h_diss_factor = env.scaleeinthdiss if any(is_h_dissociation_flags) else 1.0
    sum_reaction_energy_prior_ev, ker_average_ev_previous_steps, ea0_formation_ev, de0_formation_ev, nat0_ts_formation, nvib0_ts_formation = 0.0, 0.0, 0.0, 0.0, 0, 0
    if fragmentation_level > 1: 
        # Warning: Logic for reading these prior terms needs to ensure precursor_dir_path.parent correctly points to the reaction directory that formed this precursor.
        # This is a known simplification point if product directories are flattened.
        print("  INFO: Reading prior energy terms for multi-level fragmentation. Pathing is critical here.")
        formation_reaction_dir = precursor_dir_path.parent # Assumed to be the reaction dir that formed current precursor
        sum_reaction_energy_prior_ev = iomod.rdshort_real_py(formation_reaction_dir / f"sumreac_{env.tslevel}", default=0.0)
        ea0_formation_ev = iomod.rdshort_real_py(formation_reaction_dir / f"barrier_{env.tslevel}", default=0.0)
        de0_formation_ev = iomod.rdshort_real_py(formation_reaction_dir / f"de_{env.tslevel}", default=0.0)
        ts_geom_for_parent = formation_reaction_dir / "ts" / "ts.xyz"
        if ts_geom_for_parent.exists():
            nat0_ts_formation, _, _ = utility.get_atomic_numbers_and_coords_py(ts_geom_for_parent)
            if nat0_ts_formation > 0: nvib0_ts_formation = max(1, 3*nat0_ts_formation - (5 if utility.is_linear_molecule_py(str(ts_geom_for_parent)) else 6))
        ker_average_ev_previous_steps = iomod.rdshort_real_py(precursor_dir_path / "keravold", default=0.0) # keravold IS in the current precursor dir
    current_eiee_ev = [max(0.0, e - sum_reaction_energy_prior_ev) for e in current_eiee_ev]
    alldgs_kjmol, all_nvib_ts_channel, all_delta_g_reaction_ev = [[0.0]*env.nthermosteps for _ in range(npairs)], [num_vib_dof_precursor]*npairs, [[0.0]*env.nthermosteps for _ in range(npairs)]
    if env.eyring or env.eyzpve:
        num_e_increments = env.nthermosteps if env.nthermosteps > 0 else 1
        max_eiee_for_thermo_ev = max(current_eiee_ev) if current_eiee_ev else env.eimp0 
        if env.mode.lower() == "cid": max_eiee_for_thermo_ev += env.cid_elab
        rrho_start_hartree_list = _xtb_thermo_py(env, num_e_increments, num_vib_dof_precursor, max_eiee_for_thermo_ev, precursor_full_path, True)
        for i in range(npairs):
            pair_id_str, item1_dir_rel, item2_dir_rel = fragdirs[i]
            reaction_dir, ts_xyz_path = precursor_dir_path / pair_id_str, precursor_dir_path / pair_id_str / "ts" / "ts.xyz"
            nat_ts, _, _ = utility.get_atomic_numbers_and_coords_py(ts_xyz_path); nat_ts = nat_ts or num_atoms_precursor
            all_nvib_ts_channel[i] = max(1, 3*nat_ts - (5 if utility.is_linear_molecule_py(str(ts_xyz_path)) else 6))
            rrho_ts_hartree_list = _xtb_thermo_py(env, num_e_increments, all_nvib_ts_channel[i], max_eiee_for_thermo_ev, ts_xyz_path, True, 0.0)
            prod1_xyz_abs_path = Path(env.startdir)/item1_dir_rel/("isomer.xyz" if not item2_dir_rel else "fragment.xyz")
            nat_f1, _, _ = utility.get_atomic_numbers_and_coords_py(prod1_xyz_abs_path); nat_f1 = nat_f1 or (num_atoms_precursor//2 if item2_dir_rel else num_atoms_precursor)
            nvib_f1 = max(1, 3*nat_f1 - (5 if utility.is_linear_molecule_py(str(prod1_xyz_abs_path)) else 6))
            rrho_f1_hartree_list = _xtb_thermo_py(env, num_e_increments, nvib_f1, max_eiee_for_thermo_ev, prod1_xyz_abs_path, True)
            rrho_f2_hartree_list = None
            if item2_dir_rel:
                prod2_xyz_abs_path = Path(env.startdir)/item2_dir_rel/"fragment.xyz"
                nat_f2, _, _ = utility.get_atomic_numbers_and_coords_py(prod2_xyz_abs_path); nat_f2 = nat_f2 or (num_atoms_precursor - nat_f1 if 0 < nat_f1 < num_atoms_precursor else num_atoms_precursor//2)
                nvib_f2 = max(1, 3*nat_f2 - (5 if utility.is_linear_molecule_py(str(prod2_xyz_abs_path)) else 6))
                rrho_f2_hartree_list = _xtb_thermo_py(env, num_e_increments, nvib_f2, max_eiee_for_thermo_ev, prod2_xyz_abs_path, True)
            for k_inc in range(num_e_increments):
                alldgs_kjmol[i][k_inc] = (barriers_ev[i] + (rrho_ts_hartree_list[k_inc] - rrho_start_hartree_list[k_inc]) * HARTREE_TO_EV) * EVTOKCAL * 4.184
                all_delta_g_reaction_ev[i][k_inc] = delta_e_ev[i] + ((rrho_f1_hartree_list[k_inc] + (rrho_f2_hartree_list[k_inc] if rrho_f2_hartree_list else 0.0)) - rrho_start_hartree_list[k_inc]) * HARTREE_TO_EV
    sumdekin_val, sumdeint_val, xtrav_val = 0.0, 0.0, 0.0 # For CID outputs
    if env.mode.lower() == "cid": print("    CID activation call placeholder.") # Actual call to cid.simulate_cid_activation_py would update these.
    sum_ktot_reciprocal, ker_total_weighted_sum, num_fragmentation_events = 0.0, 0.0, 0
    k_uni, p_react_norm = [0.0]*npairs, [0.0]*npairs
    for sample_idx in range(env.nsamples):
        if not current_eiee_ev : continue
        internal_e_ev = random.choices(current_eiee_ev, weights=current_piee, k=1)[0]; internal_e_ev = max(0.0, internal_e_ev)
        ktot_s_minus_1 = 0.0
        num_e_inc = num_e_increments if (env.eyring or env.eyzpve) else 1
        e_bin_idx = utility.get_index_from_energy_py(max_eiee_for_thermo_ev if (env.eyring or env.eyzpve) else internal_e_ev, internal_e_ev, num_e_inc)
        for j in range(npairs):
            e_rate_calc = internal_e_ev * scale_iee_for_h_diss_factor if is_h_dissociation_flags[j] else internal_e_ev
            barrier_j = alldgs_kjmol[j][e_bin_idx-1]/(EVTOKCAL*4.184) if (env.eyring or env.eyzpve) else barriers_ev[j]
            if e_rate_calc < barrier_j and barrier_j > 0: k_uni[j] = 0.0
            elif env.eyring or env.eyzpve: k_uni[j] = utility.calc_eyring_py(e_rate_calc, barrier_j, num_vib_dof_precursor)
            else: k_uni[j] = utility.calc_rrkm_py(e_rate_calc, barrier_j, num_vib_dof_precursor, irc_cm_minus_1[j])
            ktot_s_minus_1 += k_uni[j]
        if ktot_s_minus_1 > 1e-10: sum_ktot_reciprocal += (1.0/ktot_s_minus_1); p_react_norm = [k/ktot_s_minus_1 for k in k_uni]
        else: sum_ktot_reciprocal += (1.0/1e-10); p_react_norm = [0.0]*npairs
        rand_f, cum_prob, reacted_ch_idx = random.random(), 0.0, -1
        for j in range(npairs):
            cum_prob += p_react_norm[j]
            if rand_f <= cum_prob: reacted_ch_idx = j; break
        current_sample_prob = current_piee[sample_idx] if sample_idx < len(current_piee) else (sum(current_piee)/len(current_piee) if current_piee else 0) # Handle potential mismatch
        if reacted_ch_idx != -1: 
            pfrag_abs[reacted_ch_idx][0] += current_sample_prob; num_fragmentation_events +=1 
            ea_eff = alldgs_kjmol[reacted_ch_idx][e_bin_idx-1]/(EVTOKCAL*4.184) if (env.eyring or env.eyzpve) else barriers_ev[reacted_ch_idx]
            de_eff = all_delta_g_reaction_ev[reacted_ch_idx][e_bin_idx-1] if (env.eyring or env.eyzpve) else delta_e_ev[reacted_ch_idx]
            nvib_ts_ch = all_nvib_ts_channel[reacted_ch_idx]
            ker_eex = (internal_e_ev - ea_eff)/(0.44*nvib_ts_ch) if internal_e_ev > ea_eff and nvib_ts_ch > 0 else 0.0
            ker_ear = 0.33*(ea_eff - de_eff) if ea_eff > de_eff else 0.0
            ker_sample = (max(0.0, ker_eex) + max(0.0, ker_ear)) * env.scaleker
            ker_total_weighted_sum += ker_sample * current_sample_prob
        else: pfrag_abs[npairs][0] += current_sample_prob
    total_piee_s = sum(current_piee); pnorm = total_piee_s if total_piee_s > 1e-9 else 1.0
    for j in range(npairs + 1): pfrag_abs[j][0] /= pnorm
    tav_s = sum_ktot_reciprocal / env.nsamples if env.nsamples > 0 else float('inf')
    ker_avg_ev = ker_total_weighted_sum / pnorm if pnorm > 0 and num_fragmentation_events > 0 else 0.0
    for i in range(npairs):
        pair_id_str, item1_dir_rel, item2_dir_rel = fragdirs[i]
        path_f1_abs = Path(env.startdir) / item1_dir_rel; path_f1_abs.mkdir(exist_ok=True)
        iomod.wrshort_real_py(path_f1_abs / "pfrag", pfrag_abs[i][0] * p0_intensity_fraction, '')
        if item2_dir_rel: 
            path_f2_abs = Path(env.startdir) / item2_dir_rel; path_f2_abs.mkdir(exist_ok=True)
            iomod.wrshort_real_py(path_f2_abs / "pfrag", pfrag_abs[i][0] * p0_intensity_fraction, '')
            f1_xyz = path_f1_abs/("isomer.xyz" if not item2_dir_rel else "fragment.xyz"); f2_xyz = path_f2_abs/"fragment.xyz"
            qmass = 1.0
            if f1_xyz.exists() and f2_xyz.exists():
                ats1, ats2 = utility.get_atomic_numbers_from_xyz_py(str(f1_xyz)), utility.get_atomic_numbers_from_xyz_py(str(f2_xyz))
                if ats1 and ats2: m1,m2 = utility.get_average_mol_mass_py(ats1), utility.get_average_mol_mass_py(ats2); qmass = m2/m1 if m1 > 1e-6 else 1.0
            iomod.wrshort_real_py(precursor_dir_path / pair_id_str / "qmass", qmass, '')
            iomod.wrshort_real_py(precursor_dir_path / pair_id_str / "kerav", ker_avg_ev if pfrag_abs[i][0] > 1e-9 else 0.0, '')
        else: iomod.wrshort_real_py(precursor_dir_path / pair_id_str / "kerav", 0.0, '') 
    iomod.wrshort_real_py(precursor_dir_path / "tav", tav_s, ''); iomod.wrshort_real_py(precursor_dir_path / "pfs", pfrag_abs[npairs][0] * p0_intensity_fraction, '')
    iomod.wrshort_real_py(precursor_dir_path / "keravold", ker_average_ev_previous_steps + ker_avg_ev, '')
    if env.mode.lower() == "cid" and env.cid_elab > 0.0:
        iomod.wrshort_real_py(precursor_dir_path / "sumdekin", sumdekin_val, ''); iomod.wrshort_real_py(precursor_dir_path / "sumdeint", sumdeint_val, ''); iomod.wrshort_real_py(precursor_dir_path / "x_trav", xtrav_val, '')
    return True

def mcesim_py(env: RunTypeData, eiee_dist_energies: List[WP], eiee_dist_probs: List[WP]) -> None:
    print("\n--- Starting Multi-Generational Monte Carlo Simulation (mcesim_py) ---")
    original_main_run_dir = Path.cwd() 
    if not (hasattr(env, 'startdir') and env.startdir and env.startdir.is_dir()):
        print("Warning: env.startdir not properly set or accessible. Using current directory as base for mcesim.")
        env.startdir = Path.cwd()
    os.chdir(env.startdir)
    iomod.cleanup_qcxms2_for_restart_py(env)
    allfrags_file = env.startdir / "allfrags.dat"
    if not allfrags_file.exists():
        print(f"Error: {allfrags_file} not found."); os.chdir(original_main_run_dir); return
    total_frags_count = iomod.rdshort_int_py(allfrags_file, default=0)
    all_frag_paths: List[str] = []
    try:
        with open(allfrags_file, 'r') as f: lines = f.readlines()
        if len(lines) > 1: all_frag_paths = [line.strip() for line in lines[1:] if line.strip()]
    except Exception as e: print(f"Error reading {allfrags_file}: {e}"); os.chdir(original_main_run_dir); return
    if total_frags_count != len(all_frag_paths): print(f"Warning: Mismatch in frag count in {allfrags_file}.")
    restart_frags: List[str] = []
    for idx, frag_path_rel_str in enumerate(all_frag_paths):
        current_frag_dir_abs = (env.startdir / frag_path_rel_str).resolve()
        print(f"\n  Processing fragment {idx+1}/{len(all_frag_paths)}: {current_frag_dir_abs}")
        if not current_frag_dir_abs.is_dir(): print(f"    Warning: Dir not found. Skipping."); continue
        os.chdir(current_frag_dir_abs)
        frag_level = current_frag_dir_abs.name.count('p') # Simplified level calculation
        frags_file = Path("fragments")
        if not frags_file.exists():
            pfrag_val = iomod.rdshort_real_py(Path("pfrag"), default=0.0) # Intensity here
            if pfrag_val >= env.pthr : restart_frags.append(str(current_frag_dir_abs))
            else: print(f"    'fragments' file not found and pfrag ({pfrag_val:.4f}) < pthr ({env.pthr}). Skipping.")
            os.chdir(env.startdir); continue
        npairs_curr_frag = iomod.rdshort_int_py(Path("npairs2"), default=0)
        if npairs_curr_frag == 0: print(f"    No reaction pairs. Skipping MC."); os.chdir(env.startdir); continue
        fragdirs_mc: List[List[str]] = []
        try:
            with open(frags_file, 'r') as f_frgs:
                for line in f_frgs:
                    parts = line.strip().split()
                    if len(parts) == 2: fragdirs_mc.append([parts[0], parts[1], ""])
                    elif len(parts) == 3: fragdirs_mc.append(list(parts))
        except Exception as e: print(f"    Error reading {frags_file}: {e}"); os.chdir(env.startdir); continue
        if not fragdirs_mc and npairs_curr_frag > 0: print(f"    Warning: npairs2 > 0 but no entries in {frags_file}."); os.chdir(env.startdir); continue
        p0_intensity_curr_frag = iomod.rdshort_real_py(Path("pfrag"), default=0.0)
        print(f"    Fragment level: {frag_level}, Max level: {env.nfrag}, Intensity: {p0_intensity_curr_frag:.4f}, Threshold: {env.pthr}")
        if npairs_curr_frag > 0 and p0_intensity_curr_frag >= env.pthr and frag_level <= env.nfrag + 1 :
            print(f"    Proceeding with Monte Carlo for {current_frag_dir_abs.name}")
            # Pathing for prior energy terms in montecarlo_py (reading from parent reaction dir) is critical.
            # The CWD for montecarlo_py is current_frag_dir_abs. Its parent is current_frag_dir_abs.parent.
            # This parent must be the reaction directory that formed the current fragment.
            montecarlo_py(env, npairs_curr_frag, fragdirs_mc, eiee_dist_energies, eiee_dist_probs, frag_level, p0_intensity_curr_frag)
        else: print("    Skipping MC based on conditions (npairs, intensity, or level).")
        os.chdir(env.startdir)
    if restart_frags: print("\nWarning: Fragments needing restart/check:"); [print(f"  {p}") for p in restart_frags]
    if env.spectrum == 1: print("\nWarning: allspec.dat generation (env.spectrum=1) not implemented.")
    os.chdir(original_main_run_dir)
    print("\n--- Multi-Generational Monte Carlo Simulation (mcesim_py) Finished ---")

if __name__ == '__main__':
    print("Testing mcsimu.py (montecarlo_py structure and mcesim_py placeholder)")
    env_test = RunTypeData()
    env_test.geolevel = "gfn2"; env_test.tslevel = "gfn2"; env_test.sthr = 100; env_test.ithr = 50
    env_test.temp = 298.15; env_test.nsamples = 10; env_test.eyring = True; env_test.nthermosteps = 2
    env_test.mode = "ei"; env_test.infile = "dummy_mol.xyz"; env_test.chrg = 0
    env_test.startdir = Path("mcsimu_test_run_main").resolve()
    env_test.pthr = 0.01 # Probability threshold
    env_test.nfrag = 1 # Max fragmentation generations for mcesim

    precursor_rel_path = Path("gen0_precursor")
    current_precursor_dir = env_test.startdir / precursor_rel_path
    current_precursor_dir.mkdir(parents=True, exist_ok=True)
    precursor_xyz_content = "3\nPrecursor\nO 0 0 0\nH 0 0 1\nH 0 1 0\n"
    (current_precursor_dir / "fragment.xyz").write_text(precursor_xyz_content)
    (env_test.startdir / env_test.infile).write_text(precursor_xyz_content) # For IEE get_molecular_specs
    (current_precursor_dir / ".CHRG").write_text(f"{env_test.chrg}\n") # For xtb_thermo via _get_chrg_uhf
    (current_precursor_dir / ".UHF").write_text("0\n")

    pair1_react_subdir = current_precursor_dir / "p0"
    pair1_react_subdir.mkdir(exist_ok=True)
    (pair1_react_subdir / f"de_{env_test.tslevel}").write_text("0.5\n")
    (pair1_react_subdir / f"barrier_{env_test.tslevel}").write_text("1.0\n")
    (pair1_react_subdir / f"ircmode_{env_test.geolevel}").write_text("-200.0\n")
    (pair1_react_subdir / "ts").mkdir(exist_ok=True)
    (pair1_react_subdir / "ts" / "ts.xyz").write_text(precursor_xyz_content)
    (pair1_react_subdir / "ts" / ".CHRG").write_text(f"{env_test.chrg}\n")
    (pair1_react_subdir / "ts" / ".UHF").write_text("0\n")

    prod1_dir_rel_start = Path("p0f1_product") # Product of gen0_precursor via p0 reaction
    prod1_dir_abs = env_test.startdir / prod1_dir_rel_start
    prod1_dir_abs.mkdir(exist_ok=True)
    (prod1_dir_abs / "fragment.xyz").write_text("1\nH\nH 0 0 0\n")
    (prod1_dir_abs / ".CHRG").write_text(f"{env_test.chrg}\n")
    (prod1_dir_abs / ".UHF").write_text("0\n")
    # This product (p0f1_product) would have its own "pfrag" file from the montecarlo_py run on gen0_precursor
    (prod1_dir_abs / "pfrag").write_text(f"{0.8*1.0}\n") # Assume 80% yield from precursor of intensity 1.0
    # It might also have an npairs2 and fragments file if it can fragment further
    (prod1_dir_abs / "npairs2").write_text("0\n") # Say this one doesn't fragment further for the test
    (prod1_dir_abs / "fragments").write_text("") # Empty fragments file


    fragdirs_test_mc = [["p0", str(prod1_dir_rel_start), ""]]
    dummy_p0_intensity_for_mc = 1.0 
    eiee_test, piee_test = [1.0, 2.0, 3.0], [0.3, 0.4, 0.3]
    
    # Mocking utilities
    utility.get_atomic_numbers_and_coords_py = lambda fname: (3, [8,1,1], [[0,0,0],[0,0,1],[0,1,0]]) if "fragment.xyz" in str(fname) or "mol.xyz" in str(fname) or "ts.xyz" in str(fname) else ((1, [1], [[0,0,0]]) if "p0f1_product" in str(fname) else (0,[],[]))
    utility.get_atomic_numbers_from_xyz_py = lambda fname: [8,1,1] if "fragment.xyz" in str(fname) else ([1] if "p0f1_product" in str(fname) else [])
    utility.is_linear_molecule_py = lambda fname: False
    utility.get_average_mol_mass_py = lambda ats: sum(ats) # Simplified mass
    utility.calctemp_py = lambda nvib, energy_ev: 298.15 + energy_ev * 50
    utility.get_index_from_energy_py = lambda max_e, e, ninc: min(max(1, int(e/ (max_e/ninc if max_e > 0 and ninc > 0 else 1.0) ) +1 if max_e > 0 and ninc > 0 else 1), ninc if ninc > 0 else 1)
    utility.calc_eyring_py = lambda E, Ea, Nvib: 1e12 * math.exp(-Ea*EVTOKCAL*1000 / (8.314 * (298.15 + E*10))) if (E > Ea and (298.15+E*10)>0) else 0.0
    utility.calc_rrkm_py = lambda E, Ea, Nvib, nu: 1e12 * ( (E-Ea)/E )**(Nvib-1) if E > Ea else 0.0
    qmmod._get_chrg_uhf = lambda env, fname, chrg_in=None, uhf_in=None: (env.chrg, 0)
    iomod.rdshort_real_py = iomod.rdshort_real_py # Use actual ones if available in iomod.py
    iomod.wrshort_real_py = iomod.wrshort_real_py
    iomod.rdshort_int_py = iomod.rdshort_int_py
    iomod.cleanup_qcxms2_for_restart_py = lambda env_ignored: None


    # Test montecarlo_py
    print(f"--- Testing montecarlo_py directly in {current_precursor_dir.resolve()} ---")
    original_cwd = Path.cwd()
    os.chdir(current_precursor_dir) 
    montecarlo_py(env_test, 1, fragdirs_test_mc, eiee_test, piee_test, 0, dummy_p0_intensity_for_mc)
    os.chdir(original_cwd)
    print(f"--- Finished testing montecarlo_py ---")

    # Test mcesim_py
    # Create allfrags.dat for mcesim test
    (env_test.startdir / "allfrags.dat").write_text(f"2\n{str(precursor_rel_path)}\n{str(prod1_dir_rel_start)}\n")
    # Precursor (gen0_precursor) needs npairs2 and fragments to be processed by mcesim
    (current_precursor_dir / "npairs2").write_text("1\n")
    with open(current_precursor_dir / "fragments", "w") as f:
        f.write(f"{fragdirs_test_mc[0][0]} {fragdirs_test_mc[0][1]} {fragdirs_test_mc[0][2]}\n")
    (current_precursor_dir / "pfrag").write_text(f"{dummy_p0_intensity_for_mc}\n") # Intensity of this precursor

    print(f"\n--- Testing mcesim_py in {env_test.startdir.resolve()} ---")
    # mcesim_py itself changes CWD to env.startdir
    mcesim_py(env_test, eiee_test, piee_test)
    print(f"--- Finished testing mcesim_py ---")

    shutil.rmtree(env_test.startdir)
    print("\nMCsimu test structure execution finished.")

```
