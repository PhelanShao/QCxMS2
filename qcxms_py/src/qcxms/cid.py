import math
import random
from pathlib import Path
from typing import Tuple, Optional, List, Dict

try:
    from .data import RunTypeData, WP
    from . import iomod
    from . import utility
    from .boxmuller import vary_energies_py
    from .constants import PI, AMU_TO_KG, AUTOAA, EV_TO_JOULE, KB_J_K # Bohr to Angstrom for AUTOAA
except ImportError:
    # Fallbacks for standalone/testing
    print("Attempting to import dummy/mock modules for cid.py standalone run.")
    from data import RunTypeData, WP # type: ignore
    import iomod_mock as iomod # type: ignore
    import utility_mock as utility # type: ignore
    from boxmuller import vary_energies_py # type: ignore
    PI = math.pi
    AMU_TO_KG = 1.66053906660e-27
    AUTOAA = 0.529177210903 # Bohr to Angstrom
    EV_TO_JOULE = 1.602176634e-19
    KB_J_K = 1.38064852E-23


# Atomic Radii in Angstroms, used for collision cross-section calculation.
# Indexed by atomic number Z. Rad[0] is dummy.
# Source: Hardcoded in QCxMS cid.f90, likely from a standard set of covalent or vdW radii.
# These are converted to Bohr by multiplying by AUTOAA in the Fortran code where Rad is defined,
# but the values themselves seem to be Angstrom-like.
# The Fortran code has `Rad(118) = aatoau * [ ... angstrom values ... ]`
# This means Rad(Z) in Fortran is already in Bohr.
# For Python, let's store them in Angstrom and convert when used if needed, or store directly in Bohr.
# Storing in Bohr as per Fortran `Rad` array's effective units after multiplication by aatoau.
# However, the comment says "Radius used in QCxMS (in au)".
# Let's assume the Fortran `Rad` array values are already effectively in Bohr.
# The values are: 0.32_wp,0.37_wp for H, He then 1.30_wp * 0.529... for Li etc.
# No, the Fortran line `real(wp), parameter :: Rad(118) = aatoau * [ ... ]` implies the list
# in square brackets are Angstrom values that are then converted to Bohr.
# So, let's store the Angstrom values and convert to Bohr where needed, or store Bohr directly.
# For simplicity, storing Bohr directly by pre-applying AUTOAA (if these are indeed Angstroms).
# The Fortran code snippet is: Rad(118) = aatoau * [0.32_wp, 0.37_wp, ...].
# This means Rad(1) = aatoau * 0.32, Rad(2) = aatoau * 0.37.
# So the list should be the Angstrom values.

_ATOMIC_RADII_ANGSTROM: List[WP] = [
    0.0,  # Dummy for Z=0
    0.32, 0.37,  # H, He
    1.30, 0.99, 0.84, 0.75, 0.71, 0.64, 0.60, 0.62,  # Li-Ne
    1.60, 1.40, 1.24, 1.14, 1.09, 1.04, 1.00, 1.01,  # Na-Ar
    2.00, 1.74,  # K, Ca
    1.59, 1.48, 1.44, 1.30, 1.29, 1.24, 1.18, 1.17, 1.22, 1.20,  # Sc-Zn
    1.23, 1.20, 1.20, 1.18, 1.17, 1.16,  # Ga-Kr
    2.15, 1.90,  # Rb, Sr
    1.76, 1.64, 1.56, 1.46, 1.38, 1.36, 1.34, 1.30, 1.36, 1.40,  # Y-Cd
    1.42, 1.40, 1.40, 1.37, 1.36, 1.36,  # In-Xe
    2.38, 2.06,  # Cs, Ba
    1.94, 1.84, 1.90, 1.88, 1.86, 1.85, 1.83, 1.82, 1.81, 1.80, 1.79, 1.77, 1.77, 1.78,  # La-Lu (Lu is 71)
    1.74, 1.64, 1.58, 1.50, 1.41, 1.36, 1.32, 1.30, 1.30, 1.32,  # Hf-Hg (Hg is 80)
    1.44, 1.45, 1.50, 1.42, 1.48, 1.46,  # Tl-Rn (Rn is 86)
    # Fortran data continues with more values (likely estimated/less common radii for heavier elements)
    # Copied from Fortran `Rad` array values, assuming these are Angstrom for direct use.
    # The Fortran code had `Rad(118) = aatoau * [ ... angstrom values ... ]`
    # So, if we want Rad(iat(i)) to be in Bohr directly, we should pre-multiply.
    # For clarity, let's keep them as Angstroms here and convert where used.
    2.42, 2.11,  # Fr, Ra (87, 88)
    2.01, 1.90, 1.84, 1.83, 1.80, 1.80, # Ac-Pu (89-94)
    1.49, # Am (95) - Note: Fortran data seems to have values for heavier elements.
    1.49, 1.51, 1.51, 1.48, 1.50, 1.56, 1.58, # Cm-No (96-102)
    1.45, 1.41, 1.34, 1.29, 1.27, # Lr- (103-107)
    1.21, 1.16, 1.15, 1.09, 1.22, # - Cn (108-112)
    1.36, 1.43, 1.46, 1.58, 1.48, 1.57 # Nh-Og (113-118)
]

def get_atomic_radius_bohr(atomic_number: int) -> WP:
    """Returns atomic radius in Bohr."""
    if 0 < atomic_number < len(_ATOMIC_RADII_ANGSTROM):
        return _ATOMIC_RADII_ANGSTROM[atomic_number] / AUTOAA # Convert Angstrom to Bohr
    print(f"Warning: Atomic radius for Z={atomic_number} not found. Using default 1.5 Bohr.")
    return 1.5 # Default radius in Bohr if not found


def _collision_setup_py(
    env: RunTypeData, 
    num_atoms: int, 
    atomic_numbers: List[int], 
    coordinates_angstrom: List[List[WP]] # Assuming Nx3 list of lists for coords in Angstrom
) -> Tuple[WP, WP]:
    """
    Calculates collision parameters: mean free path (mfpath) and 
    average number of collisions (calc_coll) for an ion.
    
    Args:
        env: RunTypeData object.
        num_atoms: Number of atoms in the ion.
        atomic_numbers: List of atomic numbers.
        coordinates_angstrom: List of [x,y,z] coordinates in Angstrom.
        
    Returns:
        Tuple (mean_free_path_meters, avg_num_collisions).
    """
    mol_radius_total_bohr: WP = 0.0
    
    # Convert coordinates to Bohr for consistency if needed, or ensure Rad is also Angstrom
    # The Fortran code calculates center_of_geometry on xyz (presumably Bohr from xtb)
    # then uses Rad(iat(i)) which is already in Bohr.
    # Let's assume coordinates_angstrom are indeed Angstrom and convert them to Bohr.
    coordinates_bohr: List[List[WP]] = [[c / AUTOAA for c in atom_coord] for atom_coord in coordinates_angstrom]

    center_of_geom_bohr = utility.calculate_center_of_geometry_py(num_atoms, coordinates_bohr)

    # Determine the molecular radius by checking the largest distance from CoG to atom's vdW surface
    max_extent_bohr: WP = 0.0
    for i in range(num_atoms):
        atom_rad_bohr = get_atomic_radius_bohr(atomic_numbers[i])
        dist_sq: WP = sum([(coordinates_bohr[i][j] - center_of_geom_bohr[j])**2 for j in range(3)])
        # distance from CoG to atom center + atom's own radius
        current_extent_bohr = math.sqrt(dist_sq) + atom_rad_bohr 
        if current_extent_bohr > max_extent_bohr:
            max_extent_bohr = current_extent_bohr
    
    mol_radius_total_bohr = max_extent_bohr

    # Radii in meters for cross-section calculation
    # env.cid_rgas is already in Bohr from Fortran context
    gas_radius_meters = (env.cid_rgas * AUTOAA) * 1e-10  # Bohr to Angstrom, then Angstrom to meters
    ion_radius_meters = (mol_radius_total_bohr * AUTOAA) * 1e-10 # Bohr to Angstrom, then Angstrom to meters

    # Calculate collisional cross section and mean free path
    collisional_cross_section_m2 = PI * ((ion_radius_meters + gas_radius_meters)**2)

    # Mean Free Path (mfpath) = (kB * T_gas) / (cross_section * P_gas)
    # kB_J_K is Boltzmann const in J/K
    # env.cid_Tgas in Kelvin
    # env.cid_Pgas in Pascal (N/m^2)
    if collisional_cross_section_m2 == 0 or env.cid_Pgas == 0:
        mean_free_path_meters = float('inf') # Avoid division by zero
    else:
        mean_free_path_meters = (KB_J_K * env.cid_Tgas) / (collisional_cross_section_m2 * env.cid_Pgas)
    
    print(f"  Ion radius: {ion_radius_meters*1e10:.2f} A, Gas radius: {gas_radius_meters*1e10:.2f} A")
    print(f"  Collisional cross section: {collisional_cross_section_m2:.2e} m^2")
    print(f"  Mean free path: {mean_free_path_meters:.2e} m")

    # Average number of collisions expected in the chamber length
    avg_num_collisions: WP
    if mean_free_path_meters == 0 or mean_free_path_meters == float('inf'): # Avoid issues
        avg_num_collisions = 0.0 if mean_free_path_meters == float('inf') else float('inf')
    else:
        avg_num_collisions = env.cid_lchamb / mean_free_path_meters
    
    print(f"  Calculated average number of collisions in chamber (L={env.cid_lchamb}m): {avg_num_collisions:.2f}")
    
    return mean_free_path_meters, avg_num_collisions


def _simulate_single_collision_py(
    env: RunTypeData, 
    current_eiee_dist: List[WP], # Will be modified
    ion_mass_amu: WP, 
    current_ion_ekin_ev: WP
) -> Tuple[WP, WP]:
    """
    Simulates a single collision, updating the internal energy distribution (eiee)
    and calculating kinetic energy loss.
    
    Modifies current_eiee_dist in place.
    Returns (delta_ekin_ev, average_delta_eint_ev).
    """
    eta_inelasticity = 0.43 # Hardcoded inelasticity factor from Fortran
    beta_mass_factor = env.cid_mgas / (env.cid_mgas + ion_mass_amu)

    # Center-of-mass collision energy
    e_com_ev = beta_mass_factor * current_ion_ekin_ev

    # Average internal energy deposited per collision
    avg_delta_eint_ev = eta_inelasticity * e_com_ev
    
    # Vary deposited internal energy for each sample in the distribution
    # and add to existing internal energy.
    # distribution_width for varying dEint is env.cid_collw
    for i in range(len(current_eiee_dist)):
        varied_delta_eint = vary_energies_py(avg_delta_eint_ev, env.cid_collw)
        if varied_delta_eint < 0: varied_delta_eint = 0.0 # Cannot deposit negative energy
        current_eiee_dist[i] += varied_delta_eint
        # eiee should not be negative (though adding positive dEint won't make it so if it wasn't already)
        if current_eiee_dist[i] < 0: current_eiee_dist[i] = 0.0


    # Calculate kinetic energy change (loss) of the ion
    # gamma_factor = (2 * ion_mass_amu + eta_inelasticity * env.cid_mgas) / (ion_mass_amu + env.cid_mgas)
    # The Fortran formula was: gamma = (2*mp + (eta*env%cid_mgas))/(mp + env%cid_mgas)
    # This gamma seems to be > 1. dEkin = -gamma * Ecom means Ekin loss is > Ecom.
    # This specific formula might be from a particular collision model.
    # Let's recheck common CID energy transfer models. A simpler model:
    # dEkin (loss for ion) = - dEint (gain for ion) - dE_gas (energy gained by collision gas)
    # If we assume dEint is what's converted from Ecom to internal:
    # The change in ion's kinetic energy in COM frame is -dEint.
    # Then transform back to lab frame.
    # However, to match Fortran:
    gamma_factor = (2.0 * ion_mass_amu + eta_inelasticity * env.cid_mgas) / (ion_mass_amu + env.cid_mgas)
    
    # Ensure beta * gamma is not >= 1, as per Fortran comment implies it's a check on physical possibility
    # Though the Fortran code only printed error if gamma*beta >=1, it didn't change behavior.
    # This check is not in the simcoll routine but implied by the overall energy transfer model.
    # The formulas are from J. Am. Soc. Mass Spectrom. 2013, 24, 7, 1064–1071
    
    delta_ekin_ev = -gamma_factor * e_com_ev # This is the change (loss)
    
    return delta_ekin_ev, avg_delta_eint_ev


def simulate_cid_activation_py(
    env: RunTypeData,
    tf_max_reaction_time: WP, # Max time for reaction (from mcsimu)
    max_iee_for_indexing: WP, # Max IEE for binning barriers (from mcsimu)
    eiee_list_in: List[WP],   # Input internal energy distribution
    piee_list_in: List[WP],   # Corresponding probabilities/weights
    num_samples: int,
    all_reaction_barriers_dgs: List[List[WP]], # Barriers [reaction_idx][energy_bin_idx]
    num_energy_increments: int, # nincr for barrier indexing
    current_frag_level: int,
    ker_average_from_previous_step: WP,
    delta_e0_formation: WP, # Reaction energy for formation of current ion
    ea0_formation: WP,      # Activation energy for formation of current ion
    nat0_ts_formation: int, # Num atoms in TS forming current ion
    nvib0_ts_formation: int,# Num vib modes in TS forming current ion
    num_possible_reactions: int,
    is_h_dissociation_flags: List[bool],
    scale_eint_for_h_diss_factor: WP,
    is_rearrangement_flags: List[bool],
    consider_only_fragmentations: bool # noisos flag from mcsimu
) -> Tuple[List[WP], List[WP]]:
    """
    Simulates Collision-Induced Dissociation by modifying the internal energy distribution.
    This is the Python equivalent of Fortran `simcid`.
    Returns the modified (eiee_list, piee_list). piee is often kept uniform.
    """
    
    print(f"CID Mode Activation (Level {current_frag_level}): Modifying internal energy distribution.")
    if env.cid_mode == 2: # Temprun mode, no collisions simulated
        print("  CID Temprun mode: No collisions performed. Using initial ESI distribution.")
        return list(eiee_list_in), list(piee_list_in) # Return copies

    if not env.eyring: # Current CID logic in Fortran depends on Eyring rates
        raise ValueError("CID mode currently only works with Eyring rate calculations enabled.")

    # Get properties of the current ion
    # Current precursor's geometry file (must be present in CWD)
    current_xyz_fname = "isomer.xyz" if Path("isomer.xyz").exists() else "fragment.xyz"
    if not Path(current_xyz_fname).exists():
        raise FileNotFoundError(f"Precursor XYZ file {current_xyz_fname} not found for CID simulation.")
    
    # nat, atomic_numbers, coords_angstrom = utility.read_full_xyz_py(current_xyz_fname)
    # Mocking this read for now as utility.read_full_xyz_py is not defined yet.
    num_atoms_current_ion = iomod.rdshort_int(current_xyz_fname)
    atomic_numbers_current_ion = [1]*num_atoms_current_ion # Dummy atom types
    coords_angstrom_current_ion = [[0.0,0.0,0.0]]*num_atoms_current_ion # Dummy coords
    if num_atoms_current_ion == 0:
        raise ValueError(f"Could not read nat from {current_xyz_fname} for CID.")


    ion_mass_amu = utility.get_average_mol_mass_py(atomic_numbers_current_ion)
    if ion_mass_amu <= 0:
        raise ValueError("Could not determine mass of ion for CID.")
    
    num_vib_modes_current_ion = 3 * num_atoms_current_ion - 6
    if num_vib_modes_current_ion <= 0: num_vib_modes_current_ion = 1


    # Initialize kinetic energy and travel parameters
    ekin_current_ev: WP
    sum_delta_eint_ev: WP = 0.0
    sum_delta_ekin_ev: WP = 0.0
    distance_traveled_m: WP = 0.0
    time_traveled_s: WP = 0.0

    # Read parameters from previous fragmentation step if not the first step
    # `noisos` from mcsimu means the current precursor is a direct fragmentation product, not an isomer of initial mol.
    # Fortran: if (nfragl == 1 .and. .not. noisos)
    is_first_fragmentation_step_not_isomer = (current_frag_level == 1 and not consider_only_fragmentations)

    if is_first_fragmentation_step_not_isomer:
        ekin_current_ev = env.cid_elab
    else:
        sum_delta_ekin_ev = iomod.rdshort_real("sumdekin", default=0.0)
        sum_delta_eint_ev = iomod.rdshort_real("sumdeint", default=0.0)
        distance_traveled_m = iomod.rdshort_real("x_trav", default=0.0)
        time_traveled_s = iomod.rdshort_real("t_trav", default=0.0)
        mass_ratio_q = iomod.rdshort_real("qmass", default=0.0) # m_other_frag / m_current_ion
        
        # Ekin = (Elab_initial + Sum_dEkin_precursor + KER_formation_current_ion) * partitioning_factor
        partitioning_factor = 1.0 / (1.0 + mass_ratio_q) if (1.0 + mass_ratio_q) != 0 else 1.0
        ekin_current_ev = (env.cid_elab + sum_delta_ekin_ev + ker_average_from_previous_step) * partitioning_factor
        if ekin_current_ev < 0: ekin_current_ev = 0.0
    
    print(f"  Initial E_kin for current ion: {ekin_current_ev:.2f} eV")

    # Collision setup
    mean_free_path_m, avg_total_collisions = _collision_setup_py(
        env, num_atoms_current_ion, atomic_numbers_current_ion, coords_angstrom_current_ion
    )
    
    num_collisions_performed = 0
    # Use average total collisions as max, but loop can break earlier
    max_collisions_to_simulate = int(round(avg_total_collisions)) 
    if env.cid_maxcoll > 0 : # User can override maxcoll
        max_collisions_to_simulate = min(max_collisions_to_simulate, env.cid_maxcoll)


    # Internal energy distribution for this step (start with input, potentially modified by prior sum_dEint)
    current_eiee_dist = [e + (sum_delta_eint_ev if not is_first_fragmentation_step_not_isomer else 0.0) for e in eiee_list_in]

    # Collision Loop
    for coll_idx in range(max_collisions_to_simulate):
        if distance_traveled_m >= env.cid_lchamb:
            print(f"  Ion exited collision chamber after {num_collisions_performed} collisions ({distance_traveled_m:.2e} m).")
            break

        # Average internal energy of the population
        avg_internal_e = sum(e * p for e, p in zip(current_eiee_dist, piee_list_in)) / sum(piee_list_in) if sum(piee_list_in) > 0 else 0.0

        # KER calculation for formation of current ion (if applicable)
        # This KER is released, some part might go to ion's Ekin
        # Fortran logic reduces Eavg by KER before rate calc, which is unusual.
        # For now, replicating, but this needs review.
        ker_formation_ev = 0.0
        if current_frag_level > 1 and nat0_ts_formation > 0 and env.calcKER: # Formed from a previous reaction
            # This KER is from the reaction that FORMED the current ion
            ker_eex = (avg_internal_e - ea0_formation) / (0.44 * nvib0_ts_formation) if nvib0_ts_formation > 0 else 0.0
            ker_ear = 0.33 * (ea0_formation - delta_e0_formation)
            if ker_eex < 0: ker_eex = 0.0
            if ker_ear < 0: ker_ear = 0.0
            ker_formation_ev = ker_eex + ker_ear + ker_average_from_previous_step # Total KER available from formation
            if env.scaleker != 1.0: ker_formation_ev *= env.scaleker
        
        avg_internal_e_for_rate = avg_internal_e - ker_formation_ev # Fortran logic
        
        # Find fastest unimolecular reaction rate k_max and its barrier
        energy_bin_idx = utility.get_index_from_energy_py(max_iee_for_indexing, avg_internal_e_for_rate, num_energy_increments)
        
        k_max = 0.0
        fastest_reaction_idx = -1
        min_barrier_for_fastest_eV = float('inf')

        for r_idx in range(num_possible_reactions):
            # Apply filters from mcsimu (noisos, isrearr)
            if consider_only_fragmentations and is_rearrangement_flags[r_idx]: continue
            if not consider_only_fragmentations and not is_rearrangement_flags[r_idx]: continue # If only isomers and this is fragmentation

            e_for_rate = avg_internal_e_for_rate * scale_eint_for_h_diss_factor if is_h_dissociation_flags[r_idx] else avg_internal_e_for_rate
            current_barrier_ev = all_reaction_barriers_dgs[r_idx][energy_bin_idx -1] # -1 for 0-indexed Python
            
            # Assuming utility.calc_eyring_py takes energy, barrier (both eV), and nvib
            k_val = utility.calc_eyring_py(e_for_rate, current_barrier_ev, num_vib_modes_current_ion)
            if k_val > k_max:
                k_max = k_val
                fastest_reaction_idx = r_idx
                min_barrier_for_fastest_eV = current_barrier_ev
        
        # Time for unimolecular reaction (half-life) vs. time to next collision
        t_unimolecular_min_s = math.log(2.0) / k_max if k_max > 0 else float('inf')
        
        time_to_next_collision_s: WP
        if ekin_current_ev <= 0 or mean_free_path_m == float('inf'): # Ion not moving or infinite mfp
            time_to_next_collision_s = float('inf')
        else:
            ion_velocity_m_s = math.sqrt(2.0 * ekin_current_ev * EV_TO_JOULE / (ion_mass_amu * AMU_TO_KG))
            if ion_velocity_m_s == 0: time_to_next_collision_s = float('inf')
            else: time_to_next_collision_s = mean_free_path_m / ion_velocity_m_s
        
        print(f"  Collision {coll_idx+1}: E_kin={ekin_current_ev:.2f} eV, E_int_avg={avg_internal_e:.2f} eV")
        print(f"    t_unimolecular={t_unimolecular_min_s:.2e} s (for barrier {min_barrier_for_fastest_eV:.2f} eV), t_collision={time_to_next_collision_s:.2e} s")

        if t_unimolecular_min_s > time_to_next_collision_s and time_to_next_collision_s != float('inf'):
            # Collision happens
            num_collisions_performed += 1
            delta_ekin_ev, avg_delta_eint_ev = _simulate_single_collision_py(
                env, current_eiee_dist, ion_mass_amu, ekin_current_ev
            )
            print(f"    Collision occurs: dE_int_avg={avg_delta_eint_ev:.2f} eV, dE_kin={delta_ekin_ev:.2f} eV")

            sum_delta_eint_ev += avg_delta_eint_ev # Track total average internal energy gained
            sum_delta_ekin_ev += delta_ekin_ev   # Track total kinetic energy change (usually negative)
            
            ekin_current_ev += delta_ekin_ev * env.cid_scool # Apply cooling factor to loss
            if ekin_current_ev < 0: ekin_current_ev = 0.0
            
            distance_traveled_m += mean_free_path_m 
            time_traveled_s += time_to_next_collision_s

            # Optionally log eiee distribution after this collision
            if env.printlevel >=3 : # Or a specific flag
                 coll_fname = f"ecoll_{'iso' if not consider_only_fragmentations else 'frag'}_{coll_idx+1}.dat"
                 iomod.printiee(num_samples, current_eiee_dist, piee_list_in, coll_fname)

        else: # Reaction is faster or no more collisions possible
            print(f"  Ion reacts/fragments before next collision (or no more collisions). Total collisions this step: {num_collisions_performed}")
            break 
    
    # Write out final state for next fragmentation step if any
    iomod.wrshort_real("sumdekin", sum_delta_ekin_ev + (ker_formation_ev if ker_formation_ev > 0 else 0.0)) # Add formation KER to dEkin for next step's Ekin calc
    iomod.wrshort_real("sumdeint", sum_delta_eint_ev)
    iomod.wrshort_real("x_trav", distance_traveled_m)
    iomod.wrshort_real("t_trav", time_traveled_s)
    # qmass is specific to a pair, written by mcsimu/reaction, not this routine

    return current_eiee_dist, piee_list_in # piee is usually not modified here

if __name__ == '__main__':
    print("Testing cid.py functions...")
    # Dummy env
    env_test_cid = RunTypeData()
    env_test_cid.cid_mgas = 39.948 # Argon
    env_test_cid.cid_rgas = 3.55266638 # Argon radius in Bohr
    env_test_cid.cid_Tgas = 300 # K
    env_test_cid.cid_Pgas = 0.132 # Pa
    env_test_cid.cid_lchamb = 0.25 # m
    env_test_cid.cid_collw = 0.5
    env_test_cid.cid_scool = 1.0
    env_test_cid.eyring = True # Important for simcid

    # Test _collision_setup_py
    # Dummy molecule (e.g. Benzene C6H6)
    nat_test = 12
    ats_test = [6]*6 + [1]*6
    coords_test_angstrom = [ # Simplified planar coords
        [0.0, 1.39, 0.0], [1.20, 0.695, 0.0], [1.20, -0.695, 0.0],
        [0.0, -1.39, 0.0], [-1.20, -0.695, 0.0], [-1.20, 0.695, 0.0],
        [0.0, 2.47, 0.0], [2.14, 1.235, 0.0], [2.14, -1.235, 0.0],
        [0.0, -2.47, 0.0], [-2.14, -1.235, 0.0], [-2.14, 1.235, 0.0]
    ]
    print("\n--- Testing _collision_setup_py ---")
    mfp, n_coll_avg = _collision_setup_py(env_test_cid, nat_test, ats_test, coords_test_angstrom)
    print(f"Calculated MFP: {mfp:.3e} m, Avg Collisions: {n_coll_avg:.2f}")

    # Test _simulate_single_collision_py
    print("\n--- Testing _simulate_single_collision_py ---")
    test_eiee = [1.0, 1.5, 2.0] # eV
    ion_mass_test = utility.get_average_mol_mass_py(ats_test) # Benzene ~78 amu
    ion_ekin_test = 30.0 # eV
    
    print(f"Initial EIEE: {test_eiee}, Ion Mass: {ion_mass_test:.2f} amu, Ion Ekin: {ion_ekin_test:.2f} eV")
    dEkin, dEint_avg = _simulate_single_collision_py(env_test_cid, test_eiee, ion_mass_test, ion_ekin_test)
    print(f"After collision: dEkin = {dEkin:.2f} eV, dEint_avg = {dEint_avg:.2f} eV")
    print(f"Modified EIEE: {test_eiee}")

    # Full simulate_cid_activation_py is complex to test standalone without extensive mocking of inputs
    # like alldgs, and other precursor data files.
```
