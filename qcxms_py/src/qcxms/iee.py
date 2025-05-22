import math
import random
from pathlib import Path
from typing import Tuple, Optional, List

try:
    from .data import RunTypeData, WP
    from .boxmuller import vary_energies_py
    from . import iomod
    from . import utility
except ImportError:
    # Fallbacks for standalone/testing
    print("Attempting to import dummy/mock modules for iee.py standalone run.")
    from data import RunTypeData, WP # type: ignore
    from boxmuller import vary_energies_py # type: ignore
    import iomod_mock as iomod # type: ignore
    import utility_mock as utility # type: ignore

# --- Distribution Functions ---

def _calculate_gaussian_like_prob_py(iee_a: WP, iee_b: WP, num_valence_electrons: WP, energy_x: WP) -> WP:
    """
    Calculates probability P(E) using a Gaussian-like function.
    p = exp(-iee_a * (x - num_valence_electrons*iee_b)**2 / num_valence_electrons)
    """
    if num_valence_electrons == 0: return 0.0 # Avoid division by zero
    exponent = -iee_a * (energy_x - num_valence_electrons * iee_b)**2 / num_valence_electrons
    try:
        return math.exp(exponent)
    except OverflowError:
        return 0.0 # Result of exponent too small

def _calculate_poisson_like_prob_py(iee_a: WP, iee_b: WP, num_valence_electrons: WP, energy_x: WP) -> WP:
    """
    Calculates probability P(E) using a Poisson-like function from Angew. paper.
    """
    if iee_a == 0 or num_valence_electrons == 0 or energy_x == 0: # Avoid division by zero or log(0)
        return 0.0

    z = iee_b
    k = 1.0 / iee_a
    
    t2 = k / num_valence_electrons  # This is 'c' in some notations

    # Argument of log: z/k * num_valence_electrons / x
    log_arg = (z / k * num_valence_electrons) / energy_x
    if log_arg <= 0: return 0.0 # Log undefined or will lead to issues

    t8_log_term = math.log(log_arg)
    
    # Exponent: t2*x*(1 + t8_log_term) - z
    exponent = t2 * energy_x * (1.0 + t8_log_term) - z
    
    try:
        t14_exp_term = math.exp(exponent)
    except OverflowError: # If exponent is too large or too small
        return 0.0 

    # Denominator for t17: (t2*x + 1)
    t17_denom_sqrt_arg = t2 * energy_x + 1.0
    if t17_denom_sqrt_arg <= 0: return 0.0 # Sqrt undefined or complex

    t17_factor = t17_denom_sqrt_arg**(-0.5)
    
    t18_prob = t14_exp_term * t17_factor
    return t18_prob

# --- Helper Functions for Parameter Determination ---

def _calculate_dist_properties_py(
    iee_a: WP, iee_b: WP, num_valence_electrons: WP, 
    is_poisson: bool, max_scan_energy: WP
) -> Tuple[WP, WP, WP]:
    """
    Calculates p_max (max probability), iee_max (energy at p_max), and e_avg (average energy)
    for a given distribution by scanning.
    """
    x = 0.001
    step = 0.01
    p_max = -1.0
    iee_max = 0.0
    e_avg_numerator = 0.0
    sum_probabilities = 0.0

    while x < max_scan_energy:
        prob_val: WP
        if is_poisson:
            prob_val = _calculate_poisson_like_prob_py(iee_a, iee_b, num_valence_electrons, x)
        else:
            prob_val = _calculate_gaussian_like_prob_py(iee_a, iee_b, num_valence_electrons, x)

        if prob_val > p_max:
            p_max = prob_val
            iee_max = x
        
        e_avg_numerator += prob_val * x
        sum_probabilities += prob_val
        x += step
        if x >= max_scan_energy and p_max < 0 : # Ensure p_max is non-negative if loop finishes
             p_max = 0.0


    e_avg = e_avg_numerator / sum_probabilities if sum_probabilities > 0 else 0.0
    if p_max < 0.0: p_max = 0.0 # Ensure p_max is not negative if no peak found
    
    return p_max, iee_max, e_avg


def _determine_iee_distribution_params_py(
    env: RunTypeData, 
    num_valence_electrons: WP, 
    is_poisson: bool, 
    max_excess_energy: WP, # exc from Fortran
    num_atoms: int
) -> Tuple[WP, WP, WP]:
    """
    Determines parameters iee_a and iee_b for the distribution by fitting
    to match env.ieeatm (average energy per atom).
    Also returns p_max.
    """
    step_param = 0.005
    current_iee_a = 0.0
    current_iee_b = 0.0
    
    # Target average energy per atom
    target_avg_energy_per_atom = env.ieeatm

    for k_iter in range(1, 10001): # Max 10000 iterations
        current_iee_a = min(current_iee_a + step_param, 0.3) # Cap iee_a at 0.3
        current_iee_b += step_param * 7.0

        p_max, _, e_avg = _calculate_dist_properties_py(
            current_iee_a, current_iee_b, num_valence_electrons, 
            is_poisson, max_excess_energy
        )
        
        current_avg_energy_per_atom = e_avg / num_atoms if num_atoms > 0 else 0.0

        if current_avg_energy_per_atom >= target_avg_energy_per_atom:
            return current_iee_a, current_iee_b, p_max
        
        if k_iter == 10000:
            print("Error: _determine_iee_distribution_params_py failed to converge within 10000 iterations.")
            # Fallback or raise error - Fortran used STOP
            # For now, return current params which might be inaccurate
            return current_iee_a, current_iee_b, p_max
            
    return current_iee_a, current_iee_b, 0.0 # Should be unreachable if loop runs


def _get_molecular_specs_py(env: RunTypeData) -> Tuple[int, WP, WP]:
    """
    Reads molecule, gets number of atoms, number of valence electrons (neutral), and E_HOMO (from env).
    """
    # This relies on utility.py having a robust XYZ reader
    atomic_numbers = utility.get_atomic_numbers_from_xyz_py(env.infile)
    if not atomic_numbers:
        raise ValueError(f"Could not read atomic composition from {env.infile}")
    
    num_atoms = len(atomic_numbers)
    
    # Calculate valence electrons for the (potentially charged) input species
    # utility.velectrons_amount_py should take list of Z values
    n_valence_charged = utility.velectrons_amount_py(num_atoms, atomic_numbers, env.chrg)
    
    # Number of valence electrons for the neutral species
    n_valence_neutral = float(n_valence_charged + env.chrg)
    
    e_homo_ev = env.ehomo # This is likely IP in eV from input or prior calc
    
    return num_atoms, n_valence_neutral, e_homo_ev


# --- Main IEE Distribution Generation Functions ---

def _get_iee_distribution_py(
    env: RunTypeData, 
    nsamples: int, 
    plotlog: bool
) -> Tuple[List[WP], List[WP]]:
    """
    Generates IEE distribution for EI mode using specified parameters and functions.
    """
    eiee_list: List[WP] = [0.0] * nsamples
    piee_list: List[WP] = [0.0] * nsamples

    num_atoms, num_valence_electrons, e_homo_ev = _get_molecular_specs_py(env)
    
    target_avg_energy_per_atom = env.ieeatm
    electron_impact_energy = env.eimp0 # e.g., 70 eV
    impact_energy_width = env.eimpw   # e.g., 0.1 (fractional width for vary_energies)

    # Max excess energy available after ionization by electron_impact_energy
    max_excess_energy_after_ionization = electron_impact_energy - e_homo_ev
    if max_excess_energy_after_ionization <= 0:
        print(f"Warning: E_HOMO ({e_homo_ev} eV) is >= electron impact energy ({electron_impact_energy} eV). No excess energy for IEE.")
        return eiee_list, piee_list # All zeros

    is_poisson_dist = (env.edist.lower() == 'poisson')
    dist_type_str = "Poisson" if is_poisson_dist else "Gaussian-like"
    print(f"Generating IEE using {dist_type_str} distribution.")

    iee_a, iee_b, p_max_dist = _determine_iee_distribution_params_py(
        env, num_valence_electrons, is_poisson_dist, 
        max_excess_energy_after_ionization, num_atoms
    )
    if p_max_dist == 0: # Could happen if determination failed or dist is flat
        print("Warning: p_max for IEE distribution is zero. Probabilities might be incorrect.")
        # Avoid division by zero later if p_max_dist is used for normalization.
        # The Fortran code normalizes by p_max for rejection sampling.

    for i in range(nsamples):
        # 1. Generate a trial electron energy, ensuring it's above HOMO energy
        trial_electron_e: WP
        while True:
            trial_electron_e = vary_energies_py(electron_impact_energy, impact_energy_width)
            if trial_electron_e >= e_homo_ev:
                break
        
        # 2. Calculate max IEE for this trial electron
        current_max_iee = trial_electron_e - e_homo_ev
        if current_max_iee <= 0: # Should be rare if above loop works
            eiee_list[i] = 0.0
            piee_list[i] = 0.0 # Or some minimal probability
            continue

        # 3. Rejection sampling for the actual IEE value
        accepted_energy: WP
        accepted_prob_unnormalized: WP # Probability from gauss0/poiss0 before dividing by p_max
        while True:
            random_fraction = random.random()
            trial_iee = random_fraction * current_max_iee # Uniformly sample energy up to current_max_iee
            
            prob_at_trial_iee: WP
            if is_poisson_dist:
                prob_at_trial_iee = _calculate_poisson_like_prob_py(iee_a, iee_b, num_valence_electrons, trial_iee)
            else:
                prob_at_trial_iee = _calculate_gaussian_like_prob_py(iee_a, iee_b, num_valence_electrons, trial_iee)

            # Normalize for rejection: prob_to_compare = prob_at_trial_iee / p_max_dist
            # (if p_max_dist is truly the max of the function up to max_excess_energy_after_ionization)
            # The sampling range for trial_iee (up to current_max_iee) can be smaller than the range
            # for which p_max_dist was determined. This could affect rejection efficiency
            # but the principle should hold if p_max_dist is a global max for relevant energies.
            
            if p_max_dist > 0: # Avoid division by zero
                 prob_threshold_for_acceptance = prob_at_trial_iee / p_max_dist
            else: # If p_max is zero, accept if prob_at_trial_iee is also zero (flat zero dist), else reject
                 prob_threshold_for_acceptance = 1.0 if prob_at_trial_iee == 0 else 0.0


            if random.random() <= prob_threshold_for_acceptance:
                accepted_energy = trial_iee
                accepted_prob_unnormalized = prob_at_trial_iee 
                break
        
        eiee_list[i] = accepted_energy
        piee_list[i] = accepted_prob_unnormalized # Storing the P(E) value from dist func

        if env.fixe != 0.0: # Override for testing with fixed energy
            eiee_list = [env.fixe] * nsamples
            # piee might need adjustment too if fixe is used; Fortran didn't specify piee for fixe.
            # Assume uniform probability if fixed energy is used, or 1.0 for the single energy.
            # For now, let it be overwritten by the loop if fixe is set after loop.
            # Fortran sets eiee = env%fixe *inside* the loop, but piee(i) is still prob.
            # This means if fixe is on, all eiee are same, but piee can vary.
            # A more logical approach for fixe would be eiee = [fixed_E], piee = [1.0]
            # For now, matching Fortran:
            eiee_list[i] = env.fixe


    if env.fixe != 0.0 and nsamples > 0 : # If fixe was used, ensure all energies are set
        # This is redundant if done inside loop, but ensures final state
        # piee would still contain varied probabilities from the last accepted sample before override.
        # This part of logic might need refinement based on intended use of piee with fixe.
        # For now, piee will just be what the last sample before fixe override had.
        # A common use of fixe is to have a single energy with probability 1.0.
        # If so: eiee_list = [env.fixe]; piee_list = [1.0]; nsamples_actual = 1
        pass


    if plotlog:
        try:
            with open('eimp.dat', 'w') as f:
                for i in range(nsamples):
                    f.write(f"{eiee_list[i]} {piee_list[i]}\n")
        except IOError:
            print("Error writing eimp.dat")
            
    return eiee_list, piee_list


def _get_esi_scaled_distribution_py(
    env: RunTypeData, 
    nsamples: int, 
    plotlog: bool
) -> Tuple[List[WP], List[WP]]:
    """Generates energy distribution for CID mode (scaled ESI energies)."""
    eiee_list: List[WP] = [0.0] * nsamples
    # For CID, piee is often treated as uniform as these are starting internal energies
    # before collision, not probabilities of formation.
    piee_list: List[WP] = [1.0 / nsamples if nsamples > 0 else 0.0] * nsamples

    if env.cid_mode == 3: # Forced collisions only, no initial internal energy
        # eiee_list is already all zeros
        return eiee_list, piee_list

    # Read number of atoms from input file
    num_atoms = iomod.rdshort_int(env.infile)
    if num_atoms == 0 and Path(env.infile).exists(): # If rdshort_int failed but file is there
        print(f"Warning: Could not read atom count from {env.infile} for ESI scaling. Assuming 0 atoms.")
    elif not Path(env.infile).exists():
        print(f"Error: Input file {env.infile} not found for ESI scaling.")
        return eiee_list, piee_list # All zeros

    cid_esiatom = env.cid_esiatom
    if cid_esiatom == 0: # Default if not set by user
        cid_esiatom = 0.4 if env.cid_mode == 2 else 0.1 # temprun vs auto/other

    # Mean internal energy for ESI-generated ions
    mean_escale = cid_esiatom * num_atoms + env.cid_esi
    
    print(f"ESI internal energy scaling: {cid_esiatom:.2f} eV/atom, shift: {env.cid_esi:.2f} eV.")
    print(f"Mean internal energy (ESCALE) for ESI: {mean_escale:.2f} eV for {num_atoms} atoms.")

    # Width of the distribution for ESI internal energies
    esi_energy_width = env.cid_esiw 

    for i in range(nsamples):
        eiee_list[i] = vary_energies_py(mean_escale, esi_energy_width)
        # piee_list[i] remains 1.0/nsamples (uniform probability for this starting energy)
    
    if plotlog:
        try:
            with open('eimp.dat', 'w') as f:
                for i in range(nsamples):
                    f.write(f"{eiee_list[i]} {piee_list[i]}\n")
        except IOError:
            print("Error writing eimp.dat for ESI distribution.")
            
    return eiee_list, piee_list


# --- Main Public Function ---
def get_energy_distribution_py(
    env: RunTypeData, 
    nsamples: int, 
    plotlog: bool = False
) -> Tuple[List[WP], List[WP]]:
    """
    Public interface to get energy distribution based on mode (EI or CID).
    Returns (list_of_energies, list_of_probabilities_or_weights).
    """
    if env.mode.lower() == "ei":
        return _get_iee_distribution_py(env, nsamples, plotlog)
    elif env.mode.lower() == "cid":
        return _get_esi_scaled_distribution_py(env, nsamples, plotlog)
    else:
        print(f"Warning: Unknown mode '{env.mode}' in get_energy_distribution. Returning empty.")
        return [], []


if __name__ == '__main__':
    print("Testing iee.py functions...")
    
    # Dummy env for EI
    env_ei_test = RunTypeData()
    env_ei_test.mode = "ei"
    env_ei_test.infile = "test_ch4.xyz" # Needs a mock utility for get_atomic_numbers_from_xyz_py
    env_ei_test.ehomo = 12.6 # Example IP for CH4 in eV
    env_ei_test.eimp0 = 70.0
    env_ei_test.eimpw = 0.1
    env_ei_test.ieeatm = 0.8 
    env_ei_test.edist = "poisson" # or "gaussian"
    env_ei_test.fixe = 0.0

    # Mock utility functions for standalone test
    _original_get_atoms = utility.get_atomic_numbers_from_xyz_py
    _original_velectrons = utility.velectrons_amount_py
    def mock_get_atoms_ch4(fname): return [6, 1, 1, 1, 1] # CH4
    def mock_velectrons(nat, iats, chrg): return 8 # For CH4 neutral
    utility.get_atomic_numbers_from_xyz_py = mock_get_atoms_ch4 # type: ignore
    utility.velectrons_amount_py = mock_velectrons # type: ignore


    print("\n--- Testing EI mode (get_iee_distribution_py) ---")
    nsamples_test = 100 # Small number for testing
    energies_ei, probs_ei = get_energy_distribution_py(env_ei_test, nsamples_test, plotlog=True)
    
    if energies_ei:
        avg_e_ei = sum(e * p for e, p in zip(energies_ei, probs_ei)) / sum(probs_ei) if sum(probs_ei) > 0 else 0
        print(f"Generated {len(energies_ei)} EI samples. Example first 5:")
        for k in range(min(5, len(energies_ei))):
            print(f"  E: {energies_ei[k]:.2f} eV, P: {probs_ei[k]:.4f}")
        print(f"  Average energy (weighted by P(E) from dist func): {avg_e_ei:.2f} eV")
        if Path("eimp.dat").exists(): print("  EI eimp.dat created.")


    # Dummy env for CID
    env_cid_test = RunTypeData()
    env_cid_test.mode = "cid"
    env_cid_test.infile = "test_protein_fragment.xyz" # Needs mock utility and iomod
    env_cid_test.cid_mode = 1 # auto
    env_cid_test.cid_esiatom = 0.05 # eV/atom
    env_cid_test.cid_esi = 1.0 # eV base
    env_cid_test.cid_esiw = 0.2 # fractional width
    
    _original_rdshort_int = iomod.rdshort_int
    def mock_rdshort_int_cid(fname, default=0): return 50 # Say 50 atoms for protein fragment
    iomod.rdshort_int = mock_rdshort_int_cid # type: ignore


    print("\n--- Testing CID mode (get_esi_scaled_distribution_py) ---")
    energies_cid, weights_cid = get_energy_distribution_py(env_cid_test, nsamples_test, plotlog=True)
    if energies_cid:
        avg_e_cid = sum(energies_cid) / len(energies_cid) if energies_cid else 0
        print(f"Generated {len(energies_cid)} CID samples. Example first 5:")
        for k in range(min(5, len(energies_cid))):
            print(f"  E: {energies_cid[k]:.2f} eV, Weight: {weights_cid[k]:.4f}")
        print(f"  Average energy (simple avg of samples): {avg_e_cid:.2f} eV")
        if Path("eimp.dat").exists(): print("  CID eimp.dat created.")

    # Restore original mocks/functions if they were changed
    utility.get_atomic_numbers_from_xyz_py = _original_get_atoms # type: ignore
    utility.velectrons_amount_py = _original_velectrons # type: ignore
    iomod.rdshort_int = _original_rdshort_int # type: ignore
    
    # Clean up dummy files
    Path("eimp.dat").unlink(missing_ok=True)
```
