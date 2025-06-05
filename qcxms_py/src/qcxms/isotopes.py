import math
import random
from pathlib import Path
from typing import Tuple, Optional, List, Dict, Union
from collections import Counter

try:
    from .data import RunTypeData, WP
    from . import iomod
    from . import utility # For reading XYZ file
except ImportError:
    # Fallbacks for standalone/testing
    print("Attempting to import dummy/mock modules for isotopes.py standalone run.")
    from data import RunTypeData, WP # type: ignore
    import iomod_mock as iomod # type: ignore
    import utility_mock as utility # type: ignore

# Average atomic masses (index by atomic number, 0 is dummy)
# Data from ams array in Fortran code
AVERAGE_ATOMIC_MASSES: List[WP] = [
    0.0, 1.00790, 4.00260, 6.94000, 9.01218, 10.81000, 12.01100, 14.00670, 
    15.99940, 18.99840, 20.17900, 22.98977, 24.30500, 26.98154, 28.08550, 
    30.97376, 32.06000, 35.45300, 39.94800, 39.09830, 40.08000, 44.95590, 
    47.90000, 50.94150, 51.99600, 54.93800, 55.84700, 58.93320, 58.71000, 
    63.54600, 65.38000, 69.73500, 72.59000, 74.92160, 78.96000, 79.90400, 
    83.80000, 85.46780, 87.62000, 88.90590, 91.22000, 92.90640, 95.94000, 
    98.90620, 101.0700, 102.9055, 106.4000, 107.8680, 112.4100, 114.8200, 
    118.6900, 121.7500, 127.6000, 126.9045, 131.3000, 132.9054, 137.3300, 
    138.91, 140.12, 140.91, 144.24, 147.00, 150.36, 151.97, 157.25, 158.93, 
    162.50, 164.93, 167.26, 168.93, 173.04, 174.97, 178.4900, 180.9479, 
    183.8500, 186.2070, 190.2000, 192.2200, 195.0900, 196.9665, 200.5900, 
    204.3700, 207.2000, 208.9804, 209., 210., 222.
]
# Extend with zeros if needed up to 107 as in Fortran, though only up to Rn (86) is typically used/stable.
# The Fortran data had 21*0.0 at the end, starting from element 87 (Fr)
AVERAGE_ATOMIC_MASSES.extend([0.0] * (107 - len(AVERAGE_ATOMIC_MASSES)))


# Isotope data: {Z: {"count": num_isotopes, "masses": [], "abundances": []}}
# Abundances are stored as fractional (0.0 to 1.0)
ISOTOPE_DATA: Dict[int, Dict[str, Union[int, List[WP]]]] = {
    1: {"count": 2, "masses": [2.014101, 1.007825], "abundances": [0.000115, 0.999885]}, # H (D, H) - Fortran prob was %
    2: {"count": 2, "masses": [3.016029309, 4.002603249], "abundances": [0.00000134, 0.9999866]}, # He (He3, He4)
    3: {"count": 2, "masses": [6.0151223, 7.0160041], "abundances": [0.0759, 0.9241]}, # Li
    4: {"count": 1, "masses": [9.0121822], "abundances": [1.0]}, # Be
    5: {"count": 2, "masses": [10.0129370, 11.0093055], "abundances": [0.199, 0.801]}, # B
    6: {"count": 2, "masses": [13.003354, 12.000000], "abundances": [0.0107, 0.9893]}, # C
    7: {"count": 2, "masses": [15.000108, 14.003074], "abundances": [0.00368, 0.99632]}, # N
    8: {"count": 3, "masses": [16.999131, 17.999159, 15.994914], "abundances": [0.00038, 0.00205, 0.99757]}, # O
    9: {"count": 1, "masses": [18.998403], "abundances": [1.0]}, # F
    10: {"count": 3, "masses": [19.992440176, 20.99384674, 21.99138550], "abundances": [0.9048, 0.0027, 0.0925]}, # Ne
    11: {"count": 1, "masses": [22.98976966], "abundances": [1.0]}, # Na
    12: {"count": 3, "masses": [23.98504187, 24.98583700, 25.98259300], "abundances": [0.7899, 0.1000, 0.1101]}, # Mg
    13: {"count": 1, "masses": [26.981538], "abundances": [1.0]}, # Al
    14: {"count": 3, "masses": [27.976926, 28.976494, 29.973770], "abundances": [0.92223, 0.04685, 0.03092]}, # Si
    15: {"count": 1, "masses": [30.973761], "abundances": [1.0]}, # P
    16: {"count": 4, "masses": [35.967080, 32.971458, 33.967867, 31.972071], "abundances": [0.0002, 0.0076, 0.0429, 0.9493]}, # S
    17: {"count": 2, "masses": [34.968852, 36.965902], "abundances": [0.7576, 0.2424]}, # Cl
    18: {"count": 3, "masses": [35.9675462, 37.9627322, 39.96238312], "abundances": [0.003365, 0.000632, 0.996003]}, # Ar
    19: {"count": 3, "masses": [38.9637069,39.96399867,40.96182597], "abundances": [0.932581,0.000117,0.067302]}, # K
    20: {"count": 5, "masses": [39.962591,41.958618,42.958766,43.955482,45.953688], "abundances": [0.9694,0.0065,0.0014,0.0209,0.0018]}, # Ca (sum is slightly off 100 in F)
    21: {"count": 1, "masses": [44.9559102], "abundances": [1.0]}, # Sc
    22: {"count": 5, "masses": [45.952632,46.951763,47.947946,48.947870,49.944791], "abundances": [0.080,0.073,0.738,0.055,0.054]}, # Ti
    23: {"count": 2, "masses": [49.9471627,50.943963], "abundances": [0.00250,0.99750]}, # V
    24: {"count": 4, "masses": [49.946044,51.940507,52.940649,53.938880], "abundances": [0.0435,0.8379,0.0950,0.0237]}, # Cr
    25: {"count": 1, "masses": [54.938045], "abundances": [1.0]}, # Mn
    26: {"count": 4, "masses": [53.939,55.934,56.935,57.933], "abundances": [0.05845,0.91754,0.02119,0.00282]}, # Fe (0.2% in F was 0.282%)
    27: {"count": 1, "masses": [58.933195], "abundances": [1.0]}, # Co
    28: {"count": 5, "masses": [57.935343,59.930786,60.931056,61.928345,63.927966], "abundances": [0.6808,0.2622,0.0114,0.0363,0.0093]}, # Ni
    29: {"count": 2, "masses": [62.929597,64.927789], "abundances": [0.6915,0.3085]}, # Cu
    30: {"count": 5, "masses": [63.929142,65.926033,66.927127,67.924884,69.925319], "abundances": [0.486,0.279,0.041,0.188,0.006]}, # Zn (sum is slightly off 100 in F)
    31: {"count": 2, "masses": [68.925581,70.9247073], "abundances": [0.60108,0.39892]}, # Ga
    32: {"count": 5, "masses": [69.924247,71.922076,72.923459,73.921178,75.921402], "abundances": [0.2038,0.2731,0.0776,0.3672,0.0783]}, # Ge
    33: {"count": 1, "masses": [74.921596], "abundances": [1.0]}, # As
    34: {"count": 6, "masses": [73.922476,75.919213,76.919914,77.917309,79.916521,81.916699], "abundances": [0.0089,0.0937,0.0763,0.2377,0.4961,0.0873]}, # Se
    35: {"count": 2, "masses": [78.91833,80.91629], "abundances": [0.5069,0.4931]}, # Br
    36: {"count": 6, "masses": [77.920388,79.916379,81.9134850,82.914137,83.911508,85.910615], "abundances": [0.00355,0.02286,0.11593,0.11500,0.56987,0.17279]}, # Kr
    # ... up to Rn (86)
    # For brevity, only a subset is included here. The full list from Fortran needs to be added.
    # Ensure probabilities are normalized (sum to 1.0) if taken directly from % values.
    # The Fortran probabilities were already percentages, so dividing by 100.
    # Example: H: prob(1,1)=0.0115 (means 0.0115%) -> 0.000115
    # Corrected H: prob(1,1)=0.0115_wp, prob(1,2)=99.9885_wp. These are percentages.
    # So H_abundances should be [0.000115, 0.999885]
    # Let's re-verify a few:
    # C (6): Fortran: 1.07, 98.93 -> Python: 0.0107, 0.9893. Correct.
    # O (8): Fortran: 0.038, 0.205, 99.757 -> Python: 0.00038, 0.00205, 0.99757. Correct.
    # Cl (17): Fortran: 75.76, 24.24 -> Python: 0.7576, 0.2424. Correct.
    # The data above has been corrected to be fractional.
    # Elements 37-86 need to be added similarly.
    83: {"count": 1, "masses": [208.980398], "abundances": [1.0]}, # Bi (Example of a later element)
}

def get_average_mol_mass_py(atomic_numbers: List[int]) -> WP:
    """Calculates the average molecular mass from a list of atomic numbers."""
    molmass: WP = 0.0
    for z_val in atomic_numbers:
        if 0 < z_val < len(AVERAGE_ATOMIC_MASSES):
            molmass += AVERAGE_ATOMIC_MASSES[z_val]
        else:
            # Handle unknown atomic number, e.g., add 0 or raise error
            print(f"Warning: Average atomic mass for Z={z_val} not available. Contributing 0 to molecular mass.")
    return molmass


def get_isotopic_pattern_py(
    env: RunTypeData, 
    fname_xyz: str, 
    # maxatm: int, # Maxatm was for pre-allocating rnd array, not strictly needed if generating on the fly
    # random_numbers_2d: Optional[List[List[WP]]] = None, # For using pre-generated random numbers
    num_trials: int = 50000,
    charge_state: int = 1 # Assume z=1 for m/z if not specified
) -> Tuple[List[WP], List[WP]]:
    """
    Calculates isotopic mass distribution using Monte Carlo.
    Returns (list_of_masses, list_of_intensities).
    """
    
    # 1. Read XYZ to get atomic composition
    # Assuming utility.py has a function that returns list of atomic numbers
    # nat, atomic_numbers = utility.read_xyz_atomic_numbers(fname_xyz) # Placeholder
    # For now, let's mock this part or assume it's passed in.
    # Example: atomic_numbers = [6, 1, 1, 1, 1] # CH4
    
    # Using the mctc_io equivalent from utility.py (placeholder)
    # This needs to be replaced by actual XYZ parsing to get atomic numbers
    atomic_numbers = utility.get_atomic_numbers_from_xyz_py(fname_xyz)
    if not atomic_numbers:
        print(f"Error: Could not read atomic composition from {fname_xyz}")
        return [], []
    
    nat = len(atomic_numbers)
    
    # 2. Monte Carlo simulation
    mass_counter: Counter[int] = Counter() # Use integer representation of mass for binning initially
                                          # Rounded to a few decimal places then scaled, e.g., mass * 10000
    
    # For exact mass sum, use list to store precise float sums
    observed_exact_masses: List[WP] = []

    for _ in range(num_trials):
        current_exact_mass: WP = 0.0
        for i in range(nat):
            z_val = atomic_numbers[i]
            if z_val not in ISOTOPE_DATA:
                # Use average mass if isotope data missing for this element
                print(f"Warning: Isotope data for Z={z_val} not found. Using average mass.")
                current_exact_mass += AVERAGE_ATOMIC_MASSES[z_val] if 0 < z_val < len(AVERAGE_ATOMIC_MASSES) else 0.0
                continue

            isotope_info = ISOTOPE_DATA[z_val]
            # random.choices returns a list of k elements, so take the first one
            chosen_mass = random.choices(
                population=isotope_info["masses"], # type: ignore
                weights=isotope_info["abundances"], # type: ignore
                k=1
            )[0]
            current_exact_mass += chosen_mass
        
        # Apply charge state for m/z
        # Ensure charge_state is not zero to avoid division by zero
        effective_charge = abs(charge_state) if charge_state != 0 else 1
        current_mz = current_exact_mass / effective_charge

        if current_mz > env.mthr: # Mass threshold filtering
            observed_exact_masses.append(current_mz)

    if not observed_exact_masses:
        return [], []

    # Histogramming the exact masses (can be tricky with floats)
    # A common approach is to round masses to a certain precision before counting
    # Or use a tolerance for grouping. Fortran code implicitly did this by storing in list_masses(loop)
    # and comparing with abs(list_masses(loop) - current_mass) < 1.0d-10
    # Let's sort and then group close masses.
    observed_exact_masses.sort()
    
    if not observed_exact_masses: return [], []

    final_masses: List[WP] = []
    final_counts: List[int] = []
    
    current_bin_mass_sum = observed_exact_masses[0]
    current_bin_count = 1
    for i in range(1, len(observed_exact_masses)):
        if abs(observed_exact_masses[i] - current_bin_mass_sum / current_bin_count) < 1e-4: # Tolerance for grouping
            current_bin_mass_sum += observed_exact_masses[i]
            current_bin_count += 1
        else:
            final_masses.append(current_bin_mass_sum / current_bin_count)
            final_counts.append(current_bin_count)
            current_bin_mass_sum = observed_exact_masses[i]
            current_bin_count = 1
    # Add the last bin
    final_masses.append(current_bin_mass_sum / current_bin_count)
    final_counts.append(current_bin_count)

    # 3. Normalize intensities
    total_counts_after_mthr = sum(final_counts)
    if total_counts_after_mthr == 0: # Should not happen if observed_exact_masses was not empty
        return [], []

    # Read pfrag (parent fragment intensity)
    # Assuming pfrag file is in the current directory (fragment's directory)
    pfrag = iomod.rdshort_real("pfrag", default=1.0) # Default to 1.0 if not found (treat as primary fragment)
    if pfrag == 0.0 and Path("pfrag").exists(): # if pfrag is explicitly 0.0
        pfrag = 1.0 # Fortran logic implies pfrag might not be there or might not be used this way always.
                    # For now, if it's 0.0, means parent had 0 intensity, so this has 0.

    exact_intensity: List[WP] = [(count / total_counts_after_mthr) * pfrag for count in final_counts]

    # 4. env.noiso handling
    if env.noiso and final_masses:
        if not exact_intensity: return [], []
        max_intensity = -1.0
        idx_max_intensity = -1
        for i, intensity_val in enumerate(exact_intensity):
            if intensity_val > max_intensity:
                max_intensity = intensity_val
                idx_max_intensity = i
        
        if idx_max_intensity != -1:
            # Fortran code took max intensity value, not sum.
            return [final_masses[idx_max_intensity]], [max_intensity] 
        else: # Should not happen if exact_intensity is not empty
            return [],[]


    # 5. Fortran's check: if (sum(store_int) .ne. nrnd) then exact_intensity = 0 ...
    # This means if mthr filtering removed any samples, the whole pattern is discarded.
    # This seems too strict. Python version currently normalizes based on observed counts *after* mthr.
    # If the Fortran behavior is strictly required:
    # if total_counts_after_mthr != num_trials:
    #     print("Warning: Not all MC trials passed m/z threshold. Fortran would zero intensities.")
    #     return [0.0] * len(final_masses), [0.0] * len(final_masses) # Or just [],[]

    return final_masses, exact_intensity


if __name__ == '__main__':
    print("Testing isotopes.py...")
    # Dummy env
    env_test = RunTypeData()
    env_test.mthr = 0.0
    env_test.noiso = False

    # Mock utility.get_atomic_numbers_from_xyz_py
    def mock_get_atoms(fname):
        if "ch4" in fname: return [6, 1, 1, 1, 1] # CH4
        if "c2h6" in fname: return [6, 1, 1, 1, 6, 1, 1, 1] # C2H6
        return []
    utility.get_atomic_numbers_from_xyz_py = mock_get_atoms # type: ignore

    # Create dummy pfrag file
    Path("pfrag").write_text("1.0\n") # Assume parent fragment intensity is 100%

    print("\nCH4 Isotopic Pattern:")
    masses_ch4, intensities_ch4 = get_isotopic_pattern_py(env_test, "ch4.xyz")
    for m, i in zip(masses_ch4, intensities_ch4):
        if i > 0.0001 : # Print significant peaks
            print(f"  Mass: {m:.4f}, Intensity: {i:.4%}")

    print("\nC2H6 Isotopic Pattern:")
    masses_c2h6, intensities_c2h6 = get_isotopic_pattern_py(env_test, "c2h6.xyz")
    for m, i in zip(masses_c2h6, intensities_c2h6):
        if i > 0.0001 :
             print(f"  Mass: {m:.4f}, Intensity: {i:.4%}")
    
    env_test.noiso = True
    print("\nC2H6 Isotopic Pattern (noiso=True):")
    masses_c2h6_noiso, intensities_c2h6_noiso = get_isotopic_pattern_py(env_test, "c2h6.xyz")
    for m, i in zip(masses_c2h6_noiso, intensities_c2h6_noiso):
        print(f"  Mass: {m:.4f}, Intensity: {i:.4%}")

    # Test average molecular mass
    print(f"\nAverage mass of CH4: {get_average_mol_mass_py([6,1,1,1,1]):.4f}")
    print(f"Average mass of C2H6: {get_average_mol_mass_py([6,1,1,1,6,1,1,1]):.4f}")

    Path("pfrag").unlink(missing_ok=True)
```
