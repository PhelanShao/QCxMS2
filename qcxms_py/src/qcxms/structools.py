import math
from pathlib import Path
from typing import Tuple, Optional, List, Dict, Union

try:
    from .data import RunTypeData, WP
    from . import iomod
    from . import utility
    from . import qmmod # For potential XTB calls if wbo file is missing
    from .isotopes import AVERAGE_ATOMIC_MASSES # For CMAv
    from .constants import AUTOAA # Bohr to Angstrom
except ImportError:
    # Fallbacks for standalone/testing
    print("Attempting to import dummy/mock modules for structools.py standalone run.")
    from data import RunTypeData, WP # type: ignore
    import iomod_mock as iomod # type: ignore
    import utility_mock as utility # type: ignore
    import qmmod_mock as qmmod # type: ignore
    from isotopes import AVERAGE_ATOMIC_MASSES # type: ignore
    AUTOAA = 0.529177210903


# --- Data from Fortran module ---

# Pauling electronegativities (used for cn_damp_d4, though that was commented out)
PAULING_EN: List[WP] = [
    0.0, # Dummy Z=0
    2.20, 3.00,  # H,He (He EN is often undefined or very high, 3.00 is unusual, often given as N/A or ~4.5-5.5 if forced)
    0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98, 4.50,  # Li-Ne
    0.93, 1.31, 1.61, 1.90, 2.19, 2.58, 3.16, 3.50,  # Na-Ar
    0.82, 1.00,  # K,Ca
    1.36, 1.54, 1.63, 1.66, 1.55, 1.83, 1.88, 1.91, 1.90, 1.65,  # Sc-Zn
    1.81, 2.01, 2.18, 2.55, 2.96, 3.00,  # Ga-Kr
    0.82, 0.95,  # Rb,Sr
    1.22, 1.33, 1.60, 2.16, 1.90, 2.20, 2.28, 2.20, 1.93, 1.69,  # Y-Cd
    1.78, 1.96, 2.05, 2.10, 2.66, 2.60,  # In-Xe
    0.79, 0.89,  # Cs,Ba
    1.10, 1.12, 1.13, 1.14, 1.15, 1.17, 1.18, 1.20, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26,  # La-Lu
    1.27, 1.30, 1.50, 2.36, 1.90, 2.20, 2.20, 2.28, 2.54, 2.00,  # Hf-Hg
    1.62, 2.33, 2.02, 2.00, 2.20, 2.20,  # Tl-Rn
] + [1.50]*32 # Dummy values for Z > 86 as in Fortran

# Covalent radii from Pyykko and Atsumi, scaled as in D3 definition (R0 = R_pyykko * const)
# In Fortran: RCOV(i) = R_pyykko_angstrom(i) * AUTOAA * 4.0/3.0. This gives R0 in Bohr.
# We'll store R_pyykko_angstrom and apply scaling factor when used.
_RCOV_PYYKKO_ANGSTROM: List[WP] = [
    0.0, # Dummy Z=0
    0.32, 0.46,  # H,He
    1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67,  # Li-Ne
    1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96,  # Na-Ar
    1.76, 1.54,  # K,Ca
    1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09,  # Sc-Zn
    1.12, 1.09, 1.15, 1.10, 1.14, 1.17,  # Ga-Kr
    1.89, 1.67,  # Rb,Sr
    1.47, 1.39, 1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23,  # Y-Cd
    1.28, 1.26, 1.26, 1.23, 1.32, 1.31,  # In-Xe
    2.09, 1.76,  # Cs,Ba
    1.62, 1.47, 1.58, 1.57, 1.56, 1.55, 1.51, 1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53,  # La-Lu
    1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32,  # Hf-Hg
    1.30, 1.30, 1.36, 1.31, 1.38, 1.42,  # Tl-Rn
] + [1.50]*32 # Placeholder for Z > 86, as in Fortran's EN array

_RCOV_SCALING_FACTOR = AUTOAA * (4.0 / 3.0) # To get R0_ij in Bohr for damping functions

def get_scaled_covalent_radius_bohr(atomic_number: int) -> WP:
    if 0 < atomic_number < len(_RCOV_PYYKKO_ANGSTROM):
        return _RCOV_PYYKKO_ANGSTROM[atomic_number] * _RCOV_SCALING_FACTOR
    print(f"Warning: Scaled covalent radius for Z={atomic_number} not found. Using default 1.5 Bohr.")
    return 1.5 * _RCOV_SCALING_FACTOR # Should be a Bohr value

# Radii for fragment_structure (these are R0_ij / rcut, where rcut is typically 1.2-1.3 for VdW or 1.0 for covalent)
# The Fortran Rad was: Rad(118) = aatoau * [ ... angstrom values ... ]
# This implies the list was raw radii in Angstrom, then converted to Bohr.
# fragment_structure uses rcov = rcut * 0.5 * (Rad(oz(i)) + Rad(oz(j)))
# This means Rad(oz(i)) should be sum of covalent radii / rcut.
# For now, using _ATOMIC_RADII_ANGSTROM directly assuming they are appropriate covalent radii in Angstrom for this.
# And rcut will be applied. This matches the structure of `fragment_structure`.
_RADII_FOR_CONNECTIVITY_ANGSTROM: List[WP] = _ATOMIC_RADII_ANGSTROM # From cid.py (originally from structools.f90)


# --- CN Damping Functions ---
def _cn_damp_exp_py(r_sum_cov_bohr: WP, r_ij_bohr: WP) -> WP:
    """Classic exponential CN damping factor."""
    kcn = 4.0  # Fortran: kcn = 16.0_wp, then changed to 4.0_wp in comments
    k2 = 1.5
    # Check for r_ij_bohr being zero or very small
    if r_ij_bohr < 1e-6: return 1.0 # Max damping if atoms are on top of each other
    val = k2 * r_sum_cov_bohr / r_ij_bohr - 1.0
    try:
        return 1.0 / (1.0 + math.exp(-kcn * val))
    except OverflowError: # exp(-kcn*val) too large or too small
        return 0.0 if -kcn * val < 0 else 1.0


def _cn_damp_exp_lrange_py(r_sum_cov_bohr: WP, r_ij_bohr: WP) -> WP:
    """Exponential CN damping factor with different kcn, for longer range."""
    kcn = 8.0
    if r_ij_bohr < 1e-6: return 1.0
    val = r_sum_cov_bohr / r_ij_bohr - 1.0
    try:
        return 1.0 / (1.0 + math.exp(-kcn * val))
    except OverflowError:
        return 0.0 if -kcn * val < 0 else 1.0

# --- Core CN and Bond Matrix Calculation ---
def get_coordination_numbers_and_bond_matrix_py(
    fname_xyz: str, 
    use_long_range_damping: bool = False
) -> Tuple[Optional[List[WP]], Optional[List[List[WP]]]]:
    """
    Calculates coordination numbers (CN) and a bond order-like matrix.
    Reads XYZ from fname_xyz. Coordinates are assumed/converted to Bohr.
    Returns (list_of_CNs, bond_matrix) or (None, None) on failure.
    """
    # Read structure using utility function (assuming it returns nat, atomic_numbers, coords_bohr)
    # nat, atomic_numbers, coords_bohr = utility.read_xyz_structure_py(fname_xyz) # Placeholder
    # For now, use mctc_io equivalent (assuming it's in utility or iomod)
    
    # This needs a robust XYZ reader. Using a placeholder.
    # nat, atomic_numbers, coords_angstrom = iomod.read_xyz_detailed(fname_xyz)
    # if nat == 0: return None, None
    # coords_bohr = [[c * (1.0/AUTOAA) for c in atom_coord] for atom_coord in coords_angstrom]
    
    # Using utility.read_structure_py which should provide coordinates in Angstrom
    # This function is not yet defined in utility.py based on previous analysis.
    # Let's assume a simplified direct read here for now for testing.
    try:
        with open(fname_xyz, 'r') as f:
            nat = int(f.readline().strip())
            f.readline() # Comment
            atomic_symbols: List[str] = []
            coords_angstrom: List[List[WP]] = []
            for _ in range(nat):
                parts = f.readline().split()
                atomic_symbols.append(parts[0])
                coords_angstrom.append([float(c) for c in parts[1:4]])
        atomic_numbers = [utility.element_symbol_to_atomic_number_py(sym) for sym in atomic_symbols]
        coords_bohr = [[c / AUTOAA for c in atom_coord] for atom_coord in coords_angstrom]

    except Exception as e:
        print(f"Error reading XYZ file {fname_xyz} in get_CN: {e}")
        return None, None

    cn = [0.0] * nat
    bond_matrix = [[0.0] * nat for _ in range(nat)]

    for i in range(nat):
        z_i = atomic_numbers[i]
        r_cov_i_bohr = get_scaled_covalent_radius_bohr(z_i) # R0_i in Bohr

        for j in range(i + 1, nat): # Loop j from i+1 to nat-1 (0-indexed)
            z_j = atomic_numbers[j]
            r_cov_j_bohr = get_scaled_covalent_radius_bohr(z_j) # R0_j in Bohr
            
            r_sum_cov_bohr = r_cov_i_bohr + r_cov_j_bohr # R0_ij in Bohr

            diff_x = coords_bohr[i][0] - coords_bohr[j][0]
            diff_y = coords_bohr[i][1] - coords_bohr[j][1]
            diff_z = coords_bohr[i][2] - coords_bohr[j][2]
            r_ij_sq_bohr = diff_x**2 + diff_y**2 + diff_z**2
            
            # Fortran: if (r2 > cn_thr) cycle; cn_thr = 625 (Bohr^2) (i.e., 25 Bohr cutoff)
            if r_ij_sq_bohr > 625.0: 
                continue
            
            r_ij_bohr = math.sqrt(r_ij_sq_bohr)
            if r_ij_bohr < 1e-6 : # Atoms are essentially on top of each other
                damp = 1.0 # Treat as fully bonded/connected
            elif use_long_range_damping:
                damp = _cn_damp_exp_lrange_py(r_sum_cov_bohr, r_ij_bohr)
            else:
                damp = _cn_damp_exp_py(r_sum_cov_bohr, r_ij_bohr)
            
            bond_matrix[i][j] = damp
            bond_matrix[j][i] = damp # Symmetric
            cn[i] += damp
            cn[j] += damp
            
    # if env.printlevel >=3: # Assuming env is available or passed
    #     for i in range(nat): print(f"CN of atom {i+1} ({atomic_symbols[i]}): {cn[i]:.3f}")
            
    return cn, bond_matrix


# --- Further functions to be translated: ---
# readrhbond -> read_rh_bond_partners_py
# compare_rhbond -> compare_rh_bond_partners_py
# findhdiss -> find_h_dissociation_sumformula_py
# findhdisswbo -> find_h_dissociation_wbo_py
# get_active_cn -> get_active_atom_average_cn_py
# compare_bond -> _compare_bond_matrices_py (internal helper)
# get_active_cnstring -> get_active_atom_string_for_orca_py
# isdissociated -> is_structure_dissociated_py
# fragment_structure -> _get_connected_fragments_py (core connectivity algorithm)
# wrxyz_file, wrxyz_file_mask, wrxyz -> utility or iomod
# CMAv -> utility.calculate_center_of_mass_py
# center_of_geometry -> utility.calculate_center_of_geometry_py (already in utility.py)
# rdxmol, coordline -> utility or iomod (XYZ parsing)

if __name__ == '__main__':
    print("Testing structools.py (primarily CN calculation for now)")
    # Create a dummy RunTypeData if needed by any functions directly
    # env_test_struct = RunTypeData()

    # Create a dummy water.xyz file
    water_xyz_content = """3
Water molecule
O  0.000000  0.000000  0.117300
H  0.000000  0.757200 -0.469200
H  0.000000 -0.757200 -0.469200
"""
    Path("water.xyz").write_text(water_xyz_content)

    # Mock utility functions if not fully implemented in a utility_mock.py
    if not hasattr(utility, 'element_symbol_to_atomic_number_py'):
        def mock_e2z(sym):
            return {'O': 8, 'H': 1}.get(sym, 0)
        utility.element_symbol_to_atomic_number_py = mock_e2z # type: ignore

    cn_water, bond_matrix_water = get_coordination_numbers_and_bond_matrix_py("water.xyz")
    if cn_water and bond_matrix_water:
        print("\nWater (standard damping):")
        for i, c_val in enumerate(cn_water):
            print(f"  Atom {i} CN: {c_val:.3f}")
        # print("  Bond Matrix:")
        # for row in bond_matrix_water:
        #     print(f"    [{', '.join(f'{x:.3f}' for x in row)}]")
    
    cn_water_lr, _ = get_coordination_numbers_and_bond_matrix_py("water.xyz", use_long_range_damping=True)
    if cn_water_lr:
        print("\nWater (long-range damping):")
        for i, c_val in enumerate(cn_water_lr):
            print(f"  Atom {i} CN: {c_val:.3f}")

    Path("water.xyz").unlink(missing_ok=True)
```
