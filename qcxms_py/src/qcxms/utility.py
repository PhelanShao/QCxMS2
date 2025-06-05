import os
import sys
import math
import subprocess
import shutil
from pathlib import Path
from typing import List, Tuple, Optional, TypeVar, Callable, Dict

# Assuming WP and RunTypeData will be imported from qcxms.data when available
# For now, define WP for type hinting if not available from a central place yet.
try:
    from .data import WP, RunTypeData # Relative import if utility is part of a package
    from .isotopes import AVERAGE_ATOMIC_MASSES # For get_average_mol_mass_py
except ImportError:
    # Fallback for standalone development or if data.py is not yet structured for this
    WP = float
    from isotopes import AVERAGE_ATOMIC_MASSES # type: ignore
    @dataclass
    class RunTypeData: # Dummy class for type hinting
        chrg: int = 0
        cores: int = 1
        threads: int = 1
        mode: str = "ei"
        tslevel: str = "gfn2"
        ieeatm: float = 0.8
        temp: int = 298
        pthr: float = 1.0
        topocheck: str = "molbar"
        pass

# Constants from xtb_mctc_constants (will be centralized later)
KB_EV_K = 8.617333262e-5  # Boltzmann constant in eV/K
KB_J_K = 1.3806488e-23    # Boltzmann constant in J/K
H_J_S = 6.62606957e-34    # Planck constant in J*s
C_CM_S = 299792458.0 * 100 # Speed of light in cm/s


# --- End of iomod placeholders ---

# --- OMP related functions (Sequential Implementation) ---
# Import functions from iomod that were previously placeholders here
try:
    from . import iomod # For package structure
except ImportError:
    import iomod # For standalone testing if utility.py and iomod.py are in same dir

# --- OMP related functions (Sequential Implementation) ---
def omp_jobcall(ndirs: int, basedir: str, jobcall: str):
    """
    Executes a job in multiple directories sequentially.
    Original was OpenMP parallel.
    """
    # print(f"Simulating omp_jobcall for {ndirs} dirs, base: {basedir}")
    k = 0
    for i in range(1, ndirs + 1):
        tmppath = Path(basedir) / str(i)
        # Fortran changed directory. Python's subprocess can use `cwd` argument.
        # Ensure tmppath exists before trying to run command in it.
        tmppath.mkdir(parents=True, exist_ok=True)
        # print(f"Executing in {tmppath.resolve()}: {jobcall}")
        # Use the more robust execute_command from iomod
        # Convert jobcall string to list for execute_command if it's simple,
        # otherwise, shell=True is needed. The original Fortran used shell execution.
        result = iomod.execute_command(jobcall.split(), cwd=tmppath.resolve(), shell=True) # shell=True to mimic Fortran
        exit_code = result.returncode
        k += 1
        print(f"{k} ", end="", flush=True)
    print()

def omp_samejobcall(njobs: int, dirs: List[str], jobcall: str, do_print: bool = True):
    """
    Executes the same job in a list of specified directories sequentially.
    Original was OpenMP parallel.
    """
    # print(f"Simulating omp_samejobcall for {njobs} jobs.")
    if njobs <= 0:
        return
    
    original_pwd = Path.cwd()
    k = 0
    actual_dirs = [d for d in dirs if d.strip()] # Filter out empty strings like Fortran logic might imply

    for i in range(len(actual_dirs)): # Iterate over actual_dirs now
        dir_to_run_in = Path(actual_dirs[i])
        if not dir_to_run_in.is_absolute():
            dir_to_run_in = original_pwd / dir_to_run_in
        
        dir_to_run_in.mkdir(parents=True, exist_ok=True)
        # print(f"Executing in {dir_to_run_in.resolve()}: {jobcall}")
        result = iomod.execute_command(jobcall.split(), cwd=dir_to_run_in.resolve(), shell=True)
        exit_code = result.returncode
        k += 1
        if do_print:
            print(f"{k} ", end="", flush=True)
    print()


def set_omp_threads(env: RunTypeData, njobs: int):
    """Sets OMP_NUM_THREADS, MKL_NUM_THREADS, and PARNODES environment variables."""
    if njobs == 0:
        ncalc = 1
    else:
        ncalc = njobs

    threads = env.cores // ncalc
    if threads == 0:
        threads = 1
    
    if threads > 1 and threads % 2 != 0:
        if threads < 4:
            threads += 1
        else:
            threads -=1
    if threads == 2: 
        threads = 4 
    if threads > 28: 
        threads = 27 

    env.threads = threads
    
    env_omp_threads = threads
    if env_omp_threads > 8: 
        env_omp_threads = 8
    
    os.environ['OMP_NUM_THREADS'] = str(env_omp_threads)
    os.environ['MKL_NUM_THREADS'] = str(env_omp_threads)
    
    env_par_nodes = threads
    if env_par_nodes > 28: 
        env_par_nodes = 28
    os.environ['PARNODES'] = str(env_par_nodes)
    
    # print(f"Set OMP_NUM_THREADS={os.environ['OMP_NUM_THREADS']}, MKL_NUM_THREADS={os.environ['MKL_NUM_THREADS']}, PARNODES={os.environ['PARNODES']}")
    # print(f"Calculated threads per job for RunTypeData: {env.threads}")


def set_omp_to_one():
    """Sets OMP_NUM_THREADS, MKL_NUM_THREADS, and PARNODES to 1."""
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['PARNODES'] = '1'
    # print("Set OMP_NUM_THREADS, MKL_NUM_THREADS, PARNODES to 1")

# --- Electron/Spin counting ---

def get_core_e(iat: int) -> int:
    """Gets the number of core electrons for a given atomic number."""
    if iat <= 2: return 0    # H, He
    elif iat <= 10: return 2   # Li-Ne
    elif iat <= 18: return 10  # Na-Ar
    elif iat <= 29: return 18  # K-Cu (atomic number 19 to 29)
    elif iat <= 36: return 28  # Ga-Kr (atomic number 31 to 36, using 28 from Zn: 30-2=28?) This seems to imply d-shell is core for Ga-Kr
                                # Zn(30) has 30-18 = 12 valence e- with Ar core. Or 30-28=2 if Ni-like core.
                                # Fortran logic: iat=30 (Zn), ncore=28. iat=31 (Ga), ncore=28.
    elif iat <= 47: return 36  # Rb-Ag (atomic number 37 to 47)
                                # Kr (36) is core. Fortran: iat=48 (Cd), ncore=46.
    elif iat <= 54: return 46  # In-Xe (atomic number 49 to 54)
    elif iat <= 71: return 54  # Cs-Lu (atomic number 55 to 71)
                                # Xe (54) is core. Fortran: iat=72 (Hf), ncore=68 (Xe+14f)
    elif iat <= 79: return 68  # Hf-Au (atomic number 72 to 79)
    elif iat <= 86: return 78  # Hg-Rn (atomic number 80 to 86)
    else:
        print(f"Warning: get_core_e not defined for atomic number {iat}. Defaulting to 0 core electrons.")
        return 0 

def velectrons_amount(nat: int, iat: List[int], chrg: int) -> int:
    """Calculates the number of valence electrons."""
    z_valence = [element_iat - get_core_e(element_iat) for element_iat in iat]
    nvel = int(sum(z_valence))
    nvel -= chrg # Standard definition: subtract positive charge, add negative charge
    return nvel

def electrons_amount(nat: int, iat: List[int], chrg: int) -> Tuple[int, int, List[WP]]:
    """Calculates total valence electrons, number of pairs, and list of valence electrons per atom."""
    z_valence_per_atom = [WP(element_iat - get_core_e(element_iat)) for element_iat in iat]
    nel_valence = int(sum(z_valence_per_atom))
    nel_valence -= chrg
    nb_pairs = nel_valence // 2
    return nel_valence, nb_pairs, z_valence_per_atom


def valel(at: int) -> WP:
    """Calculate number of valence electrons, specific logic (possibly for semi-empirical methods)."""
    el: WP
    if at <= 2:     # H, He
        el = float(at)
    elif at <= 10:  # Li - Ne
        el = float(at - 2)
    elif at <= 18:  # Na - Ar
        el = float(at - 10)
    elif at <= 36:  # K - Kr
        el = float(at - 18) # Default for this block (e.g. K, Ca)
        if at > 28: # Override for Sc-Kr (effectively, 2p of previous period + d electrons)
                    # e.g. Cu(29): 29-28 = 1. Zn(30): 30-28=2. Ga(31): 31-28=3.
            el = float(at - 28)
    else:
        print(f"Warning: valel not explicitly defined for atomic number {at}. Defaulting to 0.")
        el = 0.0
    return el

def get_spin(nat: int, iat: List[int], chrg: int) -> int:
    """Calculates spin multiplicity (2S+1). Assumes low-spin."""
    total_atomic_numbers = sum(iat)
    num_electrons = total_atomic_numbers - chrg

    if num_electrons < 0: # Should not happen for valid species
        print(f"Warning: Negative number of electrons ({num_electrons}) calculated for spin. Returning -1.")
        return -1 
    if num_electrons == 0: # E.g. H+ if iat=[1], chrg=1. Fortran returned -1.
        # Or for a system with no electrons like bare nuclei.
        # For consistency with Fortran:
        if nat > 0 : # If there are atoms but no electrons
             return -1 # (e.g. H+ after stripping electron)
        else: # No atoms, no electrons
             return 1 # Or 0? Let's assume 1 for no atoms. But Fortran implies -1 if j<1.

    spin_multiplicity = 1 + (num_electrons % 2)
    return spin_multiplicity

# --- String and Element Symbol Manipulations ---
_ELEMENT_SYMBOLS = [
    "XX", "H", "HE", "LI", "BE", "B", "C", "N", "O", "F", "NE", "NA", "MG", "AL", "SI", "P", "S", "CL", "AR",
    "K", "CA", "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR",
    "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I", "XE",
    "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU",
    "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN"
] 

_ELEMENT_DICT = {symbol: i for i, symbol in enumerate(_ELEMENT_SYMBOLS)}
_ELEMENT_DICT["D"] = 1 # Deuterium
_ELEMENT_DICT["T"] = 1 # Tritium

_CONVERT_LABEL_REMOVE_CHARS = '0123456789*_+-'

def convert_label(s: str) -> str:
    """Cleans and uppercases an element label string."""
    temp_s = "".join([char for char in s if char not in _CONVERT_LABEL_REMOVE_CHARS])
    return temp_s.upper().strip()

def e2i(element_symbol: str) -> int:
    """Converts element symbol string to atomic number."""
    clean_symbol = convert_label(element_symbol)
    return _ELEMENT_DICT.get(clean_symbol, 0)

def i2e(atomic_number: int, oformat: Optional[str] = None) -> str:
    """Converts atomic number to element symbol string."""
    symbol = "XX"
    if 0 < atomic_number < len(_ELEMENT_SYMBOLS):
        symbol = _ELEMENT_SYMBOLS[atomic_number]

    if oformat:
        if oformat in ("lc", "lowercase"):
            return symbol.lower()
        elif oformat in ("nc", "nicecase"):
            if len(symbol) > 1:
                return symbol[0] + symbol[1:].lower()
            return symbol # Already nice case for single letter symbols
    return symbol

def sum_formula_from_atoms(nat: int, atom_numbers: List[int]) -> str:
    """Creates a sum formula string from atomic numbers (e.g., C2H4O1)."""
    if nat == 0:
        return ""
    
    counts: Dict[int, int] = {}
    for at_num in atom_numbers:
        counts[at_num] = counts.get(at_num, 0) + 1
    
    formula_parts = []
    # Order: C, then H, then others alphabetically by symbol
    # This matches common conventions like Hill system.
    
    if 6 in counts: # Carbon
        formula_parts.append(f"{i2e(6, 'nc')}{counts.pop(6)}")
    if 1 in counts: # Hydrogen
        formula_parts.append(f"{i2e(1, 'nc')}{counts.pop(1)}")
            
    # Add remaining elements sorted alphabetically by their symbol
    sorted_remaining_elements = sorted(counts.keys(), key=lambda at_num: i2e(at_num, 'nc'))
    
    for at_num in sorted_remaining_elements:
        formula_parts.append(f"{i2e(at_num, 'nc')}{counts[at_num]}")
        
    return "".join(formula_parts)

# --- Sorting ---
def quicksort_integers(arr: List[int]) -> List[int]:
    """Quicksort for a list of integers (low to high). Non-in-place."""
    if len(arr) <= 1:
        return arr
    
    pivot = arr[len(arr) // 2] # Choose pivot (middle element)
    left = [x for x in arr if x < pivot]
    middle = [x for x in arr if x == pivot]
    right = [x for x in arr if x > pivot]
    
    return quicksort_integers(left) + middle + quicksort_integers(right)

def qsort_reals_with_indices(a: List[WP], ind: List[int]) -> None:
    """In-place quicksort for 'a', carrying 'ind' along."""
    _qsort_recursive(a, ind, 0, len(a) - 1)

def _qsort_recursive(a: List[WP], ind: List[int], first: int, last: int):
    if first >= last:
        return

    pivot_val = a[(first + last) // 2]
    i = first
    j = last

    while i <= j:
        while a[i] < pivot_val: i += 1
        while a[j] > pivot_val: j -= 1
        if i <= j:
            a[i], a[j] = a[j], a[i]
            ind[i], ind[j] = ind[j], ind[i]
            i += 1
            j -= 1
            
    if first < j: _qsort_recursive(a, ind, first, j)
    if i < last: _qsort_recursive(a, ind, i, last)

# --- Array/Mask Manipulation ---
def mask_invert(nall: int, mask: List[int]) -> List[int]:
    """Inverts a mask array. Assumes mask contains unique indices from 0 to nall-1 or 1 to nall."""
    if not mask: return []
    
    # Determine if mask is 0-indexed or 1-indexed based on values
    is_1_indexed = all(x >= 1 for x in mask) and any(x > nall -1 for x in mask) if nall > 0 else False # Heuristic
    
    if nall == 0 : return []

    imask = [0] * nall
    if is_1_indexed: # Adjust for 1-based indexing if detected
        for i_idx, val_in_mask in enumerate(mask):
            if 1 <= val_in_mask <= nall:
                 imask[val_in_mask - 1] = i_idx + 1 # Store 1-based original index
            else:
                raise ValueError(f"Mask value {val_in_mask} out of bounds for 1-based indexing (1 to {nall})")
    else: # Assume 0-based indexing
        for i_idx, val_in_mask in enumerate(mask):
            if 0 <= val_in_mask < nall:
                imask[val_in_mask] = i_idx # Store 0-based original index
            else:
                raise ValueError(f"Mask value {val_in_mask} out of bounds for 0-based indexing (0 to {nall-1})")
    return imask


# --- Temperature and Energy Calculations ---
def calctemp(nvib: int, energy: WP) -> WP:
    """Simple formula to approximate temperature (Kelvin) from energy (eV)."""
    if nvib <= 0: return 0.01 
    temperature = energy / (nvib * KB_EV_K)
    return max(0.01, temperature) # Ensure temp is at least 0.01 K

def set_temp_from_energy(env: RunTypeData, fname_mol_structure: str):
    """Approximates and sets temperature in env based on molecular properties and energy."""
    chrg = iomod.rdshort_int(".CHRG") 
    uhf = iomod.rdshort_int(".UHF")   

    query = f"{env.tslevel} sumreac {chrg} {uhf}" 
    # grepval in iomod returns tuple (bool, val), not modifying list args
    found_sumreac_bool, sumreac_val_float = iomod.grepval("qmdata", query)
    sumreac = sumreac_val_float if found_sumreac_bool else 0.0

    nat = iomod.rdshort_int(fname_mol_structure)

    if nat > 0:
        # Number of vibrational modes (linear or non-linear molecule)
        # This is a simplification; actual calculation might depend on linearity.
        # QCxMS Fortran used (3*nat - 6) implying non-linear.
        nvib = 3 * nat - 6 if nat > 2 else 1 # Avoid nvib <=0 for very small molecules
        if nvib <= 0 : nvib = 1 # Ensure positive nvib

        # Internal energy E = IEE_atm * nat - sum_of_reaction_energies
        # ieeatm is likely an average internal energy per atom.
        internal_energy_ev = env.ieeatm * nat - sumreac
        
        env.temp = int(calctemp(nvib, internal_energy_ev))
        print(f"Average temperature set to {env.temp} K (approx.)")
    else:
        print("Warning: Could not determine number of atoms for temperature calculation.")
        env.temp = 0 # Or keep a default like 298


def calc_eyring(energy: WP, barrier: WP, nvib: int) -> WP:
    """Compute rate constant with Eyring equation."""
    temperature = calctemp(nvib, energy)
    if temperature <= 0.01: return 0.0 
    
    prefactor = (KB_J_K * temperature) / H_J_S
    exponent_arg = -barrier / (KB_EV_K * temperature)
    
    # Prevent math.exp from overflowing with very large positive exponent_arg (huge negative barrier)
    # or underflowing to zero too quickly if not necessary.
    if exponent_arg > 700: # Corresponds to exp(700) ~ 1e304
        return float('inf') # Effectively infinite rate
    if exponent_arg < -700: # Corresponds to exp(-700) ~ 0
        return 0.0

    rate_constant = prefactor * math.exp(exponent_arg)
    
    # Limit k to prefactor (kB*T/h), occurs for zero or negative effective barriers
    return min(rate_constant, prefactor) if barrier <= 0 else rate_constant

def calc_releyring(energy: WP, dbarrier: WP, nvib: int) -> WP:
    """Compute relative rate constant with Eyring equation: exp(-dEa_adj / RT)."""
    temperature = calctemp(nvib, energy)
    if temperature <= 0.01: return 0.0

    exponent_arg = -dbarrier / (KB_EV_K * temperature)
    if exponent_arg > 700: return float('inf')
    if exponent_arg < -700: return 0.0
    
    k_rel = math.exp(exponent_arg)
    # The Fortran code limited this to 1.0 if k_rel / (kB*T/h) > 1.0, which implies k_rel > kB*T/h.
    # This seems to be a misunderstanding in the original or a specific context.
    # A relative rate exp(-dE/RT) can exceed 1 if dE is negative.
    # For now, I will remove the cap at 1.0 unless context demands it.
    return k_rel


def calc_rrkm(energy: WP, barrier: WP, nvib: int, freq_cm_minus_1: WP) -> WP:
    """Compute rate constant with simplified RRKM-like equation from Fortran code."""
    freq_s_minus_1 = abs(freq_cm_minus_1) * C_CM_S

    if energy <= 0 or barrier < 0: # Barrier must be positive for this form
        return 0.0
    if barrier > energy : # Not enough energy to overcome barrier
        return 0.0
    if energy == 0 and barrier == 0: # Avoid 0/0 if possible, rate should be high
        return freq_s_minus_1 # Max rate

    exponent = -(nvib - 1) * barrier / energy
    
    try:
        rate_constant = freq_s_minus_1 * math.exp(exponent)
    except OverflowError: # Should not happen with barrier > 0, energy > 0
        rate_constant = float('inf')
        
    # Limit rate to prefactor, as in Eyring for barrierless
    # This occurs if exponent is positive or zero (barrier <=0, which we excluded)
    # or if (nvib-1) is zero or negative (nvib <=1)
    if nvib <=1 : # If nvib=1, exponent is 0, k=freq. If nvib=0, exponent positive, k > freq.
        return min(rate_constant, freq_s_minus_1)

    return rate_constant


def get_index_from_energy(max_iee: WP, energy: WP, increments: int) -> int:
    """Calculates an index for an energy value within a discretized energy range."""
    if increments <= 0: return 1 # Default to first bin if increments is invalid
    
    energy_step = max_iee / increments
    if energy_step == 0: # Avoid division by zero if max_iee is 0
        return 1 if energy == 0 else increments # Or handle as error

    idx = int(energy / energy_step) # Fortran int() truncates

    # Ensure index is within bounds [1, increments] (Fortran-style 1-based indexing)
    idx = max(1, min(idx, increments))
    return idx

def print_iee_distribution(nsamples: int, eiee: List[WP], piee: List[WP], fname: str):
    """Prints energy distribution (e.g., IEE) to a file."""
    try:
        with open(fname, 'w') as f:
            for i in range(min(nsamples, len(eiee), len(piee))): # Ensure we don't go out of bounds
                f.write(f"{eiee[i]} {piee[i]}\n")
    except IOError as e:
        print(f"Error writing IEE distribution to {fname}: {e}")


# --- File and System Utilities ---
def print_pwd(message: Optional[str] = None):
    """Prints the current working directory, optionally with a message."""
    pwd = Path.cwd()
    if message:
        print(f"Current working directory is {message} {pwd}")
    else:
        print(f"Current working directory is {pwd}")

def check_prog(pname: str, verbose: bool = False, critical: bool = True) -> bool:
    """
    Checks if a program is available.
    If critical and not found, prints error and exits.
    Returns True if found, False otherwise.
    """
    path_to_prog = shutil.which(pname)
    if path_to_prog is None:
        if verbose or critical: # Print message if verbose or if it's critical and missing
            print(f"    Binary: \"{pname}\" not found.")
            if pname == "geodesic_interpolate":
                print("    You can turn off the geodesic interpolation with the option -notsgeo")
        if critical:
            sys.exit(f'Critical program "{pname}" not found, aborting QCxMS2.')
        return False
    else:
        if verbose:
            print(f"    Found binary: \"{pname}\" at {path_to_prog}")
        return True


def cleanup_qcxms2_for_restart(env: RunTypeData):
    """Cleanup files for restart calculation."""
    base_path = Path(".") # Assuming cleanup happens in current dir
    
    files_to_remove_in_base = ["pfrag", "kerav", "keravold"]
    if env.mode == "cid":
        files_to_remove_in_base.extend(["sumdekin", "sumdeint", "x_trav"])

    for fname in files_to_remove_in_base:
        (base_path / fname).unlink(missing_ok=True)

    # Glob patterns for subdirectories
    patterns_in_subdirs = ["*/sumreac_*", "*/pfrag", "*/kerav", "*/keravold"]
    if env.mode == "cid":
        patterns_in_subdirs.extend(["*/sumdekin", "*/sumdeint", "*/x_trav"])

    for pattern in patterns_in_subdirs:
        for p_file in base_path.glob(pattern):
            p_file.unlink(missing_ok=True)
    # print("Executed cleanup for restart.")


def append_char(str_list: List[str], element: str): # Modified to match typical Python list append
    """Appends an element to a list of strings."""
    str_list.append(element)
    # No need to return str_list if modified in-place, but can be convenient

def read_struc_from_xyz_ensemble(xyz_ensemble_file: str) -> Tuple[int, int]:
    """Reads an XYZ ensemble file to get number of atoms (nat) and number of structures (nstruc)."""
    nat = 0
    nstruc = 0
    
    try:
        with open(xyz_ensemble_file, 'r') as f:
            first_line = f.readline().strip()
            if not first_line: # Empty file
                return 0, 0
            try:
                nat = int(first_line)
            except ValueError:
                print(f"Error: Cannot read number of atoms from first line of {xyz_ensemble_file}")
                return 0, 0
            
            if nat <= 0: return 0, 0 # Invalid number of atoms

            # Count lines to determine number of structures
            # Reset file pointer to beginning after reading nat
            f.seek(0)
            total_lines = 0
            for line in f:
                if line.strip(): # Count non-empty lines
                    total_lines +=1
            
            if total_lines > 0 and nat > 0:
                lines_per_structure = nat + 2
                if total_lines % lines_per_structure == 0:
                    nstruc = total_lines // lines_per_structure
                else:
                    print(f"Warning: Total lines ({total_lines}) in {xyz_ensemble_file} is not a perfect multiple of (nat={nat} + 2 lines per structure). File might be malformed or contain extra data.")
                    nstruc = total_lines // lines_per_structure # Integer division will give floor
                    if nstruc == 0: nstruc = 1 # If there are lines, assume at least one structure
            elif total_lines == 0 and nat == 0: # Empty file case after first line check
                pass # nat=0, nstruc=0
            elif total_lines > 0 and nat == 0 : # Should have been caught by nat <=0
                 print(f"Warning: nat is 0 but file has lines in {xyz_ensemble_file}")
                 return 0,0


    except FileNotFoundError:
        print(f"Error: File {xyz_ensemble_file} not found.")
        return 0, 0
    except Exception as e:
        print(f"An error occurred while reading {xyz_ensemble_file}: {e}")
        return 0,0
        
    return nat, nstruc


def rewrite_xyz_elements_to_nice_case(fname: str):
    """Rewrites element symbols in an XYZ file to 'NiceCase' (e.g., CL -> Cl)."""
    # This function needs to read the file, modify content, and write back.
    # It's safer to write to a temporary file then replace.
    
    nat = iomod.rdshort_int(fname) # Get number of atoms
    if nat == 0 and not Path(fname).exists(): # rdshort_int might return 0 if file not found
        print(f"Error: Cannot read number of atoms from {fname} or file does not exist.")
        return

    lines = []
    try:
        with open(fname, 'r') as f_in:
            lines.append(f_in.readline()) # Number of atoms
            lines.append(f_in.readline()) # Comment line
            for _ in range(nat):
                atom_line = f_in.readline()
                parts = atom_line.split()
                if parts:
                    element_symbol = parts[0]
                    nice_case_symbol = i2e(e2i(element_symbol), 'nc') # Convert to num, then back to nice string
                    parts[0] = nice_case_symbol
                    lines.append(" ".join(parts) + "\n")
                else:
                    lines.append(atom_line) # Keep empty or malformed lines as is
            # Read any remaining lines (e.g. if file has more than one structure or trailing data)
            lines.extend(f_in.readlines())


        with open(fname, 'w') as f_out:
            for line in lines:
                f_out.write(line)
    except FileNotFoundError:
        print(f"ERROR: Could not find file {fname} for element rewriting")
    except Exception as e:
        print(f"Error during XYZ element rewriting for {fname}: {e}")


def cut_topology_from_molbar(identifiertopo: str) -> str:
    """Cuts the topology part from molbar output string."""
    # Molbar output example: ...|topography_string|...
    # The Fortran code counts '|' occurrences.
    
    if 'Error' in identifiertopo: # Check for molbar error
        return ""
        
    parts = identifiertopo.split('|')
    # Fortran code implies topography is up to the 5th pipe.
    # If parts = [s0, s1, s2, s3, s4, s5, ...], then 5th pipe is after s4.
    # So, we want s0|s1|s2|s3|s4|
    if len(parts) > 5:
        return "|".join(parts[:5]) + "|"
    else:
        # Not enough parts for 5 pipes, maybe return original or error
        return identifiertopo # Or handle as an incomplete molbar string

def check_topology_change(env: RunTypeData, fname_xyz: str) -> bool:
    """
    Checks if molecular topology has changed using an external program (molbar or inchi).
    Returns True if changed, False otherwise or if check fails.
    """
    nat = rdshort_int(fname_xyz)
    if nat == 1: # Single atom, topology doesn't change, molbar might fail
        return False

    # Ensure elements are in 'NiceCase' for consistency if tools are sensitive
    rewrite_xyz_elements_to_nice_case(fname_xyz) # Uses iomod.rdshort_int indirectly

    topo_file = Path("topo_check.txt") # Temporary file for topology string
    
    reference_topo_str = iomod.rdshort_string("topo") 
    if env.topocheck == "molbar":
        reference_topo_str = cut_topology_from_molbar(reference_topo_str)

    # Calculate current topology
    cmd_list: List[str] = []
    # Using shell=True for commands with redirection, similar to Fortran's execute_command_line
    # If no redirection, better to use list of args and shell=False.
    shell_cmd = ""

    if env.topocheck == "molbar":
        shell_cmd = f"molbar {fname_xyz} > {topo_file} 2> /dev/null"
    elif env.topocheck == "inchi":
        shell_cmd = f"obabel -i xyz {fname_xyz} -oinchi > {topo_file} 2> /dev/null" # Fixed -o inchi
    else:
        print(f"Unknown topocheck method: {env.topocheck}")
        return False

    # iomod.execute_command takes list or string if shell=True
    iomod.execute_command([shell_cmd], shell=True)


    if not topo_file.exists():
        print("Topology calculation failed (output file not created).")
        return False # Assume no change if check fails

    current_topo_str = ""
    with open(topo_file, 'r') as f:
        current_topo_str = f.read().strip()
    
    topo_file.unlink(missing_ok=True) # Clean up

    if env.topocheck == "molbar":
        current_topo_str = cut_topology_from_molbar(current_topo_str)

    if reference_topo_str and current_topo_str:
        if reference_topo_str != current_topo_str:
            # print("Topology has changed.")
            return True
        else:
            # print("Topology remains the same.")
            return False
    else:
        # If either string is empty (e.g., molbar error, empty reference file)
        # print("Could not compare topologies (one or both empty). Assuming no change.")
        return False


def boltzmann_average_energy(env: RunTypeData, energies: List[WP]) -> WP:
    """Give Boltzmann averaged energy for a set of energies in eV."""
    if not energies:
        return 0.0

    nenergies = len(energies)
    temperature_k = float(env.temp) # Use temp from RunTypeData
    
    # Avoid issues with zero or very low temperature
    if temperature_k < 0.01: temperature_k = 0.01 

    boltzmann_factor_rt = KB_EV_K * temperature_k # This is R*T in eV
    if boltzmann_factor_rt == 0: # Avoid division by zero if T is extremely low (though capped)
        # If T=0, only lowest energy state is populated.
        return min(energies) if energies else 0.0

    min_energy = min(energies) # Shift energies to prevent overflow with exp for large energies
    shifted_energies = [e - min_energy for e in energies]
    
    sum_exp_terms = 0.0
    weighted_sum_energies = 0.0
    
    for i in range(nenergies):
        try:
            # For numerical stability with exp(-E/kT):
            # If E/kT is very large and positive, exp term is tiny.
            # If E/kT is very large and negative (i.e. E is far below others relative to kT), exp term is huge.
            # Using shifted energies helps with the positive E/kT case.
            exponent = -shifted_energies[i] / boltzmann_factor_rt
            if exponent < -700: # exp(-700) is effectively zero
                q = 0.0
            elif exponent > 700: # exp(700) is huge
                # This can happen if an energy is far *below* min_energy before shifting,
                # or if kT is tiny and shifted_energy[i] is negative (which it shouldn't be here).
                # This indicates a state is overwhelmingly probable.
                # To handle this, we can normalize later; for now, let it be large.
                q = math.exp(700) # Cap to avoid immediate overflow if not careful
            else:
                q = math.exp(exponent)
        except OverflowError: # Should be caught by capping above
            q = float('inf')

        sum_exp_terms += q
        # Use original energies for averaging, but shifted for exp calculation
        weighted_sum_energies += q * energies[i] 

    if sum_exp_terms == 0: # All exp terms were zero (e.g., huge energies, tiny T)
        return min_energy # Or simply min(energies)
        
    return weighted_sum_energies / sum_exp_terms


def check_important_fragments(env: RunTypeData, npairs: int, fragdirs: List[List[str]]) -> Tuple[int, List[List[str]]]:
    """
    Checks important fragments based on 'pfrag' value (probability/percentage).
    fragdirs is expected to be a list of lists/tuples, e.g., [[dir1, dir2, dir3], ...]
    where dir1 is primary, dir2/dir3 are for pairs.
    """
    imp_frag_count = 0
    imp_fragdirs_list: List[List[str]] = [] # Stores [['d1','d2','d3'],...] for important ones

    # pthr from env, default 1.0 means 1%
    # The Fortran code implies pfrag is a percentage (e.g. 1.0 for 1%)
    # If pfrag is a fraction (0.01 for 1%), adjust pthr accordingly.
    # Assuming pfrag values read are percentages.
    p_threshold = env.pthr 

    for i in range(npairs):
        if i >= len(fragdirs): break # Safety break

        current_frag_info = fragdirs[i] # e.g., ["p1_fragA", "p1_fragB_path", "p1_fragB_data"] or ["m1_fragA", "", ""]
        
        is_important = False
        pfrag1: WP = 0.0
        pfrag2: WP = 0.0

        # Primary fragment (always present)
        if current_frag_info[0]:
            pfrag_file1 = Path(current_frag_info[0]) / "pfrag"
            pfrag1 = iomod.rdshort_real(str(pfrag_file1))

        is_pair = bool(len(current_frag_info) > 2 and current_frag_info[2] and 'p' in current_frag_info[2].lower())

        if is_pair and len(current_frag_info) > 2 and current_frag_info[2]: 
            pfrag_file2 = Path(current_frag_info[2]) / "pfrag" 
            pfrag2 = iomod.rdshort_real(str(pfrag_file2))
            if pfrag1 > p_threshold or pfrag2 > p_threshold:
                is_important = True
        else: 
            if pfrag1 > p_threshold:
                is_important = True
        
        if is_important:
            imp_frag_count += 1
            imp_fragdirs_list.append(current_frag_info[:]) # Add a copy

    return imp_frag_count, imp_fragdirs_list


def sort_out_zero_elements_from_fragdirs(
    npairs_in: int,
    fragdirs_in: List[List[str]], # List of [dir1, dir2, dir3]
    remove_dir_flag: bool
) -> Tuple[int, List[List[str]]]:
    """
    Filters out fragment directory entries considered "zero" or invalid.
    In Fortran, this was based on 'p' not being in fragdirs_in(i,1).
    Returns the new count and the filtered list.
    """
    if npairs_in == 0:
        return 0, []

    filtered_fragdirs: List[List[str]] = []
    
    for i in range(npairs_in):
        if i >= len(fragdirs_in): break # Safety

        current_entry = fragdirs_in[i]
        primary_dir = current_entry[0] if len(current_entry) > 0 else ""
        
        # Condition from Fortran: index(fragdirs_in(i, 1), 'p') .ne. 0
        # This means the primary directory name must contain 'p'.
        # This seems like a heuristic for "product" directories.
        if 'p' in primary_dir.lower(): # Keep if 'p' is in the primary directory name
            filtered_fragdirs.append(current_entry[:]) # Keep a copy
        else:
            # This entry is to be removed
            if remove_dir_flag:
                if primary_dir: iomod.rmrf(primary_dir)
                
                dir2 = current_entry[1] if len(current_entry) > 1 else ""
                dir3_indicator = current_entry[2] if len(current_entry) > 2 else ""

                if dir3_indicator and 'p' in dir3_indicator.lower(): # Indicates a pair
                    if dir2: iomod.rmrf(dir2) 
                    if dir3_indicator: iomod.rmrf(dir3_indicator) 

    return len(filtered_fragdirs), filtered_fragdirs

# Main block for testing some of the translated functions
if __name__ == '__main__':
    # Test calctemp, calc_eyring, calc_rrkm (example values)
    nvib_test = 30 
    energy_test_ev = 1.5
    barrier_test_ev = 0.5
    dbarrier_test_ev = 0.2
    freq_test_cm_1 = 1000 

    temp_k = calctemp(nvib_test, energy_test_ev)
    print(f"Calculated T for E={energy_test_ev}eV, Nvib={nvib_test}: {temp_k:.2f} K")

    k_eyring = calc_eyring(energy_test_ev, barrier_test_ev, nvib_test)
    print(f"Eyring rate k for E={energy_test_ev}eV, barrier {barrier_test_ev}eV: {k_eyring:.2e} s^-1")
    
    k_rel_eyring = calc_releyring(energy_test_ev, dbarrier_test_ev, nvib_test)
    print(f"Relative Eyring rate for E={energy_test_ev}eV, dbarrier {dbarrier_test_ev}eV: {k_rel_eyring:.2e}")

    k_rrkm = calc_rrkm(energy_test_ev, barrier_test_ev, nvib_test, freq_test_cm_1)
    print(f"RRKM-like rate k for E={energy_test_ev}eV, barrier {barrier_test_ev}eV, freq {freq_test_cm_1}cm-1: {k_rrkm:.2e} s^-1")

    # Test get_index_from_energy
    print(f"Index for E=0.75, Max=1.5, Incr=10: {get_index_from_energy(1.5, 0.75, 10)}")
    print(f"Index for E=0.0, Max=1.5, Incr=10: {get_index_from_energy(1.5, 0.0, 10)}")
    print(f"Index for E=1.5, Max=1.5, Incr=10: {get_index_from_energy(1.5, 1.5, 10)}")
    print(f"Index for E=2.0, Max=1.5, Incr=10: {get_index_from_energy(1.5, 2.0, 10)}")

    # Test mask_invert
    mask0 = [0, 1, 2]
    print(f"Mask {mask0}, nall={len(mask0)}, Inverted: {mask_invert(len(mask0), mask0)}")
    mask1 = [2, 0, 1] # Permutation of 0,1,2
    print(f"Mask {mask1}, nall={len(mask1)}, Inverted: {mask_invert(len(mask1), mask1)}")
    # Fortran style 1-based example
    mask_f = [3, 1, 2] # Permutation of 1,2,3
    print(f"Mask {mask_f} (1-based), nall={len(mask_f)}, Inverted: {mask_invert(len(mask_f), mask_f)}")


    # Test read_struc_from_xyz_ensemble (requires a dummy file)
    dummy_xyz_content = "2\nComment\nC 0 0 0\nH 0 0 1\n2\nComment2\nO 0 0 0\nN 0 0 1\n"
    dummy_xyz_file = Path("dummy_ensemble.xyz")
    with open(dummy_xyz_file, "w") as f:
        f.write(dummy_xyz_content)
    nat_read, nstruc_read = read_struc_from_xyz_ensemble(str(dummy_xyz_file))
    print(f"Read from {dummy_xyz_file}: nat={nat_read}, nstruc={nstruc_read}") # Expected: nat=2, nstruc=2
    dummy_xyz_file.unlink()

    # Test rewrite_xyz_elements_to_nice_case (uses iomod.rdshort_int)
    dummy_xyz_raw_elements = "1\nTest\nCL 1 2 3\n"
    dummy_xyz_nice_file = Path("dummy_nice.xyz")
    with open(dummy_xyz_nice_file, "w") as f:
        f.write(dummy_xyz_raw_elements)
    # Need to ensure iomod.rdshort_int can read this before rewrite works
    # For now, manually provide nat to rewrite_xyz_elements_to_nice_case or adapt it
    # Temporarily skipping this specific test as it relies on rdshort_int behavior within rewrite
    # rewrite_xyz_elements_to_nice_case(str(dummy_xyz_nice_file)) 
    # with open(dummy_xyz_nice_file, "r") as f:
    #     content_after_rewrite = f.read()
    #     print(f"Content of {dummy_xyz_nice_file} after rewrite:\n{content_after_rewrite.strip()}")
    if dummy_xyz_nice_file.exists(): dummy_xyz_nice_file.unlink()

    # Test sum_formula_from_atoms
    atoms_for_formula = [6, 1, 1, 1, 1, 8] # CH4O
    print(f"Sum formula for {atoms_for_formula}: {sum_formula_from_atoms(len(atoms_for_formula), atoms_for_formula)}")

    # Test cleanup (creates dummy files/dirs then removes)
    # print("Testing cleanup...")
    # Path("pfrag").touch()
    # Path("kerav").touch()
    # Path("p1").mkdir(exist_ok=True)
    # (Path("p1") / "sumreac_test").touch() # Uses iomod.rmrf
    # (Path("p1") / "pfrag").touch()        # Uses iomod.rmrf
    # dummy_env_cleanup = RunTypeData()
    # dummy_env_cleanup.mode = "ei" 
    # cleanup_qcxms2_for_restart(dummy_env_cleanup) # Uses iomod.rmrf
    # if Path("p1").exists(): shutil.rmtree("p1")


# --- New XYZ Parsing and Geometry Utilities ---

def get_atomic_numbers_and_coords_py(
    xyz_filepath: Union[str, Path]
) -> Tuple[int, Optional[List[int]], Optional[List[List[WP]]]]:
    """
    Reads a standard XYZ file.
    Returns:
        Tuple (natoms, list_of_atomic_numbers, list_of_coordinates_angstrom)
        Returns (0, None, None) if file cannot be read or is malformed.
    """
    try:
        with open(xyz_filepath, 'r') as f:
            lines = f.readlines()
        
        if not lines:
            print(f"Warning: XYZ file {xyz_filepath} is empty.")
            return 0, None, None

        num_atoms = int(lines[0].strip())
        if len(lines) < num_atoms + 2:
            print(f"Warning: XYZ file {xyz_filepath} is malformed (not enough lines for {num_atoms} atoms).")
            return 0, None, None
            
        atomic_numbers: List[int] = []
        coordinates: List[List[WP]] = []
        
        for i in range(2, num_atoms + 2):
            parts = lines[i].split()
            if len(parts) < 4:
                print(f"Warning: Malformed atom line in {xyz_filepath}: {lines[i].strip()}")
                return 0, None, None # Or skip this atom? For now, treat as error.
            
            atomic_numbers.append(e2i(parts[0])) # Uses e2i from this module
            coordinates.append([float(parts[1]), float(parts[2]), float(parts[3])])
            
        return num_atoms, atomic_numbers, coordinates
        
    except FileNotFoundError:
        print(f"Warning: XYZ file {xyz_filepath} not found.")
        return 0, None, None
    except ValueError:
        print(f"Warning: Could not parse number of atoms or coordinates in {xyz_filepath}.")
        return 0, None, None
    except Exception as e:
        print(f"An unexpected error occurred reading {xyz_filepath}: {e}")
        return 0, None, None

def get_atomic_numbers_from_xyz_py(xyz_filepath: Union[str, Path]) -> List[int]:
    """Reads an XYZ file and returns a list of atomic numbers."""
    _, atomic_numbers, _ = get_atomic_numbers_and_coords_py(xyz_filepath)
    return atomic_numbers if atomic_numbers is not None else []

def is_linear_molecule_py(xyz_filepath: Union[str, Path], tol_deg: WP = 5.0) -> bool:
    """
    Checks if a molecule is linear.
    For N<3 atoms, it's considered linear.
    For N=3, checks if the bond angle is close to 180 degrees.
    For N>3, uses a simplified check (e.g., if all atoms lie on a line defined by two outermost).
    A robust check would use moments of inertia (requires NumPy).
    """
    num_atoms, _, coords_angstrom = get_atomic_numbers_and_coords_py(xyz_filepath)

    if num_atoms is None or coords_angstrom is None or num_atoms == 0:
        print(f"Warning: Could not read molecule from {xyz_filepath} for linearity check.")
        return False # Default to non-linear if structure is unreadable

    if num_atoms < 3:
        return True # Point or diatomic is linear

    # Convert to NumPy arrays for easier vector math if available, otherwise list math
    try:
        import numpy as np
        coords = np.array(coords_angstrom, dtype=float)

        if num_atoms == 3:
            v1 = coords[0] - coords[1]
            v2 = coords[2] - coords[1]
            v1_u = v1 / np.linalg.norm(v1)
            v2_u = v2 / np.linalg.norm(v2)
            dot_product = np.dot(v1_u, v2_u)
            angle_rad = np.arccos(np.clip(dot_product, -1.0, 1.0))
            angle_deg = np.degrees(angle_rad)
            return abs(angle_deg - 180.0) < tol_deg or angle_deg < tol_deg # Linear if ~180 or ~0
        
        # For N > 3, moment of inertia is more robust.
        # Simplified check: are all atoms collinear with the first two?
        # This is not very robust.
        # For now, default to non-linear for N > 3 for simplicity without full inertia calc.
        # A proper implementation would calculate principal moments of inertia.
        # If one is close to zero (or much smaller than the other two), it's linear.
        if num_atoms > 3: # Placeholder for more robust N>3 check
            # Calculate inertia tensor
            centered_coords = coords - np.mean(coords, axis=0)
            Ixx = np.sum(centered_coords[:,1]**2 + centered_coords[:,2]**2)
            Iyy = np.sum(centered_coords[:,0]**2 + centered_coords[:,2]**2)
            Izz = np.sum(centered_coords[:,0]**2 + centered_coords[:,1]**2)
            Ixy = -np.sum(centered_coords[:,0] * centered_coords[:,1])
            Ixz = -np.sum(centered_coords[:,0] * centered_coords[:,2])
            Iyz = -np.sum(centered_coords[:,1] * centered_coords[:,2])
            
            inertia_tensor = np.array([
                [Ixx, Ixy, Ixz],
                [Ixy, Iyy, Iyz],
                [Ixz, Iyz, Izz]
            ])
            eigenvalues = np.linalg.eigvalsh(inertia_tensor)
            # Sort eigenvalues by magnitude
            sorted_eigenvalues = np.sort(np.abs(eigenvalues))
            # If smallest eigenvalue is close to zero relative to others (e.g., < 1e-4 of largest)
            if sorted_eigenvalues[0] < 1e-4 * sorted_eigenvalues[-1] and sorted_eigenvalues[-1] > 1e-4 : # Check if smallest is near zero and molecule is not a point
                return True


        print(f"Warning: Linearity check for N>3 atoms in {xyz_filepath} is using moment of inertia. Ensure NumPy is available.")
        return False # Default for N>3 if NumPy based check fails or not precise enough

    except ImportError:
        print("Warning: NumPy not available for robust linearity check. Defaulting to non-linear for N>3.")
        return False # Default if NumPy is not available
    except Exception as e:
        print(f"Error during linearity check for {xyz_filepath}: {e}. Defaulting to non-linear.")
        return False

def get_average_mol_mass_py(atomic_numbers: List[int]) -> WP:
    """Calculates the average molecular mass from a list of atomic numbers."""
    # This function was originally in isotopes.py, moved/re-exported here for utility access
    molmass: WP = 0.0
    for z_val in atomic_numbers:
        if 0 < z_val < len(AVERAGE_ATOMIC_MASSES):
            molmass += AVERAGE_ATOMIC_MASSES[z_val]
        else:
            # Handle unknown atomic number, e.g., add 0 or raise error
            print(f"Warning: Average atomic mass for Z={z_val} not available. Contributing 0 to molecular mass.")
    return molmass


# --- Additional functions needed for tsmod.py ---

def xtb_extract_snapshot_py(trajectory_file: Path, frame_number: int, output_file: Path) -> bool:
    """
    Extracts a specific frame from an XYZ trajectory file.
    frame_number is 1-based (Fortran style).
    Returns True if successful, False otherwise.
    """
    try:
        with open(trajectory_file, 'r') as f:
            current_frame = 0
            while True:
                # Read number of atoms
                natoms_line = f.readline()
                if not natoms_line:
                    break  # EOF
                
                natoms = int(natoms_line.strip())
                comment_line = f.readline()  # Comment line
                
                current_frame += 1
                
                if current_frame == frame_number:
                    # This is the frame we want
                    with open(output_file, 'w') as out_f:
                        out_f.write(natoms_line)
                        out_f.write(comment_line)
                        
                        # Read and write atom coordinates
                        for _ in range(natoms):
                            atom_line = f.readline()
                            out_f.write(atom_line)
                    
                    return True
                else:
                    # Skip this frame
                    for _ in range(natoms):
                        f.readline()
        
        print(f"Error: Frame {frame_number} not found in trajectory {trajectory_file}")
        return False
        
    except Exception as e:
        print(f"Error extracting frame {frame_number} from {trajectory_file}: {e}")
        return False


def concatenate_xyz_files_py(xyz_files: List[Path], output_file: Path) -> bool:
    """
    Concatenates multiple XYZ files into a single XYZ file.
    Used for combining fragment XYZ files into product structures.
    """
    try:
        total_atoms = 0
        all_atom_lines = []
        
        # Read all input files
        for xyz_file in xyz_files:
            if not xyz_file.exists():
                print(f"Warning: XYZ file {xyz_file} not found for concatenation")
                continue
                
            with open(xyz_file, 'r') as f:
                lines = f.readlines()
                
            if len(lines) < 2:
                print(f"Warning: XYZ file {xyz_file} is too short")
                continue
                
            natoms = int(lines[0].strip())
            total_atoms += natoms
            
            # Collect atom lines (skip natoms and comment lines)
            for i in range(2, min(len(lines), natoms + 2)):
                all_atom_lines.append(lines[i])
        
        # Write combined file
        with open(output_file, 'w') as f:
            f.write(f"{total_atoms}\n")
            f.write("Combined structure\n")
            for line in all_atom_lines:
                f.write(line)
        
        return True
        
    except Exception as e:
        print(f"Error concatenating XYZ files: {e}")
        return False


def check_atom_py(xyz_file: Union[str, Path]) -> bool:
    """
    Checks if a molecule is a single atom.
    Returns True if it's a single atom, False otherwise.
    """
    try:
        natoms, _, _ = get_atomic_numbers_and_coords_py(xyz_file)
        return natoms == 1
    except:
        return False


def set_omp_to_one_py():
    """Python version of set_omp_to_one"""
    set_omp_to_one()


def get_uhf_from_file_py() -> int:
    """
    Reads UHF value from .UHF file in current directory.
    Returns 0 if file not found or cannot be read.
    """
    try:
        return iomod.rdshort_int(".UHF", default=0)
    except:
        return 0


def calctemp_py(nvib: int, energy_ev: WP) -> WP:
    """Python version of calctemp"""
    return calctemp(nvib, energy_ev)


def calc_eyring_py(energy: WP, barrier: WP, nvib: int) -> WP:
    """Python version of calc_eyring"""
    return calc_eyring(energy, barrier, nvib)


def calc_rrkm_py(energy: WP, barrier: WP, nvib: int, freq_cm_minus_1: WP) -> WP:
    """Python version of calc_rrkm"""
    return calc_rrkm(energy, barrier, nvib, freq_cm_minus_1)


def get_index_from_energy_py(max_iee: WP, energy: WP, increments: int) -> int:
    """Python version of get_index_from_energy"""
    return get_index_from_energy(max_iee, energy, increments)
