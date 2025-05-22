import os
import re
import shutil
from pathlib import Path
from typing import Tuple, Optional, List, Dict, Union

try:
    from .data import RunTypeData, WP
    from . import iomod
    from . import utility
    # Placeholder for structools if/when translated
    # from . import structools 
except ImportError:
    # Fallbacks for standalone/testing
    # Ensure these mock/dummy modules exist in your test environment
    print("Attempting to import dummy/mock modules for qmmod.py standalone run.")
    from data import RunTypeData, WP # type: ignore
    import iomod_mock as iomod # type: ignore
    import utility_mock as utility # type: ignore
    class StructoolsMock: # Mock structools for standalone testing
        @staticmethod
        def get_active_cnstring(env, filename_xyz): 
            print(f"Warning: Using MOCK structools.get_active_cnstring for {filename_xyz}")
            return "1 2 3" # Dummy active atoms string
    structools = StructoolsMock()


# --- Constants ---
QMDATA_FILE = "qmdata" # Name of the file storing QM calculation results

# --- Helper Functions (already defined in previous turn, ensure they are complete) ---

def _get_chrg_uhf(env: RunTypeData, fname_xyz: str, 
                  chrg_in: Optional[int] = None, 
                  uhf_in: Optional[int] = None) -> Tuple[int, int]:
    """
    Determines charge and UHF value.
    Reads from .CHRG/.UHF files or input, calculates spin if needed.
    Writes .CHRG and .UHF files.
    Assumes fname_xyz is relative to current working directory if used for spin calculation.
    """
    chrg: int
    uhf: int

    current_dir = Path(".") # Operations are in the CWD

    if chrg_in is not None:
        chrg = chrg_in
        iomod.wrshort_int(current_dir / ".CHRG", chrg)
        iomod.remove_file(current_dir / ".UHF") 
    else:
        chrg_file = current_dir / ".CHRG"
        if chrg_file.exists():
            chrg = iomod.rdshort_int(chrg_file)
        else:
            chrg = env.chrg
            iomod.wrshort_int(chrg_file, chrg)
            iomod.remove_file(current_dir / ".UHF")

    if uhf_in is not None:
        uhf = uhf_in
    else:
        uhf_file = current_dir / ".UHF"
        if uhf_file.exists():
            uhf = iomod.rdshort_int(uhf_file)
        else:
            # This part needs a robust way to get nat and iat from fname_xyz
            # For now, using a placeholder utility.get_nat_iat_from_xyz
            try:
                # Assuming utility.py will have a function to read XYZ for nat/iat
                # nat_val, iat_val = utility.read_xyz_for_elements(Path(fname_xyz))
                # spin_multiplicity = utility.get_spin(nat_val, iat_val, chrg)
                # For now, using a simplified placeholder as in the original conversion attempt:
                spin_multiplicity = utility.print_spin_py(env, fname_xyz, chrg) # Assumes this utility exists and is adapted

            except Exception as e:
                print(f"Warning: Could not determine nat/iat for spin calculation from {fname_xyz}: {e}. Using default spin.")
                spin_multiplicity = 1 if chrg == 0 else 2 # Fallback default

            uhf = spin_multiplicity - 1
            if uhf == -2: uhf = 0 
    
    iomod.wrshort_int(current_dir / ".UHF", uhf)
    return chrg, uhf


def _construct_qmdata_query(level: str, job: str, chrg: int, uhf: int) -> str:
    return f"{level.strip()} {job.strip()} {chrg} {uhf}"

def _check_qmdata_cache(query: str) -> Tuple[bool, Optional[WP]]:
    qmdata_path = Path(QMDATA_FILE)
    if qmdata_path.exists():
        found, value = iomod.grepval(qmdata_path, query)
        if found and abs(value) > 1e-5:
            return True, value
    return False, None

def _write_to_qmdata(query: str, value: WP, filepath: Union[str, Path] = QMDATA_FILE):
    entry = f"{query}   {value:.10f}" 
    
    path_obj = Path(filepath)
    # Prepend to make sure the latest value is found first by grepval (if it reads top-down)
    content_to_write = entry + "\n"
    if path_obj.exists():
        existing_content = path_obj.read_text()
        content_to_write += existing_content
    
    try:
        path_obj.write_text(content_to_write)
    except IOError as e:
        print(f"Error writing to {filepath}: {e}")


def _remove_entries_from_qmdata(queries_to_remove: List[str], filepath: Union[str, Path] = QMDATA_FILE):
    """Removes lines starting with any of the specified queries from qmdata file."""
    path_obj = Path(filepath)
    if not path_obj.exists():
        return

    lines_to_keep = []
    try:
        with open(path_obj, 'r') as f:
            for line in f:
                stripped_line = line.strip()
                if not any(stripped_line.startswith(query) for query in queries_to_remove):
                    lines_to_keep.append(line)
        
        with open(path_obj, 'w') as f:
            f.writelines(lines_to_keep)
            
    except IOError as e:
        print(f"Error updating {filepath}: {e}")


# --- QM Program Specific Preparers (from previous turn, ensure complete) ---

def _prep_xtb_py(env: RunTypeData, fname_xyz: str, level: str, job: str, 
                 chrg: int, uhf: int, restart: bool
                 ) -> Tuple[List[str], str, str, str]:
    xtb_cmd_parts = ["xtb", "xtbin.xyz"] 
    level_arg = f"--{level.replace('spinpol', '')}"
    xtb_cmd_parts.append(level_arg)
    if "spinpol" in level:
        xtb_cmd_parts.extend(["--tblite", "--spinpol"])
    
    xtb_cmd_parts.extend([f"--{job}", f"--chrg", str(chrg), f"--uhf", str(uhf)])
    
    etemp = 300
    if restart: etemp = 5000
    
    xtb_inp_path = Path("xtb.inp")
    if job == 'bhess' and not env.notemp:
        xtb_inp_path.write_text(f"$thermo\n  temp={env.temp}\n  sthr={env.sthr}\n  ithr={env.ithr}\nend\n")
        xtb_cmd_parts.extend(["--input", str(xtb_inp_path)]) # XTB needs explicit input file for $thermo
    
    xtb_cmd_parts.extend(["--etemp", str(etemp)])
    if env.solv: xtb_cmd_parts.extend(["--alpb", "water"])
        
    iomod.copy_file(fname_xyz, "xtbin.xyz")
    output_file = "xtb.out"
    jobcall_str_list = xtb_cmd_parts # Will be joined by execute_command if shell=True

    success_pattern = ""
    cleanup_parts = ["xtbrestart", "xtb.out", "wbo", "charges", "xtbin.xyz", str(xtb_inp_path)]
    if job == 'sp':
        success_pattern = "| TOTAL ENERGY"
    elif job == 'opt':
        success_pattern = "| TOTAL ENERGY"
        # jobcall_str += " && cp xtbopt.xyz opt.xyz" # Handled after successful run
        cleanup_parts.append("xtbopt.log")
    elif job in ['hess', 'bhess', 'ohess']:
        success_pattern = " G(RRHO) contrib." if not env.notemp else "zero point energy"
        cleanup_parts.extend(["xtbopt.log", "xtbhess.xyz", "hessian"])
        if job == 'ohess': cleanup_parts.append("xtbopt.xyz") # opt.xyz is final name

    if level == "gff": cleanup_parts.extend(["gfnff_charges", "gfnff_topo"])
    cleanup_parts.extend(["xtbtopo.mol", ".xtboptok"])
    cleanup_command = "rm -f " + " ".join(p for p in cleanup_parts if Path(p).exists() or "*" in p) + " >/dev/null 2>&1"
    
    return jobcall_str_list, output_file, success_pattern, cleanup_command


def _get_orca_level_keyword(level: str) -> Tuple[str, int, bool]:
    default_etemp = 5000 
    needs_fermi_for_xtb = False # For XTBn methods, fermi is effectively on or not used.
    
    # Simplified mapping, expand as needed from previous version
    if level == 'gfn2': return 'XTB2', 300, True
    if level == 'gfn1': return 'XTB1', 300, True
    if level == 'gfn2spinpol': return 'XTB2', 300, True # XTB2 + spinpol keyword
    if level == 'r2scan3c': return 'R2SCAN-3c', 5000, False
    if level == 'wb97x3c': return 'wB97X-3c', 15000, False
    # ... add all other mappings from previous implementation
    return level, default_etemp, needs_fermi_for_xtb


def _prep_orca_py(env: RunTypeData, fname_xyz: str, level: str, job: str, 
                  chrg: int, uhf: int, restart: bool
                  ) -> Tuple[List[str], str, str, str]:
    orca_inp_lines = []
    level_keyword, default_etemp, xtb_level_no_explicit_fermi = _get_orca_level_keyword(level)
    
    etemp = default_etemp
    apply_fermi_smearing = env.fermi
    if restart:
        orca_inp_lines.append("! SlowConv")
        if not xtb_level_no_explicit_fermi : etemp = 5000 
        apply_fermi_smearing = True 

    if level != "ccsdt": orca_inp_lines.append("! DEFGRID2") # Default unless CCSD(T)
    orca_inp_lines.append(f"! {level_keyword}")

    if env.solv:
        solv_keyword = "! ALPB(water)" if "gfn" in level or "xtb" in level.lower() else "! CPCM(WATER)"
        orca_inp_lines.append(solv_keyword)

    job_keyword_map = {'sp': 'SP', 'opt': 'LOOSEOPT', 'hess': 'FREQ', 
                       'ohess': 'OPT FREQ', 'optts': 'OptTS LOOSEOPT', 'tddft': 'SP'}
    orca_job_keyword = job_keyword_map.get(job, job.upper())
    if job == 'optts' and env.tsoptgmf: orca_job_keyword = 'OptTS(GMF) LOOSEOPT'
    if job == 'hess' and apply_fermi_smearing and not xtb_level_no_explicit_fermi:
        apply_fermi_smearing = False
        print("Warning: Fermi smearing disabled for ORCA Hessian calculation with DFT.")
        
    orca_inp_lines.append(f"! {orca_job_keyword}")

    if level == 'gfn2spinpol': # XTB2 is level_keyword, need to add spinpol via %xtb
        orca_inp_lines.extend(['%xtb', 'XTBINPUTSTRING2 "--tblite --spinpol"', 'end'])
    
    if job == 'optts':
        nmode = iomod.rdshort_int('nmode', default=1) # Default to mode 1 if file missing
        orca_inp_lines.append('%geom')
        orca_inp_lines.append(f'  TS_Mode {{ M {nmode -1} }}')
        orca_inp_lines.append('  inhess read')
        orca_inp_lines.append('  inhessname "orca.hess"')
        orca_inp_lines.append('  maxiter 100')
        if env.tsoptact:
            active_atoms_str = structools.get_active_cnstring(env, fname_xyz)
            orca_inp_lines.append(f'  TS_Active_Atoms {{ {active_atoms_str} }}')
            orca_inp_lines.append(f'  TS_Active_Atoms_Factor 1.5')
        orca_inp_lines.append('end') # end %geom
        if env.geolevel == 'gxtb':
            xtb_driver_opts = f"--driver 'gxtb -c orca.xtbdriver.xyz -symthr 0.0{'' if not apply_fermi_smearing else ' -tel 15000'}'"
            orca_inp_lines.extend(['%xtb', f'XTBINPUTSTRING2 "{xtb_driver_opts}"', 'end'])
            Path('.GRAD').touch(exist_ok=True)

    orca_inp_lines.append(f"%pal nprocs {env.threads} end")
    orca_inp_lines.append(f"%maxcore {env.threads * 1000}") # Adjusted to use env.threads

    if apply_fermi_smearing and not xtb_level_no_explicit_fermi:
        orca_inp_lines.append(f"%scf SmearTemp {etemp} end")

    multiplicity = uhf + 1
    orca_inp_lines.append(f"*xyzfile {chrg} {multiplicity} {fname_xyz}")

    if job == 'tddft' and env.exstates > 0:
        orca_inp_lines.append(f"%TDDFT NROOTS {env.exstates} END")
    
    if level not in ["gfn2", "gfn1", "gxtb", "ccsdt"] and job == "optts":
         if Path(fname_xyz).name == "ts.xyz": 
            orca_inp_lines.append("! UKS")

    Path("orca.inp").write_text("\n".join(orca_inp_lines) + "\n")
    
    orca_executable = shutil.which("orca")
    if not orca_executable:
        return ["echo", "ORCA_NOT_FOUND"], "orca.out", "FATAL_ERROR", "echo 'cleanup'"
    job_command_list = [orca_executable, "orca.inp"]
    
    output_file = "orca.out"
    cleanup_parts = ["orca.gbw", "orca.prop", "orca.inp", "orca.bibtex", 
                     "orca.property.txt", "orca.engrad", "orca.xtbrestart", "orca_main_removeme.xyz", "orca_opt_removeme.xyz"] # Added common temp files
    success_pattern = ""

    if job == 'sp': success_pattern = 'FINAL SINGLE POINT ENERGY'
    elif job == 'opt': success_pattern = 'THE OPTIMIZATION HAS CONVERGED'; cleanup_parts.extend(["orca.opt", "orca.gu.tmp"])
    elif job == 'hess':
        success_pattern = "Total VIBRATIONAL ENTHALPY" # More specific
        if env.notemp: success_pattern = "Zero point energy" 
        else: success_pattern = "Total THERMAL CORRECTION TO ENTHALPY" # From ORCA output
        cleanup_parts.extend(["orca.opt", "orca.gu.tmp", "orca.hess"])
    elif job == 'ohess': success_pattern = 'THE OPTIMIZATION HAS CONVERGED'; cleanup_parts.extend(["orca.opt", "orca.gu.tmp", "orca.hess"])
    elif job == 'optts': success_pattern = 'THE OPTIMIZATION HAS CONVERGED'; cleanup_parts.extend(["orca.opt", "orca.gu.tmp", "orca.xyz.opt"])
    elif job == 'tddft': success_pattern = 'FINAL SINGLE POINT ENERGY'
        
    cleanup_command = "rm -f " + " ".join(p for p in cleanup_parts if Path(p).exists() or "*" in p) + " >/dev/null 2>&1"
    return job_command_list, output_file, success_pattern, cleanup_command


def _prep_gxtb_py(env: RunTypeData, fname_xyz: str, level: str, job: str, 
                 chrg: int, uhf: int, restart: bool
                 ) -> Tuple[List[str], str, str, str]:
    iomod.copy_file(fname_xyz, "xtbin.xyz")
    iomod.remove_file("gxtbrestart"); iomod.remove_file("ceh.charges")

    job_command_list: List[str] = []
    output_file = "gxtb.out"
    success_pattern = ""
    cleanup_parts = ["gxtbrestart", "energy", "ceh.charges", "gxtb.out", "xtbin.xyz", "errorfile"] # keep errorfile if printlevel >=2

    gxtb_base_cmd = ["gxtb", "-c", "xtbin.xyz"]
    if env.fermi or restart: gxtb_base_cmd.extend(["-tel", "15000"])

    if job == 'sp':
        iomod.remove_file(".GRAD"); iomod.remove_file(".HESS")
        job_command_list = gxtb_base_cmd
        success_pattern = " GXTBSPENERGYFILE " 
    elif job == 'opt':
        Path(".GRAD").touch(exist_ok=True); iomod.remove_file("coord")
        driver_cmd_str = f"gxtb -c xtbdriver.xyz -symthr 0.0"
        if env.fermi or restart: driver_cmd_str += " -tel 15000"
        job_command_list = ["xtb", "xtbin.xyz", "--opt", "--driver", driver_cmd_str]
        output_file = "opt.out"
        success_pattern = "| TOTAL ENERGY"
        cleanup_parts.extend(["xtbdriver.xyz", "opt.out", ".GRAD", "xtbopt.xyz"])
    elif job == 'hess':
        Path(".HESS").touch(exist_ok=True)
        job_command_list = gxtb_base_cmd
        success_pattern = " GXTBSPENERGYFILE "
        cleanup_parts.extend([".GRAD", ".HESS"])
    else:
        return ["echo", f"Job_type_{job}_not_supported_for_gXTB"], "error.out", "FATAL_ERROR", "echo 'cleanup'"
    
    cleanup_command = "rm -f " + " ".join(p for p in cleanup_parts if Path(p).exists() or "*" in p)
    if env.printlevel >= 2 and "errorfile" in cleanup_parts : cleanup_parts.remove("errorfile")
    cleanup_command += " >/dev/null 2>&1"
        
    return job_command_list, output_file, success_pattern, cleanup_command

# --- Main QM Interface Functions ---

PrepResult = Tuple[Optional[List[str]], Optional[str], Optional[str], Optional[str], bool, Optional[WP]]

def prepqm_py(env: RunTypeData, fname_xyz: str, level: str, job: str, 
              chrg_in: Optional[int] = None, uhf_in: Optional[int] = None, 
              restart: bool = False) -> PrepResult:
    """
    Prepares QM calculation: determines chrg/uhf, checks cache, calls specific preparer.
    Returns: (job_command_list, output_file, success_pattern, cleanup_command, already_there_bool, cached_energy_float)
    """
    fname_path = Path(fname_xyz)
    if not fname_path.exists():
        print(f"Error: Input file {fname_xyz} not found for prepqm.")
        return None, None, None, None, False, None

    chrg, uhf = _get_chrg_uhf(env, fname_xyz, chrg_in, uhf_in)
    
    query = _construct_qmdata_query(level, job, chrg, uhf)
    is_cached, cached_value = _check_qmdata_cache(query)
    if is_cached and cached_value is not None:
        return None, None, None, None, True, cached_value

    job_command_list: Optional[List[str]] = None
    output_file: Optional[str] = None
    success_pattern: Optional[str] = None
    cleanup_command: Optional[str] = None

    # Dispatch to specific preparers
    if job == 'bhess': # Always XTB GFN2 for bhess
        prep_level = 'gfn2spinpol' if env.geolevel == "gfn2spinpol" else 'gfn2'
        job_command_list, output_file, success_pattern, cleanup_command = _prep_xtb_py(env, fname_xyz, prep_level, job, chrg, uhf, restart)
    elif job == 'optts' and level in ['gfn2', 'gfn2spinpol', 'gfn1', 'gxtb']:
        job_command_list, output_file, success_pattern, cleanup_command = _prep_orca_py(env, fname_xyz, level, job, chrg, uhf, restart)
    elif job == 'tddft':
        job_command_list, output_file, success_pattern, cleanup_command = _prep_orca_py(env, fname_xyz, level, job, chrg, uhf, restart)
    elif job == 'hess' and level not in ['gfn2', 'gfn1', 'gfn2spinpol']: # If not GFNx, use ORCA GFN2 for hess
        prep_level = 'gfn2spinpol' if env.geolevel == "gfn2spinpol" else 'gfn2'
        job_command_list, output_file, success_pattern, cleanup_command = _prep_orca_py(env, fname_xyz, prep_level, job, chrg, uhf, restart)
    else: # General dispatcher
        if level in ['gfn2', 'gfn2spinpol', 'gfn1', 'pm6', 'dxtb', 'gff']:
            job_command_list, output_file, success_pattern, cleanup_command = _prep_xtb_py(env, fname_xyz, level, job, chrg, uhf, restart)
        elif level in ['pbe', 'b973c', 'r2scan3c', 'pbeh3c', 'wb97x3c', 'pbe0', 'ccsdt', 
                       'kpr2scan50d4', 'pbe0tzvpd', 'wb97xd4tz', 'wb97xd4matzvp']: # Add all supported ORCA levels
            job_command_list, output_file, success_pattern, cleanup_command = _prep_orca_py(env, fname_xyz, level, job, chrg, uhf, restart)
        elif level == 'gxtb':
            if job in ['sp', 'opt', 'hess']: # gXTB supports these
                 job_command_list, output_file, success_pattern, cleanup_command = _prepgxtb_py(env, fname_xyz, level, job, chrg, uhf, restart)
            else:
                print(f"Warning! Job '{job}' not directly supported for gXTB via main dispatcher. Check logic.")
                return None, None, None, None, False, None
        else:
            print(f"Warning! QM level '{level}' is not supported by prepqm.")
            return None, None, None, None, False, None

    if env.printlevel == 3: # Debug mode, don't clean up
        cleanup_command = "echo 'QM files kept due to printlevel 3'"
        
    # If opt job is being (re)run, invalidate old sp/tddft/opt entries
    if job == 'opt' and not is_cached:
        queries_to_remove = [_construct_qmdata_query(level, "sp", chrg, uhf), 
                             _construct_qmdata_query(level, "opt", chrg, uhf)]
        if env.exstates > 0:
            queries_to_remove.append(_construct_qmdata_query(level, "tddft", chrg, uhf))
        _remove_entries_from_qmdata(queries_to_remove)

    return job_command_list, output_file, success_pattern, cleanup_command, False, None


def _readout_gxtb_sp_py(output_file: str) -> Tuple[bool, Optional[WP]]:
    """Reads energy from gXTB 'energy' file (specific format)."""
    # Fortran: reads 2nd value on 2nd line of 'energy' file
    energy_file = Path("energy") # gXTB SP writes to 'energy'
    if not energy_file.exists():
        print(f"Warning: gXTB 'energy' file not found after SP calculation via {output_file}.")
        return False, None
    try:
        with open(energy_file, 'r') as f:
            f.readline() # Skip first line
            second_line = f.readline().strip()
            parts = second_line.split()
            if len(parts) >= 2:
                return True, float(parts[1])
            else:
                print(f"Warning: Not enough values on second line of gXTB 'energy' file: {second_line}")
                return False, None
    except Exception as e:
        print(f"Error reading gXTB 'energy' file: {e}")
        return False, None

def _readout_tddft_states_py(env: RunTypeData, orca_out_file: str) -> Tuple[Optional[WP], bool]:
    """Parses ORCA TDDFT output for excited state energies, returns Boltzmann weighted energy (in Hartrees)."""
    # Fortran version read up to nstates=10, energies in eV, then converted bwenergie to Hartrees.
    # This needs robust parsing of ORCA output for "STATE X: E=" lines.
    
    energies_ev: List[WP] = [0.0] # Ground state reference is 0 eV relative for excited states
    n_actual_states_found = 0
    max_states_to_read = env.exstates if env.exstates > 0 else 10 # Default to 10 if not specified

    try:
        with open(orca_out_file, 'r') as f:
            for line in f:
                if "STATE " in line and "eV" in line : # Example: "STATE  1:  E=   2.751412 eV"
                    match = re.search(r"STATE\s+\d+:\s+E=\s*([-\d.]+)\s*eV", line)
                    if match:
                        energies_ev.append(float(match.group(1)))
                        n_actual_states_found += 1
                        if n_actual_states_found >= max_states_to_read:
                            break
        if not energies_ev or n_actual_states_found == 0: # Only GS or no states found
             return 0.0, False # No contribution from excited states, not failed

        # Call utility.boltzmann_average_energy (adapted from Fortran boltzmannweights)
        # It expects energies relative to some reference. Here, they are excitation energies.
        # The Fortran code's `boltzmannweights` calculated `bwenergie = sum(pop(i)*energies(i))`.
        # If energies(i) are excitation energies, then bwenergie is the average excitation energy.
        # This average excitation energy (in eV) is then converted to Hartrees and added.
        avg_excitation_ev = utility.boltzmann_average_energy(env, energies_ev) 
        # Convert average excitation energy from eV to Hartrees
        ev_to_au = 1.0 / 27.211386245988 # Standard conversion factor
        bwenergie_au = avg_excitation_ev * ev_to_au
        return bwenergie_au, False # False means not failed

    except FileNotFoundError:
        print(f"Error: {orca_out_file} not found for TDDFT state readout.")
        return None, True
    except Exception as e:
        print(f"Error parsing TDDFT states from {orca_out_file}: {e}")
        return None, True


ReadResult = Tuple[Optional[WP], bool]

def readoutqm_py(env: RunTypeData, fname_orig_xyz: str, level: str, job: str, 
                 qm_out_file: str, success_pattern: str) -> ReadResult:
    """
    Reads QM output, extracts value, writes to qmdata, updates geometry if opt.
    Returns (value_float, failed_bool)
    """
    # Ensure fname_orig_xyz is a Path object for consistency
    fname_path = Path(fname_orig_xyz)

    chrg, uhf = _get_chrg_uhf(env, str(fname_path), None, None) # Read current .CHRG/.UHF
    query = _construct_qmdata_query(level, job, chrg, uhf)

    # Re-check cache: This might be redundant if prepqm already returned cached value.
    # However, useful if readoutqm is called independently or if prepqm's cache check was bypassed.
    is_cached, cached_value = _check_qmdata_cache(query)
    if is_cached and cached_value is not None:
        return cached_value, False

    # Handle H+ special case from Fortran
    sumformula = utility.sum_formula_from_xyz(str(fname_path)) # Assumes utility has this
    if sumformula == "H1" and chrg == 1 and level not in ['gfn1', 'gfn2', 'gfn2spinpol', 'gxtb']:
        _write_to_qmdata(query, 0.0)
        return 0.0, False

    # Actual parsing
    value: Optional[WP] = None
    failed = False
    
    if success_pattern == " GXTBSPENERGYFILE ": # Special flag for gXTB SP
        found_gxtb, value_gxtb = _readout_gxtb_sp_py(qm_out_file)
        if not found_gxtb:
            print(f"Warning: gXTB SP energy pattern not found or error reading 'energy' file from {qm_out_file}.")
            failed = True
        else:
            value = value_gxtb
    elif job == 'optts' and Path(qm_out_file).name == "orca.out": # ORCA optts specific
        # Check for convergence first
        converged = iomod.minigrep(qm_out_file, success_pattern) # success_pattern is "THE OPTIMIZATION HAS CONVERGED"
        if not converged:
            print(f"Warning: ORCA OptTS did not converge (pattern '{success_pattern}' not found in {qm_out_file}).")
            failed = True
        else: # If converged, the energy is usually not the primary result of OptTS, but the structure.
              # Fortran code didn't parse energy for OptTS, it only moved opt.xyz.
              # For consistency, we might parse the final energy if available.
            found_sp, sp_val = iomod.grepval(qm_out_file, "FINAL SINGLE POINT ENERGY")
            if found_sp: value = sp_val
            else: value = 0.0 # Or mark as failed to get energy
    else: # General case
        found_pattern, parsed_value = iomod.grepval(qm_out_file, success_pattern)
        if not found_pattern:
            print(f"Warning: Success pattern '{success_pattern}' not found in {qm_out_file}.")
            failed = True
        else:
            value = parsed_value

    if failed:
        iomod.print_pwd(f"QM calculation failed in directory containing {qm_out_file}")
        return None, True

    # Handle TDDFT post-processing for ORCA
    if job == 'tddft' and level not in ['gfn1', 'gfn2', 'gfn2spinpol', 'gxtb'] and value is not None: # Assuming value is ground state SCF+Disp
        bw_contrib_au, tddft_failed = _readout_tddft_states_py(env, qm_out_file)
        if tddft_failed:
            failed = True
        elif bw_contrib_au is not None:
            value += bw_contrib_au # Add Boltzmann weighted excitation energy
        # Write ground state SP energy also to qmdata for reference
        sp_query = _construct_qmdata_query(level, "sp", chrg, uhf)
        _write_to_qmdata(sp_query, value - (bw_contrib_au or 0.0)) # Subtract to get original SP

    if value is not None:
        _write_to_qmdata(query, value)
    else: # Should have been caught by failed=True earlier
        print(f"Warning: Value is None after parsing {qm_out_file} for job {job}, level {level}.")
        return None, True

    # Geometry update for optimizations
    if job in ['opt', 'optts', 'ohess'] and not failed:
        opt_geom_file_src: Optional[Path] = None
        if Path("opt.xyz").exists(): # General name used by some wrappers
            opt_geom_file_src = Path("opt.xyz")
        elif Path("orca.xyz").exists() and "orca" in qm_out_file.lower() : # ORCA output opt geom
             opt_geom_file_src = Path("orca.xyz")
        elif Path("xtbopt.xyz").exists() and "xtb" in qm_out_file.lower(): # XTB output opt geom
             opt_geom_file_src = Path("xtbopt.xyz")
        
        if opt_geom_file_src and opt_geom_file_src.exists():
            try:
                shutil.move(str(opt_geom_file_src), str(fname_path))
                print(f"Updated {fname_path} with optimized geometry from {opt_geom_file_src}.")
            except Exception as e:
                print(f"Error moving optimized geometry: {e}")
        else:
            print(f"Warning: Optimized geometry file (e.g. opt.xyz, orca.xyz, xtbopt.xyz) not found after {job} for {fname_path}")
            # This might not be a failure of the energy calculation, but good to note.

    return value, failed

# --- Placeholder for utility functions that would be in utility.py ---
# These are simplified versions or direct calls that assume utility.py is structured accordingly.

if not hasattr(utility, 'print_spin_py'):
    def _print_spin_py_placeholder(env, fname_xyz, chrg_val):
        print(f"Placeholder: Estimating spin for {fname_xyz}, charge {chrg_val}")
        # Simplified logic, replace with actual XYZ parsing and electron counting
        if chrg_val == 0: return 1 # Singlet for neutral
        return 2 # Doublet for charged (common default)
    utility.print_spin_py = _print_spin_py_placeholder

if not hasattr(utility, 'sum_formula_from_xyz'):
    def _sum_formula_from_xyz_placeholder(fname_xyz_str):
        print(f"Placeholder: Reading sum formula for {fname_xyz_str}")
        # Extremely simplified, assumes H1 for testing H+ case
        if Path(fname_xyz_str).name == "hplus.xyz": return "H1" 
        return "Unknown"
    utility.sum_formula_from_xyz = _sum_formula_from_xyz_placeholder

if not hasattr(utility, 'boltzmann_average_energy'):
    def _boltzmann_average_energy_placeholder(env, energies_list_ev):
        print(f"Placeholder: Calculating Boltzmann average for {energies_list_ev} at T={env.temp}K")
        if not energies_list_ev: return 0.0
        # Simplified: just returns the first energy, not a real Boltzmann average
        return energies_list_ev[0] 
    utility.boltzmann_average_energy = _boltzmann_average_energy_placeholder


if __name__ == '__main__':
    print("QMMod main execution for testing.")
    # Create dummy env and files for testing
    test_env = RunTypeData(cores=1, threads=1) # Minimal env
    test_env.printlevel = 3 # Keep files for debugging this test
    
    # Create a dummy XYZ file
    xyz_content = "1\nH atom\nH 0.0 0.0 0.0\n"
    Path("h_atom.xyz").write_text(xyz_content)
    Path("hplus.xyz").write_text(xyz_content) # For H+ test

    # Test 1: XTB SP calculation for H atom (neutral, doublet)
    print("\n--- Test 1: XTB SP for H atom (neutral, doublet) ---")
    cmd_list, out_f, patt, clean_cmd, cached, val = prepqm_py(test_env, "h_atom.xyz", "gfn2", "sp", chrg_in=0, uhf_in=1)
    if not cached and cmd_list:
        print(f"Running XTB: {' '.join(cmd_list)} > {out_f} 2>/dev/null")
        result = iomod.execute_command(cmd_list, shell=False) # Use shell=False for list of args
        if result.returncode == 0:
            print(f"XTB finished. Output in {out_f}")
            energy, failed = readoutqm_py(test_env, "h_atom.xyz", "gfn2", "sp", out_f, patt)
            print(f"XTB SP H atom: Energy={energy}, Failed={failed}")
        else:
            print(f"XTB failed with code {result.returncode}. Stderr:\n{result.stderr}")
        if clean_cmd: iomod.execute_command(clean_cmd.split(), shell=True) # Cleanup uses shell
    elif cached:
        print(f"XTB SP H atom: Cached value={val}")

    # Test 2: ORCA Opt for H2O (example, needs ORCA and proper setup)
    # This test will likely fail if ORCA is not installed or configured.
    # h2o_xyz = "3\nWater\nO 0.000000 0.000000 0.117300\nH 0.000000 0.757200 -0.469200\nH 0.000000 -0.757200 -0.469200\n"
    # Path("h2o.xyz").write_text(h2o_xyz)
    # print("\n--- Test 2: ORCA Opt for H2O (neutral, singlet) ---")
    # test_env_orca = RunTypeData(cores=1, threads=1, qmprog='orca', geolevel='r2scan3c', tslevel='r2scan3c')
    # test_env_orca.printlevel = 3
    # cmd_list_orca, out_f_orca, patt_orca, clean_cmd_orca, cached_orca, val_orca = \
    #     prepqm_py(test_env_orca, "h2o.xyz", "r2scan3c", "opt", chrg_in=0, uhf_in=0)
    
    # if not cached_orca and cmd_list_orca and cmd_list_orca[0] != "echo": # Check if ORCA found
    #     print(f"Running ORCA: {' '.join(cmd_list_orca)} > {out_f_orca} 2>&1") # Redirect to see output for test
    #     # For testing, might want to see ORCA output directly:
    #     # result_orca = iomod.execute_command(cmd_list_orca, capture_output=False, text=True)
    #     with open(out_f_orca, "w") as orca_out_handle, open("orca_err.txt", "w") as orca_err_handle:
    #          result_orca = subprocess.run(cmd_list_orca, stdout=orca_out_handle, stderr=orca_err_handle, text=True)

    #     if result_orca.returncode == 0:
    #         print(f"ORCA finished. Output in {out_f_orca}")
    #         energy_orca, failed_orca = readoutqm_py(test_env_orca, "h2o.xyz", "r2scan3c", "opt", out_f_orca, patt_orca)
    #         print(f"ORCA Opt H2O: Energy={energy_orca}, Failed={failed_orca}")
    #         if Path("h2o.xyz").read_text() != h2o_xyz: print("h2o.xyz was updated.")
    #     else:
    #         print(f"ORCA failed with code {result_orca.returncode}. Output in {out_f_orca}, errors in orca_err.txt")
    #     if clean_cmd_orca: iomod.execute_command(clean_cmd_orca.split(), shell=True)
    # elif cached_orca:
    #     print(f"ORCA Opt H2O: Cached value={val_orca}")
    # else:
    #     print("ORCA command not generated (possibly ORCA not found or other prep error).")

    # Cleanup test files
    for f_name in ["h_atom.xyz", "hplus.xyz", ".CHRG", ".UHF", "qmdata", "xtbin.xyz", "xtb.out", "orca.inp", "orca.out", "orca_err.txt", "opt.xyz", "nmode", "energy"]:
        Path(f_name).unlink(missing_ok=True)
    for d_name in ["gxtb_files"]: # Example if gxtb creates specific dirs
        if Path(d_name).is_dir(): shutil.rmtree(d_name, ignore_errors=True)

```
