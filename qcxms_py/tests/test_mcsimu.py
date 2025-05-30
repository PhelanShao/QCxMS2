import pytest
import shutil
import os
from pathlib import Path
from unittest.mock import patch # For more complex/targeted mocking if needed via @patch decorator

# Assuming qcxms is installed or PYTHONPATH is set up for imports
from qcxms.data import RunTypeData, WP
from qcxms.mcsimu import montecarlo_py # Target function
from qcxms import iomod
from qcxms import utility # To be mocked
from qcxms import cid # To be mocked
from qcxms.constants import HARTREE_TO_EV, EVTOKCAL, AUTOEV, KB_EV_K, PI # Ensure all needed are here

# Minimal mock for qmmod if _get_chrg_uhf is indirectly called by a non-mocked utility
# For this test, it seems not directly needed if utility functions using it are mocked.

@pytest.fixture
def simple_env(tmp_path: Path) -> RunTypeData:
    """Pytest fixture for a simple RunTypeData instance."""
    env = RunTypeData()
    env.startdir = str(tmp_path) # Set startdir to tmp_path for isolated testing
    # env.path will be set by chdir in the test or by functions if they manage it
    
    env.tslevel = "gfn2" # Example, used for file naming
    env.geolevel = "gfn2" # Example
    env.iplevel = "gfn2"  # Example

    env.tf = 50.0  # Time of flight (us)
    env.tfscale = 0.0 # No scaling of ToF with fragmentation level for simplicity
    
    env.nsamples = 100 # Number of MC samples (keep low for speed, but >0)
    env.eyring = False # Test RRKM by default
    env.eyzpve = False
    env.nots = False # Enable TS processing (barriers will be read)
    env.bhess = False # No hessian/RRHO for this simple test (mock _xtb_thermo_py)
    env.nthermosteps = 1 # Simplifies thermo mocking
    
    # KER calculation settings (can be simplified if not asserting KER values)
    env.scaleker = 1.0

    env.mode = "ei" 
    env.chrg = 1 # Example charge
    env.pthr = 0.0 # Probability threshold, set low to not interfere with small N samples
    
    env.scaleeinthdiss = 1.0 # No scaling for H-dissociation IEE
    env.cneintscale = False # No CNEINT scaling

    env.cid_elab = 0.0 # Not a CID run

    # For _xtb_thermo_py if it were called (it will be mocked)
    env.temp = 298.15
    env.sthr = 100.0
    env.ithr = 0.1
    env.geolevel_for_xtb_thermo = "gfn2"


    return env

def test_montecarlo_simple_competitive_fragmentation(simple_env: RunTypeData, tmp_path: Path, mocker):
    """
    Test montecarlo_py with two competitive reaction channels:
    1. Isomerization (p0)
    2. Fragmentation (p1) into f1 and f2.
    Channel p1 has a lower barrier and should be favored.
    Thermal corrections are mocked to zero.
    """
    precursor_dir = tmp_path / "precursor_mol"
    precursor_dir.mkdir()

    # Product directories (relative to env.startdir which is tmp_path)
    # These are where pfrag files for products will be written.
    prod_p0_isomer_dir_rel = Path("product_isomer_p0")
    prod_p1f1_dir_rel = Path("product_frag_p1f1")
    prod_p1f2_dir_rel = Path("product_frag_p1f2")

    (tmp_path / prod_p0_isomer_dir_rel).mkdir(exist_ok=True)
    (tmp_path / prod_p1f1_dir_rel).mkdir(exist_ok=True)
    (tmp_path / prod_p1f2_dir_rel).mkdir(exist_ok=True)

    # --- Create mock input files ---
    # Precursor directory files
    (precursor_dir / "fragment.xyz").write_text("3\nPrecursor\nC 0 0 0\nH 1 0 0\nH -1 0 0\n")
    # For `montecarlo_py`, `p0_intensity_fraction` is passed as param. No need for pfrag in precursor_dir itself.
    # `keravold` is read if fragmentation_level > 1. For level 0, it's not.
    # `tav0` is read if fragmentation_level > 1. For level 0, it's not.

    # Channel p0 (Isomerization) - reaction subdir of precursor_dir
    channel_p0_dir = precursor_dir / "p0"
    channel_p0_dir.mkdir()
    # (channel_p0_dir / "isomer.xyz").write_text("3\nIsomer P0\nC 0 0 0\nH 0 1 0\nH 0 -1 0\n") # Not read by mcsimu
    (channel_p0_dir / f"barrier_{simple_env.tslevel}").write_text("0.5\n") # eV
    (channel_p0_dir / f"de_{simple_env.tslevel}").write_text("-0.1\n")      # eV
    (channel_p0_dir / f"ircmode_{simple_env.geolevel}").write_text("500.0\n")# cm-1
    # TS XYZ for thermo (mocked, but file existence might be checked)
    (channel_p0_dir / "ts").mkdir()
    (channel_p0_dir / "ts" / "ts.xyz").write_text("3\nTS for P0\nC 0 0.1 0\nH 1 0.1 0\nH -1 0.1 0\n")


    # Channel p1 (Fragmentation) - reaction subdir of precursor_dir
    channel_p1_dir = precursor_dir / "p1"
    channel_p1_dir.mkdir()
    # (channel_p1_dir / "pair.xyz").write_text("3\nFragment Pair P1\n...\n") # Not read by mcsimu
    (channel_p1_dir / f"barrier_{simple_env.tslevel}").write_text("0.3\n") # eV (Lower barrier)
    (channel_p1_dir / f"de_{simple_env.tslevel}").write_text("0.2\n")       # eV
    (channel_p1_dir / f"ircmode_{simple_env.geolevel}").write_text("400.0\n")# cm-1
    (channel_p1_dir / "ts").mkdir()
    (channel_p1_dir / "ts" / "ts.xyz").write_text("3\nTS for P1\nC 0 0.2 0\nH 1 0.2 0\nH -1 0.2 0\n")

    # Product structure files needed for qmass calculation and _xtb_thermo (though mocked)
    # These go into the product directories defined relative to env.startdir
    (tmp_path / prod_p0_isomer_dir_rel / "isomer.xyz").write_text("3\nIsomer P0\nC 0 0 0\nH 0 1 0\nH 0 -1 0\n")
    (tmp_path / prod_p1f1_dir_rel / "fragment.xyz").write_text("1\nFrag P1F1\nH 0 0 0\n")
    (tmp_path / prod_p1f2_dir_rel / "fragment.xyz").write_text("2\nFrag P1F2\nC 0 0 0\nH 1 0 0\n")


    # --- IEE Distribution ---
    eiee_dist_energies = [1.0]  # Single energy point (eV)
    eiee_dist_probs = [1.0]   # Probability for this point

    # --- Fragdirs data ---
    # fragdirs[i][0] = reaction channel subdir name (e.g. "p0")
    # fragdirs[i][1] = product1 directory path (relative to env.startdir) for pfrag writing & qmass
    # fragdirs[i][2] = product2 directory path (relative to env.startdir) for pfrag writing & qmass
    fragdirs_data = [
        ["p0", str(prod_p0_isomer_dir_rel), ""],  # Isomerization: prod1 is the isomer dir, no prod2
        ["p1", str(prod_p1f1_dir_rel), str(prod_p1f2_dir_rel)]   # Fragmentation
    ]

    # --- Mocking external dependencies ---
    # Mock _xtb_thermo_py: returns list of zeros (no thermal correction)
    mocker.patch("qcxms.mcsimu._xtb_thermo_py", return_value=[0.0] * simple_env.nthermosteps)

    # Mock utility functions for nat, nvib, mass
    mocker.patch("qcxms.utility.get_atomic_numbers_and_coords_py", side_effect=[
        (3, [6,1,1], []), # Precursor (fragment.xyz)
        (3, [6,1,1], []), # TS for p0
        (3, [6,1,1], []), # Isomer product for p0 (for thermo, though G is 0)
        (3, [6,1,1], []), # TS for p1
        (1, [1], []),     # Product p1f1 (for thermo)
        (2, [6,1], []),   # Product p1f2 (for thermo)
    ])
    mocker.patch("qcxms.utility.is_linear_molecule_py", return_value=False) # Assume all non-linear
    
    # Mock for qmass calculation (get_atomic_numbers_from_xyz_py is different from get_atomic_numbers_and_coords_py)
    def mock_get_atomic_numbers(filepath_str):
        if str(prod_p0_isomer_dir_rel / "isomer.xyz") in filepath_str : return [6,1,1] # Isomer
        if str(prod_p1f1_dir_rel / "fragment.xyz") in filepath_str : return [1]   # p1f1
        if str(prod_p1f2_dir_rel / "fragment.xyz") in filepath_str : return [6,1] # p1f2
        return []
    mocker.patch("qcxms.utility.get_atomic_numbers_from_xyz_py", side_effect=mock_get_atomic_numbers)
    
    def mock_get_avg_mass(atomic_numbers_list):
        mass = 0
        for at_num in atomic_numbers_list:
            if at_num == 6: mass += 12
            elif at_num == 1: mass += 1
            else: mass += at_num # fallback
        return float(mass)
    mocker.patch("qcxms.utility.get_average_mol_mass_py", side_effect=mock_get_avg_mass)


    # Mock CID simulation (not used in EI mode, but good practice)
    mocker.patch("qcxms.cid.simulate_cid_activation_py", 
                 return_value=(eiee_dist_energies, eiee_dist_probs, 0.0, 0.0, 0.0))


    # --- Call montecarlo_py ---
    # It operates in the current working directory, which should be the precursor's directory.
    os.chdir(precursor_dir)
    
    success = montecarlo_py(
        simple_env,
        npairs=2,
        fragdirs=fragdirs_data,
        eiee_dist_energies=eiee_dist_energies,
        eiee_dist_probs=eiee_dist_probs,
        fragmentation_level=0, # First level of MC from a precursor
        p0_intensity_fraction=1.0
    )
    os.chdir(tmp_path) # Change back for cleanup/other tests

    assert success, "montecarlo_py reported failure"

    # --- Assertions ---
    # Check pfrag files in product directories
    pfrag_p0_val = iomod.rdshort_real_py(tmp_path / prod_p0_isomer_dir_rel / "pfrag", default=-1.0)
    # For fragmentation, pfrag is written to *both* fragment product dirs with the same value (channel probability)
    pfrag_p1f1_val = iomod.rdshort_real_py(tmp_path / prod_p1f1_dir_rel / "pfrag", default=-1.0)
    pfrag_p1f2_val = iomod.rdshort_real_py(tmp_path / prod_p1f2_dir_rel / "pfrag", default=-1.0)
    
    # Check pfs (surviving precursor fraction) in precursor_dir
    pfs_precursor_val = iomod.rdshort_real_py(precursor_dir / "pfs", default=-1.0)

    assert pfrag_p0_val >= 0.0, "pfrag for channel p0 not found or invalid."
    assert pfrag_p1f1_val >= 0.0, "pfrag for channel p1 (f1) not found or invalid."
    assert pfrag_p1f2_val == pfrag_p1f1_val, "pfrag for p1f1 and p1f2 should be identical (channel prob)."
    assert pfs_precursor_val >= 0.0, "pfs for surviving precursor not found or invalid."

    # Sum of probabilities should be close to 1.0 (p0_intensity_fraction)
    # pfrag_abs[j][0] is normalized by total_piee_sum (which is 1.0 here).
    # Then written value is pfrag_abs[j][0] * p0_intensity_fraction.
    # So, sum of pfrag_p0_val, pfrag_p1f1_val (representing channel p1), and pfs_precursor_val should be 1.0.
    print(f"pfrags: p0={pfrag_p0_val}, p1={pfrag_p1f1_val}, precursor_stable={pfs_precursor_val}")
    assert abs((pfrag_p0_val + pfrag_p1f1_val + pfs_precursor_val) - 1.0) < 1e-6, "Probabilities do not sum to 1.0"

    # Channel p1 (barrier 0.3 eV) should be favored over p0 (barrier 0.5 eV) at E_int=1.0 eV
    assert pfrag_p1f1_val > pfrag_p0_val, \
        f"Channel p1 (lower barrier) was not favored: p1_frag={pfrag_p1f1_val}, p0_iso={pfrag_p0_val}"

    # Check qmass for channel p1
    qmass_p1 = iomod.rdshort_real_py(precursor_dir / "p1" / "qmass", default=-1.0)
    # mass_f1 (H) = 1, mass_f2 (CH) = 13. Expected qmass = 13/1 = 13.0
    assert abs(qmass_p1 - (12.0+1.0)/1.0) < 1e-6, f"qmass for p1 is incorrect: {qmass_p1}"
    
    # Check kerav for channel p1 (more complex to predict exact value without deep dive, ensure non-negative)
    kerav_p1 = iomod.rdshort_real_py(precursor_dir / "p1" / "kerav", default=-1.0)
    assert kerav_p1 >= 0.0, f"kerav for p1 is negative: {kerav_p1}"

    # Check keravold in precursor dir
    keravold_precursor = iomod.rdshort_real_py(precursor_dir / "keravold", default=-1.0)
    assert keravold_precursor >= 0.0, f"keravold for precursor is negative: {keravold_precursor}"

    # Check tav in precursor dir
    tav_precursor = iomod.rdshort_real_py(precursor_dir / "tav", default=-1.0)
    assert tav_precursor >= 0.0, f"tav for precursor is negative: {tav_precursor}"

```
