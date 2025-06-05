import pytest
from pathlib import Path
import shutil

# Adjust import path based on how qcxms_py is structured and installed for testing
# If qcxms_py is installed in editable mode (-e .), this should work from the root.
# If tests are run as part of the package, relative imports might be needed.
from qcxms.data import RunTypeData, WP
from qcxms.tsmod import _pick_ts_from_path_py # Target function
from qcxms import utility # For mocking

# Define dummy content for trajectory and snapshot files
DUMMY_TRAJ_CONTENT_GOOD_PATH = """2
Snapshot 0 Energy -10.0
C 0 0 0
H 0 0 1.0
2
Snapshot 1 Energy -5.0
C 0 0 0
H 0 0 1.1
2
Snapshot 2 Energy -2.0 (TS Candidate)
C 0 0 0
H 0 0 1.2
2
Snapshot 3 Energy -8.0
C 0 0 0
H 0 0 1.15
2
Snapshot 4 Energy -12.0
C 0 0 0
H 0 0 0.9
"""

DUMMY_MEP_ENERGY_CONTENT_GOOD_PATH = """0 -10.0
1 -5.0
2 -2.0
3 -8.0
4 -12.0
"""

# This is what xtb_extract_snapshot_py would write to ts.xyz for frame index 3 (1-based)
# which corresponds to python index 2 (0-based for energy list, then +1 for frame index)
TS_CANDIDATE_XYZ_CONTENT = """2
energy: -2.00000000 kcal/mol, EMAX
C   0.000000   0.000000   0.000000
H   0.000000   0.000000   1.200000
"""

# Content for a path that should trigger cascade detection (TS is an endpoint image)
DUMMY_TRAJ_CONTENT_CASCADE_ENDPOINT = """2
Snapshot 0 Energy -10.0
C 0 0 0
H 0 0 1.0
2
Snapshot 1 Energy -2.0 (TS Candidate at endpoint)
C 0 0 0
H 0 0 1.1
2
Snapshot 2 Energy -12.0
C 0 0 0
H 0 0 0.9
"""
DUMMY_MEP_ENERGY_CONTENT_CASCADE_ENDPOINT = """0 -10.0
1 -2.0
2 -12.0
"""

# Content for a path with multiple maxima
DUMMY_TRAJ_CONTENT_MULTI_MAXIMA = """2
Snapshot 0 Energy -10.0
C 0 0 0
H 0 0 1.0
2
Snapshot 1 Energy -5.0 (Max 1)
C 0 0 0
H 0 0 1.1
2
Snapshot 2 Energy -8.0
C 0 0 0
H 0 0 1.05
2
Snapshot 3 Energy -4.0 (Max 2, Global Max)
C 0 0 0
H 0 0 1.2
2
Snapshot 4 Energy -12.0
C 0 0 0
H 0 0 0.9
"""
DUMMY_MEP_ENERGY_CONTENT_MULTI_MAXIMA = """0 -10.0
1 -5.0
2 -8.0
3 -4.0
4 -12.0
"""


@pytest.fixture
def mock_utility_functions(monkeypatch):
    """Mocks utility functions used by _pick_ts_from_path_py."""
    def mock_xtb_extract_snapshot(traj_file: Path, frame_idx_1_based: int, out_file: Path):
        # Simulate xtb_extract_snapshot:
        # In a real scenario, it would parse traj_file for the specific frame.
        # Here, we just write predefined TS content if the index matches the expected TS.
        # This mock assumes the test will always try to extract the correct TS frame.
        # For DUMMY_TRAJ_CONTENT_GOOD_PATH, TS is frame 3 (index 2 in energies[1:-1], so overall index 2+1=3)
        if frame_idx_1_based == 3 and "good_path" in str(traj_file): # Check for specific test context
             out_file.write_text(TS_CANDIDATE_XYZ_CONTENT)
        elif frame_idx_1_based == 2 and "cascade_endpoint" in str(traj_file): # For cascade_endpoint test
             out_file.write_text("2\nTS at endpoint\nC 0 0 0\nH 0 0 1.1")
        elif frame_idx_1_based == 4 and "multi_maxima" in str(traj_file): # For multi_maxima test
             out_file.write_text("2\nTS at multi_maxima\nC 0 0 0\nH 0 0 1.2")
        else:
            # Generic placeholder if specific content isn't crucial for a particular test branch
            out_file.write_text("2\nMocked Snapshot\nA 0 0 0\nB 0 0 1\n")
        # print(f"MOCK: xtb_extract_snapshot called: traj={traj_file}, frame={frame_idx_1_based}, out={out_file}")

    monkeypatch.setattr(utility, "xtb_extract_snapshot_py", mock_xtb_extract_snapshot)
    # If _pick_ts_from_path_py uses other utility functions that need mocking (e.g. RMSD)
    # they would be added here. For now, focusing on file creation.


def test_pick_ts_simple_case(tmp_path, mock_utility_functions):
    """Tests _pick_ts_from_path_py for a simple, valid NEB path."""
    reaction_dir = tmp_path / "p0"
    reaction_dir.mkdir()

    (reaction_dir / "orca_MEP_energy.dat").write_text(DUMMY_MEP_ENERGY_CONTENT_GOOD_PATH)
    # Add a context hint to the filename for the mock
    (reaction_dir / "orca_MEP_trj_good_path.xyz").rename(reaction_dir / "orca_MEP_trj.xyz")
    (reaction_dir / "orca_MEP_trj.xyz").write_text(DUMMY_TRAJ_CONTENT_GOOD_PATH)


    env = RunTypeData()
    env.sortoutcascade = False # Test without cascade detection first

    found_ts = _pick_ts_from_path_py(env, reaction_dir)

    assert found_ts is True
    ts_xyz_file = reaction_dir / "ts.xyz"
    assert ts_xyz_file.exists()
    
    content = ts_xyz_file.read_text()
    # Check some key parts of the expected content
    assert "2" in content.splitlines()[0] # Atom count
    assert "energy: -2.00000000 kcal/mol, EMAX" in content.splitlines()[1] # Comment line
    assert "H   0.000000   0.000000   1.200000" in content # Coordinates of TS

    assert not (reaction_dir / "cascade").exists()
    assert not (reaction_dir / "endeqts").exists()


def test_pick_ts_cascade_detection_endpoint(tmp_path, mock_utility_functions):
    """Tests cascade detection when TS is an endpoint image."""
    reaction_dir = tmp_path / "p1_cascade_endpoint"
    reaction_dir.mkdir()

    (reaction_dir / "orca_MEP_energy.dat").write_text(DUMMY_MEP_ENERGY_CONTENT_CASCADE_ENDPOINT)
    (reaction_dir / "orca_MEP_trj_cascade_endpoint.xyz").rename(reaction_dir / "orca_MEP_trj.xyz")
    (reaction_dir / "orca_MEP_trj.xyz").write_text(DUMMY_TRAJ_CONTENT_CASCADE_ENDPOINT)
    
    env = RunTypeData()
    env.sortoutcascade = True # Enable cascade detection

    found_ts = _pick_ts_from_path_py(env, reaction_dir)

    assert found_ts is False # Should detect cascade (TS is endpoint)
    # Depending on implementation, a "cascade" file might be created.
    # The current _pick_ts_from_path_py doesn't create it, just returns False.
    # If it were to create it: assert (reaction_dir / "cascade").exists()


def test_pick_ts_cascade_detection_multi_maxima(tmp_path, mock_utility_functions):
    """Tests cascade detection when path has multiple maxima."""
    reaction_dir = tmp_path / "p2_multi_maxima"
    reaction_dir.mkdir()

    (reaction_dir / "orca_MEP_energy.dat").write_text(DUMMY_MEP_ENERGY_CONTENT_MULTI_MAXIMA)
    (reaction_dir / "orca_MEP_trj_multi_maxima.xyz").rename(reaction_dir / "orca_MEP_trj.xyz")
    (reaction_dir / "orca_MEP_trj.xyz").write_text(DUMMY_TRAJ_CONTENT_MULTI_MAXIMA)

    env = RunTypeData()
    env.sortoutcascade = True # Enable cascade detection

    found_ts = _pick_ts_from_path_py(env, reaction_dir)

    assert found_ts is False # Should detect cascade (multiple maxima)
    # If it were to create it: assert (reaction_dir / "cascade").exists()

def test_pick_ts_no_intermediate_images(tmp_path, mock_utility_functions):
    """Tests behavior when there are no intermediate images (only R and P)."""
    reaction_dir = tmp_path / "p3_no_intermediate"
    reaction_dir.mkdir()

    (reaction_dir / "orca_MEP_energy.dat").write_text("0 -10.0\n1 -12.0\n") # Reactant and Product only
    (reaction_dir / "orca_MEP_trj.xyz").write_text("2\nR\nC 0 0 0\nH 0 0 1\n2\nP\nC 0 0 0\nH 0 0 2\n")
    
    env = RunTypeData()
    env.sortoutcascade = False

    found_ts = _pick_ts_from_path_py(env, reaction_dir)
    assert found_ts is False
    assert not (reaction_dir / "ts.xyz").exists()

```
