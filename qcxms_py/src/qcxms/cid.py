import os
import math
import random
from pathlib import Path
from typing import List, Tuple, Optional, Dict
import numpy as np

try:
    from .data import RunTypeData, WP
    from . import iomod
    from . import utility
    from .constants import AUTOEV, EVTOKCAL, KB_EV_K, BOHR
except ImportError:
    # Fallbacks for standalone/testing
    print("Attempting to import dummy/mock modules for cid.py standalone run.")
    from data import RunTypeData, WP # type: ignore
    import iomod_mock as iomod # type: ignore
    import utility_mock as utility # type: ignore
    AUTOEV = 27.211386245988
    EVTOKCAL = 23.060547830618307
    KB_EV_K = 8.6173332621415e-5
    BOHR = 0.52917721067


def simulate_cid_activation_py(
    env: RunTypeData,
    precursor_xyz_file: str,
    collision_energy_lab_ev: float,
    target_gas: str = "Ar"
) -> Tuple[float, float, float]:
    """
    Simulates Collision-Induced Dissociation (CID) activation.
    
    This function models the collision between a precursor ion and a neutral target gas,
    calculating the energy transfer and internal energy deposition.
    
    Args:
        env: Runtime environment data
        precursor_xyz_file: Path to precursor molecule XYZ file
        collision_energy_lab_ev: Laboratory collision energy in eV
        target_gas: Target gas species (default: "Ar")
    
    Returns:
        Tuple of (internal_energy_deposited_ev, kinetic_energy_loss_ev, scattering_angle_deg)
    """
    print(f"  Simulating CID activation: E_lab = {collision_energy_lab_ev:.2f} eV, target = {target_gas}")
    
    try:
        # Get molecular properties
        natoms, atomic_numbers, coords = utility.get_atomic_numbers_and_coords_py(precursor_xyz_file)
        if natoms == 0 or atomic_numbers is None:
            print(f"    Error: Could not read precursor structure from {precursor_xyz_file}")
            return 0.0, 0.0, 0.0
        
        # Calculate molecular mass
        precursor_mass_amu = utility.get_average_mol_mass_py(atomic_numbers)
        
        # Target gas properties
        target_masses = {
            "He": 4.003, "Ne": 20.18, "Ar": 39.948, "Kr": 83.798, "Xe": 131.29,
            "N2": 28.014, "O2": 31.998, "CO2": 44.01, "H2": 2.016
        }
        target_mass_amu = target_masses.get(target_gas, 39.948)  # Default to Ar
        
        # Convert to center-of-mass collision energy
        reduced_mass = (precursor_mass_amu * target_mass_amu) / (precursor_mass_amu + target_mass_amu)
        collision_energy_cm_ev = collision_energy_lab_ev * (target_mass_amu / (precursor_mass_amu + target_mass_amu))
        
        print(f"    Precursor mass: {precursor_mass_amu:.1f} amu, Target mass: {target_mass_amu:.1f} amu")
        print(f"    CM collision energy: {collision_energy_cm_ev:.2f} eV")
        
        # Sample scattering angle using realistic distribution
        scattering_angle_rad = _sample_scattering_angle_py(collision_energy_cm_ev, precursor_mass_amu, target_mass_amu)
        scattering_angle_deg = math.degrees(scattering_angle_rad)
        
        # Calculate energy transfer using classical mechanics
        # Energy transfer depends on scattering angle and masses
        max_energy_transfer = 4 * reduced_mass * collision_energy_cm_ev / precursor_mass_amu
        energy_transfer_factor = (1 - math.cos(scattering_angle_rad)) / 2
        kinetic_energy_transfer_ev = max_energy_transfer * energy_transfer_factor
        
        # Not all kinetic energy transfer becomes internal energy
        # Some fraction goes to rotational/translational modes
        internal_efficiency = _calculate_internal_efficiency_py(
            env, natoms, collision_energy_cm_ev, scattering_angle_deg
        )
        
        internal_energy_deposited_ev = kinetic_energy_transfer_ev * internal_efficiency
        
        print(f"    Scattering angle: {scattering_angle_deg:.1f}°")
        print(f"    Energy transfer: {kinetic_energy_transfer_ev:.3f} eV")
        print(f"    Internal efficiency: {internal_efficiency:.3f}")
        print(f"    Internal energy deposited: {internal_energy_deposited_ev:.3f} eV")
        
        return internal_energy_deposited_ev, kinetic_energy_transfer_ev, scattering_angle_deg
        
    except Exception as e:
        print(f"    Error in CID simulation: {e}")
        return 0.0, 0.0, 0.0


def _sample_scattering_angle_py(collision_energy_cm_ev: float, precursor_mass: float, target_mass: float) -> float:
    """
    Sample scattering angle from a realistic distribution.
    Uses a combination of hard sphere and attractive potential models.
    """
    try:
        # Convert energy to velocity for classical calculation
        # E = 1/2 * mu * v^2, where mu is reduced mass
        reduced_mass_amu = (precursor_mass * target_mass) / (precursor_mass + target_mass)
        
        # Sample impact parameter (0 to b_max)
        # b_max is roughly the sum of molecular radii
        molecular_radius_bohr = _estimate_molecular_radius_py(precursor_mass)
        target_radius_bohr = _estimate_molecular_radius_py(target_mass)
        b_max_bohr = molecular_radius_bohr + target_radius_bohr
        
        # Random impact parameter
        b_bohr = random.uniform(0, b_max_bohr)
        
        # Classical scattering angle for hard sphere collision
        if b_bohr >= b_max_bohr:
            # No collision
            return 0.0
        else:
            # Hard sphere scattering
            sin_half_theta = b_bohr / b_max_bohr
            half_theta = math.asin(min(1.0, sin_half_theta))
            theta_rad = 2 * half_theta
            
            # Add some quantum/thermal broadening
            thermal_broadening = random.gauss(0, 0.1)  # Small random component
            theta_rad = max(0, min(math.pi, theta_rad + thermal_broadening))
            
            return theta_rad
            
    except Exception as e:
        print(f"    Error sampling scattering angle: {e}")
        # Return a reasonable default
        return random.uniform(0, math.pi/2)


def _estimate_molecular_radius_py(mass_amu: float) -> float:
    """
    Estimate molecular radius in Bohr based on mass.
    Uses empirical scaling relationship.
    """
    # Rough scaling: radius ~ mass^(1/3) for spherical molecules
    # Calibrated to give reasonable values for common molecules
    radius_bohr = 3.0 * (mass_amu / 40.0)**(1/3)  # Ar as reference
    return max(2.0, radius_bohr)  # Minimum radius


def _calculate_internal_efficiency_py(
    env: RunTypeData, 
    natoms: int, 
    collision_energy_ev: float, 
    scattering_angle_deg: float
) -> float:
    """
    Calculate the efficiency of converting kinetic energy to internal energy.
    This depends on molecular size, collision energy, and scattering dynamics.
    """
    try:
        # Base efficiency depends on molecular size
        # Larger molecules have more vibrational modes to accept energy
        nvib = max(1, 3 * natoms - 6)  # Number of vibrational modes
        
        # Base efficiency scales with number of vibrational modes
        base_efficiency = min(0.8, 0.1 + 0.02 * nvib)
        
        # Efficiency increases with collision energy (more violent collisions)
        energy_factor = 1.0 + 0.1 * math.log(1 + collision_energy_ev / 5.0)
        
        # Efficiency depends on scattering angle
        # Head-on collisions (large angles) are more efficient
        angle_factor = 0.5 + 0.5 * (scattering_angle_deg / 180.0)
        
        # Temperature effect - higher temperatures reduce efficiency
        temp_factor = 1.0 / (1.0 + env.temp / 1000.0)
        
        efficiency = base_efficiency * energy_factor * angle_factor * temp_factor
        
        # Clamp to reasonable range
        return max(0.05, min(0.95, efficiency))
        
    except Exception as e:
        print(f"    Error calculating internal efficiency: {e}")
        return 0.3  # Default reasonable value


def calculate_cid_energy_distribution_py(
    env: RunTypeData,
    precursor_xyz_file: str,
    num_collisions: int = 1000
) -> Tuple[List[float], List[float]]:
    """
    Calculate the distribution of internal energies from CID activation.
    
    Args:
        env: Runtime environment data
        precursor_xyz_file: Path to precursor molecule XYZ file
        num_collisions: Number of collision events to simulate
    
    Returns:
        Tuple of (energy_values_ev, probabilities)
    """
    print(f"  Calculating CID energy distribution with {num_collisions} collision events")
    
    try:
        internal_energies = []
        
        # Simulate multiple collision events
        for i in range(num_collisions):
            # Sample collision energy from experimental distribution
            collision_energy_lab = _sample_collision_energy_py(env)
            
            # Simulate single collision
            internal_energy, _, _ = simulate_cid_activation_py(
                env, precursor_xyz_file, collision_energy_lab, env.cid_target_gas
            )
            
            if internal_energy > 0:
                internal_energies.append(internal_energy)
        
        if not internal_energies:
            print("    Warning: No successful collision events simulated")
            return [0.0], [1.0]
        
        # Create energy distribution histogram
        energy_bins = 50
        min_energy = min(internal_energies)
        max_energy = max(internal_energies)
        
        if max_energy <= min_energy:
            return [min_energy], [1.0]
        
        bin_width = (max_energy - min_energy) / energy_bins
        energy_values = []
        probabilities = []
        
        for i in range(energy_bins):
            bin_center = min_energy + (i + 0.5) * bin_width
            bin_min = min_energy + i * bin_width
            bin_max = min_energy + (i + 1) * bin_width
            
            count = sum(1 for e in internal_energies if bin_min <= e < bin_max)
            probability = count / len(internal_energies)
            
            if probability > 0:
                energy_values.append(bin_center)
                probabilities.append(probability)
        
        # Normalize probabilities
        total_prob = sum(probabilities)
        if total_prob > 0:
            probabilities = [p / total_prob for p in probabilities]
        
        print(f"    Generated distribution with {len(energy_values)} energy points")
        print(f"    Energy range: {min(energy_values):.2f} - {max(energy_values):.2f} eV")
        
        return energy_values, probabilities
        
    except Exception as e:
        print(f"    Error calculating CID energy distribution: {e}")
        return [env.cid_elab], [1.0]  # Fallback to single energy


def _sample_collision_energy_py(env: RunTypeData) -> float:
    """
    Sample collision energy from experimental distribution.
    In real CID experiments, there's a distribution around the nominal energy.
    """
    try:
        # Base energy from environment
        base_energy = env.cid_elab
        
        # Add experimental broadening (typically 10-20% FWHM)
        sigma = base_energy * 0.1  # 10% standard deviation
        
        # Sample from Gaussian distribution
        collision_energy = random.gauss(base_energy, sigma)
        
        # Ensure positive energy
        return max(0.1, collision_energy)
        
    except Exception as e:
        print(f"    Error sampling collision energy: {e}")
        return env.cid_elab


def write_cid_results_py(
    env: RunTypeData,
    output_dir: Path,
    energy_values: List[float],
    probabilities: List[float],
    total_internal_energy: float,
    total_kinetic_loss: float,
    average_scattering_angle: float
) -> None:
    """
    Write CID simulation results to files.
    """
    try:
        output_dir.mkdir(exist_ok=True)
        
        # Write energy distribution
        with open(output_dir / "cid_energy_distribution.dat", 'w') as f:
            f.write("# CID Energy Distribution\n")
            f.write("# Energy(eV)  Probability\n")
            for energy, prob in zip(energy_values, probabilities):
                f.write(f"{energy:.6f}  {prob:.6f}\n")
        
        # Write summary statistics
        with open(output_dir / "cid_summary.dat", 'w') as f:
            f.write("# CID Simulation Summary\n")
            f.write(f"Laboratory collision energy: {env.cid_elab:.3f} eV\n")
            f.write(f"Target gas: {getattr(env, 'cid_target_gas', 'Ar')}\n")
            f.write(f"Total internal energy deposited: {total_internal_energy:.3f} eV\n")
            f.write(f"Total kinetic energy loss: {total_kinetic_loss:.3f} eV\n")
            f.write(f"Average scattering angle: {average_scattering_angle:.1f} degrees\n")
            f.write(f"Energy transfer efficiency: {total_internal_energy/total_kinetic_loss:.3f}\n")
        
        # Write individual energy values for Monte Carlo sampling
        iomod.wrshort_real(output_dir / "sumdekin", total_kinetic_loss)
        iomod.wrshort_real(output_dir / "sumdeint", total_internal_energy)
        iomod.wrshort_real(output_dir / "x_trav", average_scattering_angle)
        
        print(f"    CID results written to {output_dir}")
        
    except Exception as e:
        print(f"    Error writing CID results: {e}")


def run_cid_simulation_py(env: RunTypeData, precursor_xyz_file: str) -> Tuple[List[float], List[float]]:
    """
    Main function to run complete CID simulation.
    
    Args:
        env: Runtime environment data
        precursor_xyz_file: Path to precursor molecule XYZ file
    
    Returns:
        Tuple of (energy_values_ev, probabilities) for use in Monte Carlo simulation
    """
    print(f"\n--- CID Simulation for {precursor_xyz_file} ---")
    
    try:
        # Check if CID is enabled
        if not hasattr(env, 'cid_elab') or env.cid_elab <= 0:
            print("  CID simulation disabled (no collision energy specified)")
            return [0.0], [1.0]
        
        # Set default target gas if not specified
        if not hasattr(env, 'cid_target_gas'):
            env.cid_target_gas = "Ar"
        
        # Calculate energy distribution
        energy_values, probabilities = calculate_cid_energy_distribution_py(
            env, precursor_xyz_file, num_collisions=getattr(env, 'cid_nsamples', 1000)
        )
        
        # Calculate summary statistics
        total_internal_energy = sum(e * p for e, p in zip(energy_values, probabilities))
        
        # Single collision for kinetic loss and scattering angle
        _, kinetic_loss, scattering_angle = simulate_cid_activation_py(
            env, precursor_xyz_file, env.cid_elab, env.cid_target_gas
        )
        
        # Write results
        output_dir = Path.cwd() / "cid_results"
        write_cid_results_py(
            env, output_dir, energy_values, probabilities,
            total_internal_energy, kinetic_loss, scattering_angle
        )
        
        print(f"--- CID Simulation Complete ---")
        print(f"  Average internal energy: {total_internal_energy:.3f} eV")
        print(f"  Energy distribution points: {len(energy_values)}")
        
        return energy_values, probabilities
        
    except Exception as e:
        print(f"  Error in CID simulation: {e}")
        return [env.cid_elab * 0.3], [1.0]  # Fallback: 30% efficiency


if __name__ == '__main__':
    print("Testing CID simulation module...")
    
    # Create test environment
    env_test = RunTypeData()
    env_test.cid_elab = 25.0  # 25 eV collision energy
    env_test.cid_target_gas = "Ar"
    env_test.cid_nsamples = 100
    env_test.temp = 298
    
    # Create test molecule
    test_xyz_content = """3
Test molecule
C 0.0 0.0 0.0
H 1.0 0.0 0.0
H -1.0 0.0 0.0
"""
    
    test_xyz_file = Path("test_molecule.xyz")
    with open(test_xyz_file, 'w') as f:
        f.write(test_xyz_content)
    
    # Mock utility functions for testing
    utility.get_atomic_numbers_and_coords_py = lambda fname: (3, [6, 1, 1], [[0,0,0], [1,0,0], [-1,0,0]])
    utility.get_average_mol_mass_py = lambda ats: sum([12.01 if a==6 else 1.008 for a in ats])
    
    try:
        # Test single collision
        print("\n=== Testing single collision ===")
        internal_e, kinetic_e, angle = simulate_cid_activation_py(
            env_test, str(test_xyz_file), 25.0, "Ar"
        )
        print(f"Results: Internal={internal_e:.3f} eV, Kinetic={kinetic_e:.3f} eV, Angle={angle:.1f}°")
        
        # Test energy distribution
        print("\n=== Testing energy distribution ===")
        energies, probs = calculate_cid_energy_distribution_py(env_test, str(test_xyz_file), 50)
        print(f"Distribution: {len(energies)} points, range {min(energies):.2f}-{max(energies):.2f} eV")
        
        # Test full simulation
        print("\n=== Testing full CID simulation ===")
        final_energies, final_probs = run_cid_simulation_py(env_test, str(test_xyz_file))
        print(f"Final distribution: {len(final_energies)} points")
        
    finally:
        # Cleanup
        test_xyz_file.unlink(missing_ok=True)
        if Path("cid_results").exists():
            import shutil
            shutil.rmtree("cid_results")
    
    print("CID simulation testing completed.")
