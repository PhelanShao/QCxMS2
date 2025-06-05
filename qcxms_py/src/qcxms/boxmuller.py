import math
import random
from typing import Tuple

try:
    from .data import WP # Assuming WP is float or similar
except ImportError:
    WP = float # Fallback for standalone execution

PI = math.pi

def _box_muller_transform() -> Tuple[WP, WP]:
    """
    Helper function to generate two independent standard normal deviates (mean 0, stddev 1)
    using the Box-Muller transform.
    """
    # Ensure u1 is not zero for log, random.random() is in [0.0, 1.0)
    # Some implementations use (0.0, 1.0] or handle u1=0 explicitly.
    # random.random() is [0.0, 1.0), so u1 can be 0.0.
    # math.log(0) is undefined.
    u1 = 0.0
    while u1 == 0.0: # Ensure u1 is not zero
        u1 = random.random()
    
    u2 = random.random()

    magnitude = math.sqrt(-2.0 * math.log(u1))
    z0 = magnitude * math.cos(2.0 * PI * u2)
    z1 = magnitude * math.sin(2.0 * PI * u2)
    
    return z0, z1

def vary_collisions_py(mean_collisions: WP, distribution_param: WP) -> int:
    """
    Generates a normally distributed random integer number of collisions.
    
    Args:
        mean_collisions: Mean number of collisions.
        distribution_param: A parameter to calculate standard deviation (sigma = mean_collisions * distribution_param).
        
    Returns:
        A normally distributed integer number of collisions (non-negative).
    """
    if mean_collisions == 0: # If mean is 0, result should be 0, avoid sigma=0 if not intended
        return 0
    if distribution_param == 0: # No variation
        return int(round(mean_collisions))

    sigma = mean_collisions * distribution_param
    if sigma == 0: # Avoid issues if sigma becomes zero due to inputs
        return int(round(mean_collisions))

    z0, _ = _box_muller_transform() # Only z0 is used for the result in Fortran
    
    # Scale and shift z0 to get the normally distributed number
    varied_value_float = z0 * sigma + mean_collisions
    
    # Convert to nearest integer and ensure non-negativity
    box_coll = int(round(varied_value_float))
    
    return max(0, box_coll)


def vary_energies_py(mean_energy: WP, distribution_width: WP) -> WP:
    """
    Generates a normally distributed random energy.

    Args:
        mean_energy: Mean energy.
        distribution_width: Fractional width to calculate standard deviation (sigma = mean_energy * distribution_width).

    Returns:
        A normally distributed random energy.
    """
    if distribution_width == 0: # No variation
        return mean_energy

    sigma = mean_energy * distribution_width
    if sigma == 0: # Avoid issues if sigma becomes zero (e.g. mean_energy is 0)
        return mean_energy

    # Store u1 for the decision, as it's used in Box-Muller
    u1_for_decision = 0.0
    while u1_for_decision == 0.0:
        u1_for_decision = random.random()
    u2_for_bm = random.random()

    magnitude = math.sqrt(-2.0 * math.log(u1_for_decision))
    z0 = magnitude * math.cos(2.0 * PI * u2_for_bm)
    z1 = magnitude * math.sin(2.0 * PI * u2_for_bm)
    
    energy_val1 = z0 * sigma + mean_energy
    energy_val2 = z1 * sigma + mean_energy
    
    # The Fortran code used the first random number (u1, here u1_for_decision)
    # to decide which of the two generated normal deviates to return.
    if u1_for_decision > 0.5: # Note: Fortran `dum` was used for both log and this condition.
                              # Using a fresh random number or the same u1 for decision is a choice.
                              # Fortran used the *same* `dum` that went into log for this check.
        return energy_val1
    else:
        return energy_val2

# Pythonic alternative using random.gauss:
def vary_collisions_gauss(mean_collisions: WP, distribution_param: WP) -> int:
    if mean_collisions == 0: return 0
    if distribution_param == 0: return int(round(mean_collisions))
    sigma = mean_collisions * distribution_param
    if sigma == 0: return int(round(mean_collisions))
    
    varied_value = random.gauss(mean_collisions, sigma)
    return max(0, int(round(varied_value)))

def vary_energies_gauss(mean_energy: WP, distribution_width: WP) -> WP:
    if distribution_width == 0: return mean_energy
    sigma = mean_energy * distribution_width
    if sigma == 0: return mean_energy
    return random.gauss(mean_energy, sigma)


if __name__ == '__main__':
    print("Testing Box-Muller implementations...")

    # Test vary_collisions_py
    mean_coll = 10.0
    dist_param = 0.2
    print(f"\nTesting vary_collisions_py (mean={mean_coll}, dist_param={dist_param}):")
    for _ in range(5):
        print(f"  Generated collisions: {vary_collisions_py(mean_coll, dist_param)}")

    # Test vary_energies_py
    mean_eng = 70.0
    dist_width = 0.1
    print(f"\nTesting vary_energies_py (mean={mean_eng}, dist_width={dist_width}):")
    for _ in range(5):
        print(f"  Generated energy: {vary_energies_py(mean_eng, dist_width):.2f} eV")

    # Test Pythonic alternatives
    print(f"\nTesting vary_collisions_gauss (mean={mean_coll}, dist_param={dist_param}):")
    for _ in range(5):
        print(f"  Generated collisions (gauss): {vary_collisions_gauss(mean_coll, dist_param)}")

    print(f"\nTesting vary_energies_gauss (mean={mean_eng}, dist_width={dist_width}):")
    for _ in range(5):
        print(f"  Generated energy (gauss): {vary_energies_gauss(mean_eng, dist_width):.2f} eV")
    
    # Test edge case for vary_collisions_py: mean_collisions = 0
    print(f"\nTesting vary_collisions_py (mean=0, dist_param=0.1): {vary_collisions_py(0.0, 0.1)}")
    # Test edge case for vary_collisions_py: distribution_param = 0
    print(f"Testing vary_collisions_py (mean=10, dist_param=0): {vary_collisions_py(10.0, 0.0)}")

    # Test edge case for vary_energies_py: mean_energy = 0
    print(f"Testing vary_energies_py (mean=0, dist_width=0.1): {vary_energies_py(0.0, 0.1):.2f}")
    # Test edge case for vary_energies_py: distribution_width = 0
    print(f"Testing vary_energies_py (mean=70, dist_width=0): {vary_energies_py(70.0, 0.0):.2f}")
```
