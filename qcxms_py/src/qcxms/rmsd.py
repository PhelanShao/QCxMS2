import math
from typing import Tuple, Optional, List
import numpy as np

try:
    from .data import WP # Assuming WP is float or similar
except ImportError:
    WP = float # Fallback for standalone execution

DP = WP # Use WP consistent with other modules for precision type

def _pythag_py(a: DP, b: DP) -> DP:
    """
    Calculates sqrt(a^2 + b^2) robustly.
    Equivalent to math.hypot(a,b) available in Python 3.8+
    """
    return math.hypot(a,b)

def _rotation_matrix_from_quaternion_py(q: np.ndarray) -> np.ndarray:
    """
    Constructs rotation matrix U from quaternion q.
    q is assumed to be a NumPy array of shape (4,).
    U will be a 3x3 NumPy array.
    """
    q0, q1, q2, q3 = q[0], q[1], q[2], q[3]

    # Coefficients from Fortran code
    # b0, b1, b2, b3 = 2.0*q0, 2.0*q1, 2.0*q2, 2.0*q3
    # q00 = b0*q0 - 1.0; q01 = b0*q1; q02 = b0*q2; q03 = b0*q3
    # q11 = b1*q1; q12 = b1*q2; q13 = b1*q3
    # q22 = b2*q2; q23 = b2*q3
    # q33 = b3*q3

    # More standard formulation (e.g., https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Quaternion-derived_rotation_matrix)
    # Assuming q is normalized (q0^2 + q1^2 + q2^2 + q3^2 = 1)
    # If not, it should be normalized before this step.
    # The Kabsch algorithm provides a normalized quaternion eigenvector.

    U = np.zeros((3, 3), dtype=DP)

    U[0, 0] = q0**2 + q1**2 - q2**2 - q3**2 # Fortran: q00+q11 if q00 = 2q0^2-1, q11=2q1^2 -> 2(q0^2+q1^2)-1
                                         # Wikipedia: 1 - 2*(q2^2 + q3^2) or q0^2+q1^2-q2^2-q3^2 (if q0 is scalar part)
                                         # The Fortran S matrix construction implies q0 is the scalar part.
    U[0, 1] = 2.0 * (q1*q2 - q0*q3)       # Fortran: q12 - q03
    U[0, 2] = 2.0 * (q1*q3 + q0*q2)       # Fortran: q13 + q02

    U[1, 0] = 2.0 * (q1*q2 + q0*q3)       # Fortran: q12 + q03
    U[1, 1] = q0**2 - q1**2 + q2**2 - q3**2 # Fortran: q00 + q22
    U[1, 2] = 2.0 * (q2*q3 - q0*q1)       # Fortran: q23 - q01

    U[2, 0] = 2.0 * (q1*q3 - q0*q2)       # Fortran: q13 - q02
    U[2, 1] = 2.0 * (q2*q3 + q0*q1)       # Fortran: q23 + q01
    U[2, 2] = q0**2 - q1**2 - q2**2 + q3**2 # Fortran: q00 + q33
    
    # The Fortran formulas for U seem to map to a common convention if b_i = 2q_i and q00 = 2q0^2-1.
    # U(1,1) = 2q0^2-1 + 2q1^2 = 2(q0^2+q1^2)-1. This matches if sum(q_i^2)=1 -> q2^2+q3^2 = 1 - (q0^2+q1^2)
    # then 1 - 2(q2^2+q3^2) = 1 - 2(1 - q0^2 - q1^2) = 1 - 2 + 2(q0^2+q1^2) = 2(q0^2+q1^2)-1. Yes.
    
    return U


def get_rmsd_for_coords_py(
    coord1: Union[List[List[DP]], np.ndarray], 
    coord2: Union[List[List[DP]], np.ndarray], 
    mask: Optional[Union[List[bool], np.ndarray]] = None,
    calculate_gradient: bool = False 
) -> Tuple[DP, Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Calculates the least-square RMSD between two coordinate sets using Kabsch algorithm (quaternion based).

    Args:
        coord1: Reference coordinates (Nx3 or 3xN list/array).
        coord2: Coordinates to compare (Nx3 or 3xN list/array).
        mask: Optional boolean mask to select atoms for RMSD calculation.
        calculate_gradient: If True, also computes and returns the gradient.

    Returns:
        A tuple containing:
        - rmsd (float): The calculated RMSD.
        - gradient (np.ndarray, optional): Gradient of RMSD w.r.t. coord1 (if requested). Shape (N,3).
        - transformation_matrix (np.ndarray, optional): 3x3 rotation matrix to align coord2 to coord1 (if requested or gradient requested).
    """
    
    x_orig = np.asarray(coord1, dtype=DP)
    y_orig = np.asarray(coord2, dtype=DP)

    # Ensure coordinates are (num_atoms, 3)
    if x_orig.shape[0] == 3 and x_orig.shape[1] != 3: x_orig = x_orig.T
    if y_orig.shape[0] == 3 and y_orig.shape[1] != 3: y_orig = y_orig.T
    
    if x_orig.shape[1] != 3 or y_orig.shape[1] != 3:
        raise ValueError("Coordinate arrays must have 3 columns (x,y,z).")
    
    num_atoms_orig = x_orig.shape[0]
    if y_orig.shape[0] != num_atoms_orig:
        raise ValueError("Coordinate arrays must have the same number of atoms.")

    x_working: np.ndarray
    y_working: np.ndarray

    if mask is not None:
        mask_np = np.asarray(mask, dtype=bool)
        if mask_np.shape[0] != num_atoms_orig:
            raise ValueError("Mask must have the same length as the number of atoms.")
        x_working = x_orig[mask_np]
        y_working = y_orig[mask_np]
        if x_working.shape[0] == 0: # No atoms selected by mask
            return 0.0, (np.zeros_like(x_orig) if calculate_gradient else None), (np.eye(3, dtype=DP) if calculate_gradient else None)
    else:
        x_working = np.copy(x_orig) # Use copies to avoid modifying originals
        y_working = np.copy(y_orig)

    num_atoms_active = x_working.shape[0]

    # 1. Calculate centroids
    x_center = np.mean(x_working, axis=0)
    y_center = np.mean(y_working, axis=0)

    # 2. Center coordinates
    x_c = x_working - x_center
    y_c = y_working - y_center

    # Calculate norms of centered coordinates (E0 in some notations)
    x_norm_sq = np.sum(x_c**2) # sum of squares of elements
    y_norm_sq = np.sum(y_c**2)

    # 3. Covariance matrix R (often denoted C or H)
    # R_ij = sum_k (x_ki * y_kj). In NumPy: R = x_c.T @ y_c
    R_matrix = np.dot(x_c.T, y_c) # (3,N) @ (N,3) -> (3,3)

    # 4. Construct the 4x4 matrix S (related to quaternion)
    S = np.zeros((4, 4), dtype=DP)
    S[0, 0] = R_matrix[0, 0] + R_matrix[1, 1] + R_matrix[2, 2]
    S[0, 1] = R_matrix[1, 2] - R_matrix[2, 1]
    S[0, 2] = R_matrix[2, 0] - R_matrix[0, 2]
    S[0, 3] = R_matrix[0, 1] - R_matrix[1, 0]

    S[1, 0] = S[0, 1]
    S[1, 1] = R_matrix[0, 0] - R_matrix[1, 1] - R_matrix[2, 2]
    S[1, 2] = R_matrix[0, 1] + R_matrix[1, 0]
    S[1, 3] = R_matrix[0, 2] + R_matrix[2, 0]

    S[2, 0] = S[0, 2]
    S[2, 1] = S[1, 2]
    S[2, 2] = -R_matrix[0, 0] + R_matrix[1, 1] - R_matrix[2, 2]
    S[2, 3] = R_matrix[1, 2] + R_matrix[2, 1]

    S[3, 0] = S[0, 3]
    S[3, 1] = S[1, 3]
    S[3, 2] = S[2, 3]
    S[3, 3] = -R_matrix[0, 0] - R_matrix[1, 1] + R_matrix[2, 2]

    # 5. Eigenvalue problem for S
    # Fortran dstmev finds the largest eigenvalue. np.linalg.eigh sorts them.
    eigenvalues, eigenvectors = np.linalg.eigh(S)
    lambda_max = eigenvalues[-1] # Largest eigenvalue
    q_optimal = eigenvectors[:, -1] # Corresponding eigenvector (quaternion)
    
    # Ensure quaternion is normalized (should be by eigh for symmetric matrix)
    # q_optimal /= np.linalg.norm(q_optimal) # Usually not needed

    # 6. Calculate RMSD
    # RMSD^2 = (x_norm_sq + y_norm_sq - 2.0 * lambda_max) / num_atoms_active
    # Add small epsilon to prevent sqrt of tiny negative due to precision
    rmsd_sq_val = (x_norm_sq + y_norm_sq - 2.0 * lambda_max) / num_atoms_active
    rmsd_val = math.sqrt(max(0.0, rmsd_sq_val))

    # 7. Calculate rotation matrix and gradient if requested
    rot_matrix: Optional[np.ndarray] = None
    gradient: Optional[np.ndarray] = None

    if calculate_gradient: # Also implies rot_matrix is needed
        rot_matrix = _rotation_matrix_from_quaternion_py(q_optimal)
        
        grad_full = np.zeros_like(x_orig, dtype=DP) # Initialize gradient for all atoms
        
        # y_c_rotated = y_c @ rot_matrix # y_c is (N,3), rot_matrix is (3,3) -> (N,3)
                                         # Fortran: tmp(:) = matmul(transpose(Rmatrix), y(:, i))
                                         # This means Rmatrix.T @ y_c[i,:].T which is (y_c[i,:] @ Rmatrix).T
                                         # So, y_c_rotated_i = y_c[i,:] @ rot_matrix
        
        # Denominator for gradient, ensure not zero
        denom = max(np.finfo(DP).eps, rmsd_val * num_atoms_active)
        
        # Apply to active atoms
        grad_active = (x_c - y_c @ rot_matrix) / denom
        
        if mask is not None:
            grad_full[mask_np, :] = grad_active
        else:
            grad_full = grad_active
        gradient = grad_full
        
    return rmsd_val, gradient, rot_matrix


if __name__ == '__main__':
    print("Testing RMSD calculations...")

    # Example coordinates (N, 3)
    coords1 = np.array([[0.0, 0.0, 0.0], 
                        [1.0, 0.0, 0.0], 
                        [0.0, 1.0, 0.0]], dtype=DP)
    
    coords2 = np.array([[0.0, 0.0, 0.0], 
                        [1.0, 0.0, 0.0], 
                        [0.0, 1.0, 0.0]], dtype=DP) # Identical

    coords3 = np.array([[0.0, 0.0, 0.0],  # Translated version of coords1
                        [1.0, 0.0, 1.0], 
                        [0.0, 1.0, 1.0]], dtype=DP)
    
    # Rotated version of coords1 (90 deg around z-axis)
    # x' = x*cos(a) - y*sin(a)
    # y' = x*sin(a) + y*cos(a)
    # For 90 deg: x' = -y, y' = x
    coords4 = np.array([[0.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                        [-1.0, 0.0, 0.0]], dtype=DP)

    rmsd12, grad12, rot12 = get_rmsd_for_coords_py(coords1, coords2, calculate_gradient=True)
    print(f"\nRMSD Coords1 vs Coords2 (identical): {rmsd12:.6f}") # Expected: 0.0
    if grad12 is not None: print(f"Gradient12:\n{grad12}")
    if rot12 is not None: print(f"Rotation12:\n{rot12}") # Expected: Identity matrix (approx)

    rmsd13, grad13, rot13 = get_rmsd_for_coords_py(coords1, coords3, calculate_gradient=True)
    print(f"\nRMSD Coords1 vs Coords3 (translated): {rmsd13:.6f}") # Expected: 0.0 (after centering)
    if grad13 is not None: print(f"Gradient13:\n{grad13}")
    if rot13 is not None: print(f"Rotation13:\n{rot13}") # Expected: Identity matrix (approx)

    rmsd14, grad14, rot14 = get_rmsd_for_coords_py(coords1, coords4, calculate_gradient=True)
    print(f"\nRMSD Coords1 vs Coords4 (rotated): {rmsd14:.6f}") # Expected: 0.0 (after rotation)
    if grad14 is not None: print(f"Gradient14:\n{grad14}")
    if rot14 is not None: print(f"Rotation14:\n{rot14}") # Expected: Rotation matrix for -90 deg around z or inverse
    
    # Test with mask
    mask_test = np.array([True, True, False], dtype=bool) # Use only first two atoms
    rmsd_mask, grad_mask, rot_mask = get_rmsd_for_coords_py(coords1, coords4, mask=mask_test, calculate_gradient=True)
    print(f"\nRMSD Coords1 vs Coords4 (masked, first 2 atoms): {rmsd_mask:.6f}")
    if grad_mask is not None: print(f"Gradient_mask:\n{grad_mask}") # Gradient should be zero for masked atoms
    if rot_mask is not None: print(f"Rotation_mask:\n{rot_mask}")

    # Test slightly perturbed coordinates
    coords5 = coords1 + np.random.rand(*coords1.shape) * 0.1
    rmsd15, grad15, rot15 = get_rmsd_for_coords_py(coords1, coords5, calculate_gradient=True)
    print(f"\nRMSD Coords1 vs Coords5 (perturbed): {rmsd15:.6f}")
    if grad15 is not None: print(f"Gradient15:\n{grad15}")
    if rot15 is not None: print(f"Rotation15:\n{rot15}")

    # Verify gradient (simple check: RMSD should decrease if moving slightly along -gradient)
    # This is a more involved test for numerical gradients.
    # For now, just checking if it runs.

```
