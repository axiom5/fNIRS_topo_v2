# registration.py
import numpy as np
from scipy.spatial.transform import Rotation as R # For registration
import sys # For stderr

def calculate_registration(source_points, target_points):
    """
    Calculates the rigid transformation (rotation + translation) to align
    source_points to target_points using Kabsch algorithm (via SVD).

    Args:
        source_points (np.ndarray): Array of source points (Nx3). Usually digitizer fiducials.
        target_points (np.ndarray): Array of target points (Nx3). Usually MNI fiducials.

    Returns:
        np.ndarray or None: 4x4 transformation matrix (source -> target) or None if failed.
    """
    source_points = np.asarray(source_points)
    target_points = np.asarray(target_points)

    if source_points.shape != target_points.shape or source_points.shape[0] < 3 or source_points.shape[1] != 3:
        print("Error: Input point arrays must have the same shape (Nx3, N>=3).", file=sys.stderr)
        return None

    # Calculate centroids
    source_centroid = np.mean(source_points, axis=0)
    target_centroid = np.mean(target_points, axis=0)

    # Center the points
    source_centered = source_points - source_centroid
    target_centered = target_points - target_centroid

    # --- SciPy's Rotation.align_vectors implementation ---
    # This uses SVD internally and handles reflections.
    try:
        rotation, rmsd = R.align_vectors(source_centered, target_centered) # Note: Order matters! (source, target)
        print(f"  Registration RMSD: {rmsd:.4f} mm") # Assuming input units are mm

        # Construct the 4x4 transformation matrix
        T = np.identity(4)
        T[:3, :3] = rotation.as_matrix()
        # Translation: target_centroid = R * source_centroid + t
        # => t = target_centroid - R * source_centroid
        T[:3, 3] = target_centroid - rotation.apply(source_centroid)
        return T

    except Exception as e:
        print(f"Error during registration calculation: {e}", file=sys.stderr)
        return None


def apply_transform(points_3d, T):
    """
    Applies a 4x4 transformation matrix to 3D points.

    Args:
        points_3d (np.ndarray): Array of 3D points (Nx3 or 1x3).
        T (np.ndarray or None): 4x4 transformation matrix. If None, returns original points.

    Returns:
        np.ndarray: Transformed points (Nx3 or 1x3). Returns original if T is None.
    """
    if T is None:
        # print("Warning: Transformation matrix is None. Returning original points.", file=sys.stderr) # Optional warning
        return np.asarray(points_3d) # Ensure output is numpy array

    points_3d = np.asarray(points_3d)
    is_single = points_3d.ndim == 1
    if is_single:
        points_3d = points_3d.reshape(1, 3)

    if points_3d.shape[1] != 3:
        raise ValueError("Input points must be 3D (Nx3).")

    # Convert to homogeneous coordinates
    points_homogeneous = np.hstack((points_3d, np.ones((points_3d.shape[0], 1))))

    # Apply transformation: T @ points_homogeneous.T
    # Transpose points for multiplication, then transpose result back
    points_transformed_homogeneous = (T @ points_homogeneous.T).T

    # Extract 3D coordinates
    transformed_3d = points_transformed_homogeneous[:, :3]

    return transformed_3d[0] if is_single else transformed_3d

def extract_fiducials(required_names, all_ch_pos, case_sensitive=False):
    """
    Extracts specific fiducial coordinates from the loaded digitizer positions.

    Args:
        required_names (list or set): List/set of fiducial names to extract (e.g., ['LPA', 'Nz', 'RPA']).
        all_ch_pos (dict): Dictionary of all loaded digitizer positions {name: [x,y,z]}.
        case_sensitive (bool): Whether the matching of names should be case-sensitive. Default is False.

    Returns:
        tuple: (fiducials_dict, missing_names)
               fiducials_dict (dict): {name: [x,y,z]} for found fiducials.
               missing_names (list): List of names from required_names not found in all_ch_pos.
    """
    found_fiducials = {}
    missing = []
    required_set = set(required_names) # Ensure unique names

    # Create a mapping from lowercase/uppercase name to original name if needed
    lookup_map = {}
    if not case_sensitive:
        for name in all_ch_pos.keys():
            lookup_map[name.lower()] = name
    else:
        lookup_map = {name: name for name in all_ch_pos.keys()}

    for req_name in required_set:
        lookup_key = req_name.lower() if not case_sensitive else req_name
        original_name = lookup_map.get(lookup_key)

        if original_name and original_name in all_ch_pos:
            found_fiducials[req_name] = all_ch_pos[original_name] # Store with the requested name as key
            # print(f"  Found digitizer point for {req_name}: {found_fiducials[req_name]}")
        else:
            missing.append(req_name)
            print(f"  Warning: Required fiducial '{req_name}' not found in digitizer positions.")

    if missing:
        print(f"Error: Missing required fiducials in digitizer data: {missing}", file=sys.stderr)

    return found_fiducials, missing