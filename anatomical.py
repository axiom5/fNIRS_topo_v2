# anatomical.py
import os
import re
import sys
import time
import numpy as np
import nibabel as nib # For loading NIfTI files
from scipy.spatial import cKDTree # For nearest neighbor search

# --- Constants ---
ATLAS_FOLDER = "Brodmann_MM" # Folder containing atlas files (relative path assumed for now)
ATLAS_FILENAME = "Brodmann_Mai_Matajnik.nii.gz"
LABEL_FILENAME = "balabels_lr.txt"
ATLAS_FILEPATH = os.path.join(ATLAS_FOLDER, ATLAS_FILENAME)
LABEL_FILEPATH = os.path.join(ATLAS_FOLDER, LABEL_FILENAME)

# Define your coarse Broadmann Area groupings
# ADJUST THIS MAPPING AS NEEDED FOR YOUR ANALYSIS
BA_TO_COARSE_REGION = {
    # PFC Subregions (Unique Assignments)
    9: 'DLPFC', 46: 'DLPFC',             # Dorsolateral
    44: 'VLPFC/IFG', 45: 'VLPFC/IFG', 47: 'VLPFC/IFG', # Ventrolateral / IFG
    11: 'OFC', 12: 'OFC',                 # Orbitofrontal (BA12 less common)
    25: 'VMPFC',                          # Ventromedial (includes subgenual)
    10: 'Frontopolar',                             # Frontopolar
    8: 'Dorsal PFC',                      # Dorsal PFC (near premotor/DLPFC)
    # Cingulate
    24: 'ACC', 32: 'ACC', 33: 'ACC',      # Anterior Cingulate
    23: 'Post Cingulate', 31: 'Post Cingulate', # Posterior Cingulate
    # Motor/Premotor
    4: 'Motor', # 'Premotor', # Changed to Premotor based on original script
    6: 'Premotor',
    # Somatosensory
    1: 'Somatosensory', 2: 'Somatosensory', 3: 'Somatosensory', # Assuming 3a/3b combined
    5: 'Somatosensory',                   # Somatosensory Association
    # Parietal (excluding primary somatosensory)
    7: 'Parietal', 39: 'Parietal', 40: 'Parietal',
    # Temporal
    20: 'Temporal', 21: 'Temporal', 22: 'Temporal', 37: 'Temporal', 38: 'Temporal',
    41: 'Temporal (Aud)', 42: 'Temporal (Aud)', # Primary/Secondary Auditory focus
    # Occipital
    17: 'Occipital', 18: 'Occipital', 19: 'Occipital',
    # Other Limbic/Insula/Subcortical related (using BA approx)
    26: 'Limbic', 29: 'Limbic', 30: 'Limbic', # Retrosplenial/Limbic
    34: 'Limbic', 35: 'Limbic', 36: 'Limbic', # Entorhinal/Perirhinal
    52: 'Temporal (Ins)',                 # Parainsular (often grouped w/ Temporal/Auditory/Insula)
    # Add BA13 if needed based on your atlas labels:
    # 13: 'OFC/Insula', # Example grouping if BA13 is present
}

UNKNOWN_REGION_NAME = "Unknown"
OUTSIDE_BRAIN_NAME = "Outside" # Used if point falls outside KDTree range, though lookup uses nearest valid

# --- Atlas and Label Loading ---

def parse_label_file(filepath):
    """
    Parses the AFNI label file (like balabels_lr.txt).

    Args:
        filepath (str): Path to the label file.

    Returns:
        dict: Mapping {nifti_value (int): descriptive_label (str)}, or None on error.
    """
    print(f"Parsing label file: {filepath}...")
    label_map = {}
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'): continue # Skip empty/comment lines
                parts = line.split(None, 1) # Split on first whitespace
                if len(parts) == 2:
                    try:
                        value = int(parts[0])
                        label = parts[1].strip()
                        label_map[value] = label
                    except ValueError:
                        print(f"  Warning: Could not parse integer value on line: '{line}'. Skipping.")
                else:
                     print(f"  Warning: Unexpected format on line: '{line}'. Skipping.")
        if not label_map:
            print(f"  Error: No valid labels found in {filepath}.", file=sys.stderr)
            return None
        print(f"  Successfully parsed {len(label_map)} labels.")
        return label_map
    except FileNotFoundError:
        print(f"Error: Label file not found at {filepath}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred parsing label file {filepath}: {e}", file=sys.stderr)
        return None

def create_value_to_region_map(nifti_value_to_desc_map, ba_to_coarse_map):
    """
    Creates a direct map from NIfTI integer value to coarse region name.

    Args:
        nifti_value_to_desc_map (dict): Map from {nifti_value: descriptive_label}.
        ba_to_coarse_map (dict): Map from {ba_number: coarse_region_name}.

    Returns:
        dict: Mapping {nifti_value (int): coarse_region_name (str)}. Returns empty dict on critical error.
    """
    print("Creating NIfTI value to coarse region map...")
    value_to_region = {}
    if nifti_value_to_desc_map is None:
        print("Error: Cannot create map, input nifti_value_to_desc_map is None.", file=sys.stderr)
        return {}

    ba_pattern = re.compile(r'BA(\d+)', re.IGNORECASE) # Pattern to find BA number
    found_ba_count = 0
    mapped_to_coarse = 0

    for value, desc_label in nifti_value_to_desc_map.items():
        match = ba_pattern.search(desc_label)
        region_assigned = UNKNOWN_REGION_NAME # Default
        if match:
            try:
                ba_num = int(match.group(1))
                found_ba_count += 1
                region_assigned = ba_to_coarse_map.get(ba_num, UNKNOWN_REGION_NAME)
                if region_assigned != UNKNOWN_REGION_NAME:
                    mapped_to_coarse += 1
            except ValueError:
                 print(f"  Warning: Could not convert extracted BA '{match.group(1)}' to int for label '{desc_label}'. Assigning Unknown.")
        # else: # Label doesn't contain 'BA##', keep region_assigned as Unknown
              # Add specific non-BA label handling here if needed
              # e.g., if 'Entorhinal' in desc_label: region_assigned = 'Limbic'

        value_to_region[value] = region_assigned

    print(f"  Processed {len(nifti_value_to_desc_map)} descriptive labels.")
    print(f"  Found {found_ba_count} labels containing 'BA##'.")
    print(f"  Mapped {mapped_to_coarse} BA labels to a coarse region.")
    print(f"  Final map contains {len(value_to_region)} NIfTI value -> coarse region entries.")
    return value_to_region

def load_atlas(filepath):
    """
    Loads the NIfTI atlas image data and affine matrices.

    Args:
        filepath (str): Path to the NIfTI atlas file.

    Returns:
        tuple: (atlas_data, atlas_affine, inverse_atlas_affine) or (None, None, None) on error.
               atlas_data is returned as int16 for memory efficiency.
    """
    print(f"Loading NIfTI atlas: {filepath}...")
    try:
        img = nib.load(filepath)
        # Ensure data type is integer, use int16 for potential memory savings
        atlas_data = img.get_fdata(dtype=np.float32).astype(np.int16)
        atlas_affine = img.affine
        inverse_atlas_affine = np.linalg.inv(atlas_affine)
        print(f"  Atlas loaded. Shape: {atlas_data.shape}, Affine:\n{atlas_affine}")
        return atlas_data, atlas_affine, inverse_atlas_affine
    except FileNotFoundError:
        print(f"Error: Atlas file not found at {filepath}", file=sys.stderr)
        return None, None, None
    except Exception as e:
        print(f"An unexpected error occurred loading atlas {filepath}: {e}", file=sys.stderr)
        return None, None, None

# --- Coordinate Transformation ---

def mni_to_voxel(mni_coords, inverse_affine):
    """
    Converts MNI coordinates to voxel coordinates using the inverse affine.

    Args:
        mni_coords (np.ndarray): Array of MNI coordinates (Nx3 or 1x3).
        inverse_affine (np.ndarray): 4x4 inverse affine matrix of the atlas.

    Returns:
        np.ndarray: Array of voxel coordinates (Nx3, dtype=int).
    """
    mni_coords = np.asarray(mni_coords)
    is_single = mni_coords.ndim == 1
    if is_single:
        mni_coords = mni_coords.reshape(1, 3)

    # Add homogeneous coordinate (1)
    mni_homogeneous = np.hstack((mni_coords, np.ones((mni_coords.shape[0], 1))))
    # Apply inverse affine transformation
    voxel_homogeneous = (inverse_affine @ mni_homogeneous.T).T
    # Get voxel indices by rounding and converting to int
    voxel_indices = np.rint(voxel_homogeneous[:, :3]).astype(int)

    return voxel_indices[0] if is_single else voxel_indices

def voxel_to_mni(voxel_coords, affine):
    """
    Converts voxel coordinates to MNI coordinates using the forward affine.

    Args:
        voxel_coords (np.ndarray): Array of voxel coordinates (Nx3 or 1x3).
        affine (np.ndarray): 4x4 forward affine matrix of the atlas.

    Returns:
        np.ndarray: Array of MNI coordinates (Nx3).
    """
    voxel_coords = np.asarray(voxel_coords)
    is_single = voxel_coords.ndim == 1
    if is_single:
        voxel_coords = voxel_coords.reshape(1, 3)

    # Add homogeneous coordinate (1)
    voxel_homogeneous = np.hstack((voxel_coords, np.ones((voxel_coords.shape[0], 1))))
    # Apply affine transformation
    mni_homogeneous = (affine @ voxel_homogeneous.T).T
    mni_coords = mni_homogeneous[:, :3]

    return mni_coords[0] if is_single else mni_coords


# --- KDTree and Region Lookup ---

def build_atlas_kdtree(atlas_data, atlas_affine):
    """
    Builds a KDTree from the MNI coordinates of non-zero atlas voxels.

    Args:
        atlas_data (np.ndarray): Loaded 3D atlas data array.
        atlas_affine (np.ndarray): 4x4 forward affine matrix of the atlas.

    Returns:
        tuple: (kdtree, index_to_voxel_map) or (None, None) on error.
               kdtree: scipy.spatial.cKDTree object.
               index_to_voxel_map: dict mapping KDTree index to (vx, vy, vz) tuple.
    """
    print("\n--- Pre-calculating Atlas Brain Voxel MNI Coords & KDTree ---")
    start_time = time.time()
    try:
        # Find voxel indices of all non-zero labels (brain voxels)
        brain_voxel_indices = np.array(np.where(atlas_data != 0)).T
        if brain_voxel_indices.shape[0] == 0:
            print("Error: No non-zero voxels found in the atlas data!", file=sys.stderr)
            return None, None
        print(f"  Found {brain_voxel_indices.shape[0]} non-zero voxels in atlas.")

        # Convert these voxel indices to MNI coordinates
        atlas_brain_mni_coords = voxel_to_mni(brain_voxel_indices, atlas_affine)

        # Build KDTree from the MNI coordinates of brain voxels
        print("  Building KDTree from atlas brain MNI coordinates...")
        atlas_brain_kdtree = cKDTree(atlas_brain_mni_coords)

        # Create a map from the KDTree's internal index to the original voxel index tuple
        kdtree_index_to_voxel_index = {i: tuple(brain_voxel_indices[i])
                                      for i in range(brain_voxel_indices.shape[0])}

        print(f"  Atlas brain KDTree built in {time.time() - start_time:.2f} seconds.")
        return atlas_brain_kdtree, kdtree_index_to_voxel_index

    except Exception as e:
        print(f"Error building atlas KDTree: {e}", file=sys.stderr)
        return None, None


def get_regions_for_mni_coords(mni_coords_channels,
                               atlas_brain_kdtree,
                               kdtree_index_to_voxel_index,
                               atlas_data,
                               value_to_region_map):
    """
    Gets coarse region name(s) for channel MNI coordinate(s) by finding the
    nearest labeled *brain* voxel using the precomputed KDTree.

    Args:
        mni_coords_channels (np.ndarray): MNI coordinates of fNIRS channels (Nx3 or 1x3).
        atlas_brain_kdtree (scipy.spatial.cKDTree): KDTree built from MNI coordinates
                                                    of non-zero atlas voxels.
        kdtree_index_to_voxel_index (dict): Maps KDTree query result index to the
                                            original (vx, vy, vz) voxel index tuple.
        atlas_data (np.ndarray): Loaded 3D atlas data array (integers).
        value_to_region_map (dict): Map from {nifti_value: coarse_region_name}.

    Returns:
        list or str: List of coarse region names (str) corresponding to input coords,
                     or a single string if only one coordinate was input.
                     Returns UNKNOWN_REGION_NAME on error or if mapping fails.
    """
    # print("  DEBUG: Looking up regions using Nearest Brain Voxel KDTree...") # Optional verbose debug
    mni_coords_channels = np.asarray(mni_coords_channels)
    is_single_coord = mni_coords_channels.ndim == 1
    if is_single_coord:
        mni_coords_channels = mni_coords_channels.reshape(1, 3)

    if atlas_brain_kdtree is None or kdtree_index_to_voxel_index is None:
        print("Error: KDTree or index map is not available for region lookup.", file=sys.stderr)
        return [UNKNOWN_REGION_NAME] * mni_coords_channels.shape[0] if not is_single_coord else UNKNOWN_REGION_NAME

    region_names = []
    try:
        # Query the KDTree: find the index (in the KDTree's data) of the nearest brain voxel
        distances, kdtree_indices = atlas_brain_kdtree.query(mni_coords_channels, k=1)

        # --- Optional Debug Print ---
        # if kdtree_indices.shape[0] > 0:
        #     print(f"  DEBUG: Nearest brain voxel KDTree indices (sample): {kdtree_indices[:min(5, len(kdtree_indices))]}")
        #     print(f"  DEBUG: Distances to nearest brain voxel (sample): {distances[:min(5, len(distances))]}")
        # --- End Optional Debug ---

        for i, kdtree_idx in enumerate(kdtree_indices):
            # Handle potential edge case where query returns index outside map bounds
            # (Shouldn't happen with valid KDTree and query)
            if kdtree_idx >= len(kdtree_index_to_voxel_index):
                 print(f"  Error: KDTree query returned out-of-bounds index {kdtree_idx}. Assigning Unknown.")
                 region_names.append(UNKNOWN_REGION_NAME)
                 continue

            # Get the original (vx, vy, vz) voxel index tuple using the map
            voxel_idx_tuple = kdtree_index_to_voxel_index.get(kdtree_idx)

            if voxel_idx_tuple:
                vx, vy, vz = voxel_idx_tuple
                # Check bounds before accessing atlas_data (important!)
                if 0 <= vx < atlas_data.shape[0] and 0 <= vy < atlas_data.shape[1] and 0 <= vz < atlas_data.shape[2]:
                    nifti_value = atlas_data[vx, vy, vz]
                    # --- Optional Debug Print ---
                    # if i < 5: print(f"  DEBUG: Chan {i} -> MNI {mni_coords_channels[i]} -> Nearest Voxel [{vx},{vy},{vz}] (Dist: {distances[i]:.2f}) -> Atlas Value: {nifti_value}")
                    # --- End Optional Debug ---

                    if nifti_value == 0:
                        # This *shouldn't* happen if KDTree only contains non-zero voxels, but safety check
                        print(f"  Warning: Nearest brain voxel {voxel_idx_tuple} has value 0! Assigning Unknown.")
                        region_names.append(UNKNOWN_REGION_NAME)
                    else:
                        region_name = value_to_region_map.get(nifti_value, UNKNOWN_REGION_NAME)
                        region_names.append(region_name)
                else:
                    print(f"  Error: Mapped voxel index {voxel_idx_tuple} is out of atlas bounds {atlas_data.shape}. Assigning Unknown.")
                    region_names.append(UNKNOWN_REGION_NAME)
            else:
                # Should not happen if KDTree query is successful and map is correct
                print(f"  Error: Could not find voxel index for KDTree result index {kdtree_idx}! Assigning Unknown.")
                region_names.append(UNKNOWN_REGION_NAME)

    except Exception as e:
        print(f"Error during KDTree region lookup: {e}", file=sys.stderr)
        # Fallback to Unknown for all requested points
        region_names = [UNKNOWN_REGION_NAME] * mni_coords_channels.shape[0]

    return region_names[0] if is_single_coord and len(region_names) == 1 else region_names