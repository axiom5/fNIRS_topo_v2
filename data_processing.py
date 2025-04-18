# data_processing.py
import pandas as pd
import numpy as np
import re
import math
import sys
import time
from collections import defaultdict
import warnings

# Import functions from other modules created in this refactoring
from utils import project_point_azimuthal_equidistant
from registration import apply_transform
from anatomical import get_regions_for_mni_coords, UNKNOWN_REGION_NAME, OUTSIDE_BRAIN_NAME

# --- HTA Helper (Internal or could be moved to utils if more general) ---
def _perform_hta(input_data, hta_times_sec, srate, subject_id):
    """Internal helper to perform HTA averaging."""
    # print(f"    Performing HTA for {subject_id}...") # Optional debug
    processed_data = np.array([])
    processed_times = np.array([])
    try:
        if len(hta_times_sec) < 2: raise ValueError("Not enough HTA timestamps (need at least 2).")
        # Assume data starts at t=0 relative to itself for indexing
        # Calculate indices based on timestamps relative to the *first* HTA timestamp
        first_hta_time = hta_times_sec[0]
        hta_indices_relative = np.round((hta_times_sec - first_hta_time) * srate).astype(int)

        # Ensure indices are within the bounds of the input data
        max_idx_data = input_data.shape[0] - 1
        hta_indices_clipped = np.clip(hta_indices_relative, 0, max_idx_data)
        hta_indices_unique = np.unique(hta_indices_clipped)

        if len(hta_indices_unique) < 2: raise ValueError("Not enough unique HTA indices after clipping/rounding.")

        averaged_signal_list = []
        segment_midpoint_times_relative = [] # Times relative to first HTA event

        for i in range(len(hta_indices_unique) - 1):
            start_idx, end_idx = hta_indices_unique[i], hta_indices_unique[i+1]
            if start_idx >= end_idx: continue # Skip empty or invalid segments

            segment_data = input_data[start_idx:end_idx, :]
            if segment_data.shape[0] == 0: continue

            # Handle potential NaNs during averaging
            with warnings.catch_warnings(): # Suppress RuntimeWarning for empty slice mean
                warnings.simplefilter("ignore", category=RuntimeWarning)
                segment_avg = np.nanmean(segment_data, axis=0)

            if np.all(np.isnan(segment_avg)): # Skip if entire segment was NaN
                # print(f"      Skipping HTA segment {i+1} for {subject_id} (all NaNs).") # Debug
                continue

            averaged_signal_list.append(segment_avg)

            # Calculate midpoint time relative to the start of the *data* segment (first HTA time)
            # Use the original (non-clipped) HTA times for accurate midpoint calculation
            start_time_sec = hta_times_sec[i] - first_hta_time
            end_time_sec = hta_times_sec[i+1] - first_hta_time
            segment_midpoint_times_relative.append((start_time_sec + end_time_sec) / 2.0)

        if averaged_signal_list:
            processed_data = np.array(averaged_signal_list)
            processed_times = np.array(segment_midpoint_times_relative)
            # print(f"    HTA successful for {subject_id}: {processed_data.shape[0]} points.") # Debug
        else:
            raise ValueError("No valid HTA segments generated after processing.")

    except ValueError as ve:
        print(f"  Warning: HTA ValueError for {subject_id}: {ve}. Falling back.")
        return None, None # Indicate fallback needed
    except Exception as e_hta:
        print(f"  Warning: Unexpected error during HTA for {subject_id} ({e_hta}). Falling back.")
        return None, None # Indicate fallback needed

    return processed_data, processed_times


# --- Main Data Preparation Function ---

def prepare_trial_visualization_data(
    trial_data,                # Dict {subj_id: pd.DataFrame}
    trial_num,                 # int
    ch_pos,                    # Dict {ch_name: np.array([x,y,z])} from digpts
    chromophore,               # str ("HbO" or "HbR")
    srate,                     # float (original sampling rate)
    downsample_rate,           # int (used if HTA is off/fails)
    do_hta_averaging,          # bool
    hta_df,                    # pd.DataFrame or None (HTA event times)
    hta_col_rename_func,       # function (maps subj_id to HTA column name)
    dsfreq,                    # float (downsampled frequency, srate/downsample_rate)
    ss_detector_numbers,       # set (detector numbers to exclude)
    channel_subgroups_dict,    # dict {group_name: list_of_pairs}
    default_subgroup_name,     # str
    # --- Registration & Anatomical ---
    T_digi_to_template,        # np.ndarray (4x4 transformation matrix)
    value_to_region_map,       # dict {nifti_val: region_name}
    atlas_data,                # np.ndarray (loaded atlas NIfTI data)
    atlas_brain_kdtree,        # scipy.spatial.cKDTree
    kdtree_index_to_voxel_index, # dict {kdtree_idx: (vx, vy, vz)}
    # --- Options ---
    ignore_fiducials_for_center=False # bool (whether to ignore Nz,LPA,RPA etc for head center calc)
    ):
    """
    Processes data for a single trial, calculates geometry, performs anatomical mapping,
    and prepares data structure for JS visualization.
    """
    print(f"\nProcessing Trial {trial_num}, Chromophore {chromophore}...")
    start_time_trial = time.time()

    all_subjects_data_processed = {}
    potential_subject_list = sorted(trial_data.keys())
    common_channel_names = None # List of names like 'S1_D1_HbO'
    selected_channel_info = [] # List of dicts with {'name', 's', 'd', 's_num', 'd_num', 'index_orig'}
    active_sources_trial = set()
    active_detectors_trial = set()
    subgroup_indices_map = {} # {group_name: [idx1, idx2...]} (indices relative to common_channel_names)
    final_subject_list = [] # Subjects successfully processed

    # --- Prepare subgroup info ---
    if channel_subgroups_dict:
        print(f"  Using predefined subgroups: {list(channel_subgroups_dict.keys())}")
        # Convert subgroup S/D pairs to integer sets for efficient lookup
        subgroup_sets = {name: set((int(s), int(d)) for s, d in pairs if isinstance(s, (int, str)) and isinstance(d, (int, str))) # Basic check
                         for name, pairs in channel_subgroups_dict.items()}
    else:
        print("  No channel subgroups defined.")
        subgroup_sets = {}
        default_subgroup_name = "All Channels" # Force default if no subgroups

    # === Identify Common Channels (based on first valid subject) ===
    first_subject_processed = False
    for subject_id in potential_subject_list:
        df = trial_data[subject_id]
        subject_all_channel_names = df.columns.tolist()

        temp_selected_channels = []
        temp_common_names = []
        temp_active_sources = set()
        temp_active_detectors = set()

        channel_pattern = re.compile(r'(S\d+)_(D\d+)_(' + chromophore + r')$') # Compile regex once

        for i, name in enumerate(subject_all_channel_names):
            match = channel_pattern.match(name)
            if match:
                s_name, d_name, _ = match.groups()
                # Ensure corresponding optodes exist in digitizer data (ch_pos)
                if s_name in ch_pos and d_name in ch_pos:
                    try:
                        s_num = int(s_name[1:])
                        d_num = int(d_name[1:])
                    except ValueError:
                        print(f"    Warning: Could not parse S/D number from {name}. Skipping.")
                        continue

                    if d_num not in ss_detector_numbers: # Exclude Short Separation channels
                        ch_data = { 'name': name, 'index_orig': i, # Index in *this subject's* original df
                                    's': s_name, 'd': d_name,
                                    's_num': s_num, 'd_num': d_num }
                        temp_selected_channels.append(ch_data)
                        temp_common_names.append(name)
                        temp_active_sources.add(s_name)
                        temp_active_detectors.add(d_name)

        if temp_selected_channels: # Found valid channels for this subject
            selected_channel_info = temp_selected_channels
            common_channel_names = temp_common_names
            active_sources_trial = temp_active_sources
            active_detectors_trial = temp_active_detectors
            print(f"  Identified {len(common_channel_names)} standard channels based on subject '{subject_id}' and digpts.")
            first_subject_processed = True
            break # Use the first subject with valid channels as the template
        else:
            print(f"  Warning: No standard channels found for {chromophore} matching digpts for subject '{subject_id}'. Checking next subject.")

    if not first_subject_processed:
         print(f"  Error: No subjects found with standard channels for {chromophore} in Trial {trial_num} matching digitizer points. Skipping trial.")
         return None

    # --- Map subgroup names to indices (relative to common_channel_names) ---
    if subgroup_sets:
         for name, sd_set in subgroup_sets.items():
             indices = [idx for idx, ch_info in enumerate(selected_channel_info)
                        if (ch_info['s_num'], ch_info['d_num']) in sd_set]
             if indices:
                  subgroup_indices_map[name] = indices
                  # print(f"    Subgroup '{name}' maps to {len(indices)} channel indices.")
             else:
                  print(f"    Warning: No channels found for predefined subgroup '{name}'.")
    # Validate default subgroup
    if default_subgroup_name != "All Channels" and default_subgroup_name not in subgroup_indices_map:
        print(f"  Warning: Default subgroup '{default_subgroup_name}' not found or empty. Defaulting to 'All Channels'.")
        default_subgroup_name = "All Channels"


    # === Process Each Subject's Data ===
    for subject_id in potential_subject_list:
        df = trial_data[subject_id]
        subject_all_channel_names = df.columns.tolist()
        processing_method = "downsampling" # Default

        # Select data for common channels for *this* subject
        try:
            # Map common_channel_names to their indices in *this subject's* dataframe
            current_subj_indices = [subject_all_channel_names.index(name) for name in common_channel_names]
            signal_data_selected = df.iloc[:, current_subj_indices].values.astype(np.float32) # Ensure float type
        except ValueError as e:
             print(f"  Warning: Channel selection error for {subject_id} (Missing channels like: {e}). Skipping subject.")
             continue
        except IndexError as e:
             print(f"  Warning: Indexing error during channel selection for {subject_id} ({e}). Skipping subject.")
             continue

        # --- HTA or Downsampling Logic ---
        processed_data = np.array([])
        processed_times = np.array([])

        if do_hta_averaging and hta_df is not None:
            hta_col_name = hta_col_rename_func(subject_id)
            if hta_col_name in hta_df.columns:
                hta_times_sec = hta_df[hta_col_name].dropna().values
                if len(hta_times_sec) > 0:
                    processed_data_hta, processed_times_hta = _perform_hta(signal_data_selected, hta_times_sec, srate, subject_id)
                    if processed_data_hta is not None: # HTA was successful
                        processed_data = processed_data_hta
                        processed_times = processed_times_hta
                        processing_method = "HTA"
                    else: # HTA failed, fallback message printed in helper
                         processing_method = "downsampling" # Explicitly set fallback
                else:
                    print(f"  Warning: No HTA times found in column '{hta_col_name}' for {subject_id}. Skipping.")
                    #processing_method = "downsampling"
                    continue
            else:
                print(f"  Warning: HTA column '{hta_col_name}' not found for {subject_id}. Skipping.")
                #processing_method = "downsampling"
                continue

        # Perform downsampling if HTA not done or failed
        if processing_method == "downsampling":
            if signal_data_selected.shape[0] > 0:
                 processed_data = signal_data_selected[::downsample_rate, :]
                 processed_times = np.arange(processed_data.shape[0]) / dsfreq
                 # print(f"    Downsampled data for {subject_id} to {processed_data.shape[0]} points.") # Debug
            else:
                 print(f"  Warning: No data points available for downsampling for subject {subject_id}.")
                 continue # Skip subject if no data

        # Final checks and storage
        if processed_data.size == 0:
             print(f"  Warning: No processed data generated for {subject_id} after {processing_method}. Skipping subject.")
             continue

        # Replace any remaining NaNs (e.g., from nanmean if a whole HTA segment was NaN) with 0.0
        if np.isnan(processed_data).any():
            # print(f"    Note: Replacing NaNs with 0.0 in processed data for {subject_id}.") # Debug
            processed_data = np.nan_to_num(processed_data, nan=0.0)

        if processed_data.ndim != 2 or processed_data.shape[1] != len(common_channel_names):
             print(f"  ERROR: Processed data shape mismatch for {subject_id}! Expected: (*,{len(common_channel_names)}), Got: {processed_data.shape}. Skipping subject.")
             continue

        all_subjects_data_processed[subject_id] = {
            'times': processed_times.tolist(), # Convert numpy arrays to lists for JSON
            'data': processed_data.tolist()
        }
        final_subject_list.append(subject_id)
    # === End of Subject Loop ===

    if not all_subjects_data_processed:
        print(f"  Error: No subjects processed successfully for Trial {trial_num}, {chromophore}.")
        return None
    print(f"  Successfully processed data for {len(final_subject_list)} subjects: {final_subject_list}")

    # --- Geometry Calculations (Based on DIGITIZER coordinates in ch_pos) ---
    print("  Calculating geometry from digitizer coordinates...")
    points_for_center = []
    fiducials = {'Nz', 'Iz', 'LPA', 'RPA'} # Case-insensitive check later if needed

    # Collect points for calculating head center (active optodes and optionally fiducials)
    for name, pos in ch_pos.items():
        is_fiducial = name in fiducials # Simple check, refine if case needed
        is_active_optode = name in active_sources_trial or name in active_detectors_trial
        if is_active_optode or (not ignore_fiducials_for_center and is_fiducial):
             points_for_center.append(pos)

    if not points_for_center:
        print(f"  Error: No points found (active optodes/fiducials) in digitizer data to calculate head center for Trial {trial_num}. Check digpts.txt and channel matching.")
        return None
    head_center_3d_digi = np.mean(np.array(points_for_center), axis=0)
    # print(f"    Head center (Digi): {head_center_3d_digi}") # Debug

    # Calculate midpoints in DIGITIZER space for ALL common channels
    midpoints_3d_dict_digi = {} # {channel_name: [x,y,z]}
    midpoint_radii_digi = []
    final_midpoints_3d_digi_ordered = [] # For JS output, ordered like common_channel_names
    midpoints_3d_mni_ordered = []      # For anatomical lookup, ordered like common_channel_names

    for i, ch_info in enumerate(selected_channel_info):
        ch_name = ch_info['name']
        try:
            s_pos = ch_pos[ch_info['s']]
            d_pos = ch_pos[ch_info['d']]
            midpoint_3d_digi = (s_pos + d_pos) / 2.0
            midpoints_3d_dict_digi[ch_name] = midpoint_3d_digi
            final_midpoints_3d_digi_ordered.append(midpoint_3d_digi.tolist()) # Add for JS output

            # Calculate radius from center
            dist_from_center = np.linalg.norm(midpoint_3d_digi - head_center_3d_digi)
            if dist_from_center > 1e-6: # Avoid zero distances
                midpoint_radii_digi.append(dist_from_center)

        except KeyError as e:
            print(f"  Warning: Optode {e} for channel {ch_name} not found in digitizer positions. Cannot calculate midpoint. Skipping channel for geometry.")
            # Need to handle this missing channel in dependent steps
            final_midpoints_3d_digi_ordered.append([np.nan, np.nan, np.nan]) # Placeholder for JS
            continue # Skip this channel for further processing needing midpoint

    if not midpoints_3d_dict_digi:
        print(f"  Error: No channel midpoints could be calculated in digitizer space. Check optode names and digpts.txt.")
        return None

    head_radius_3d_digi = np.mean(midpoint_radii_digi) if midpoint_radii_digi else 50.0 # Default radius if needed
    # print(f"    Head radius (Digi): {head_radius_3d_digi:.2f}") # Debug


    # --- Anatomical Mapping (Transform Digi Midpoints -> MNI -> Get Regions) ---
    print("  Mapping channels to anatomical regions...")
    # Prepare the DIGITIZER midpoints that were successfully calculated
    valid_midpoint_indices = [i for i, pos in enumerate(final_midpoints_3d_digi_ordered) if not np.isnan(pos[0])]
    midpoints_digi_array_valid = np.array([final_midpoints_3d_digi_ordered[i] for i in valid_midpoint_indices])

    channel_ba_regions = [UNKNOWN_REGION_NAME] * len(common_channel_names) # Initialize with Unknown

    if T_digi_to_template is not None and midpoints_digi_array_valid.shape[0] > 0:
        # Apply transformation (Digi -> MNI/Template)
        midpoints_template_array_valid = apply_transform(midpoints_digi_array_valid, T_digi_to_template)
        # print(f"    Transformed {midpoints_template_array_valid.shape[0]} midpoints to MNI space.") # Debug

        # Get regions using the transformed MNI coordinates
        # Ensure anatomical resources are valid before calling
        if atlas_brain_kdtree and kdtree_index_to_voxel_index and atlas_data is not None and value_to_region_map:
            region_names_list_valid = get_regions_for_mni_coords(
                midpoints_template_array_valid, # Pass the valid MNI coordinates
                atlas_brain_kdtree,
                kdtree_index_to_voxel_index,
                atlas_data,
                value_to_region_map
            )

            # Map results back to the full channel list using valid_midpoint_indices
            if len(region_names_list_valid) == len(valid_midpoint_indices):
                for i, region_name in enumerate(region_names_list_valid):
                    original_index = valid_midpoint_indices[i]
                    channel_ba_regions[original_index] = region_name
            else:
                print(f"  Error: Mismatch between number of valid midpoints ({len(valid_midpoint_indices)}) and returned regions ({len(region_names_list_valid)}). Anatomical mapping incomplete.")
        else:
            print("  Warning: Anatomical resources (KDTree, atlas, map) not fully available. Skipping region lookup.")

    else:
        if T_digi_to_template is None:
            print("  Warning: Registration matrix (T_digi_to_template) is None. Cannot perform anatomical mapping.")
        if midpoints_digi_array_valid.shape[0] == 0:
             print("  Warning: No valid digitizer midpoints available to transform for anatomical mapping.")
        # All channels remain 'Unknown'

    mapped_count = len([r for r in channel_ba_regions if r not in [UNKNOWN_REGION_NAME, OUTSIDE_BRAIN_NAME]])
    print(f"  Mapped {mapped_count} of {len(common_channel_names)} channels to known coarse regions.")
    # Example: Print region for a few channels
    # for i, name in enumerate(common_channel_names[:min(5, len(common_channel_names))]):
    #      print(f"    {name}: Assigned Region '{channel_ba_regions[i]}'")


    # --- 2D Projection (Based on DIGITIZER coordinates) ---
    print("  Calculating 2D projection using digitizer coordinates...")
    # Project midpoints (ALL common channels with valid digi coords) for scale factor calculation
    projected_midpoints_2d_all_digi = {}
    valid_midpoint_names = [common_channel_names[i] for i in valid_midpoint_indices]
    for i, ch_name in enumerate(valid_midpoint_names):
         midpoint_3d_digi = midpoints_digi_array_valid[i] # Use the valid array
         projected_midpoints_2d_all_digi[ch_name] = project_point_azimuthal_equidistant(midpoint_3d_digi, head_center_3d_digi)

    # Project active optodes and fiducials for drawing lines/markers
    projected_optodes_2d_digi = {}
    all_active_optodes = active_sources_trial.union(active_detectors_trial)
    for opt_name in all_active_optodes:
        if opt_name in ch_pos:
            projected_optodes_2d_digi[opt_name] = project_point_azimuthal_equidistant(ch_pos[opt_name], head_center_3d_digi)
        else:
             print(f"  Warning: Active optode {opt_name} used by a channel not found in digitizer data for 2D projection.")

    projected_fiducials_2d_digi = {}
    for fid_name in fiducials: # Use the standard set
        if fid_name in ch_pos:
            projected_fiducials_2d_digi[fid_name] = project_point_azimuthal_equidistant(ch_pos[fid_name], head_center_3d_digi)
        # else: print(f"    Fiducial {fid_name} not found for projection.") # Debug

    # Calculate scale factor based on max projected distance (rho) from midpoints and optodes
    max_rho_midpoint = max((math.sqrt(p[0]**2 + p[1]**2) for p in projected_midpoints_2d_all_digi.values()), default=0)
    max_rho_optode = max((math.sqrt(p[0]**2 + p[1]**2) for p in projected_optodes_2d_digi.values()), default=0)
    # Consider fiducials too? Maybe not, they might be far off.
    max_rho = max(max_rho_midpoint, max_rho_optode)
    scale_factor = max_rho * 1.15 if max_rho > 1e-9 else 1.0 # Add margin, prevent division by zero

    # Normalize 2D coordinates for JS (scale to approx -1 to +1 range)
    norm_optodes_xy = {n: {'x': p[0] / scale_factor, 'y': p[1] / scale_factor}
                        for n, p in projected_optodes_2d_digi.items()}
    norm_fiducials_xy = {n: {'x': p[0] / scale_factor, 'y': p[1] / scale_factor}
                         for n, p in projected_fiducials_2d_digi.items()}

    # Specific handling for nose ('Nz')
    nose_norm = {'x': 0, 'y': 1} # Default upwards nose if Nz not found/projected
    if 'Nz' in norm_fiducials_xy:
        nose_norm = norm_fiducials_xy['Nz']
    # print(f"    Projection scale factor: {scale_factor:.3f}") # Debug


    # --- Final JS data structure ---
    js_data = {
        "trialNum": trial_num,
        "chromophore": chromophore,
        "channelType": "standard", # Or potentially add 'shortsep' later if needed
        "channels": {
            "names": common_channel_names,
            "positions3d": final_midpoints_3d_digi_ordered, # Pass ordered DIGITIZER midpoints (with NaNs if calc failed)
            "ba_regions": channel_ba_regions, # List of coarse region names per channel
        },
        "optodes": {
             # Separate sources/detectors based on names found in common channels
            "sources": {n: norm_optodes_xy.get(n, {'x':0,'y':0}) for n in active_sources_trial if n in norm_optodes_xy},
            "detectors": {n: norm_optodes_xy.get(n, {'x':0,'y':0}) for n in active_detectors_trial if n in norm_optodes_xy},
             # Pairs use names from selected_channel_info (which is ordered like common_channel_names)
            "pairs": [[ch['s'], ch['d']] for ch in selected_channel_info]
        },
        "head": { "nose": nose_norm }, # Use normalized nose position
        "geometry": { # Based on DIGITIZER coordinates
            "headCenter3d": head_center_3d_digi.tolist(),
            "headRadius3d": head_radius_3d_digi,
            "projectionScaleFactor": scale_factor
        },
        "subjects": all_subjects_data_processed, # Contains 'times' and 'data' lists per subject
        "subjectList": final_subject_list, # List of subject IDs successfully processed
        "fiducials": norm_fiducials_xy, # Normalized projected fiducial positions
        "channelSubgroups": subgroup_indices_map, # {name: [indices...]}
        "defaultSubgroup": default_subgroup_name
    }

    print(f"  Prepared data structure for JS. Processing time: {time.time() - start_time_trial:.2f} seconds.")
    # print(f"  DEBUG: Sample ba_regions: {channel_ba_regions[:min(15, len(channel_ba_regions))]}") # Print first 15 regions
    # unique_regions = set(r for r in channel_ba_regions if r) # Exclude None/empty if any
    # print(f"  DEBUG: Unique regions found in output: {unique_regions}")
    return js_data