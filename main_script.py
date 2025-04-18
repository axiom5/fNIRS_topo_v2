# main_script.py
import os
import sys
import glob
import re
import time
import pandas as pd
import numpy as np
from collections import defaultdict
import warnings # To suppress warnings if needed, e.g., during HTA nanmean

# --- Suppress specific warnings ---
# warnings.filterwarnings("ignore", category=FutureWarning, module='mne') # If using MNE directly here
warnings.filterwarnings("ignore", category=RuntimeWarning, message='Mean of empty slice') # For HTA nanmean

# --- Import from our refactored modules ---
import utils
import anatomical
import registration
import data_processing
import visualization

# --- Import from external dependency ---
try:
    # Assuming get_fiducials_MNI.py is in the same directory or Python path
    from get_fiducials_MNI import get_fsaverage_mni_fiducials_coreg
    print("Imported fiducial function from get_fiducials_MNI.py")
except ImportError as e:
    print(f"Error: Could not import from get_fiducials_MNI.py: {e}", file=sys.stderr)
    print("Ensure get_fiducials_MNI.py is in the Python path or the same directory.", file=sys.stderr)
    sys.exit(1)
except Exception as e_imp: # Catch other potential import errors
    print(f"An unexpected error occurred during import from get_fiducials_MNI.py: {e_imp}", file=sys.stderr)
    sys.exit(1)


# --- Main Execution Logic ---
if __name__ == "__main__":
    start_main_time = time.time()

    # --- Configuration ---
    # Data and File Paths
    # Use absolute paths or ensure relative paths are correct based on execution location
    INPUT_DATA_DIR = "Trialwise_txts"
    DIGPTS_FILENAME = "digpts_PC_est.txt" # Or 'digpts.txt'
    DIGPTS_FILEPATH = os.path.join(os.path.dirname(__file__), DIGPTS_FILENAME) # Assumes digpts is with scripts, adjust if elsewhere
    OUTPUT_VIS_DIR_BASE = "visualizations" # Base output dir for HTML files
    TEMPLATE_HTML_PATH = os.path.join(os.path.dirname(__file__), "template.html") # Assumes template is with scripts

    # Anatomical Files (Paths are now relative to anatomical.py location by default, override if needed)
    # anatomical.ATLAS_FOLDER = "path/to/Brodmann_MM" # Override if not in same dir as anatomical.py
    # anatomical.ATLAS_FILEPATH = os.path.join(anatomical.ATLAS_FOLDER, anatomical.ATLAS_FILENAME)
    # anatomical.LABEL_FILEPATH = os.path.join(anatomical.ATLAS_FOLDER, anatomical.LABEL_FILENAME)

    # Processing Parameters
    CHROMOPHORES_TO_PROCESS = ["HbO", "HbR"]
    DO_HTA_AVERAGING = False # True # Set to False to use downsampling instead
    HTA_FILE_PATH = r"C:\Users\Aseem\Documents\Colab Notebooks\testspace_files\HTA extraction\ETI Learning Curve\HTA_forattn_ETI_LC_trial4.xlsx" # Full path needed if not relative
    HTA_SHEET_NAME = 'Sheet1'
    HTA_TRIALNUM_FILTER = 4 # Process only this trial if HTA is enabled (set to None to process all HTA trials)
    # Function to map subject ID (from filename) to column name in HTA excel sheet
    HTACOL_RENAME_FN = lambda subj_id: subj_id.replace("_", "") # Example: maps "Subject_01" to "Subject01"

    # Channel Subgroup Definitions (Example)
    # Keys are subgroup names, values are lists of [Source_Num, Detector_Num] pairs
    ALL_CHANNEL_SUBGROUPS = {}
    # Example 1: PFC
    pfc_lpfc = [[1,1],[1,15],[2,3],[3,1],[3,3],[3,5],[4,1],[4,3],[4,15],[5,3]]
    pfc_rpfc = [[9,8],[9,15],[10,8],[10,9],[10,13],[11,9],[12,8],[12,9],[12,15],[13,9]]
    ALL_CHANNEL_SUBGROUPS['PFC'] = pfc_lpfc + pfc_rpfc
    # Example 2: Top 5 (ensure these match your actual layout)
    top5_lpfc = [[1,1],[1,15],[4,3],[5,3]] # Made up example
    top5_rpfc = [[9,8],] # Made up example
    ALL_CHANNEL_SUBGROUPS['Top5'] = top5_lpfc + top5_rpfc
    # Set a default subgroup to show initially in the visualization
    DEFAULT_SUBGROUP = "All Channels" # Or 'PFC', 'Top5' if defined above
    
    

    # Acquisition Parameters (adjust based on FOLDER_NAME or specific experiment)
    SRATE = 7.8125
    DOWNSAMPLE_RATE = 8
    SS_DETECTOR_NUMBERS = set(range(20, 29)) # Detectors considered short-separation


    DSFREQ = SRATE / DOWNSAMPLE_RATE # Calculated downsampled frequency

    # Registration Fiducials
    # Names of fiducials to use for the coregistration (must exist in digpts and MNI source)
    REGISTRATION_FIDUCIALS = ['LPA', 'Nz', 'RPA'] # Usually 3 points are sufficient and stable
    # --- End Configuration ---


    # --- Setup Output Directory ---
    output_visualization_dir = OUTPUT_VIS_DIR_BASE
    if DO_HTA_AVERAGING:
        output_visualization_dir = os.path.join(output_visualization_dir, "HTA")
        print(f"HTA Averaging is enabled. Output will be saved in: {output_visualization_dir}")
    else:
        print(f"HTA Averaging is disabled. Downsampling will be used. Output directory: {output_visualization_dir}")


    # --- Load HTA data (if enabled) ---
    hta_dataframe = None
    if DO_HTA_AVERAGING:
        try:
            hta_dataframe = pd.read_excel(HTA_FILE_PATH, sheet_name=HTA_SHEET_NAME)
            print(f"Successfully loaded HTA data from: {HTA_FILE_PATH} (Sheet: {HTA_SHEET_NAME})")
        except FileNotFoundError:
            print(f"Error: HTA file not found at {HTA_FILE_PATH}. HTA Averaging disabled.", file=sys.stderr)
            DO_HTA_AVERAGING = False # Disable if file not found
        except Exception as e:
            print(f"Error loading HTA file '{HTA_FILE_PATH}': {e}. HTA Averaging disabled.", file=sys.stderr)
            DO_HTA_AVERAGING = False # Disable on other errors


    # --- Load Digitizer Positions ---
    print("\n--- Loading Digitizer Positions ---")
    try:
        channel_positions_digi = utils.parse_digpts_file(DIGPTS_FILEPATH)
    except (FileNotFoundError, ValueError) as e:
         print(f"Fatal Error: Could not load or parse digitizer file '{DIGPTS_FILEPATH}'. {e}", file=sys.stderr)
         sys.exit(1)
    # print(f"Loaded {len(channel_positions_digi)} positions from {DIGPTS_FILEPATH}")


    # --- Load Anatomical Resources ---
    print("\n--- Loading Anatomical Resources ---")
    # 1. Load Atlas NIfTI (data, affine matrix, inverse affine)
    atlas_data, atlas_affine, inverse_atlas_affine = anatomical.load_atlas(anatomical.ATLAS_FILEPATH)
    if atlas_data is None: sys.exit(1) # Exit if atlas loading failed

    # 2. Parse Label File (maps NIfTI integer values to descriptive text labels)
    nifti_value_to_desc_map = anatomical.parse_label_file(anatomical.LABEL_FILEPATH)
    if nifti_value_to_desc_map is None: sys.exit(1) # Exit if label parsing failed

    # 3. Create Mapping from NIfTI value to Coarse Region Name (using BA_TO_COARSE_REGION)
    value_to_coarse_region_map = anatomical.create_value_to_region_map(
        nifti_value_to_desc_map, anatomical.BA_TO_COARSE_REGION
    )
    if not value_to_coarse_region_map:
        print("Warning: Could not create mapping from NIfTI values to coarse regions. Proceeding with 'Unknown'.")
        # Decide whether to proceed or exit. Proceeding is often acceptable.
        # sys.exit(1)

    # 4. Build KDTree for efficient nearest brain voxel lookup
    atlas_kdtree, kdtree_idx_to_voxel_map = anatomical.build_atlas_kdtree(atlas_data, atlas_affine)
    if atlas_kdtree is None: sys.exit(1) # Exit if KDTree building failed


    # --- Get MNI Fiducials & Calculate Registration ---
    print("\n--- Getting Standard MNI Fiducials & Calculating Registration ---")
    # 1. Get standard MNI fiducial coordinates using the external script/function
    mni_fiducials_all = get_fsaverage_mni_fiducials_coreg() # Returns dict {'LPA':..., 'Nz':..., 'RPA':..., 'Iz':...}
    if mni_fiducials_all is None:
        print("Fatal Error: Failed to retrieve MNI fiducials. Exiting.", file=sys.stderr)
        sys.exit(1)
    else:
        print("Successfully retrieved MNI Fiducials (including placeholder Iz):")
        for name, pos in mni_fiducials_all.items():
             print(f"  {name}: [{pos[0]:.2f}, {pos[1]:.2f}, {pos[2]:.2f}]")
        print("  Reminder: Verify the Inion (Iz) coordinate if used for anything critical.")

    # 2. Extract corresponding Digitizer fiducials using the specified names for registration
    print(f"\nExtracting digitizer points for registration using keys: {REGISTRATION_FIDUCIALS}")
    digi_fiducials_reg, missing_digi = registration.extract_fiducials(
        REGISTRATION_FIDUCIALS, channel_positions_digi, case_sensitive=False
    )
    if missing_digi:
        print(f"Fatal Error: Cannot proceed with registration due to missing digitizer fiducials: {missing_digi}", file=sys.stderr)
        sys.exit(1)

    # 3. Prepare points arrays for registration (ensure order matches)
    source_points = np.array([digi_fiducials_reg[key] for key in REGISTRATION_FIDUCIALS])
    target_points = np.array([mni_fiducials_all[key] for key in REGISTRATION_FIDUCIALS]) # Use MNI points matching the digi keys

    # 4. Calculate the transformation matrix (Digitizer -> MNI/Template)
    print("\nCalculating registration transformation...")
    T_digi_to_mni = registration.calculate_registration(source_points, target_points)
    if T_digi_to_mni is None:
        print("Fatal Error: Failed to calculate registration transformation. Exiting.", file=sys.stderr)
        sys.exit(1)
    else:
        print("Calculated Digi -> MNI Transformation Matrix (T):")
        np.set_printoptions(precision=4, suppress=True)
        print(T_digi_to_mni)
        np.set_printoptions() # Reset print options


    # --- Find and Group Data Files by Trial ---
    print("\n--- Finding Data Files ---")
    search_pattern = os.path.join(INPUT_DATA_DIR, "*.txt")
    all_files = glob.glob(search_pattern)
    if not all_files:
        print(f"Error: No '.txt' files found in '{INPUT_DATA_DIR}'. Check path and file extensions.", file=sys.stderr)
        sys.exit(1)

    grouped_files = defaultdict(dict) # {trial_num: {subject_id: filepath}}
    # Regex assumes SubjectID_TrialNum.txt format (adjust if different)
    # Example: "SubjectA_01.txt" -> subject_id='SubjectA', trial_num=1
    # Example: "P001_TaskX_Trial_3.txt" -> Needs adjusted regex
    filename_pattern = re.compile(r'(.+)_(\d+)\.txt') # Group 1: Subject ID, Group 2: Trial Num

    for f_path in all_files:
        filename = os.path.basename(f_path)
        match = filename_pattern.match(filename)
        if match:
            subject_id, trial_num_str = match.groups()
            try:
                trial_num = int(trial_num_str)
                grouped_files[trial_num][subject_id] = f_path
            except ValueError:
                print(f"Warning: Could not parse trial number from filename '{filename}'. Skipping.")
        else:
            print(f"Warning: Filename '{filename}' does not match expected pattern '{filename_pattern.pattern}'. Skipping.")

    if not grouped_files:
        print(f"Error: No files matching the pattern '{filename_pattern.pattern}' found and parsed in '{INPUT_DATA_DIR}'.", file=sys.stderr)
        sys.exit(1)

    print(f"Found data for {len(grouped_files)} trials: {sorted(grouped_files.keys())}")


    # --- Process Each Trial ---
    print("\n--- Processing Trials ---")
    processed_trial_count = 0
    for trial_num in sorted(grouped_files.keys()):

        # Skip trials if HTA is enabled and trial doesn't match the filter
        if DO_HTA_AVERAGING and HTA_TRIALNUM_FILTER is not None and trial_num != HTA_TRIALNUM_FILTER:
            # print(f"Skipping Trial {trial_num} (HTA processing filtered to Trial {HTA_TRIALNUM_FILTER})")
            continue

        print(f"\n===== Starting Trial {trial_num} =====")
        subject_files = grouped_files[trial_num]
        trial_data_loaded = {} # {subject_id: dataframe}

        # Load data for all subjects in this trial first
        print(f"  Loading data for {len(subject_files)} subjects...")
        load_errors = 0
        for subject_id, f_path in subject_files.items():
            try:
                # Consider dtype=np.float32 for memory efficiency if files are large
                df = pd.read_csv(f_path, sep=',', header=0, dtype=np.float32)
                # Basic validation: Check if dataframe is empty or has no columns
                if df.empty or len(df.columns) == 0:
                    print(f"  Warning: Loaded empty or headerless dataframe from {f_path}. Skipping subject {subject_id}.")
                    load_errors += 1
                    continue
                trial_data_loaded[subject_id] = df
                # print(f"    Loaded data for {subject_id}") # Verbose
            except FileNotFoundError:
                 print(f"  Error: File not found during data loading: {f_path}. Skipping subject {subject_id}.")
                 load_errors += 1
            except pd.errors.EmptyDataError:
                 print(f"  Error: No data or columns found in {f_path}. Skipping subject {subject_id}.")
                 load_errors += 1
            except Exception as e:
                print(f"  Error loading {f_path}: {e}. Skipping subject {subject_id}.")
                load_errors += 1
                continue

        if not trial_data_loaded:
            print(f"Warning: No data successfully loaded for any subject in Trial {trial_num}. Skipping trial.")
            continue
        if load_errors > 0:
             print(f"  Note: {load_errors} subjects skipped during loading for Trial {trial_num}.")


        # Generate visualization for each specified chromophore
        for chromo in CHROMOPHORES_TO_PROCESS:
            # Call the main processing function from data_processing module
            js_data_for_trial = data_processing.prepare_trial_visualization_data(
                trial_data=trial_data_loaded,
                trial_num=trial_num,
                ch_pos=channel_positions_digi, # Pass the loaded digitizer positions
                chromophore=chromo,
                srate=SRATE,
                downsample_rate=DOWNSAMPLE_RATE,
                do_hta_averaging=DO_HTA_AVERAGING,
                hta_df=hta_dataframe,          # Pass the loaded HTA dataframe (or None)
                hta_col_rename_func=HTACOL_RENAME_FN, # Pass the rename function
                dsfreq=DSFREQ,
                ss_detector_numbers=SS_DETECTOR_NUMBERS,
                channel_subgroups_dict=ALL_CHANNEL_SUBGROUPS,
                default_subgroup_name=DEFAULT_SUBGROUP,
                # --- Pass anatomical & registration resources ---
                T_digi_to_template=T_digi_to_mni, # Pass the calculated transformation
                value_to_region_map=value_to_coarse_region_map,
                atlas_data=atlas_data,
                atlas_brain_kdtree=atlas_kdtree,
                kdtree_index_to_voxel_index=kdtree_idx_to_voxel_map,
                # ignore_fiducials_for_center defaults to False in function def
            )

            # Write the output HTML file if data was generated
            if js_data_for_trial:
                 visualization.write_html_file(
                     js_data_for_trial,
                     TEMPLATE_HTML_PATH,
                     output_visualization_dir
                 )
                 processed_trial_count +=1 # Count successful processing per chromophore/trial

    print("\n===== Processing Complete =====")
    print(f"Total trials processed (incl. chromophores): {processed_trial_count}")
    print(f"Total execution time: {time.time() - start_main_time:.2f} seconds")