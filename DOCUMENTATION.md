
---

## DOCUMENTATION.md (Revised & Enhanced)

```markdown
# Documentation: Interactive fNIRS Topography Pipeline

This document provides detailed technical information for the fNIRS visualization pipeline.

## Table of Contents

1.  [Code Structure](#1-code-structure)
2.  [Configuration Details](#2-configuration-details)
3.  [Input File Formats](#3-input-file-formats)
    *   [3.1 Signal Data (`.txt`)](#31-signal-data-txt)
    *   [3.2 Digitizer Data (`digpts.txt`)](#32-digitizer-data-digptstxt)
    *   [3.3 HTA Data (`.xlsx` - Optional)](#33-hta-data-xlsx---optional)
    *   [3.4 Anatomical Atlas Files](#34-anatomical-atlas-files)
4.  [Workflow & Technical Details](#4-workflow--technical-details)
    *   [4.1 Resource Loading](#41-resource-loading)
    *   [4.2 Registration (Digitizer to MNI)](#42-registration-digitizer-to-mni)
    *   [4.3 Anatomical Mapping](#43-anatomical-mapping)
    *   [4.4 Data Processing (HTA/Downsampling)](#44-data-processing-htadownsampling)
    *   [4.5 Visualization Techniques](#45-visualization-techniques)
        *   [4.5.1 3D Interpolation (IDW)](#451-3d-interpolation-idw)
        *   [4.5.2 2D Projection (Azimuthal Equidistant)](#452-2d-projection-azimuthal-equidistant)
    *   [4.6 Visualization Data Preparation (`js_data`)](#46-visualization-data-preparation-js_data)
5.  [Customization](#5-customization)
    *   [5.1 Changing Coarse Regions](#51-changing-coarse-regions)
    *   [5.2 Adjusting HTA Logic](#52-adjusting-hta-logic)
6.  [Limitations and Assumptions](#6-limitations-and-assumptions)
7.  [Troubleshooting](#7-troubleshooting)

---

## 1. Code Structure

The pipeline is organized into several Python modules:

*   **`main_script.py`**: Entry point. Contains configuration settings and orchestrates the overall workflow.
*   **`utils.py`**: General utility functions (`parse_digpts_file`, `project_point_azimuthal_equidistant`, `parse_nv_file`).
*   **`anatomical.py`**: Handles atlas/label loading, region mapping, coordinate transforms (MNI<->Voxel), KDTree building, and region lookup (`get_regions_for_mni_coords`). Contains `BA_TO_COARSE_REGION` definition.
*   **`registration.py`**: Calculates/applies the Digitizer-to-MNI transformation (`calculate_registration`, `apply_transform`) and extracts fiducials (`extract_fiducials`).
*   **`data_processing.py`**: Core processing (`prepare_trial_visualization_data`) including subject handling, HTA/downsampling, geometry, calling anatomical mapping, and preparing the `js_data` structure for visualization.
*   **`visualization.py`**: Writes the `js_data` dictionary into the HTML template (`write_html_file`).
*   **`get_fiducials_MNI.py`**: External dependency (uses MNE) to retrieve standard MNI coordinates for LPA, Nasion (Nz), RPA, and a placeholder Iz.
*   **`template.html`**: The HTML/CSS/JavaScript (D3.js v7) file for the interactive visualization.

---

## 2. Configuration Details

All user configuration is done within the top section of `main_script.py`.

*   **`FOLDER_NAME`**: Identifier for the dataset/experiment.
*   **`BASE_DATA_DIR`**: Root directory for preprocessed data folders.
*   **`INPUT_DATA_DIR`**: Path to trial-wise `.txt` files.
*   **`DIGPTS_FILEPATH`**: Full path to the `digpts.txt` file.
*   **`OUTPUT_VIS_DIR_BASE`**: Base directory name for output HTML files (gets `/HTA` suffix if HTA enabled).
*   **`TEMPLATE_HTML_PATH`**: Path to `template.html`.
*   **`CHROMOPHORES_TO_PROCESS`**: List of chromophore strings (e.g., `["HbO", "HbR"]`) matching channel names.
*   **`DO_HTA_AVERAGING`**: `True` to enable HTA, `False` for downsampling.
*   **`HTA_FILE_PATH`**: Path to HTA Excel file (if `DO_HTA_AVERAGING` is `True`).
*   **`HTA_SHEET_NAME`**: Sheet name in HTA Excel file.
*   **`HTA_TRIALNUM_FILTER`**: Integer trial number to filter HTA processing for, or `None` to attempt HTA for all trials.
*   **`HTACOL_RENAME_FN`**: Lambda function mapping subject ID (from filename) to HTA Excel column header.
*   **`ALL_CHANNEL_SUBGROUPS`**: Dict defining subgroups: `{ 'GroupName': [[S_num, D_num], ...], ... }`.
*   **`DEFAULT_SUBGROUP`**: Name of the subgroup selected by default (e.g., `"All Channels"`, `'PFC'`).
*   **`SRATE`**: Original sampling rate (Hz) of input `.txt` data.
*   **`DOWNSAMPLE_RATE`**: Integer downsampling factor (if HTA is off/fails).
*   **`SS_DETECTOR_NUMBERS`**: `set` of detector numbers to exclude (e.g., `set(range(16, 24))`).
*   **`REGISTRATION_FIDUCIALS`**: List of names for registration (e.g., `['LPA', 'Nz', 'RPA']`). Must exist in `digpts.txt` and be retrievable via `get_fiducials_MNI.py`.

---

## 3. Input File Formats

Strict adherence to these formats is required.

### 3.1 Signal Data (`.txt`)

*   **Location:** Directory specified by `INPUT_DATA_DIR`.
*   **Naming:** `SubjectID_TrialNum.txt` (e.g., `LC01_D01_04.txt`, `Subject02_1.txt`). `SubjectID` can contain letters, numbers, underscores. `TrialNum` is numeric.
*   **Format:** Comma-Separated Values (CSV).
*   **Header:** The first row **must** be a header containing channel names.
*   **Channel Names:** Must follow the pattern `S<SourceNum>_D<DetectorNum>_<Chromophore>` (e.g., `S1_D5_HbO`, `S8_D12_HbR`).
    *   `<SourceNum>` and `<DetectorNum>` correspond to integer numbers identifying sources/detectors. The names `S<SourceNum>` and `D<DetectorNum>` must exist in `digpts.txt`.
    *   `<Chromophore>` must match one of the strings in the `CHROMOPHORES_TO_PROCESS` list in `main_script.py`.
*   **Content:** Subsequent rows contain numerical time-series data.

### 3.2 Digitizer Data (`digpts.txt`)

*   **Location:** Path specified by `DIGPTS_FILEPATH`.
*   **Format:** Plain text. Each line defines a 3D point.
*   **Line Format:** `Name: X Y Z` or `Name X Y Z` (whitespace separated).
*   **`Name`:** Identifier (e.g., `S1`, `D5`, `Nz`, `LPA`, `RPA`, `Iz`). Must match source/detector names derived from `.txt` headers (e.g., `S1`, `D5`) and fiducial names used for registration (`Nz`, `LPA`, `RPA`). Case-insensitivity is generally handled for fiducials.
*   **`X`, `Y`, `Z`:** Numerical 3D coordinates in a consistent Cartesian system.
*   **Comments:** Lines starting with `#` are ignored. Empty lines are ignored.
*   **Example:**
    ```
    # Fiducials
    Nz: 0.05 92.11 -6.34
    LPA: -85.0 0.0 0.0
    RPA: 85.0 0.0 0.0
    Iz: 0.0 -95.0 -10.0
    # Sources
    S1: 85.37 -15.01 41.53
    S2: 70.50 10.20 55.80
    # Detectors
    D1: 75.10 -5.00 48.30
    D5: 70.12 -25.88 60.90
    ```

### 3.3 HTA Data (`.xlsx` - Optional)

*   **Required only if:** `DO_HTA_AVERAGING` is `True`.
*   **Location:** Path specified by `HTA_FILE_PATH`.
*   **Format:** Excel (`.xlsx`). Sheet specified by `HTA_SHEET_NAME`.
*   **Columns:** Headers must match subject identifiers *after* mapping by `HTACOL_RENAME_FN` (e.g., if filename is `Sub_01_1.txt` and `HTACOL_RENAME_FN = lambda x: x.replace("_", "")`, look for column `Sub01`).
*   **Content:** Rows contain event timestamps (seconds) relative to the recording start for that subject. Averaging occurs *between* consecutive timestamps. `NaN` values are ignored. At least two timestamps needed per subject for averaging.

### 3.4 Anatomical Atlas Files

*   **Location:** Defaults to `Brodmann_MM/` subdirectory relative to `anatomical.py`, or configure paths directly in `anatomical.py`.
*   **NIfTI File (`.nii.gz`):** Standard NIfTI format where voxel values are integers representing anatomical labels (e.g., Brodmann areas). Must be in MNI space or a space compatible with the MNI fiducials used. Example: `Brodmann_Mai_Matajnik.nii.gz`.
*   **Label File (`.txt`):** Maps integer values from the NIfTI file to string labels.
    *   Format: `Value Label` (whitespace separated).
    *   Example:
        ```
        # AFNI Brodmann Labels
        17 BA17(L)
        18 BA18(L)
        ...
        1017 BA17(R)
        ```
    *   The script uses regex `BA(\d+)` on the `Label` part to extract the Brodmann number for mapping to coarse regions defined in `anatomical.py`.

---

## 4. Workflow & Technical Details

### 4.1 Resource Loading
*   Loads configuration from `main_script.py`.
*   Loads digitizer points (`utils.parse_digpts_file`).
*   Loads MNI fiducials (`get_fiducials_MNI.py`).
*   Loads NIfTI atlas data/affine (`anatomical.load_atlas`).
*   Loads NIfTI label file (`anatomical.parse_label_file`).
*   Builds NIfTI value -> Coarse Region map (`anatomical.create_value_to_region_map`).
*   Builds KDTree from atlas MNI coordinates (`anatomical.build_atlas_kdtree`).
*   Loads HTA data if enabled.
*   Finds and groups trial `.txt` files.

### 4.2 Registration (Digitizer to MNI)
*   **Goal:** Find the transformation matrix (`T_digi_to_mni`) that aligns the 3D coordinates from the digitizer (`digpts.txt`) to the standard MNI coordinate system. This is crucial for anatomical localization.
*   **Method:**
    1.  Select corresponding fiducial points present in both Digitizer space (from `digpts.txt`) and MNI space (from `get_fiducials_MNI.py`). Typically uses `LPA`, `Nz`, `RPA` as defined in `REGISTRATION_FIDUCIALS`.
    2.  Calculate the optimal rigid transformation (rotation + translation) using the Kabsch algorithm (implemented via SVD in `scipy.spatial.transform.Rotation.align_vectors`). This minimizes the Root Mean Square Distance (RMSD) between the corresponding fiducial sets after alignment.
    3.  The result is a 4x4 transformation matrix (`T_digi_to_mni`) stored in `main_script.py`.
*   **Function:** `registration.calculate_registration`.

### 4.3 Anatomical Mapping
*   **Goal:** Determine the likely coarse brain region underlying each fNIRS channel.
*   **Method:**
    1.  Calculate the midpoint between the source and detector for each channel in **Digitizer Space**.
    2.  Transform these midpoints into **MNI Space** using the `T_digi_to_mni` matrix (`registration.apply_transform`).
    3.  For each channel's MNI midpoint, query the pre-built KDTree (`anatomical.atlas_kdtree`) to find the **nearest voxel that belongs to the brain** (i.e., has a non-zero label in the atlas NIfTI file). This avoids mapping to empty space outside the brain if the channel midpoint falls slightly outside the atlas volume.
    4.  Retrieve the integer label value from the atlas NIfTI data at the coordinates of that nearest brain voxel.
    5.  Look up this integer value in the `value_to_coarse_region_map` (created earlier from the label file and `BA_TO_COARSE_REGION` dictionary) to get the final coarse region name (e.g., 'DLPFC', 'Motor', 'Temporal', 'Unknown').
*   **Functions:** `registration.apply_transform`, `anatomical.get_regions_for_mni_coords`.
*   **Accuracy:** Depends on digitizer accuracy, fiducial placement consistency, registration quality (low RMSD), and the chosen atlas's accuracy and resolution.

### 4.4 Data Processing (HTA/Downsampling)
*   **Goal:** Reduce data size for visualization and/or perform task-related averaging.
*   **HTA:** If enabled and data is available, averages the signal (`np.nanmean`) in the interval *between* consecutive event markers provided in the `.xlsx` file. Timestamps become the midpoint of these intervals.
*   **Downsampling:** If HTA is disabled or fails, takes every Nth sample (`N = DOWNSAMPLE_RATE`) and calculates corresponding timestamps based on `DSFREQ`.
*   **Output:** For each subject, produces lists of timestamps and corresponding channel data (`times`, `data` lists within the `subjects` dict in `js_data`).

### 4.5 Visualization Techniques

#### 4.5.1 3D Interpolation (IDW)
*   **Goal:** Estimate signal values across the scalp surface for the heatmap display.
*   **Method:** Inverse Distance Weighting (IDW) performed in **3D Digitizer Space**.
    *   A grid of points representing the visualization surface is generated (conceptually mapped from the 2D SVG view back to 3D Digitizer space using the inverse of the projection, assuming a spherical head model based on `headRadius3d_digi`).
    *   For each grid point, the signal value is estimated as a weighted average of the signal values from the *active* nearby fNIRS channels (midpoints) at the currently selected time point.
    *   The weight for each channel decreases with distance `d` to the grid point, typically as `1 / d^p`.
    *   The power `p` is controlled by the 'Smoothing' slider in the UI (`p = 2 / smoothingFactor`). Higher smoothing values result in smaller `p`, making the interpolation smoother (distant channels have more influence).
*   **Function:** `interpolateValue3D` in `template.html` (JavaScript). Uses `EEG_DATA.channels.positions3d` (Digitizer space).

#### 4.5.2 2D Projection (Azimuthal Equidistant)
*   **Goal:** Map the 3D layout of optodes and channels onto a 2D plane suitable for a top-down view in the HTML visualization.
*   **Method:** Azimuthal Equidistant Projection.
    *   Calculates a head center based on active optodes/fiducials in **3D Digitizer Space**.
    *   For each 3D point (optode, fiducial) in Digitizer Space:
        *   Calculates the vector from the head center to the point.
        *   Determines the angle (`theta`) between this vector and the Z-axis (polar angle).
        *   Determines the angle (`phi`) of the vector's projection onto the XY plane (azimuthal angle).
        *   The 2D projected coordinates `(proj_x, proj_y)` are calculated as: `proj_x = theta * cos(phi)`, `proj_y = theta * sin(phi)`.
    *   This projection preserves the true distance of points from the center point (vertex) along radial lines. It's commonly used for polar views.
    *   The resulting 2D coordinates are then normalized using a `projectionScaleFactor` (based on the maximum projected `theta`) to fit within the SVG view (approximately -1 to +1 range).
*   **Function:** `utils.project_point_azimuthal_equidistant` (Python) and `project_point_azimuthal_equidistant_js` (JavaScript). Applied to **Digitizer** coordinates.

### 4.6 Visualization Data Preparation (`js_data`)
*   The `data_processing.prepare_trial_visualization_data` function compiles all necessary information into the `js_data` dictionary, which is then serialized to JSON. Key components for the JS:
    *   `channels.positions3d`: 3D **Digitizer** coordinates of channel midpoints (used for IDW interpolation and can be used for future 3D rendering).
    *   `channels.ba_regions`: List of anatomical region strings for coloring/labeling.
    *   `optodes`, `fiducials`: Contain normalized **2D projected** coordinates for drawing the layout in the current `template.html`.
    *   `geometry`: Contains **Digitizer** space parameters (`headCenter3d`, `headRadius3d`, `projectionScaleFactor`) used for projection and interpolation context.
    *   `subjects`: Contains the processed (HTA'd or downsampled) `times` and `data` lists.
*   `visualization.write_html_file` injects this JSON into `template.html`.

---

## 5. Customization

### 5.1 Changing Coarse Regions
*   Modify the `BA_TO_COARSE_REGION` dictionary within `anatomical.py`. Keys are BA numbers (int), values are desired region name strings.

### 5.2 Adjusting HTA Logic
*   Modify the `_perform_hta` helper function in `data_processing.py` for different averaging windows or baseline methods.

---

## 6. Limitations and Assumptions

*   **Montage Flexibility:** Visualization works for any montage if accurate 3D coordinates are in `digpts.txt`. Layout depends on these coordinates.
*   **Coordinate System Consistency:** Assumes a consistent Cartesian system in `digpts.txt`. The 2D projection works best if the Z-axis roughly aligns with the scalp vertex normal (for a top-down view). Non-standard orientations may produce unusual projections. Accuracy of anatomical mapping relies on `digpts.txt` coordinates being in a system that can be meaningfully registered to MNI.
*   **Head Model Approximation:** 3D interpolation uses channel midpoints in Digitizer space and assumes an approximate spherical surface for visualization interpolation, not a detailed individual anatomical model.
*   **Strict Naming Conventions:** Requires exact adherence to file (`SubjID_TrialNum.txt`) and channel (`S#_D#_Chromo`) naming conventions.
*   **HTA Data Format:** Excel structure, sheet name, and column headers (post-`HTACOL_RENAME_FN`) must match expectations. Averaging occurs between consecutive timestamps.
*   **Short-Separation Channels:** Detectors listed in `SS_DETECTOR_NUMBERS` are excluded by default. Modify `data_processing.py` if your SS naming differs or if you want to include them.
*   **Data Preprocessing:** Assumes input `.txt` data is appropriately preprocessed (filtering, motion correction, concentration conversion). This script does not perform these steps.
*   **Performance:** Very large datasets (long trials with downsampling, many channels/subjects) might affect browser performance. HTA significantly reduces data size.
*   **MNI Iz Coordinate:** The default Iz coordinate from `get_fiducials_MNI.py` is a placeholder and **must be verified** if used for critical alignment or analysis. Registration primarily uses LPA, Nz, RPA.

---

## 7. Troubleshooting

*   **`FileNotFoundError`**: Check paths in `main_script.py` (input/output dirs, digpts, HTA, template) and atlas file location.
*   **Missing Fiducials**: Ensure `REGISTRATION_FIDUCIALS` names exist in `digpts.txt` (case-insensitive check) and are provided by `get_fiducials_MNI.py`.
*   **No Files Found/Matched**: Check `INPUT_DATA_DIR` and ensure `.txt` filenames match the regex pattern in `main_script.py`.
*   **HTA Column Not Found**: Verify `HTACOL_RENAME_FN` maps subject IDs correctly to Excel column headers. Check for typos.
*   **Anatomical Regions All 'Unknown'**: Verify atlas/label loading, check registration RMSD (should be low, e.g., < 5-10mm), ensure `BA_TO_COARSE_REGION` includes relevant BAs. Check console output.
*   **Dependency Errors**: Ensure `numpy, pandas, scipy, nibabel, mne, openpyxl` are installed.
*   **Memory Errors**: Large datasets may require more RAM. Consider using HTA or processing fewer subjects/trials at once if memory is limited.
*   **Projection Looks Odd**: Verify the coordinate system in `digpts.txt`. Ensure Z generally points 'up' relative to the head for a standard top-down view. Check head center calculation.