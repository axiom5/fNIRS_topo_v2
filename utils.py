# utils.py
import numpy as np
import re
import math
import os
import sys
import time
import pandas as pd # Keep pandas import if parse_nv_file is used later or for potential future utils

def parse_digpts_file(filename):
    """Parse digpts.txt file to extract channel positions"""
    print(f"Parsing digitizer points file: {filename}...")
    ch_pos = {}
    try:
        with open(filename, 'r') as f:
            for line in f:
                if not line.strip() or line.strip().startswith('#'): continue # Skip empty/comment
                # Adjusted regex to be more flexible with separators and handle potential labels
                match = re.match(r'(\w+):?\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', line)
                if match:
                    ch_name, x, y, z = match.groups()
                    # Skip Cz explicitly if needed based on previous comment
                    if ch_name.upper() == 'CZ':
                         print(f"  Skipping potentially erroneous entry: {line.strip()}")
                         continue
                    try:
                        ch_pos[ch_name] = np.array([float(x), float(y), float(z)])
                    except ValueError:
                         print(f"  Warning: Could not parse coordinates for {ch_name} in line: {line.strip()}")
                # else: print(f"  Warning: Could not parse line in digpts: {line.strip()}") # Optional debug for non-matching lines
    except FileNotFoundError:
        print(f"Error: Digitizer file '{filename}' not found.", file=sys.stderr)
        raise # Re-raise the exception to be handled by the caller
    except Exception as e:
        print(f"An unexpected error occurred parsing {filename}: {e}", file=sys.stderr)
        raise # Re-raise

    if not ch_pos:
        raise ValueError(f'No valid channel positions could be parsed from {filename}.')
    print(f"  Successfully parsed {len(ch_pos)} positions.")
    return ch_pos

def project_point_azimuthal_equidistant(point_3d, center_3d):
    """Projects a 3D point onto a 2D plane using Azimuthal Equidistant projection"""
    vec = point_3d - center_3d
    dist = np.linalg.norm(vec)
    if dist < 1e-9: return 0.0, 0.0
    # Clamp the argument for acos to avoid domain errors due to floating point inaccuracies
    acos_arg = np.clip(vec[2] / dist, -1.0, 1.0)
    theta = math.acos(acos_arg)
    phi = math.atan2(vec[1], vec[0]) # Use atan2 for correct quadrant
    rho = theta # In Azimuthal Equidistant, radial distance is the angle
    proj_x = rho * math.cos(phi)
    proj_y = rho * math.sin(phi)
    return proj_x, proj_y

def parse_nv_file(filename):
    """Parses the specific .nv file format for brain mesh (vertices/faces)."""
    # Minimal version, as it wasn't actively used in the core logic shown previously
    print(f"Parsing NV file: {filename}...")
    start_time = time.time()
    try:
        with open(filename, 'r') as f:
            f.readline() # Skip Comment
            num_vertices = int(f.readline().strip())
            vertices = np.array([list(map(float, f.readline().split())) for _ in range(num_vertices)], dtype=np.float32)
            num_faces = int(f.readline().strip())
            # Read faces carefully, handling potential 0 or 1 based indexing
            faces_raw = np.array([list(map(int, f.readline().split())) for _ in range(num_faces)], dtype=np.int32)

            min_idx = np.min(faces_raw)
            if min_idx == 1:
                faces = faces_raw - 1 # Convert 1-based to 0-based
            elif min_idx == 0:
                faces = faces_raw
            else:
                 raise ValueError(f"Unexpected minimum face index {min_idx} found in {filename}")

            if np.max(faces) >= num_vertices:
                raise ValueError(f"Face index ({np.max(faces)}) out of bounds for {num_vertices} vertices in {filename}")

            print(f"  Successfully parsed {vertices.shape[0]} vertices, {faces.shape[0]} faces in {time.time() - start_time:.2f}s.")
            return vertices, faces
    except FileNotFoundError:
        print(f"Error: NV file not found at {filename}", file=sys.stderr)
        return None, None
    except Exception as e:
        print(f"Error parsing NV file {filename}: {e}", file=sys.stderr)
        return None, None

# Potentially add other general utilities if needed (e.g., safe loading wrappers)