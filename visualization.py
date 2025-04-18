# visualization.py
import os
import json
import sys
# No longer need shutil

def write_html_file(js_data, template_path="template.html", output_dir="visualizations"):
    """
    Writes the processed data into a single, portable HTML file by embedding
    CSS and JavaScript content from external files into the template.

    Args:
        js_data (dict): The dictionary containing all data for the JS visualization.
        template_path (str): Path to the HTML template file.
        output_dir (str): Directory where the output HTML file will be saved.
    """
    if js_data is None:
        print("Info: No data provided to write_html_file. Skipping HTML generation.")
        return

    trial_num = js_data.get('trialNum', 'UnknownTrial')
    chromophore = js_data.get('chromophore', 'UnknownChromo')
    template_dir = os.path.dirname(template_path)

    css_placeholder = "/* --- CSS_PLACEHOLDER --- */"
    js_placeholder = "// --- JAVASCRIPT_PLACEHOLDER ---"
    data_placeholder = "var EEG_DATA = {};" # The data injection placeholder

    css_content = ""
    js_content = ""

    try:
        # --- Read Template ---
        with open(template_path, "r", encoding='utf-8') as f:
            html_content = f.read()

        # --- Read CSS Dependency ---
        css_path = os.path.join(template_dir, 'style.css')
        try:
            with open(css_path, "r", encoding='utf-8') as f:
                css_content = f.read()
        except FileNotFoundError:
            print(f"Warning: CSS file '{css_path}' not found. Styling will be missing.", file=sys.stderr)
        except Exception as e:
            print(f"Warning: Failed to read CSS file '{css_path}': {e}", file=sys.stderr)

        # --- Read JS Dependency ---
        js_path = os.path.join(template_dir, 'visualization.js')
        try:
            with open(js_path, "r", encoding='utf-8') as f:
                js_content = f.read()
        except FileNotFoundError:
            print(f"Warning: JavaScript file '{js_path}' not found. Visualization logic will be missing.", file=sys.stderr)
        except Exception as e:
            print(f"Warning: Failed to read JavaScript file '{js_path}': {e}", file=sys.stderr)

    except FileNotFoundError:
        print(f"Error: HTML template file not found at '{template_path}'", file=sys.stderr)
        return
    except Exception as e:
        print(f"Error reading template or dependency files: {e}", file=sys.stderr)
        return

    try:
        # --- Create Output Directory ---
        os.makedirs(output_dir, exist_ok=True)
        output_filename = os.path.join(output_dir, f"fnirs_viewer_Trial{trial_num:02d}_{chromophore}.html")

        # --- Prepare Data ---
        json_string = json.dumps(js_data, separators=(',', ':')).replace("</", "<\\/") # Compact JSON

        # --- Inject CSS, JS, and Data into Template ---
        # Replace CSS placeholder
        if css_placeholder in html_content:
            html_content = html_content.replace(css_placeholder, css_content)
        else:
            print(f"Warning: CSS placeholder '{css_placeholder}' not found in template. CSS not embedded.", file=sys.stderr)

        # Replace JS placeholder
        if js_placeholder in html_content:
             # Important: Inject JS *after* the data placeholder is dealt with,
             # or ensure JS placeholder isn't accidentally modified by data injection.
             # The current order should be fine if placeholders are distinct.
            html_content = html_content.replace(js_placeholder, js_content)
        else:
            print(f"Warning: JavaScript placeholder '{js_placeholder}' not found in template. JS not embedded.", file=sys.stderr)

        # Replace Data placeholder
        if data_placeholder in html_content:
            html_content = html_content.replace(data_placeholder, f"var EEG_DATA = {json_string};")
        else:
             print(f"Error: Data placeholder '{data_placeholder}' not found in template '{template_path}'. Cannot inject data.", file=sys.stderr)
             # Decide if we should still write the file without data
             return # Probably best not to write if data injection fails

        # --- Write Final HTML ---
        with open(output_filename, "w", encoding='utf-8') as f:
            f.write(html_content)
        print(f"  Single-file visualization saved to: {output_filename}")

    except TypeError as e:
         print(f"Error: Could not serialize data to JSON: {e}. Check data types in js_data.", file=sys.stderr)
    except OSError as e:
         print(f"Error creating output directory or writing final HTML file '{output_filename}': {e}", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred during HTML file writing: {e}", file=sys.stderr)