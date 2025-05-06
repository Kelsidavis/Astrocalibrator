import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import os
import shutil
import threading
import queue
import zipfile
from datetime import datetime
from astropy.io import fits

from gui import root, log_message, output_folder_var, progress_var, progress_label_var
from gui import session_title_var, master_dark_path, master_flat_path, master_bias_path
from gui import master_dark_enabled, master_flat_enabled, master_bias_enabled
from gui import ToolTip
from solving import plate_solve_and_update_header
from settings import load_settings, save_settings, remember_file, get_remembered_file

from gui import light_files, dark_files, flat_files, bias_files, dark_flat_files
from gui import light_btn, dark_btn, flat_btn, bias_btn, darkflat_btn
from gui import reset_btn
from gui import master_dark_btn, master_flat_btn, master_bias_btn
from gui import progress_bar
from gui import object_description_var, object_distance_var
from object_info import object_info

from astropy.coordinates import SkyCoord
import astropy.units as u
import tkinter.messagebox as mb
import glob
from calibration import run_parallel_calibration, load_fits_by_filter, normalize_filter_name


# Globals for buttons
select_output_btn = None
calibrate_btn = None
current_output_label = None
save_masters_var = None

solved_frames_wcs = set()
unsolved_frames = set()

# Prevent duplicate GUI element creation
main_widgets_initialized = False

# Constants for master calibration files
MASTER_FILES = [
    "master_dark.fits",
    "master_flat.fits",
    "master_bias.fits",
    "master_dark_flat.fits"
]

# Configure a custom Progressbar style to make it visible
style = ttk.Style()
style.theme_use('default')
style.configure(
    "TProgressbar",
    thickness=20,               # Make it thicker
    troughcolor="#333333",       # Dark background
    background="#4CAF50"         # Green moving bar
)

# main.py

def has_enough_space(folder, required_bytes=500_000_000):
    """Check if the specified folder has enough free disk space."""
    try:
        total, used, free = shutil.disk_usage(folder)
        return free > required_bytes
    except Exception as e:
        log_message(f"⚠️ Could not check disk space: {e}")
        return False

def initialize_main_widgets():
    global main_widgets_initialized
    if main_widgets_initialized:
        return
    main_widgets_initialized = True

    global select_output_btn, calibrate_btn, current_output_label

    # ✅ Create the main output folder frame
    output_folder_frame = tk.Frame(root)
    output_folder_frame.pack(pady=(5, 0))

    # ✅ Button to select output folder
    select_output_btn = tk.Button(
        output_folder_frame,
        text="Select Output Folder",
        font=("Arial", 10, "bold"),
        width=18,
        height=2,
        command=select_output_directory
    )
    ToolTip(select_output_btn, "Choose where calibrated and solved files will be saved.")
    select_output_btn.pack(padx=10, pady=5)

    # ✅ Sub-frame for output folder label
    current_output_frame = tk.Frame(output_folder_frame)
    current_output_frame.pack(fill='x', padx=10)

    current_output_label = tk.Label(
        current_output_frame,
        textvariable=output_folder_var,
        font=("Arial", 8),
        anchor='w',
        justify='left',
        wraplength=600,
        fg="gray"
    )
    current_output_label.pack(fill='x')

    # ✅ Control frame
    control_frame = tk.Frame(root)
    control_frame.pack(pady=10)

    # Save Masters
    save_masters_frame = tk.Frame(control_frame)
    save_masters_frame.pack(side='left', padx=(10, 50))

    global save_masters_var
    save_masters_var = tk.BooleanVar(value=False)
    save_masters_checkbox = tk.Checkbutton(save_masters_frame, text="Save Calibration Masters", variable=save_masters_var)
    ToolTip(save_masters_checkbox, "Save generated master dark, flat, and bias calibration frames for future use.")
    save_masters_checkbox.pack()

    # Buttons
    buttons_frame = tk.Frame(control_frame)
    buttons_frame.pack(side='left', padx=10)

    calibrate_btn = tk.Button(buttons_frame, text="Calibrate Files", font=("Arial", 12, "bold"), width=20, height=2)
    ToolTip(calibrate_btn, "Plate solve light frames and apply calibration using selected masters and settings.")
    calibrate_btn.pack(side='left', padx=10)

def global_cleanup(output_folder, light_file_paths=None):
    """
    Delete leftover .wcs, .ini, and .bak files from:
    - the output folder
    - the folders containing the original light files
    """
    patterns = ["*.wcs", "*.ini", "*.bak"]

    # Clean from output folder
    for pattern in patterns:
        for file_path in glob.glob(os.path.join(output_folder, pattern)):
            try:
                os.remove(file_path)
                print(f"🧹 Deleted from output: {file_path}")
            except Exception as e:
                print(f"⚠️ Failed to delete from output: {file_path} – {e}")

    # Clean from light file source directories
    if light_file_paths:
        cleaned_dirs = set()
        for file_path in light_file_paths:
            directory = os.path.dirname(file_path)
            if directory in cleaned_dirs:
                continue
            for pattern in patterns:
                for extra_file in glob.glob(os.path.join(directory, pattern)):
                    try:
                        os.remove(extra_file)
                        print(f"🧹 Deleted from source: {extra_file}")
                    except Exception as e:
                        print(f"⚠️ Failed to delete from source: {extra_file} – {e}")
            cleaned_dirs.add(directory)

def inject_wcs_from_sidecar(fits_path):
    """Inject WCS fields into a FITS file from its .wcs sidecar."""
    wcs_sidecar_path = os.path.splitext(fits_path)[0] + '.wcs'
    if not os.path.exists(wcs_sidecar_path):
        print(f"⚠️ No WCS sidecar found for {fits_path}")
        return False

    try:
        with fits.open(fits_path, mode='update') as hdul_cal:
            header = hdul_cal[0].header

            with fits.open(wcs_sidecar_path) as hdul_wcs:
                wcs_header = hdul_wcs[0].header

                wcs_keys = ['CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2',
                            'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
                            'CTYPE1', 'CTYPE2']

                for key in wcs_keys:
                    if key in wcs_header:
                        header[key] = wcs_header[key]

                header['WCSINJ'] = (True, 'WCS injected from sidecar')

            hdul_cal.flush()

        print(f"✅ WCS injected into {os.path.basename(fits_path)}")
        return True

    except Exception as e:
        print(f"💥 Failed to inject WCS into {fits_path}: {e}")
        return False

def generate_fallback_name(header=None):
    date_stamp = datetime.now().strftime("%Y-%m-%d")

    if header:
        try:
            ra_deg = header.get('CRVAL1')
            dec_deg = header.get('CRVAL2')
            if ra_deg is not None and dec_deg is not None:
                ra_h = int(ra_deg / 15)
                ra_m = int((ra_deg / 15 - ra_h) * 60)
                dec_sign = '+' if dec_deg >= 0 else '-'
                dec_d = int(abs(dec_deg))
                dec_m = int((abs(dec_deg) - dec_d) * 60)
                return f"RA{ra_h:02d}h{ra_m:02d}m_DEC{dec_sign}{dec_d:02d}d{dec_m:02d}m_{date_stamp}"
        except Exception as e:
            print(f"⚠️ Failed to read fallback RA/DEC: {e}")

    return f"UnknownObject_{date_stamp}"

def maintain_object_info():
    """Periodically restore object info during background operations."""
    if cached_object_description:
        object_description_var.set(cached_object_description)
    if cached_object_distance:
        object_distance_var.set(cached_object_distance)

    root.after(500, maintain_object_info)  # Keep refreshing every 500ms

# --- Global collected session names ---
session_names_collected = []

# --- Cache solved object info ---
cached_object_description = None
cached_object_distance = None

def select_output_directory():
    path = filedialog.askdirectory(title="Select Output Folder")
    if path:
        output_folder_var.set(path)
        log_message(f"📂 Output folder set to: {path}")

        # Only adjust styles, widths, and enable buttons
        select_output_btn.config(font=("Arial", 10), width=18, height=1)
        calibrate_btn.config(font=("Arial", 14, "bold"), width=25, height=3)

        light_btn.config(state='normal')
        dark_btn.config(state='normal')
        flat_btn.config(state='normal')
        darkflat_btn.config(state='normal')
        bias_btn.config(state='normal')

        if calibrate_btn:
            calibrate_btn.config(state='normal')
        else:
            log_message("⚠️ Calibrate button not yet available to re-enable.")

        reset_btn.config(state='normal')
        master_dark_btn.config(state='normal')
        master_flat_btn.config(state='normal')
        master_bias_btn.config(state='normal')

def find_nearest_known_object(fits_path, catalog):
    def safe_float(val):
        try:
            return float(str(val).strip().split()[0])
        except Exception:
            return None

    try:
        # Try reading WCS headers first
        with fits.open(fits_path) as hdul:
            hdr = hdul[0].header
            ra = safe_float(hdr.get('CRVAL1'))
            dec = safe_float(hdr.get('CRVAL2'))

        if ra is not None and dec is not None:
            if -90.0 <= dec <= 90.0:
                log_message(f"🧭 FITS WCS center: RA={ra:.4f}°, Dec={dec:.4f}°")
            else:
                log_message(f"⚠️ Invalid FITS WCS: Dec={dec:.2f}°. Trying sidecar WCS...")
                ra, dec = None, None  # Force fallback to sidecar


        # If WCS headers are missing, try reading the .wcs sidecar file
        wcs_file = os.path.splitext(fits_path)[0] + '.wcs'
        if (ra is None or dec is None) and os.path.exists(wcs_file):
            try:
                with fits.open(wcs_file) as hdul:
                    hdr = hdul[0].header
                    ra = safe_float(hdr.get('CRVAL1'))
                    dec = safe_float(hdr.get('CRVAL2'))
                    if -90.0 <= dec <= 90.0:
                        log_message(f"🧭 Sidecar WCS center: RA={ra:.4f}°, Dec={dec:.4f}°")
                    else:
                        log_message(f"⚠️ Skipping invalid sidecar WCS: Dec={dec:.2f}°")
                        return None
            except Exception as e:
                log_message(f"⚠️ Failed to read sidecar WCS: {e}")

        if ra is not None and dec is not None:
            if not (-90.0 <= dec <= 90.0):
                log_message(f"⚠️ Invalid declination value: {dec:.2f}°. Skipping object match.")
                return None
            log_message(f"🔬 Attempting SkyCoord with RA={ra}, Dec={dec}")
            print(f"DEBUG RAW: ra={ra} ({type(ra)}), dec={dec} ({type(dec)})")
            solved_coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')

            closest_object = None
            min_distance = float('inf')

            for name, (desc, dist, obj_ra, obj_dec) in catalog.items():
                obj_coord = SkyCoord(ra=obj_ra*u.deg, dec=obj_dec*u.deg, frame='icrs')
                sep = solved_coord.separation(obj_coord)
                if sep.degree < min_distance:
                    min_distance = sep.degree
                    closest_object = name

            if closest_object and min_distance < 8:
                return closest_object

    except Exception as e:
        log_message(f"⚠️ WCS matching failed: {e}")

    return None
   

def _calibration_worker():
    import time
    start_time = time.time()
    progress_label_var.set("Calibrating frames...")

    output_folder = output_folder_var.get()
    if not has_enough_space(output_folder, 500_000_000):  # Require 500MB minimum
        log_message("❌ Not enough disk space. Please free at least 500MB before running calibration.")
        progress_label_var.set("Insufficient disk space")
        progress_bar.stop()
        progress_bar.config(mode="determinate")
        progress_var.set(0)
        calibrate_btn.config(state='normal')
        return

    if cached_object_description:
        object_description_var.set(cached_object_description)
    if cached_object_distance:
        object_distance_var.set(cached_object_distance)

    progress_bar.config(mode="indeterminate")
    progress_bar.start(10)

    method = 'median'
    output_folder = output_folder_var.get()

    first_light_path = next(iter(light_files), None)
    if not first_light_path:
        log_message("❌ No light frames found.")
        if calibrate_btn:
            calibrate_btn.config(state='normal')
        return

    temp_folder = os.path.join(output_folder, "temp")
    solve_temp_folder = os.path.join(output_folder, "solve_temp")

    os.makedirs(temp_folder, exist_ok=True)

    light_by_filter = load_fits_by_filter(light_files)
    dark_by_filter = load_fits_by_filter(dark_files)
    flat_by_filter = load_fits_by_filter(flat_files)
    bias_by_filter = load_fits_by_filter(bias_files)

    # 🔧 Normalize filters before calibration to ensure master matching consistency
    light_by_filter = {normalize_filter_name(k): v for k, v in light_by_filter.items()}
    dark_by_filter = {normalize_filter_name(k): v for k, v in dark_by_filter.items()}
    flat_by_filter = {normalize_filter_name(k): v for k, v in flat_by_filter.items()}
    bias_by_filter = {normalize_filter_name(k): v for k, v in bias_by_filter.items()}


    master_dark_paths, master_flat_paths, master_bias_paths = run_parallel_calibration(
        light_by_filter=light_by_filter,
        dark_by_filter=dark_by_filter,
        flat_by_filter=flat_by_filter,
        bias_by_filter=bias_by_filter,
        output_folder=output_folder,
        session_title=session_title_var.get(),
        log_callback=log_message,
        save_masters=save_masters_var.get(),
        dark_flat_file_list=dark_flat_files  # ✅ Pass raw user input directly
    )

    elapsed = time.time() - start_time
    log_message(f"✅ Calibration complete in {elapsed:.2f} seconds.")

    # 🧹 Cleanup
    try:
        if os.path.exists(temp_folder):
            shutil.rmtree(temp_folder)
            log_message("🧹 Temp folder cleaned up.")
        if os.path.exists(solve_temp_folder):
            shutil.rmtree(solve_temp_folder)
            log_message("🧹 Solve temp folder cleaned up.")
    except Exception as e:
        log_message(f"⚠️ Cleanup warning: {e}")

    try:
        session_name_cleaned = session_title_var.get().replace(' ', '_').replace(':', '').replace('/', '_')
        imaging_date = datetime.now().strftime("%Y-%m-%d")

        if first_light_path:
            with fits.open(first_light_path) as hdul:
                header = hdul[0].header
                date_obs = header.get('DATE-OBS') or header.get('DATE') or header.get('DATEOBS')
                if date_obs:
                    imaging_date = date_obs.split('T')[0]

        calibrated_folder = os.path.join(output_folder, "calibrated")
        zip_filename = os.path.join(output_folder, f"{session_name_cleaned}__{imaging_date}_calibrated.zip")

        # 🎬 Begin progress animation
        fits_files = [
            os.path.join(root_dir, file)
            for root_dir, _, files in os.walk(calibrated_folder)
            for file in files if file.endswith('.fits')
        ]
        total_files = len(fits_files)

        progress_label_var.set("📦 Zipping calibrated frames... 0 of {}".format(total_files))
        progress_bar.config(mode='indeterminate')
        progress_bar.start()
        root.update_idletasks()

        with zipfile.ZipFile(zip_filename, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for i, file_path in enumerate(fits_files, 1):
                arcname = os.path.relpath(file_path, calibrated_folder)
                zipf.write(file_path, arcname)
                progress_label_var.set(f"📦 Zipping calibrated frames... {i} of {total_files}")
                root.update_idletasks()

        # ✅ Done
        progress_bar.stop()
        progress_bar.config(mode='determinate')
        progress_var.set(100)
        progress_label_var.set("✅ ZIP archive created.")
        log_message(f"📦 Created ZIP archive: {zip_filename}")

        # 🧹 Remove calibrated folder after zipping
        if os.path.exists(calibrated_folder):
            shutil.rmtree(calibrated_folder)
            log_message("🧹 Deleted calibrated folder after creating ZIP.")

    except Exception as e:
        log_message(f"⚠️ Failed to create calibrated ZIP or delete calibrated folder: {e}")

    # 🧹 Remove leftover master calibration files if not saving masters
    try:
        if not save_masters_var.get():
            for master_name in MASTER_FILES:
                master_path = os.path.join(output_folder, master_name)
                if os.path.exists(master_path):
                    os.remove(master_path)
                    log_message(f"🧹 Deleted leftover {master_name}.")
    except Exception as e:
        log_message(f"⚠️ Failed to delete leftover master files: {e}")

    # 🧹 Delete leftover sidecar files
    try:
        for pattern in ["*.wcs", "*.ini", "*.bak"]:
            for leftover in glob.glob(os.path.join(output_folder, pattern)):
                os.remove(leftover)
                print(f"🧹 Deleted leftover file: {leftover}")
        log_message("🧹 Deleted stray .wcs, .ini, and .bak files.")
    except Exception as e:
        log_message(f"⚠️ Failed during post-calibration cleanup: {e}")

    # 🔧 Final cleanup of leftover files from output and light source folders
    global_cleanup(output_folder, light_files)

    # ✅ Restore GUI
    progress_bar.stop()
    progress_bar.config(mode="determinate")
    progress_var.set(100)
    calibrate_btn.config(state='normal')
    progress_label_var.set("Idle")
    fade_out_progress_label()

    if cached_object_description:
        object_description_var.set(cached_object_description)
    if cached_object_distance:
        object_distance_var.set(cached_object_distance)

def run_plate_solving():
    global calibrate_btn
    result_queue = queue.Queue()
    solve_temp_folder = os.path.join(output_folder_var.get(), "solve_temp")
    os.makedirs(solve_temp_folder, exist_ok=True)

    # 🛡️ Safe check
    if calibrate_btn:
        calibrate_btn.config(state='disabled')
    else:
        log_message("⚠️ Calibrate button not initialized yet. Skipping disabling.")

    progress_label_var.set("Plate solving images...")
    log_message("📅 Starting plate solving in background...")
    progress_bar.config(mode="indeterminate")
    progress_bar.start(10)
    root.update_idletasks()

    print(f"💬 light_files from GUI: {light_files}")
    light_files_to_solve = [f for f in light_files if os.path.exists(f)]
    solver_failed = False

    def solve_worker(path):
        nonlocal solver_failed
        global cached_object_description, cached_object_distance
        try:
            log_message(f"🧪 Solving: {path}")
            try:
                print(f"👣 Entering solve_worker() for path: {path}")
                print(f"🛤️ Checking if file exists: {os.path.exists(path)}")
                # 🔧 PATCHED: handle tuple return
                result = plate_solve_and_update_header(path, log_message)
                if isinstance(result, tuple):
                    session_name_from_solver, _ = result
                else:
                    session_name_from_solver = result
                session_name_nearby = find_nearest_known_object(path, object_info)

                # Pick the best session name
                if session_name_nearby:
                    session_name = session_name_nearby
                    log_message(f"🎯 Overriding solver object '{session_name_from_solver}' with nearby known object '{session_name}'")
                elif session_name_from_solver:
                    session_name = session_name_from_solver
                    log_message(f"📋 Keeping solver object: '{session_name}'")
                else:
                    session_name = "UnknownObject"
                    log_message("⚠️ No valid object name found from solver or catalog.")

                # Normalize it properly
                if session_name:
                    session_name_upper = session_name.upper().strip()
                    if session_name_upper.startswith('M') and session_name_upper[1:].isdigit():
                        session_name = f"Messier {session_name_upper[1:]}"
                    elif session_name_upper.startswith('NGC') and session_name_upper[3:].strip().isdigit():
                        session_name = f"NGC {session_name_upper[3:].strip()}"
                    elif session_name_upper.startswith('IC') and session_name_upper[2:].strip().isdigit():
                        session_name = f"IC {session_name_upper[2:].strip()}"
                    else:
                        session_name = session_name_upper

                # ✅ Now finally set the title
                if session_name:
                    session_title_var.set(session_name)

                lookup_name = session_name.strip().upper().replace("NGC ", "NGC ").replace("IC ", "IC ").replace("MESSIER ", "Messier ")

                info = object_info.get(lookup_name)
                if info:
                    description, distance, _, _ = info
                else:
                    if session_name_nearby and session_name_nearby != session_name:
                        log_message(f"🧭 Falling back to nearest known object: {session_name_nearby}")
                        session_name = session_name_nearby
                        session_title_var.set(session_name)
                        lookup_name = session_name.strip().upper().replace("NGC ", "NGC ").replace("IC ", "IC ").replace("MESSIER ", "Messier ")
                        info = object_info.get(lookup_name)
                        if info:
                            description, distance, _, _ = info
                        else:
                            description = "No description available"
                            distance = "Unknown distance"
                    else:
                        description = "No description available"
                        distance = "Unknown distance"
                        log_message(f"⚠️ Object '{lookup_name}' not found in database.")


                object_description_var.set(description)
                object_distance_var.set(f"Distance: {distance}")
                cached_object_description = description
                cached_object_distance = f"Distance: {distance}"

                result_queue.put(session_name)
                session_names_collected.append(session_name)

                # 🛠 Debug FITS header right after solving
                try:
                    with fits.open(path) as hdul:
                        hdr = hdul[0].header
                        print("🔍 FITS Header keys:", list(hdr.keys()))
                        print("🔍 FITS OBJECT field:", hdr.get('OBJECT'))
                except Exception as e:
                    print(f"⚠️ Could not read FITS header: {e}")

            except FileNotFoundError as fnf_err:
                if not solver_failed:
                    mb.showinfo("Plate Solver Not Found", "The plate solver executable could not be found. Calibration will continue without solving.")
                    solver_failed = True
                log_message(f"❌ Solver not found: {fnf_err}")
                result_queue.put(None)
            except Exception as e:
                import traceback
                log_message(f"💥 Exception in solve_worker: {e}\n{traceback.format_exc()}")
        except Exception as outer_err:
            log_message(f"💥 Outer exception in solve_worker: {outer_err}")


    for path in light_files_to_solve:
        threading.Thread(target=solve_worker, args=(path,), daemon=True).start()

    def check_solving_results(result_queue):
        global cached_object_description, cached_object_distance
        try:
            while True:
                session_name = result_queue.get_nowait()
                if session_name:
                    # Normalize name before lookup
                    def fix_catalog_name(name):
                        name = name.strip().title()
                        name = name.replace("Ngc", "NGC").replace("Ic", "IC").replace("Messier", "Messier")
                        return name

                    lookup_name = fix_catalog_name(session_name)

                    info = object_info.get(lookup_name)
                    if info:
                        object_description_var.set(info[0])
                        object_distance_var.set(f"Distance: {info[1]}")
                        cached_object_description = info[0]
                        cached_object_distance = f"Distance: {info[1]}"
                    else:
                        object_description_var.set("No description available")
                        object_distance_var.set("Unknown distance")
                        cached_object_description = "No description available"
                        cached_object_distance = "Unknown distance"

                    log_message(f"📅 Solved frame for object: {session_name}")

        except queue.Empty:
            pass

        worker_threads_alive = any(
            thread.name.startswith("Thread-") for thread in threading.enumerate() if thread.is_alive()
        )

        if worker_threads_alive:
            root.after(500, check_solving_results, result_queue)
        else:
            # When ALL solving finished:
            if calibrate_btn:
                calibrate_btn.config(state='normal')
            else:
                log_message("⚠️ Calibrate button not yet available to re-enable.")
            progress_bar.stop()
            root.update_idletasks()
            progress_bar.config(mode="determinate")
            progress_var.set(100)
            progress_label_var.set("Idle")
            fade_out_progress_label()


            log_message(f"✅ Plate solving complete.")

            # Now choose session name for main window
            if session_names_collected:
                from collections import Counter
                if session_names_collected:
                    most_common_name, _ = Counter(session_names_collected).most_common(1)[0]
                else:
                    most_common_name = "UnknownObject"
                    log_message("⚠️ No session names collected. Using fallback.")
                session_title_var.set(most_common_name)
                # Lookup description and distance from local database
            # Lookup description and distance from local database
            try:
                def fix_catalog_name(name):
                    name = name.strip().title()
                    name = name.replace("Ngc", "NGC").replace("Ic", "IC").replace("Messier", "Messier")
                    return name

                lookup_name = fix_catalog_name(most_common_name)

                info = object_info.get(lookup_name)

                if info:
                    object_description_var.set(info[0])
                    object_distance_var.set(f"Distance: {info[1]}")
                    log_message(f"🔭 {lookup_name}: {info[0]}, {info[1]} away")
                else:
                    object_description_var.set("No description available")
                    object_distance_var.set("Unknown distance")
                    log_message(f"⚠️ Object '{lookup_name}' not found in database.")

                log_message(f"📅 Final Imaging Session: {lookup_name}")
                maintain_object_info()

            except Exception as e:
                log_message(f"💥 Failed to load object info: {e}")

            # Fallback if solving failed completely
            if session_title_var.get() == "Welcome to Astrocalibrator!":
                try:
                    first_light = next(iter(light_files), None)
                    if first_light:
                        with fits.open(first_light) as hdul:
                            header = hdul[0].header
                            fallback_name = generate_fallback_name(header)
                    else:
                        fallback_name = generate_fallback_name()
                except Exception as e:
                    print(f"⚠️ Could not read FITS for fallback: {e}")
                    fallback_name = generate_fallback_name()

                session_title_var.set(fallback_name)
                log_message(f"📅 Fallback Imaging Session set to: {fallback_name}")

                try:
                    # Clean up solved and unsolved tracking after solving
                    solved_frames_wcs.clear()
                    unsolved_frames.clear()
                    
                    if os.path.exists(solve_temp_folder):
                        shutil.rmtree(solve_temp_folder)
                        log_message("🧹 Temporary solve files cleaned up.")
                    else:
                        log_message("ℹ️ Solve temp folder already cleaned.")
                except Exception as e:
                    log_message(f"⚠️ Failed to clean solve temp folder: {e}")

    root.after(500, check_solving_results, result_queue)

def solve_then_calibrate(result_queue):
    def worker():
        try:
            _calibration_worker()
        except Exception as e:
            log_message(f"💥 Exception in solve_then_calibrate: {e}")
    threading.Thread(target=worker, daemon=True).start()
    root.after(500, check_solve_and_calibrate_results, result_queue)

def check_solve_and_calibrate_results(result_queue):
    try:
        while True:
            session_name = result_queue.get_nowait()
            if session_name:
                session_title_var.set(session_name)
                print(f"🔍 Looking for session name: '{session_title_var.get()}'")
                print(f"🔍 object_info keys: {list(object_info.keys())}")
                info = object_info.get(session_name)
                if info:
                    object_description_var.set(info[0])
                    object_distance_var.set(f"Distance: {info[1]}")
                else:
                    object_description_var.set("No description available")
                    object_distance_var.set("Unknown distance")
    except queue.Empty:
        pass

    if threading.active_count() > 1:
        root.after(500, check_solve_and_calibrate_results, result_queue)
    else:
        if calibrate_btn:
            calibrate_btn.config(state='normal')
        else:
            log_message("⚠️ Calibrate button not yet available to re-enable.")
        log_message(f"✅ Solve + Calibrate complete.")
        # 🧹 Delete leftover WCS, INI, BAK files now that both solving and calibration are complete
        try:
            for pattern in ["*.wcs", "*.ini", "*.bak"]:
                for leftover in glob.glob(os.path.join(output_folder_var.get(), pattern)):
                    os.remove(leftover)
                    print(f"🧹 Deleted leftover file: {leftover}")
            log_message("🧹 Deleted leftover .wcs, .ini, and .bak files.")
        except Exception as e:
            log_message(f"⚠️ Failed to clean leftover sidecar/config files: {e}")

def start_processing():
    if not output_folder_var.get():
        mb.showwarning("No Output Folder Selected", "⚠️ Please select an output folder before processing.")
        return

    result_queue = queue.Queue()
    calibrate_btn.config(state='disabled')
    solve_then_calibrate(result_queue)

def disable_file_buttons():
    light_btn.config(state='disabled')
    dark_btn.config(state='disabled')
    flat_btn.config(state='disabled')
    darkflat_btn.config(state='disabled')
    bias_btn.config(state='disabled')
    master_dark_btn.config(state='disabled')
    master_flat_btn.config(state='disabled')
    master_bias_btn.config(state='disabled')
    calibrate_btn.config(state='disabled')
    reset_btn.config(state='disabled')

def fade_out_progress_label():
    try:
        text = progress_label_var.get()
        if text not in ["", "Idle"]:
            # Already doing something, don't fade
            return
        
        def gradual_fade(step=0):
            if step > 10:
                progress_label_var.set("")
                return
            opacity = 1.0 - (step / 10)
            progress_label_var.set("")  # Just clear the text (no fancy fade if no label widget exists)
            root.after(100, lambda: gradual_fade(step + 1))

        gradual_fade()

    except Exception as e:
        log_message(f"⚠️ Fade out failed: {e}")

if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()

    initialize_main_widgets()
    calibrate_btn.config(command=start_processing)
    root.after(100, disable_file_buttons)
    root.mainloop()



