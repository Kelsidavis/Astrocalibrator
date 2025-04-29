import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import os
import shutil
import threading
import concurrent.futures
import queue
import zipfile
from datetime import datetime
from astropy.io import fits

from gui import root, log_message, log_textbox, output_folder_var, progress_var, progress_label_var
from gui import session_title_var, master_dark_path, master_flat_path, master_bias_path
from gui import master_dark_enabled, master_flat_enabled, master_bias_enabled
from gui import ToolTip
from calibration import run_parallel_calibration, load_fits_by_filter, create_master_frame, save_master_frame, calibrate_image
from solving import plate_solve_and_update_header
from settings import load_settings, save_settings, remember_file, get_remembered_file

from gui import light_files, dark_files, flat_files, bias_files
from gui import file_frame, light_label, dark_label, flat_label
from gui import light_btn, dark_btn, flat_btn, bias_btn, darkflat_btn
from gui import reset_btn
from gui import master_dark_btn, master_flat_btn, master_bias_btn
from gui import progress_bar
from gui import object_description_var, object_distance_var
from object_info import object_info

from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
import tkinter.messagebox as mb
import glob

# Globals for buttons
select_output_btn = None
calibrate_btn = None
current_output_label = None
save_masters_var = None

# Prevent duplicate GUI element creation
main_widgets_initialized = False

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

main_widgets_initialized = False

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

def global_cleanup(output_folder):
    """Delete leftover .wcs and .ini files from output folder."""
    patterns = ["*.wcs", "*.ini"]  # file types to delete
    for pattern in patterns:
        for file_path in glob.glob(os.path.join(output_folder, pattern)):
            try:
                os.remove(file_path)
                print(f"🧹 Deleted leftover file: {file_path}")
            except Exception as e:
                print(f"⚠️ Failed to delete {file_path}: {e}")

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

# --- Globals for plate solving status ---
solved_frames_wcs = {}  # Successfully solved frames and their WCS info
unsolved_frames = []    # Frames that failed solving

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
        try:
            ra, dec = None, None

            # Try reading WCS headers first
            with fits.open(fits_path) as hdul:
                hdr = hdul[0].header
                if 'CRVAL1' in hdr and 'CRVAL2' in hdr:
                    ra = hdr['CRVAL1']
                    dec = hdr['CRVAL2']
                    log_message(f"🧭 FITS WCS center: RA={ra:.4f}°, Dec={dec:.4f}°")

            # If WCS headers are missing, try reading the .wcs sidecar file
            if ra is None or dec is None:
                wcs_file = os.path.splitext(fits_path)[0] + '.wcs'
                if os.path.exists(wcs_file):
                    try:
                        with fits.open(wcs_file) as hdul:
                            hdr = hdul[0].header
                            ra = hdr.get('CRVAL1')
                            dec = hdr.get('CRVAL2')
                        if ra is not None and dec is not None:
                            log_message(f"🧭 Sidecar WCS center (from FITS): RA={ra:.4f}°, Dec={dec:.4f}°")
                    except Exception as e:
                        log_message(f"⚠️ Failed to read sidecar WCS: {e}")

            if ra is not None and dec is not None:
                solved_coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')

                closest_object = None
                min_distance = float('inf')

                for name, (desc, dist, obj_ra, obj_dec) in catalog.items():
                    obj_coord = SkyCoord(ra=obj_ra*u.deg, dec=obj_dec*u.deg, frame='icrs')
                    sep = solved_coord.separation(obj_coord)
                    if sep.degree < min_distance:
                        min_distance = sep.degree
                        closest_object = name

                if closest_object and min_distance < 8:  # 8 degrees tolerance
                    return closest_object

        except Exception as e:
            log_message(f"⚠️ WCS matching failed: {e}")

        return None    

def _calibration_worker():
    import time
    start_time = time.time()
    progress_label_var.set("Calibrating frames...")

    # Restore object info if cached
    if cached_object_description:
        object_description_var.set(cached_object_description)
    if cached_object_distance:
        object_distance_var.set(cached_object_distance)

    progress_bar.config(mode="indeterminate")
    progress_bar.start(10)
    method = 'median'

    first_light_path = next(iter(light_files), None)
    if not first_light_path:
        log_message("❌ No light frames found.")
        if calibrate_btn:
            calibrate_btn.config(state='normal')
        return

    output_folder = output_folder_var.get()
    temp_folder = os.path.join(output_folder, "temp")
    os.makedirs(temp_folder, exist_ok=True)

    run_parallel_calibration(
        light_images=light_files,
        dark_images=dark_files,
        flat_images=flat_files,
        bias_images=bias_files,
        output_folder=output_folder,
        session_title=session_title_var.get(),
        log_callback=log_message
    )

    progress_label_var.set("Complete!")
    elapsed = time.time() - start_time
    log_message(f"✅ Calibration complete in {elapsed:.2f} seconds.")

    if not save_masters_var.get():
        try:
            shutil.rmtree(temp_folder)
            log_message("🧹 Temporary calibration files cleaned up.")
        except Exception as e:
            log_message(f"⚠️ Failed to clean temp folder: {e}")
    else:
        log_message("💾 Temporary calibration files retained (Save Masters enabled).")

    progress_bar.stop()
    progress_bar.config(mode="determinate")
    progress_var.set(100)

    calibrate_btn.config(state='normal')
    progress_label_var.set("Idle")
    fade_out_progress_label()

    # Restore cached object info
    if cached_object_description:
        object_description_var.set(cached_object_description)
    if cached_object_distance:
        object_distance_var.set(cached_object_distance)

    # --- New: Create output ZIP archive with session name and imaging date ---
    try:
        final_output_folder = output_folder_var.get()
        session_name_cleaned = session_title_var.get().replace(' ', '_').replace(':', '').replace('/', '_')

        # 🗓 Extract DATE-OBS from the first light frame
        first_light = next(iter(light_files), None)
        if first_light:
            with fits.open(first_light) as hdul:
                header = hdul[0].header
                date_obs = header.get('DATE-OBS') or header.get('DATE') or header.get('DATEOBS')
                if date_obs:
                    imaging_date = date_obs.split('T')[0]  # Only YYYY-MM-DD
                else:
                    imaging_date = datetime.now().strftime("%Y-%m-%d")
        else:
            imaging_date = datetime.now().strftime("%Y-%m-%d")

        zip_filename = os.path.join(final_output_folder, f"{session_name_cleaned}_{imaging_date}_calibrated.zip")

        with zipfile.ZipFile(zip_filename, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root_dir, _, files in os.walk(final_output_folder):
                for file in files:
                    if file.endswith('.fits'):
                        file_path = os.path.join(root_dir, file)
                        arcname = os.path.relpath(file_path, final_output_folder)
                        zipf.write(file_path, arcname)

        log_message(f"📦 Created ZIP archive: {zip_filename}")

    except Exception as e:
        log_message(f"⚠️ Failed to create ZIP archive: {e}")

def run_calibration_pipeline():

    calibrate_btn.config(state='disabled')
    threading.Thread(target=_calibration_worker, daemon=True).start()

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
        try:
            log_message(f"🧪 Solving: {path}")
            try:
                print(f"👣 Entering solve_worker() for path: {path}")
                print(f"🛤️ Checking if file exists: {os.path.exists(path)}")
                session_name_from_solver = plate_solve_and_update_header(path, log_message)
                session_name_nearby = find_nearest_known_object(path, object_info)

                # Pick the best session name
                if session_name_nearby:
                    session_name = session_name_nearby
                    log_message(f"🎯 Overriding solver object '{session_name_from_solver}' with nearby known object '{session_name}'")
                else:
                    session_name = session_name_from_solver
                    log_message(f"📋 Keeping solver object: '{session_name}'")

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

                    info = object_info.get(session_name)
                    if info:
                        object_description_var.set(info[0])
                        object_distance_var.set(f"Distance: {info[1]}")
                        cached_object_description = info[0]
                        cached_object_distance = f"Distance: {info[1]}"
                        log_message(f"🔭 Object found: {session_name} - {info[0]}, {info[1]}")
                    else:
                        object_description_var.set("No description available")
                        object_distance_var.set("Unknown distance")
                        cached_object_description = "No description available"
                        cached_object_distance = "Unknown distance"
                        log_message(f"⚠️ Object '{session_name}' not found in database.")

                result_queue.put(session_name)
                session_names_collected.append(session_name)

                # 🛠 Debug FITS header right after solving
                try:
                    from astropy.io import fits  # local import inside thread
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

            
    print(f"💬 Lights to solve: {light_files_to_solve}")
    print(f"💬 Total light frames selected: {len(light_files_to_solve)}")

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
                        # --- Estimate WCS for Unsolved Frames ---
            if unsolved_frames and solved_frames_wcs:
                reference_path, reference_wcs = next(iter(solved_frames_wcs.items()))
                log_message(f"🛠 Estimating WCS for {len(unsolved_frames)} unsolved frames using reference: {os.path.basename(reference_path)}")

                for failed_path in unsolved_frames:
                    try:
                        with fits.open(failed_path, mode='update') as hdul:
                            hdr = hdul[0].header
                            hdr['CRVAL1'] = reference_wcs['ra']
                            hdr['CRVAL2'] = reference_wcs['dec']
                            hdr['CRPIX1'] = reference_wcs['width'] / 2
                            hdr['CRPIX2'] = reference_wcs['height'] / 2
                            hdr['CD1_1'] = -reference_wcs['scale'] / 3600
                            hdr['CD1_2'] = 0.0
                            hdr['CD2_1'] = 0.0
                            hdr['CD2_2'] = reference_wcs['scale'] / 3600
                            hdr['CTYPE1'] = 'RA---TAN'
                            hdr['CTYPE2'] = 'DEC--TAN'
                            hdr['CUNIT1'] = 'deg'
                            hdr['CUNIT2'] = 'deg'
                            hdr['WCS_EST'] = (True, "Estimated WCS based on nearby solved frame")
                            hdul.flush()
                        log_message(f"✅ Estimated WCS injected into: {os.path.basename(failed_path)}")
                    except Exception as e:
                        log_message(f"⚠️ Failed to inject WCS into {failed_path}: {e}")
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
                most_common_name, _ = Counter(session_names_collected).most_common(1)[0]
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
                shutil.rmtree(solve_temp_folder)
                log_message("🧹 Temporary solve files cleaned up.")
            except Exception as e:
                log_message(f"⚠️ Failed to clean solve temp folder: {e}")

    root.after(500, check_solving_results, result_queue)

def solve_then_calibrate(result_queue):
    def worker():
        try:
            run_plate_solving()
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
        global_cleanup(output_folder_var.get())
        log_message("🧹 Global cleanup complete.")

def start_processing():
    if not output_folder_var.get():
        mb.showwarning("No Output Folder Selected", "⚠️ Please select an output folder before processing.")
        return

    result_queue = queue.Queue()
    calibrate_btn.config(state='disabled')
    solve_then_calibrate(result_queue)

def debug_widget_list():
    print("\n🧩 Widgets inside file_frame:")
    for child in file_frame.winfo_children():
        print(" -", child, "text=", getattr(child, 'cget', lambda x: 'N/A')('text'))

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

def enable_plate_solving():
    from gui import try_run_plate_solving
    try_run_plate_solving()

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
            color = f"#{int(128 * opacity + 127):02x}{int(128 * opacity + 127):02x}{int(128 * opacity + 127):02x}"
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
    root.after(500, enable_plate_solving)
    root.after(100, disable_file_buttons)
    root.mainloop()



