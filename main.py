import tkinter as tk
from tkinter import filedialog
import os
import shutil
import threading
import concurrent.futures
import queue
import zipfile
from datetime import datetime
from astropy.io import fits

from gui import root, log_message, log_textbox, output_folder_var, progress_var
from gui import session_title_var, master_dark_path, master_flat_path, master_bias_path
from gui import master_dark_enabled, master_flat_enabled, master_bias_enabled
from gui import ToolTip
from calibration import run_parallel_calibration, load_fits_by_filter, create_master_frame, save_master_frame, calibrate_image
from solving import plate_solve_and_update_header
from settings import load_settings, save_settings, remember_file, get_remembered_file

from gui import light_files, dark_files, flat_files, bias_files
from gui import file_frame, light_label, dark_label, flat_label

from gui import object_description_var, object_distance_var
from object_info import object_info

from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u

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
                    with open(wcs_file, 'r') as f:
                        lines = f.readlines()
                        for line in lines:
                            if 'RA center' in line:
                                ra = float(line.split(':')[1].strip())
                            if 'DEC center' in line:
                                dec = float(line.split(':')[1].strip())
                    if ra is not None and dec is not None:
                        log_message(f"🧭 Sidecar WCS center: RA={ra:.4f}°, Dec={dec:.4f}°")
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


control_frame = tk.Frame(root)
control_frame.pack(pady=10)

# Save Masters + Buttons grouped into frames

save_masters_frame = tk.Frame(control_frame)
save_masters_frame.pack(side='left', padx=(10, 50))

save_masters_var = tk.BooleanVar(value=False)
save_masters_checkbox = tk.Checkbutton(save_masters_frame, text="Save Calibration Masters", variable=save_masters_var)
ToolTip(save_masters_checkbox, "Save generated master dark, flat, and bias calibration frames for future use.")
save_masters_checkbox.pack()

buttons_frame = tk.Frame(control_frame)
buttons_frame.pack(side='left', padx=10)

def select_output_directory():
    path = filedialog.askdirectory(title="Select Output Folder")
    if path:
        output_folder_var.set(path)
        log_message(f"📂 Output folder set to: {path}")

select_output_btn = tk.Button(buttons_frame, text="Select Output Folder", command=select_output_directory)
ToolTip(select_output_btn, "Choose where calibrated and solved files will be saved.")
select_output_btn.pack(side='left', padx=10)

progress_bar = tk.ttk.Progressbar(root, variable=progress_var, maximum=100)
progress_bar.pack(fill='x', padx=10, pady=5)

def log_message(msg):
    print(msg)
    log_textbox.after(0, lambda: (
        log_textbox.insert('end', msg + '\n'),
        log_textbox.see('end')
    ))

def _calibration_worker():
    import time
    start_time = time.time()
    log_message(f"📅 Session: {session_title_var.get()}")
    method = 'median'

    first_light_path = next(iter(light_files), None)
    if not first_light_path:
        log_message("❌ No light frames found.")
        calibrate_btn.config(state='normal')
        solve_btn.config(state='normal')
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
        session_title=session_title_var.get()
    )

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

    calibrate_btn.config(state='normal')
    solve_btn.config(state='normal')

def run_calibration_pipeline():

    calibrate_btn.config(state='disabled')
    solve_btn.config(state='disabled')
    threading.Thread(target=_calibration_worker, daemon=True).start()

def run_plate_solving():
    solve_temp_folder = os.path.join(output_folder_var.get(), "solve_temp")
    os.makedirs(solve_temp_folder, exist_ok=True)
    solve_btn.config(state='disabled')
    calibrate_btn.config(state='disabled')
    log_message("📅 Starting plate solving in background...")

    light_files_to_solve = [f for f in light_files if os.path.exists(f)]
    solver_failed = False

    def solve_worker(path):
        nonlocal solver_failed
        try:
            log_message(f"🧪 Solving: {path}")
            try:
                session_name = plate_solve_and_update_header(path, log_message)
                if not session_name:
                    session_name = find_nearest_known_object(path, object_info)
                if session_name:
                    log_message(f"🔭 WCS matching assigned session name: {session_name}")
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
                log_message(f"💡 Returned session name: {session_name}")
                result_queue.put(session_name)
            except FileNotFoundError as fnf_err:
                if not solver_failed:
                    import tkinter.messagebox as mb
                    mb.showinfo("Plate Solver Not Found", "The plate solver executable could not be found. Calibration will continue without solving.")
                    solver_failed = True
                log_message(f"❌ Solver not found: {fnf_err}")
                session_name = None
                result_queue.put(session_name)
            except Exception as e:
                import traceback
                log_message(f"💥 Exception in solve_worker: {e}\n{traceback.format_exc()}")
        except Exception as outer_err:
            log_message(f"💥 Outer exception in solve_worker: {outer_err}")

    for path in light_files_to_solve:
        threading.Thread(target=solve_worker, args=(path,), daemon=True).start()

    def check_results():
        while not result_queue.empty():
            session_name = result_queue.get()
            if session_name:
                session_title_var.set(session_name)
                info = object_info.get(session_name)
                if info:
                    object_description_var.set(info[0])
                    object_distance_var.set(f"Distance: {info[1]}")
                else:
                    object_description_var.set("No description available")
                    object_distance_var.set("Unknown distance")
            log_message(f"📅 Updated Imaging Session: {session_name}")

        if threading.active_count() > 1:
            root.after(500, check_results)
        else:
            solve_btn.config(state='normal')
            calibrate_btn.config(state='normal')
            log_message(f"✅ Plate solving complete.")

        if threading.active_count() <= 1:
            try:
                shutil.rmtree(solve_temp_folder)
                log_message("🧹 Temporary solve files cleaned up.")
            except Exception as e:
                log_message(f"⚠️ Failed to clean solve temp folder: {e}")

    root.after(500, check_results)

def run_solve_and_calibrate():
    calibrate_btn.config(state='disabled')
    solve_btn.config(state='disabled')

    def solve_then_calibrate():
        solve_temp_folder = os.path.join(output_folder_var.get(), "solve_temp")
        os.makedirs(solve_temp_folder, exist_ok=True)
        light_files_to_solve = [f for f in light_files if os.path.exists(f)]
        solver_failed = False
        session_set = False

        for path in light_files_to_solve:
            try:
                log_message(f"🧪 Solving: {path}")
                session_name = plate_solve_and_update_header(path, log_message)

                if not session_name:
                    session_name = find_nearest_known_object(path, object_info)
                
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

                    log_message(f"💡 Returned session name: {session_name}")

                    if session_name and not session_set:
                        session_title_var.set(session_name)
                        info = object_info.get(session_name)
                        if info:
                            object_description_var.set(info[0])
                            object_distance_var.set(f"Distance: {info[1]}")
                        else:
                            object_description_var.set("No description available")
                            object_distance_var.set("Unknown distance")
                        session_set = True


            except FileNotFoundError as fnf_err:
                if not solver_failed:
                    import tkinter.messagebox as mb
                    mb.showinfo("Plate Solver Not Found", "The plate solver executable could not be found. Calibration will continue without solving.")
                    solver_failed = True
                log_message(f"❌ Solver not found: {fnf_err}")
            except Exception as e:
                import traceback
                log_message(f"💥 Exception in solve_worker: {e}\n{traceback.format_exc()}")
                session_name = None
                result_queue.put(session_name)

        log_message("⚙️ Plate solving complete. Proceeding to calibration...")
        try:
            shutil.rmtree(solve_temp_folder)
            log_message("🧹 Temporary solve files cleaned up.")
        except Exception as e:
            log_message(f"⚠️ Failed to clean solve temp folder: {e}")

        _calibration_worker()

    threading.Thread(target=solve_then_calibrate, daemon=True).start()

calibrate_btn = tk.Button(buttons_frame, text="Solve & Calibrate", command=run_solve_and_calibrate)
ToolTip(calibrate_btn, "Plate solve light frames and apply calibration using selected masters and settings.")
calibrate_btn.pack(side='left', padx=10)
solve_btn = calibrate_btn  # Alias so both names can be used

def debug_widget_list():
    print("\n🧩 Widgets inside file_frame:")
    for child in file_frame.winfo_children():
        print(" -", child, "text=", getattr(child, 'cget', lambda x: 'N/A')('text'))

if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()
    root.mainloop()
