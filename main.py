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
from gui import light_btn, dark_btn, flat_btn, bias_btn, darkflat_btn
from gui import reset_btn
from gui import master_dark_btn, master_flat_btn, master_bias_btn

from gui import object_description_var, object_distance_var
from object_info import object_info

from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
import tkinter.messagebox as mb

from datetime import datetime

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

def wiggle_button(widget):
    def move_left():
        widget.place_configure(x=widget.winfo_x() - 5)
        widget.after(50, move_right)

    def move_right():
        widget.place_configure(x=widget.winfo_x() + 10)
        widget.after(50, move_center)

    def move_center():
        widget.place_configure(x=widget.winfo_x() - 5)

    move_left()

# after all imports and small functions like wiggle_button()

def select_output_directory():
    path = filedialog.askdirectory(title="Select Output Folder")
    if path:
        output_folder_var.set(path)
        log_message(f"📂 Output folder set to: {path}")
        select_output_btn.config(font=("Arial", 10), width=18, height=1)
        calibrate_btn.config(font=("Arial", 14, "bold"), width=25, height=3)
        light_btn.config(state='normal')
        dark_btn.config(state='normal')
        flat_btn.config(state='normal')
        darkflat_btn.config(state='normal')
        bias_btn.config(state='normal')
        calibrate_btn.config(state='normal')
        reset_btn.config(state='normal')
        master_dark_btn.config(state='normal')
        master_flat_btn.config(state='normal')
        master_bias_btn.config(state='normal')

# ✅ Now create output folder frame ONCE
output_folder_frame = tk.Frame(root)
output_folder_frame.pack(pady=(5, 0))

select_output_btn = tk.Button(
    output_folder_frame,
    text="Select Output Folder",
    font=("Arial", 12, "bold"),
    width=25,
    height=2,
    command=select_output_directory  # now this works because function is defined
)
ToolTip(select_output_btn, "Choose where calibrated and solved files will be saved.")
select_output_btn.pack(padx=10, pady=5)

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

calibrate_btn = tk.Button(buttons_frame, text="Calibrate Files", font=("Arial", 12, "bold"), width=20, height=2)
ToolTip(calibrate_btn, "Plate solve light frames and apply calibration using selected masters and settings.")
calibrate_btn.pack(side='left', padx=10)
solve_btn = calibrate_btn  # Alias so both names can be used

progress_bar = tk.ttk.Progressbar(root, variable=progress_var, maximum=100)
progress_bar.pack(fill='x', padx=10, pady=5)

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
    result_queue = queue.Queue()
    solve_temp_folder = os.path.join(output_folder_var.get(), "solve_temp")
    os.makedirs(solve_temp_folder, exist_ok=True)
    solve_btn.config(state='disabled')
    calibrate_btn.config(state='disabled')
    log_message("📅 Starting plate solving in background...")

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
                session_name = plate_solve_and_update_header(path, log_message)
                
                # 🛠 Debug FITS header right after solving
                try:
                    from astropy.io import fits  # local import inside thread
                    with fits.open(path) as hdul:
                        hdr = hdul[0].header
                        print("🔍 FITS Header keys:", list(hdr.keys()))
                        print("🔍 FITS OBJECT field:", hdr.get('OBJECT'))
                except Exception as e:
                    print(f"⚠️ Could not read FITS header: {e}")
    
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
                result_queue.put(session_name)
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

    def check_solving_results():
        try:
            while True:
                session_name = result_queue.get_nowait()
                if session_name:
                    session_name = session_name.strip().title()
                    session_title_var.set(session_name)
                    info = object_info.get(session_name)
                    if info:
                        object_description_var.set(info[0])
                        object_distance_var.set(f"Distance: {info[1]}")
                    else:
                        object_description_var.set("No description available")
                        object_distance_var.set("Unknown distance")
                    log_message(f"📅 Updated Imaging Session: {session_name}")
        except queue.Empty:
            pass

        worker_threads_alive = any(
            thread.name.startswith("Thread-") for thread in threading.enumerate() if thread.is_alive()
        )

        if worker_threads_alive:
            root.after(500, check_solving_results)
        else:
            solve_btn.config(state='normal')
            calibrate_btn.config(state='normal')
            log_message(f"✅ Plate solving complete.")

            # Fallback session title if solving failed
            if session_title_var.get() == "Welcome to Astrocalibrator!":
                try:
                    from astropy.io import fits
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
                shutil.rmtree(solve_temp_folder)
                log_message("🧹 Temporary solve files cleaned up.")
            except Exception as e:
                log_message(f"⚠️ Failed to clean solve temp folder: {e}")

    root.after(500, check_solving_results)

def start_processing():
    if not output_folder_var.get():
        mb.showwarning("No Output Folder Selected", "⚠️ Please select an output folder before processing.")
        wiggle_button(select_output_btn)
        return
    result_queue = queue.Queue()

    def solve_then_calibrate(result_queue):
        try:
            solve_images(result_queue)
            calibrate_images(result_queue)
        except Exception as e:
            log_message(f"💥 Exception in solve_then_calibrate: {e}")

    def check_solve_and_calibrate_results(result_queue):
        try:
            while True:
                session_name = result_queue.get_nowait()
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
        except queue.Empty:
            pass

        if threading.active_count() > 1:
            root.after(500, check_solve_and_calibrate_results, result_queue)
        else:
            solve_btn.config(state='normal')
            calibrate_btn.config(state='normal')
            log_message(f"✅ Solve + Calibrate complete.")

    calibrate_btn.config(state='disabled')
    solve_btn.config(state='disabled')

    threading.Thread(target=solve_then_calibrate, args=(result_queue,), daemon=True).start()
    root.after(500, check_solve_and_calibrate_results, result_queue)

calibrate_btn.config(command=start_processing)

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

if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()

    root.after(100, disable_file_buttons)  # 🛠 Schedule button disabling after GUI loads
    root.mainloop()


