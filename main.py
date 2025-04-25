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
from calibration import run_parallel_calibration, load_fits_by_filter, create_master_frame, save_master_frame, calibrate_image
from solving import plate_solve_and_update_header
from settings import load_settings, save_settings, remember_file, get_remembered_file

from gui import light_files, dark_files, flat_files, bias_files



from gui import file_frame, light_label, dark_label, flat_label

control_frame = tk.Frame(root)
control_frame.pack(pady=10)

save_masters_var = tk.BooleanVar(value=False)
save_masters_checkbox = tk.Checkbutton(control_frame, text="Save Calibration Masters", variable=save_masters_var)
save_masters_checkbox.pack(side='left', padx=10)

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

    light_images = load_fits_by_filter(light_files)
    dark_images = load_fits_by_filter(dark_files)
    flat_images = load_fits_by_filter(flat_files)
    bias_images = load_fits_by_filter(bias_files)

    first_light_path = next(iter(next(iter(light_images.values()))), None)
    if not first_light_path:
        log_message("❌ No light frames found.")
        return

    parent_folder = os.path.abspath(os.path.join(os.path.dirname(first_light_path), os.pardir))
    header = fits.getheader(first_light_path)
    object_name = session_title_var.get().replace(' ', '_').replace(':', '_') or 'UnknownObject'
    date_obs = header.get('DATE-OBS', datetime.now().strftime('%Y-%m-%d')).replace(':', '-').split('T')[0]
    zip_name = f"{object_name}_{date_obs}.zip"
    temp_output_folder = os.path.join(parent_folder, f"_temp_calib_{object_name}_{date_obs}")
    os.makedirs(temp_output_folder, exist_ok=True)

    master_dark = None
    master_bias = None
    master_flats = {}

    if save_masters_var.get():
        log_message("💾 Saving master calibration frames...")
        with concurrent.futures.ThreadPoolExecutor() as executor:
            if dark_files:
                master_dark = create_master_frame(dark_files, method)
                save_master_frame(master_dark, fits.getheader(dark_files[0]), temp_output_folder, "master_dark")

            if bias_files:
                master_bias = create_master_frame(bias_files, method)
                save_master_frame(master_bias, fits.getheader(bias_files[0]), temp_output_folder, "master_bias")

            if flat_images:
                corrected_flats = {}
                for f, images in flat_images.items():
                    corrected_flats[f] = []
                    for path in images:
                        with fits.open(path) as hdul:
                            data = hdul[0].data.astype(float)
                            if master_bias is not None:
                                data -= master_bias
                            if master_dark is not None:
                                data -= master_dark
                            corrected_flats[f].append(data)

                master_flats = dict(executor.map(lambda f: (f, create_master_frame(corrected_flats[f], method)), corrected_flats))
                for f, data in master_flats.items():
                    save_master_frame(data, fits.getheader(flat_images[f][0]), temp_output_folder, f"flat_{f}")

    log_message("🔧 Calibrating light frames...")

    progress_bar.config(mode='indeterminate')
    progress_bar.start()

    for filter_name, images in light_images.items():
        for path in images:
            with fits.open(path) as hdul:
                light_data = hdul[0].data.astype(float)
                light_header = hdul[0].header

            calibrated_data = calibrate_image(
                path,
                use_master=True,
                master_dark_path=os.path.join(temp_output_folder, "master_dark.fits") if master_dark is not None else master_dark_path.get(),
                master_flat_path=os.path.join(temp_output_folder, f"flat_{filter_name}.fits") if filter_name in master_flats else master_flat_path.get(),
                master_bias_path=os.path.join(temp_output_folder, "master_bias.fits") if master_bias is not None else master_bias_path.get()
            )

            light_header['CALIB'] = (True, "Frame has been dark and flat calibrated")

            fits.writeto(os.path.join(temp_output_folder, os.path.basename(path)), calibrated_data, header=light_header, overwrite=True)

    zip_path = os.path.join(parent_folder, zip_name)
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as archive:
        for root_dir, _, files in os.walk(temp_output_folder):
            for file in files:
                full_path = os.path.join(root_dir, file)
                arcname = os.path.relpath(full_path, temp_output_folder)
                archive.write(full_path, arcname)

    shutil.rmtree(temp_output_folder)
    progress_bar.stop()
    progress_bar.config(mode='determinate')
    progress_var.set(100)
    log_message(f"📦 Archive created: {zip_path}/n")
    try:
        os.startfile(os.path.dirname(zip_path))
    except Exception as e:
        log_message(f"⚠️ Failed to open folder: {e}")
    elapsed = time.time() - start_time
    log_message(f"✅ Calibration complete in {elapsed:.2f} seconds.")
    calibrate_btn.config(state='normal')
    solve_btn.config(state='normal')
    try:
        os.startfile(os.path.dirname(zip_path))
    except Exception as e:
        log_message(f"⚠️ Failed to open folder: {e}")
    elapsed = time.time() - start_time
    log_message(f"✅ Calibration complete in {elapsed:.2f} seconds.")
    calibrate_btn.config(state='normal')
    solve_btn.config(state='normal')

result_queue = queue.Queue()

def run_calibration_pipeline():
    calibrate_btn.config(state='disabled')
    solve_btn.config(state='disabled')
    threading.Thread(target=_calibration_worker, daemon=True).start()

def run_plate_solving():
    solve_btn.config(state='disabled')
    calibrate_btn.config(state='disabled')
    log_message("📅 Starting plate solving in background...")

    light_files_to_solve = [f for f in light_files if os.path.exists(f)]

    def solve_worker(path):
        try:
            log_message(f"🧪 Solving: {path}")
            try:
                session_name = plate_solve_and_update_header(path, log_message)
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
            session_title_var.set(session_name)
            log_message(f"📅 Updated Imaging Session: {session_name}")

        if threading.active_count() > 1:
            root.after(500, check_results)
        else:
            solve_btn.config(state='normal')
            calibrate_btn.config(state='normal')
            log_message(f"✅ Plate solving complete.")

    root.after(500, check_results)

def run_solve_and_calibrate():
    calibrate_btn.config(state='disabled')
    solve_btn.config(state='disabled')

    def solve_then_calibrate():
        light_files_to_solve = [f for f in light_files if os.path.exists(f)]
        solver_failed = False
        session_set = False

        for path in light_files_to_solve:
            try:
                log_message(f"🧪 Solving: {path}")
                session_name = plate_solve_and_update_header(path, log_message)
                log_message(f"💡 Returned session name: {session_name}")
                if session_name and not session_set:
                    session_title_var.set(session_name)
                    session_set = True
            except FileNotFoundError as fnf_err:
                if not solver_failed:
                    import tkinter.messagebox as mb
                    mb.showinfo("Plate Solver Not Found", "The plate solver executable could not be found. Calibration will continue without solving.")
                    solver_failed = True
                log_message(f"❌ Solver not found: {fnf_err}")
            except Exception as e:
                import traceback
                log_message(f"💥 Plate solving failed for {path}: {e}/n{traceback.format_exc()}")

        log_message("⚙️ Plate solving complete. Proceeding to calibration...")
        _calibration_worker()

    threading.Thread(target=solve_then_calibrate, daemon=True).start()

calibrate_btn = tk.Button(control_frame, text="Solve & Calibrate", command=run_solve_and_calibrate)
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
