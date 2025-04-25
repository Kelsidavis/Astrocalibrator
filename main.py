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

from gui import root, log_message, log_textbox, log_message, output_folder_var, progress_var
from gui import session_title_var
from calibration import run_parallel_calibration, load_fits_by_filter, create_master_frame, save_master_frame
from solving import plate_solve_and_update_header
from settings import load_settings, save_settings, remember_file, get_remembered_file

light_files, dark_files, flat_files, dark_flat_files = [], [], [], []

def select_files(file_list, label):
    files = filedialog.askopenfilenames(filetypes=[("FITS files", "*.fits")])
    if files:
        file_list.clear()
        file_list.extend(files)
        label.config(text=f"{len(files)} files selected")

file_frame = tk.Frame(root)
file_frame.pack(pady=10)

btns = [
    ("Select Light Frames", light_files),
    ("Select Dark Frames", dark_files),
    ("Select Flat Frames", flat_files),
    ("Select Dark Flats", dark_flat_files),
]
labels = []
for i, (label_text, store) in enumerate(btns):
    var = tk.Label(file_frame, text="No files selected")
    var.grid(row=i, column=1, padx=5, pady=2, sticky='w')
    btn = tk.Button(file_frame, text=label_text, command=lambda s=store, l=var: select_files(s, l))
    btn.grid(row=i, column=0, padx=5, pady=2, sticky='w')
    labels.append(var)

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
    dark_flat_images = load_fits_by_filter(dark_flat_files)

    first_light_path = next(iter(next(iter(light_images.values()))), None)

    if not first_light_path:
        log_message("❌ No light frames found.")
        return

    parent_folder = os.path.abspath(os.path.join(os.path.dirname(first_light_path), os.pardir))
    header = fits.getheader(first_light_path)
    object_name = session_title_var.get().replace(' ', '_') or 'UnknownObject'
    date_obs = header.get('DATE-OBS', datetime.now().strftime('%Y-%m-%d')).replace(':', '-').split('T')[0]
    zip_name = f"{object_name}_{date_obs}.zip"
    temp_output_folder = os.path.join(parent_folder, f"_temp_calib_{object_name}_{date_obs}")
    os.makedirs(temp_output_folder, exist_ok=True)

    if save_masters_var.get():
        log_message("💾 Saving master calibration frames...")
        with concurrent.futures.ThreadPoolExecutor() as executor:
            master_darks = dict(executor.map(lambda f: (f, create_master_frame(dark_images[f], method)), dark_images))
            master_flats = dict(executor.map(lambda f: (f, create_master_frame(flat_images[f], method)), flat_images))
            master_dark_flats = dict(executor.map(lambda f: (f, create_master_frame(dark_flat_images[f], method)), dark_flat_images))

        for f, data in master_darks.items():
            save_master_frame(data, fits.getheader(dark_images[f][0]), temp_output_folder, f"dark_{f}")
        for f, data in master_flats.items():
            save_master_frame(data, fits.getheader(flat_images[f][0]), temp_output_folder, f"flat_{f}")
        for f, data in master_dark_flats.items():
            save_master_frame(data, fits.getheader(dark_flat_images[f][0]), temp_output_folder, f"darkflat_{f}")
    if not first_light_path:
        log_message("📦 Creating archive. This may take a moment...\n")
    progress_bar.config(mode='indeterminate')
    progress_bar.start()

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
            session_name = plate_solve_and_update_header(path, log_message)
            log_message(f"💡 Returned session name: {session_name}")
            result_queue.put(session_name)
        except Exception as e:
            import traceback
            log_message(f"💥 Exception in solve_worker: {e}\n{traceback.format_exc()}")

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

        session_set = False
        for path in light_files_to_solve:
            try:
                log_message(f"🧪 Solving: {path}")
                session_name = plate_solve_and_update_header(path, log_message)
                log_message(f"💡 Returned session name: {session_name}")
                if session_name and not session_set:
                    session_title_var.set(session_name)
                    session_set = True
            except Exception as e:
                import traceback
                log_message(f"💥 Plate solving failed for {path}: {e} {traceback.format_exc()}")

        log_message("⚙️ Plate solving complete. Proceeding to calibration...")
        _calibration_worker()

    threading.Thread(target=solve_then_calibrate, daemon=True).start()

calibrate_btn = tk.Button(control_frame, text="Solve & Calibrate", command=run_solve_and_calibrate)
calibrate_btn.pack(side='left', padx=10)
solve_btn = calibrate_btn  # Alias so both names can be used

if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()
    root.mainloop()
