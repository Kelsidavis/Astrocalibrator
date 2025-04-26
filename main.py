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

control_frame = tk.Frame(root)
control_frame.pack(pady=10)

# Setup Menu
menubar = tk.Menu(root)
root.config(menu=menubar)

# Settings Menu
settings_menu = tk.Menu(menubar, tearoff=0)
menubar.add_cascade(label="Settings", menu=settings_menu)

def set_astap_location():
    path = filedialog.askopenfilename(title="Select ASTAP Executable")
    if path:
        remember_file('astap_path', path)
        log_message(f"🔧 ASTAP location set: {path}")

settings_menu.add_command(label="Set ASTAP Location", command=set_astap_location)

# Help Menu
help_menu = tk.Menu(menubar, tearoff=0)
menubar.add_cascade(label="Help", menu=help_menu)

def open_readme():
    import webbrowser
    webbrowser.open("https://github.com/Kelsidavis/astrocalibrator#readme")

help_menu.add_command(label="View README", command=open_readme)
help_menu.add_separator()

def open_about():
    about_window = tk.Toplevel(root)
    about_window.title("About Astrocalibrator")
    about_window.resizable(False, False)

    # Position about window near main root window
    x = root.winfo_x()
    y = root.winfo_y()
    about_window.geometry(f"260x360+{x+50}+{y+50}")

    if os.path.exists("icon.png"):
        icon_img = tk.PhotoImage(file="icon.png")
        about_window.iconphoto(False, icon_img)
        small_img = icon_img.subsample(max(1, int(icon_img.width()/128)), max(1, int(icon_img.height()/128)))
        icon_label = tk.Label(about_window, image=small_img)
        icon_label.image = small_img
        icon_label.pack(pady=10)

    text_label = tk.Label(
        about_window,
        text="Astrocalibrator v1.0 \n Calibrate and solve astronomical images. \n Created by Kelsi Davis.",
        justify="center",
        wraplength=250
    )
    text_label.pack(pady=10)

    def open_website():
        import webbrowser
        webbrowser.open("https://geekastro.dev")

    link_button = tk.Button(about_window, text="Official Website", command=open_website)
    link_button.pack(pady=5)

    close_button = tk.Button(about_window, text="Close", command=about_window.destroy)
    close_button.pack(pady=10)

help_menu.add_command(label="About Astrocalibrator", command=open_about)

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
                if session_name and session_name.startswith('M') and session_name[1:].isdigit():
                    session_name = f"Messier {session_name[1:]}"
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
                info = object_info.get(session_name, ("", ""))
                object_description_var.set(info[0])
                object_distance_var.set(f"Distance: {info[1]}" if info[1] else "")
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
                if session_name and session_name.startswith('M') and session_name[1:].isdigit():
                    session_name = f"Messier {session_name[1:]}"
                log_message(f"💡 Returned session name: {session_name}")
                if session_name and not session_set:
                    session_title_var.set(session_name)
                    info = object_info.get(session_name, ("", ""))
                    object_description_var.set(info[0])
                    object_distance_var.set(f"Distance: {info[1]}" if info[1] else "")
                    session_set = True
            except FileNotFoundError as fnf_err:
                if not solver_failed:
                    import tkinter.messagebox as mb
                    mb.showinfo("Plate Solver Not Found", "The plate solver executable could not be found. Calibration will continue without solving.")
                    solver_failed = True
                log_message(f"❌ Solver not found: {fnf_err}")
            except Exception as e:
                import traceback
                log_message(f"💥 Plate solving failed for {path}: {e} \n {traceback.format_exc()}")

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
