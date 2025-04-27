import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import os
import webbrowser
from urllib.request import urlopen
from PIL import Image, ImageTk
import io

# Tooltip helper
class ToolTip:
    def __init__(self, widget, text):
        self.widget = widget
        self.text = text
        self.tipwindow = None
        self.id = None
        self.widget.bind("<Enter>", self.schedule_showtip)
        self.widget.bind("<Leave>", self.hidetip)

    def schedule_showtip(self, event=None):
        self.id = self.widget.after(500, self.showtip)  # Delay tooltip popup

    def showtip(self, event=None):
        if self.tipwindow or not self.text:
            return
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 25
        y = y + self.widget.winfo_rooty() + 20
        self.tipwindow = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(True)
        tw.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(tw, text=self.text, justify='left',
                         background="#ffffe0", relief='solid', borderwidth=1,
                         font=("tahoma", "8", "normal"))
        label.pack(ipadx=1)

    def hidetip(self, event=None):
        if self.id:
            self.widget.after_cancel(self.id)
            self.id = None
        if self.tipwindow:
            self.tipwindow.destroy()
            self.tipwindow = None

def save_settings(settings):
    try:
        with open('settings.json', 'w') as f:
            json.dump(settings, f, indent=4)
    except Exception as e:
        print(f"⚠️ Failed to save settings: {e}")


def log_message(msg):
    print(msg)
    log_textbox.after(0, lambda: (
        log_textbox.insert('end', msg + '\n'),
        log_textbox.see('end')
    ))

root = tk.Tk()

# Set custom icon from local file
try:
    icon_image = Image.open("icon.png").resize((32, 32))
    icon_photo = ImageTk.PhotoImage(icon_image)
    root.iconphoto(False, icon_photo)
except Exception as e:
    print(f"Could not load window icon: {e}")

def show_about():
    def start_move(event):
        about_window.x = event.x
        about_window.y = event.y

    def do_move(event):
        x = about_window.winfo_pointerx() - about_window.x
        y = about_window.winfo_pointery() - about_window.y
        about_window.geometry(f"+{x}+{y}")

    about_window = tk.Toplevel(root)
    root_x = root.winfo_x()
    root_y = root.winfo_y()
    about_window.geometry(f"360x240+{root_x}+{root_y}")
    about_window.bind('<Button-1>', start_move)
    about_window.bind('<B1-Motion>', do_move)
    about_window.title("About Astrocalibrator")
    about_window.resizable(False, False)

    try:
        image = Image.open("icon.png").resize((100, 100))
        photo = ImageTk.PhotoImage(image)
        img_label = tk.Label(about_window, image=photo)
        img_label.image = photo
        img_label.pack(pady=5)
    except Exception as e:
        tk.Label(about_window, text="[Logo could not be loaded]").pack()

    tk.Label(about_window, text="Astrocalibrator", font=("Arial", 14, "bold")).pack()
    tk.Label(about_window, text="By Kelsi Davis").pack()
    tk.Label(about_window, text="https://geekastro.dev", fg="blue", cursor="hand2").pack()

    def open_site(e):
        webbrowser.open("https://geekastro.dev")

    about_window.bind_all("<Button-1>", open_site)

root.title("Astrocalibrator")

# Title Frame with Astrocalibrator Icon and Dynamic Session Title
title_frame = tk.Frame(root)
title_frame.pack(pady=(10, 5))

if os.path.exists("icon.png"):
    title_img = Image.open("icon.png").resize((32, 32))
    title_photo = ImageTk.PhotoImage(title_img)
    icon_label = tk.Label(title_frame, image=title_photo)
    icon_label.image = title_photo  # Keep a reference!
    icon_label.pack(side='left', padx=5)

session_title_var = tk.StringVar(value="Welcome to Astrocalibrator!")
session_title_label = tk.Label(title_frame, textvariable=session_title_var, font=("Arial", 16, "bold"))
ToolTip(session_title_label, "Displays the current imaging session name or object being calibrated.")
session_title_label.pack(side='left', padx=5)

# Container for session object description and distance
session_info_frame = tk.Frame(root)
session_info_frame.pack(pady=5)

object_description_var = tk.StringVar(value="")
object_distance_var = tk.StringVar(value="")

description_label = tk.Label(session_info_frame, textvariable=object_description_var, font=("Arial", 9))
description_label.pack()

distance_label = tk.Label(session_info_frame, textvariable=object_distance_var, font=("Arial", 9))
distance_label.pack()


output_folder_var = tk.StringVar()
max_threads_var = tk.IntVar(value=os.cpu_count())
progress_var = tk.DoubleVar()

master_dark_path = tk.StringVar()
master_flat_path = tk.StringVar()
master_bias_path = tk.StringVar()

astap_path_var = tk.StringVar()

# Load settings.json at startup
try:
    import json
    with open('settings.json', 'r') as f:
        user_settings = json.load(f)
except (FileNotFoundError, json.JSONDecodeError):
    user_settings = {}

    # Enforce default ASTAP path if missing
if 'astap_path' not in user_settings:
    user_settings['astap_path'] = "C:/Program Files/ASTAP"
    save_settings(user_settings)

if 'output_folder' not in user_settings:
    user_settings['output_folder'] = ""
    save_settings(user_settings)

# Preload saved ASTAP path if available
if 'astap_path' in user_settings:
    astap_path_var.set(user_settings['astap_path'])


master_dark_enabled = tk.BooleanVar()
master_flat_enabled = tk.BooleanVar()
master_bias_enabled = tk.BooleanVar()

# Frame selector section
light_files, dark_files, flat_files, dark_flat_files, bias_files = [], [], [], [], []

from astropy.io import fits

def select_file(path_var):
    filename = filedialog.askopenfilename(
        title="Select File",
        filetypes=[("FITS files", "*.fits"), ("All files", "*.*")]
    )
    if filename:
        path_var.set(filename)

import threading

def select_files(file_list, label, expected_type=None):
    try:
        title = f"Select {expected_type.title()}" if expected_type else "Select Files"
        files = filedialog.askopenfilenames(title=title, filetypes=[("FITS files", "*.fits")])
        if files:
            root.config(cursor="watch")

            def process_files():
                validated_files = list(files)
                file_list.clear()
                file_list.extend(validated_files)

                def update_ui():
                    label.config(text=f"{len(validated_files)} {expected_type.title()} Selected" if expected_type else f"{len(validated_files)} files selected")

                    if expected_type == "DARK":
                        master_dark_enabled.set(False)
                    elif expected_type == "FLAT":
                        master_flat_enabled.set(False)
                    elif expected_type == "BIAS":
                        master_bias_enabled.set(False)
                    elif expected_type == "DARKFLAT":
                        pass

                    if validated_files:
                        messagebox.showinfo(
                            "Selection Complete",
                            f"{len(validated_files)} {expected_type.title()} selected successfully."
                        )

                    root.config(cursor="")  # ✅ Restore normal cursor

                    # ✅ ✅ ONLY after successfully updating, then start solving
                    if expected_type == "LIGHT":
                        if output_folder_var.get():
                            from main import run_plate_solving  # <-- ADD THIS LINE *HERE* INSIDE THE IF BLOCK
                            log_message("🚀 Lights selected. Preparing plate solving...")
                            run_plate_solving()
                        else:
                            log_message("⚠️ Please select output folder before plate solving.")
                            messagebox.showwarning("Output Folder Needed", "⚠️ Please select an output folder before plate solving.")

                root.after(0, update_ui)

            root.after(10, process_files)

    except Exception as e:
        root.config(cursor="")
        raise e

    
def update_ui():
    light_label.config(text=f"{len(light_files)} lights selected")
    dark_label.config(text=f"{len(dark_files)} darks selected")
    flat_label.config(text=f"{len(flat_files)} flats selected")
    bias_label.config(text=f"{len(flat_files)} flats selected")

frame_select_container = tk.LabelFrame(root, text="Select calibration input frames", font=("Arial", 9))
frame_select_container.pack(pady=10, padx=10, fill='x')
frame_select_info = tk.Label(frame_select_container, text="Choose the raw light, dark, flat, and dark flat frames to be used for this calibration session.", font=("Arial", 8), wraplength=500, justify='center')
frame_select_info.pack(pady=(0, 5))
file_frame = tk.Frame(frame_select_container)
file_frame.pack()

light_label = tk.Label(file_frame, text="No files")
dark_label = tk.Label(file_frame, text="No files")
flat_label = tk.Label(file_frame, text="No files")
darkflat_label = tk.Label(file_frame, text="No files")

light_btn = tk.Button(file_frame, text="Select Lights", command=lambda: select_files(light_files, light_label, expected_type="LIGHT"))
ToolTip(light_btn, "Select your primary astrophotography images (LIGHT frames).")
dark_btn = tk.Button(file_frame, text="Select Darks", command=lambda: select_files(dark_files, dark_label, expected_type="DARK"))
ToolTip(dark_btn, "Select DARK frames: Same exposure as lights but with no light entering the camera.")
flat_btn = tk.Button(file_frame, text="Select Flats", command=lambda: select_files(flat_files, flat_label, expected_type="FLAT"))
ToolTip(flat_btn, "Select FLAT frames: Images of uniform light to correct optical system artifacts.")
darkflat_btn = tk.Button(file_frame, text="Select Dark Flats", command=lambda: select_files(dark_flat_files, darkflat_label, expected_type="DARKFLAT"))
ToolTip(darkflat_btn, "Select DARK FLAT frames: Dark exposures matching flat frame settings.")

bias_label = tk.Label(file_frame, text="No files")

bias_btn = tk.Button(file_frame, text="Select Bias Frames", command=lambda: select_files(bias_files, bias_label, expected_type="BIAS"))
ToolTip(bias_btn, "Select BIAS frames: Fastest exposures to calibrate camera electronics.")

for btn, label in zip([light_btn, dark_btn, flat_btn, darkflat_btn, bias_btn], [light_label, dark_label, flat_label, darkflat_label, bias_label]):
    row = [light_btn, dark_btn, flat_btn, darkflat_btn, bias_btn].index(btn)
    btn.grid(row=row, column=0, padx=5, pady=2, sticky='w')
    label.grid(row=row, column=1, sticky='w')

# Master calibration section
master_frame_container = tk.LabelFrame(root, text="Or use existing calibration masters", font=("Arial", 9))
master_frame_container.pack(pady=10, padx=10, fill='x')
master_frame_info = tk.Label(master_frame_container, text="Enable one or more of the options below to apply pre-generated master calibration files.", font=("Arial", 8), wraplength=500, justify='center')
master_frame_info.pack(pady=(0, 5))
master_frame = tk.Frame(master_frame_container)
master_frame.pack()

master_dark_btn = tk.Button(master_frame, text="Select Master Dark", command=lambda: browse_file(master_dark_path))
ToolTip(master_dark_btn, "Choose a pre-generated Master Dark frame to subtract thermal noise from your lights.")
master_dark_label = tk.Label(master_frame, text="No master selected", wraplength=300)
master_flat_btn = tk.Button(master_frame, text="Select Master Flat", command=lambda: browse_file(master_flat_path))
ToolTip(master_flat_btn, "Choose a Master Flat frame to correct uneven illumination (vignetting/dust).")
master_flat_label = tk.Label(master_frame, text="No master selected", wraplength=300)
master_bias_btn = tk.Button(master_frame, text="Select Master Bias", command=lambda: browse_file(master_bias_path))
ToolTip(master_bias_btn, "Choose a Master Bias frame to correct for electronic readout noise.")
master_bias_label = tk.Label(master_frame, text="No master selected", wraplength=300)

master_dark_btn.grid(row=0, column=0, sticky='w', padx=5, pady=2)
master_dark_label.grid(row=0, column=1, sticky='w')

master_flat_btn.grid(row=1, column=0, sticky='w', padx=5, pady=2)
master_flat_label.grid(row=1, column=1, sticky='w')

master_bias_btn.grid(row=2, column=0, sticky='w', padx=5, pady=2)
master_bias_label.grid(row=2, column=1, sticky='w')

def toggle_input_state():
    light_btn.config(state=tk.NORMAL)

    if not master_dark_path.get():
        master_dark_enabled.set(False)
    if not master_flat_path.get():
        master_flat_enabled.set(False)
    if not master_bias_path.get():
        master_bias_enabled.set(False)

    if master_dark_enabled.get() and master_dark_path.get():
        dark_btn.config(state=tk.DISABLED)
        dark_label.config(text="Master Selected")
    else:
        dark_btn.config(state=tk.NORMAL)
        dark_label.config(text="No files")

    if master_flat_enabled.get() and master_flat_path.get():
        flat_btn.config(state=tk.DISABLED)
        flat_label.config(text="Master Selected")
    else:
        flat_btn.config(state=tk.NORMAL)
        flat_label.config(text="No files")

    darkflat_btn.config(state=tk.NORMAL)

    if master_bias_enabled.get() and master_bias_path.get():
        bias_btn.config(state=tk.DISABLED)
        bias_label.config(text="Master Selected")
    else:
        bias_btn.config(state=tk.NORMAL)
        bias_label.config(text="No files")

def update_master_inputs():
    # Reset visible labels
    light_label.config(text="No files")
    dark_label.config(text="No files")
    flat_label.config(text="No files")
    darkflat_label.config(text="No files")
    bias_label.config(text="No files")
    object_description_var.set("")
    object_distance_var.set("")
    master_dark_label.config(text="No master selected")
    master_flat_label.config(text="No master selected")
    master_bias_label.config(text="No master selected")

    # Clear file lists
    light_files.clear()
    dark_files.clear()
    flat_files.clear()
    dark_flat_files.clear()
    bias_files.clear()

    # Reset internal paths
    master_dark_path.set("")
    master_flat_path.set("")
    master_bias_path.set("")
    output_folder_var.set("")
    session_title_var.set("Welcome to Astrocalibrator!")

    # Reset master enabled states
    master_dark_enabled.set(False)
    master_flat_enabled.set(False)
    master_bias_enabled.set(False)

    # Clear logs
    log_textbox.delete('1.0', tk.END)

    # 🔥 Now safely toggle input state after everything is cleaned
    root.after(10, toggle_input_state)

    # 🔥 Slight delay before showing "Reset Complete" message
    root.after(100, lambda: messagebox.showinfo(
        "Reset Complete",
        "✅ All settings, selections, and logs have been reset."
    ))

reset_btn = tk.Button(root, text="Reset Options", command=lambda: update_master_inputs())
ToolTip(reset_btn, "Clear all selected frames and reset calibration settings to default.")
reset_btn.pack(pady=5)

log_frame = tk.Frame(root)
log_frame.pack(side='bottom', fill='both', expand=True, padx=10, pady=(0, 10))

log_textbox = tk.Text(log_frame, wrap='word', height=10)
log_textbox.pack(side='left', fill='both', expand=True)

scrollbar = ttk.Scrollbar(log_frame, command=log_textbox.yview)
scrollbar.pack(side='right', fill='y')
log_textbox.config(yscrollcommand=scrollbar.set)

def browse_file(var):
    path = filedialog.askopenfilename(title="Select Master Calibration Frame", filetypes=[("FITS files", "*.fits")])
    if path:
        var.set(path)
        filename = os.path.basename(path)
        if var == master_dark_path:
            master_dark_label.config(text=filename)
            master_dark_enabled.set(True)
            dark_label.config(text="Master Selected")
        elif var == master_flat_path:
            master_flat_label.config(text=filename)
            master_flat_enabled.set(True)
            flat_label.config(text="Master Selected")
        elif var == master_bias_path:
            master_bias_label.config(text=filename)
            master_bias_enabled.set(True)
            bias_label.config(text="Master Selected")
        toggle_input_state()
        root.update_idletasks()  # Force immediate redraw here!

def set_astap_location():
    existing_path = astap_path_var.get()
    if existing_path and os.path.isdir(existing_path):
        initial_folder = existing_path
    else:
        initial_folder = "C:/Program Files/ASTAP"

    path = filedialog.askopenfilename(
        title="Select ASTAP Executable",
        initialdir=initial_folder,
        filetypes=[("Executable files", "*.exe"), ("All files", "*.*")]
    )
    if path:
        astap_path_var.set(path)
        user_settings['astap_path'] = os.path.dirname(path)  # Save FOLDER not FILE
        save_settings(user_settings)
        log_message(f"🔧 ASTAP location set: {path}")
        auto_close_message("Settings Saved", "✅ ASTAP location saved successfully!", timeout=2000)

def auto_close_message(title, message, timeout=2000):
    top = tk.Toplevel(root)
    top.title(title)
    top.geometry("300x100")
    top.resizable(False, False)
    label = tk.Label(top, text=message, font=("Arial", 10))
    label.pack(expand=True, pady=20)
    top.after(timeout, top.destroy)
    # Center the popup
    top.update_idletasks()
    x = root.winfo_x() + (root.winfo_width() // 2) - (top.winfo_width() // 2)
    y = root.winfo_y() + (root.winfo_height() // 2) - (top.winfo_height() // 2)
    top.geometry(f"+{x}+{y}")

# Now call UI initialization after everything is defined

# Create top menu bar
menubar = tk.Menu(root)

# Settings Menu
settingsmenu = tk.Menu(menubar, tearoff=0)
settingsmenu.add_command(label="Set ASTAP Location", command=set_astap_location)
menubar.add_cascade(label="Settings", menu=settingsmenu)


# Help Menu
helpmenu = tk.Menu(menubar, tearoff=0)
helpmenu.add_command(label="About", command=show_about)
menubar.add_cascade(label="Help", menu=helpmenu)

# Attach menubar to root window
root.config(menu=menubar)

root.after(100, toggle_input_state)
