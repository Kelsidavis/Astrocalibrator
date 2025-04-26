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
        self.widget.bind("<Enter>", self.showtip)
        self.widget.bind("<Leave>", self.hidetip)

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
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

def log_message(msg):
    print(msg)
    log_textbox.after(0, lambda: (
        log_textbox.insert('end', msg + '\n'),
        log_textbox.see('end')
    ))

root = tk.Tk()

# Set custom telescope icon from local file
try:
    icon_image = Image.open("icon.png").resize((32, 32))
    icon_photo = ImageTk.PhotoImage(icon_image)
    root.iconphoto(False, icon_photo)
except Exception as e:
    print(f"Could not load window icon: {e}")

# Add top menu bar
menubar = tk.Menu(root)

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

helpmenu = tk.Menu(menubar, tearoff=0)
helpmenu.add_command(label="About", command=show_about)
menubar.add_cascade(label="Help", menu=helpmenu)

root.config(menu=menubar)
root.title("Astrocalibrator")

session_title_var = tk.StringVar(value="")
session_title_label = tk.Label(root, textvariable=session_title_var, font=("Arial", 14, "bold"), anchor='center')
ToolTip(session_title_label, "Displays the current imaging session name or object being calibrated.")
session_title_label.pack(pady=5)

output_folder_var = tk.StringVar()
max_threads_var = tk.IntVar(value=os.cpu_count())
progress_var = tk.DoubleVar()

master_dark_path = tk.StringVar()
master_flat_path = tk.StringVar()
master_bias_path = tk.StringVar()

master_dark_enabled = tk.BooleanVar()
master_flat_enabled = tk.BooleanVar()
master_bias_enabled = tk.BooleanVar()

# Frame selector section
light_files, dark_files, flat_files, dark_flat_files, bias_files = [], [], [], [], []

from astropy.io import fits

def select_files(file_list, label, expected_type=None):
    validated_files = []  # Always initialize validated_files
    root.update()
    try:
        type_titles = {
            "LIGHT": "Select Light Frames",
            "DARK": "Select Dark Frames",
            "FLAT": "Select Flat Frames",
            "BIAS": "Select Bias Frames",
            "DARKFLAT": "Select Dark Flat Frames"
        }
        title = type_titles.get(expected_type.upper(), "Select Files") if expected_type else "Select Files"
        files = filedialog.askopenfilenames(title=title, filetypes=[("FITS files", "*.fits")])
        root.config(cursor="watch")
        root.update()
        if files:
            validated_files = []
            for f in files:
                try:
                    with fits.open(f) as hdul:
                        header_type = hdul[0].header.get('IMAGETYP', '').upper()
                        if expected_type:
                            if expected_type == "DARKFLAT":
                                if not ("DARK FLAT" in header_type or "DARK" in header_type):
                                    continue
                            elif expected_type not in header_type:
                                continue
                    validated_files.append(f)
                except Exception as e:
                    messagebox.showerror(
                        "FITS Read Error",
                        f"""Could not open FITS file:

{os.path.basename(f)}

Error: {e}""")

            if not validated_files:
                messagebox.showwarning("Frame Validation", f"No valid {expected_type} frames found!")

            file_list.clear()
            file_list.extend(validated_files)
            label.config(text=f"{len(validated_files)} {expected_type.title()} Selected" if expected_type else f"{len(validated_files)} files selected")

            if file_list is dark_files:
                master_dark_enabled.set(False)
            elif file_list is flat_files:
                master_flat_enabled.set(False)
        master_bias_enabled.set(False)

        if validated_files:
                type_labels = {
                    "LIGHT": "Lights",
                    "DARK": "Darks",
                    "FLAT": "Flats",
                    "BIAS": "Bias Frames",
                    "DARKFLAT": "Dark Flats"
                }
                display_type = type_labels.get(expected_type.upper(), expected_type.title() if expected_type else "")
                messagebox.showinfo(
                    "Selection Complete",
                    f"{len(validated_files)} {display_type} selected successfully."
                )
    finally:
        root.config(cursor="")
        root.update()        
# toggle_input_state()  # Removed early call

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
check_frame = tk.Frame(master_frame_container)
check_frame.pack()
master_frame = tk.Frame(master_frame_container)
master_frame.pack()
master_dark_check = tk.Checkbutton(check_frame, text="Use Master Dark", variable=master_dark_enabled, command=lambda: toggle_input_state())
ToolTip(master_dark_check, "Enable using a pre-generated Master Dark frame instead of individual dark frames.")
master_dark_check.pack(side='left', padx=10)
master_flat_check = tk.Checkbutton(check_frame, text="Use Master Flat", variable=master_flat_enabled, command=lambda: toggle_input_state())
ToolTip(master_flat_check, "Enable using a Master Flat frame to correct optical artifacts instead of individual flat frames.")
master_flat_check.pack(side='left', padx=10)
master_bias_check = tk.Checkbutton(check_frame, text="Use Master Bias", variable=master_bias_enabled, command=lambda: toggle_input_state())
ToolTip(master_bias_check, "Enable using a Master Bias frame to correct readout noise instead of individual bias frames.")
master_bias_check.pack(side='left', padx=10)

master_dark_btn = tk.Button(master_frame, text="Select Master Dark", command=lambda: browse_file(master_dark_path))
ToolTip(master_dark_btn, "Choose a pre-generated Master Dark frame to subtract thermal noise from your lights.")
master_dark_label = tk.Label(master_frame, textvariable=master_dark_path, wraplength=300)
master_flat_btn = tk.Button(master_frame, text="Select Master Flat", command=lambda: browse_file(master_flat_path))
ToolTip(master_flat_btn, "Choose a Master Flat frame to correct uneven illumination (vignetting/dust).")
master_flat_label = tk.Label(master_frame, textvariable=master_flat_path, wraplength=300)
master_bias_btn = tk.Button(master_frame, text="Select Master Bias", command=lambda: browse_file(master_bias_path))
ToolTip(master_bias_btn, "Choose a Master Bias frame to correct for electronic readout noise.")
master_bias_label = tk.Label(master_frame, textvariable=master_bias_path, wraplength=300)

def toggle_input_state():
    light_btn.config(state=tk.NORMAL)

    if master_dark_enabled.get() and master_dark_path.get():
        dark_btn.config(state=tk.DISABLED)
        if dark_files:
            dark_files.clear()
            dark_label.config(text="No files")
    else:
        dark_btn.config(state=tk.NORMAL)

    if master_flat_enabled.get() and master_flat_path.get():
        flat_btn.config(state=tk.DISABLED)
        if flat_files:
            flat_files.clear()
            flat_label.config(text="No files")
    else:
        flat_btn.config(state=tk.NORMAL)

    darkflat_btn.config(state=tk.NORMAL)

    if master_bias_enabled.get() and master_bias_path.get():
        bias_btn.config(state=tk.DISABLED)
        if bias_files:
            bias_files.clear()
            bias_label.config(text="No files")
    else:
        bias_btn.config(state=tk.NORMAL)

def update_master_inputs():
    light_label.config(text="No files")
    dark_label.config(text="No files")
    flat_label.config(text="No files")
    darkflat_label.config(text="No files")
    bias_label.config(text="No files")

    master_dark_path.set("")
    master_flat_path.set("")
    master_bias_path.set("")

    master_dark_label.config(text="No files")
    master_flat_label.config(text="No files")
    master_bias_label.config(text="No files")

    light_files.clear()
    dark_files.clear()
    flat_files.clear()
    dark_flat_files.clear()
    bias_files.clear()

    master_dark_enabled.set(False)
    master_flat_enabled.set(False)
    master_bias_enabled.set(False)

    session_title_var.set("")

    toggle_input_state()

    log_textbox.delete('1.0', tk.END)

    messagebox.showinfo(
        "Reset Complete",
        "✅ All settings, selections, and logs have been reset."
    )

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
        # Update the corresponding label
        if var == master_dark_path:
            dark_label.config(text="Master Dark selected")
        elif var == master_flat_path:
            flat_label.config(text="Master Flat selected")
        elif var == master_bias_path:
            bias_label.config(text="Master Bias selected")

# Now call UI initialization after everything is defined

toggle_input_state()
