import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import os
import webbrowser

def log_message(msg):
    print(msg)
    log_textbox.after(0, lambda: (
        log_textbox.insert('end', msg + '\n'),
        log_textbox.see('end')
    ))

root = tk.Tk()

# Add top menu bar
menubar = tk.Menu(root)

# About menu
from urllib.request import urlopen
from PIL import Image, ImageTk
import io

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
    about_window.geometry(f"360x240+{root_x}+{root_y}")  # 20% larger than default 300x200
    about_window.bind('<Button-1>', start_move)
    about_window.bind('<B1-Motion>', do_move)
    about_window.title("About Astrocalibrator")
    about_window.resizable(False, False)

    try:
        with urlopen("https://geekastro.dev/Image%20Gallery_files/title.jpg") as u:
            raw_data = u.read()
        image = Image.open(io.BytesIO(raw_data))
        image = image.resize((100, 100))
        photo = ImageTk.PhotoImage(image)
        img_label = tk.Label(about_window, image=photo)
        img_label.image = photo  # Keep reference
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

session_title_var = tk.StringVar(value="Imaging Session: Unknown")
session_title_label = tk.Label(root, textvariable=session_title_var, font=("Arial", 14, "bold"), anchor='center')
session_title_label.pack(pady=5)


output_folder_var = tk.StringVar()
max_threads_var = tk.IntVar(value=os.cpu_count())
progress_var = tk.DoubleVar()

log_frame = tk.Frame(root)
log_frame.pack(side='bottom', fill='both', expand=True, padx=10, pady=(0, 10))

log_textbox = tk.Text(log_frame, wrap='word', height=10)
log_textbox.pack(side='left', fill='both', expand=True)

scrollbar = ttk.Scrollbar(log_frame, command=log_textbox.yview)
scrollbar.pack(side='right', fill='y')
log_textbox.config(yscrollcommand=scrollbar.set)

def select_output_folder():
    folder = filedialog.askdirectory()
    if folder:
        output_folder_var.set(folder)
        log_message(f"Output directory set: {folder}")
