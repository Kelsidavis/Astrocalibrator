import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import os

def log_message(msg):
    print(msg)
    log_textbox.after(0, lambda: (
    log_textbox.insert('end', msg + ''),
    log_textbox.see('end')
    ))

root = tk.Tk()
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
