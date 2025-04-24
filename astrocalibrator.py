import os
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import numpy as np
from astropy.io import fits
from collections import defaultdict
from PIL import Image, ImageTk
import matplotlib.pyplot as plt
import zipfile
import psutil
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import threading
import queue
from datetime import datetime
from multiprocessing import Pool, cpu_count
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.simbad import Simbad

log_text = None  # Will be set later in the script

# === Globals and State ===
MESSIER_CATALOG = {
    "M1": (83.633083, 22.0145),
    "M2": (323.362417, -0.823333),
    "M3": (205.548417, 28.377778),
    "M4": (245.89625, -26.525),
    "M5": (229.638083, 2.082222),
    "M6": (265.075, -32.216667),
    "M7": (266.45, -34.8),
    "M8": (270.925, -24.380278),
    "M9": (259.799167, -18.516111),
    "M10": (254.2875, -4.100556),
    "M11": (282.46, -6.268056),
    "M12": (250.42125, -1.948333),
    "M13": (250.42125, 36.461111),
    "M14": (261.869583, -3.245),
    "M15": (322.49375, 12.166667),
    "M16": (274.7, -13.806667),
    "M17": (275.083333, -16.171944),
    "M18": (275.45, -17.133889),
    "M19": (260.083333, -26.2675),
    "M20": (270.925, -23.013333),
    "M21": (271.1, -22.489722),
    "M22": (279.1, -23.904444),
    "M23": (269.15, -19.016667),
    "M24": (270.45, -18.533333),
    "M25": (276.15, -19.233333),
    "M26": (279.25, -9.383333),
    "M27": (299.90125, 22.721667),
    "M28": (276.1375, -24.869722),
    "M29": (307.4175, 38.533611),
    "M30": (326.170833, -23.179444),
    "M31": (10.684708, 41.26875),
    "M32": (10.674167, 40.865167),
    "M33": (23.462083, 30.659944),
    "M34": (40.5, 42.75),
    "M35": (92.425, 24.35),
    "M36": (83.78, 34.136667),
    "M37": (88.066667, 32.553333),
    "M38": (85.4, 35.85),
    "M39": (326.25, 48.433333),
    "M40": (183.634583, 58.083333),
    "M41": (101.5, -20.75),
    "M42": (83.822083, -5.391111),
    "M43": (83.822083, -5.269444),
    "M44": (130.025, 19.666667),
    "M45": (56.75, 24.116667),
    "M46": (116.083333, -14.8),
    "M47": (118.0, -14.483333),
    "M48": (123.05, -5.85),
    "M49": (187.444583, 8.000556),
    "M50": (107.3, -8.333333),
    "M51": (202.469583, 47.195278),
    "M52": (356.75, 61.6),
    "M53": (198.230833, 18.169722),
    "M54": (283.76375, -30.481389),
    "M55": (294.995833, -30.964167),
    "M56": (289.15, 30.184722),
    "M57": (283.39625, 33.028333),
    "M58": (187.954167, 11.818611),
    "M59": (188.0375, 11.649444),
    "M60": (189.416667, 11.551944),
    "M61": (190.535417, 4.473056),
    "M62": (255.303333, -30.111389),
    "M63": (198.229167, 42.028611),
    "M64": (193.45, 21.683333),
    "M65": (169.545417, 13.0925),
    "M66": (170.0625, 12.991111),
    "M67": (132.825, 11.8),
    "M68": (191.064583, -26.744722),
    "M69": (276.1375, -32.346944),
    "M70": (277.846667, -32.321667),
    "M71": (298.443333, 18.779722),
    "M72": (313.36625, -12.5375),
    "M73": (313.75, -12.625),
    "M74": (24.174167, 15.783056),
    "M75": (300.161667, -21.921389),
    "M76": (9.203333, 51.571111),
    "M77": (40.669583, 0.013056),
    "M78": (86.705417, -1.201944),
    "M79": (85.9875, -24.524722),
    "M80": (244.259167, -22.975833),
    "M81": (148.88875, 69.065278),
    "M82": (148.968333, 69.679722),
    "M83": (204.253333, -29.865833),
    "M84": (186.265, 12.886111),
    "M85": (186.349167, 18.191389),
    "M86": (186.54875, 12.946111),
    "M87": (187.705833, 12.391111),
    "M88": (187.97625, 14.420278),
    "M89": (188.915417, 12.556389),
    "M90": (189.016667, 13.1625),
    "M91": (186.629167, 14.496667),
    "M92": (259.28125, 43.136389),
    "M93": (115.995, -23.858611),
    "M94": (192.721667, 41.119444),
    "M95": (161.892917, 11.703333),
    "M96": (162.347083, 11.819167),
    "M97": (168.698333, 55.019444),
    "M98": (185.531667, 14.291111),
    "M99": (185.729167, 14.423333),
    "M100": (185.72875, 15.822222),
    "M101": (210.802083, 54.348972),
    "M102": (230.174583, 55.410278),
    "M103": (10.049583, 60.716667),
    "M104": (189.997917, -11.623056),
    "M105": (161.956667, 12.581944),
    "M106": (184.739583, 47.304722),
    "M107": (246.45, -13.053056),
    "M108": (165.5175, 55.674167),
    "M109": (173.25625, 53.374444),
    "M110": (10.09125, 41.685278)
}
light_files, dark_files, flat_files, dark_flat_files = [], [], [], []
log_queue = queue.Queue()

def load_fits_by_filter(file_list):
    filtered = defaultdict(list)
    for path in file_list:
        try:
            with fits.open(path) as hdul:
                hdr = hdul[0].header
                filter_name = hdr.get('FILTER', 'UNKNOWN')
                filtered[filter_name].append(path)
        except Exception as e:
            log_message(f"Failed to read {path}: {e}")
    return filtered

def create_master_frame(files, method='median'):
    stack = []
    for path in files:
        try:
            with fits.open(path) as hdul:
                data = hdul[0].data.astype(np.float32)
                stack.append(data)
        except Exception as e:
            log_message(f"Skipping {path}: {e}")
    if not stack:
        return None
    return np.median(stack, axis=0) if method == 'median' else np.mean(stack, axis=0)

def save_master_frame(data, header, output_folder, name):
    output_path = os.path.join(output_folder, f"master_{name}.fits")
    hdu = fits.PrimaryHDU(data=data, header=header)
    hdu.writeto(output_path, overwrite=True)

def calibrate_image(args):
    path, master_dark, master_flat, master_dark_flat, output_folder, filter_name = args
    try:
        with fits.open(path) as hdul:
            data = hdul[0].data.astype(np.float32)
            header = hdul[0].header
            if master_dark is not None:
                data -= master_dark
            if master_flat is not None:
                data /= master_flat
            if master_dark_flat is not None:
                data -= master_dark_flat
            data = np.clip(data, 0, None)
            output_path = os.path.join(output_folder, f"cal_{filter_name}_{os.path.basename(path)}")
            fits.writeto(output_path, data=data, header=header, overwrite=True)
            return output_path
    except Exception as e:
        return f"Error processing {path}: {e}"

def zip_output_folder(folder):
    # Determine object name
    session_text = session_title_var.get()
    object_part = session_text.split(':')[-1].strip().split(' ')[0]

    # Get latest FITS file timestamp from light frames
    try:
        latest_file = max(light_files, key=os.path.getmtime)
        latest_date = datetime.fromtimestamp(os.path.getmtime(latest_file)).strftime('%Y%m%d')
    except Exception:
        latest_date = datetime.now().strftime('%Y%m%d')

    zip_name = f"{object_part}_{latest_date}.zip"
    parent_folder = os.path.abspath(os.path.join(folder, os.pardir))
    zip_path = os.path.join(parent_folder, zip_name)

    with zipfile.ZipFile(zip_path, 'w') as zipf:
            metadata = [
                f"Session Title: {session_title_var.get()}",
                f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
                "Included Files:",
            ]
            for file in files:
                full_path = os.path.join(rootdir, file)
                if file.startswith("stacked_") and file.endswith(".fits"):
                    zipf.write(full_path, arcname=os.path.join("stacked", file))
                    metadata.append(f"  stacked/{file}")
                elif file.startswith("cal_") and file.endswith(".fits"):
                    zipf.write(full_path, arcname=file)
                    metadata.append(f"  {file}")
    readme_path = os.path.join(folder, "README.txt")
    with open(readme_path, 'w') as f:
        f.write(''.join(metadata))
    zipf.write(readme_path, arcname="README.txt")
    os.remove(readme_path)
    log_message(f"Archived results to {zip_path}")
    log_message(f"🎉 [1mCalibration complete! Archive saved to:[0m{zip_path}")
    try:
        root.bell()  # Cheerful notification
    except Exception:
        pass

def run_calibration():
    def _run():
        from time import time
        log_message("🟢 Starting calibration process...")

        log_message("🔍 Loading and grouping light frames by filter...")
        light_images = load_fits_by_filter(light_files)
        log_message("🔍 Loading and grouping dark frames by filter...")
        dark_images = load_fits_by_filter(dark_files)
        log_message("🔍 Loading and grouping flat frames by filter...")
        flat_images = load_fits_by_filter(flat_files)
        log_message("🔍 Loading and grouping dark flat frames by filter...")
        dark_flat_images = load_fits_by_filter(dark_flat_files)

        method = stacking_method.get()
        output_folder = output_folder_var.get()
        save_masters = save_masters_var.get()

        master_darks = {}
        master_flats = {}
        master_dark_flats = {}

        for f in dark_images:
            log_message(f"🌑 Creating master dark for filter: {f} using {len(dark_images[f])} frames")
            master = create_master_frame(dark_images[f], method)
            master_darks[f] = master
            if save_masters:
                save_master_frame(master, fits.getheader(dark_images[f][0]), output_folder, f"dark_{f}")

        for f in flat_images:
            log_message(f"🟦 Creating master flat for filter: {f} using {len(flat_images[f])} frames")
            master = create_master_frame(flat_images[f], method)
            master_flats[f] = master
            if save_masters:
                save_master_frame(master, fits.getheader(flat_images[f][0]), output_folder, f"flat_{f}")

        for f in dark_flat_images:
            log_message(f"🌓 Creating master dark flat for filter: {f} using {len(dark_flat_images[f])} frames")
            master = create_master_frame(dark_flat_images[f], method)
            master_dark_flats[f] = master
            if save_masters:
                save_master_frame(master, fits.getheader(dark_flat_images[f][0]), output_folder, f"darkflat_{f}")

        tasks = []
        for filter_name, files in light_images.items():
            log_message(f"📸 Preparing {len(files)} light frames for filter: {filter_name}")
            for path in files:
                tasks.append((path, master_darks.get(filter_name), master_flats.get(filter_name), master_dark_flats.get(filter_name), output_folder, filter_name))

        start_time = time()
        log_message("⚙️ Beginning parallel calibration using multiprocessing...")
        with Pool(processes=max_threads_var.get()) as pool:
            for i, result in enumerate(pool.imap_unordered(calibrate_image, tasks), 1):
                if result.startswith("Error"):
                    log_message(result)
                else:
                    log_message(f"Saved calibrated: {result}")
                elapsed = time() - start_time
                eta = (elapsed / i) * (len(tasks) - i) if i > 0 else 0
                progress_var.set((i / len(tasks)) * 100)
                log_message(f"Progress: {i}/{len(tasks)} | ETA: {int(eta)}s")
                root.update_idletasks()

        if stack_lights_var.get():
            from astropy.stats import sigma_clip
            for filter_name in set(f[5] for f in tasks):
                stacked = []
                for file in os.listdir(output_folder):
                    if file.startswith(f"cal_{filter_name}_") and file.endswith(".fits"):
                        with fits.open(os.path.join(output_folder, file)) as hdul:
                            stacked.append(hdul[0].data.astype(np.float32))
                if stacked:
                    stack_array = np.array(stacked)
                    if sigma_clip_var.get():
                        clipped = sigma_clip(stack_array, sigma=sigma_threshold.get(), maxiters=sigma_iterations.get(), axis=0)
                        combined = np.mean(clipped, axis=0)
                        log_message(f"📊 Sigma clipping used for stacking {len(stacked)} frames on {filter_name}")
                    else:
                        combined = np.median(stack_array, axis=0) if method == 'median' else np.mean(stack_array, axis=0)
                    output_path = os.path.join(output_folder, f"stacked_{filter_name}.fits")
                    fits.writeto(output_path, combined, overwrite=True)
                    log_message(f"🧱 Stacked {len(stacked)} frames for filter: {filter_name} -> {output_path}")
                    log_message(f"📈 Stack stats for {filter_name}: min={combined.min():.2f}, max={combined.max():.2f}, mean={combined.mean():.2f}, std={combined.std():.2f}")
        zip_output_folder(output_folder)
        # Remove uncompressed calibrated files
        for file in os.listdir(output_folder):
            if file.startswith("cal_") and file.endswith(".fits"):
                try:
                    os.remove(os.path.join(output_folder, file))
                except Exception as e:
                    log_message(f"Failed to delete {file}: {e}")
        log_message("✅ Calibration complete. Results zipped and cleaned up.")

    threading.Thread(target=_run).start()

def log_message(message):
    if log_text:
        log_text.insert(tk.END, message + "\n")
        log_text.see(tk.END)
    else:
        print(message)

def select_output_folder():
    folder = filedialog.askdirectory()
    if folder:
        output_folder_var.set(folder)
        output_display_label.config(text=f"Output Directory: {folder}")

def select_files(file_list, label):
    files = filedialog.askopenfilenames(filetypes=[("FITS files", "*.fits")])
    if files:
        file_list.clear()
        file_list.extend(files)
        label.config(text=f"{len(files)} files selected")
        if file_list is light_files:
            run_calibration_btn.config(state='normal')
            plate_solve_btn.config(state='normal')
            output_folder = os.path.dirname(files[0])
            output_folder_var.set(output_folder)
            output_display_label.config(text=f"Output Directory: {output_folder}")
            try:
                with fits.open(files[0]) as hdul:
                    hdr = hdul[0].header
                    if 'RA' in hdr and 'DEC' in hdr:
                        coord = SkyCoord(ra=hdr['RA'], dec=hdr['DEC'], unit=(u.hourangle, u.deg))
                        log_message(f"Querying Simbad near RA={coord.ra.deg:.6f}, DEC={coord.dec.deg:.6f}")
                        Simbad.add_votable_fields("ids", "coo_const")
                        result_table = Simbad.query_region(coord, radius='0d40m0s')
                        if result_table is not None and len(result_table) > 0:
                            object_name = result_table[0]['MAIN_ID'].decode('utf-8')
                            session_title_var.set(f"Imaging Session: {object_name}")
            except Exception:
                session_title_var.set("Imaging Session: Unknown")

def plate_solve_and_update_header(fits_path):
    try:
        import subprocess
        import shutil
        cmd = ["astap.exe", "-f", fits_path, "-wcs"]
        if shutil.which("astap.exe") is None:
            default_path = r"C:\Program Files\astap\astap.exe"
            if os.path.exists(default_path):
                cmd[0] = default_path
            else:
                log_message("Local solver 'ASTAP' not found. Please install it from https://www.hnsky.org/astap.htm and add it to your PATH.")
                return
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0 or not os.path.exists(fits_path):
            log_message(f"[ERROR] ASTAP failed with code {result.returncode}")
            if result.stderr:
                log_message(f"[STDERR] {result.stderr.strip()}")
            return
        log_message(f"[INFO] ASTAP stdout: {result.stdout.strip()}")
        log_message(f"Solved with ASTAP: {os.path.basename(fits_path)}")

        result_table = None
        object_name = None

        try:
            with fits.open(fits_path) as hdul:
                hdr = hdul[0].header
                if 'CRVAL1' in hdr and 'CRVAL2' in hdr:
                    crval1, crval2 = hdr['CRVAL1'], hdr['CRVAL2']
                    log_message(f"WCS update confirmed in FITS: CRVAL1={crval1}, CRVAL2={crval2}")
                else:
                    # Try loading from .wcs file
                    wcs_file = os.path.splitext(fits_path)[0] + '.wcs'
                    if os.path.exists(wcs_file):
                        try:
                            with fits.open(wcs_file) as hdul:
                                hdr = hdul[0].header
                                crval1 = hdr.get('CRVAL1')
                                crval2 = hdr.get('CRVAL2')
                        except Exception as e:
                            log_message(f"⚠️ Failed to read .wcs file as FITS: {e}")
                            return
                        if crval1 is not None and crval2 is not None:
                            with fits.open(fits_path, mode='update') as hdul:
                                hdul[0].header['CRVAL1'] = crval1
                                hdul[0].header['CRVAL2'] = crval2
                                hdul.flush()
                            log_message(f"Injected WCS from .wcs: CRVAL1={crval1}, CRVAL2={crval2}")
                        else:
                            log_message(f"⚠️ WCS keys missing in .wcs file for {os.path.basename(fits_path)}")
                            return
                    else:
                        log_message(f"⚠️ No WCS in FITS or .wcs file for {os.path.basename(fits_path)}")
                        return

            coord = SkyCoord(ra=crval1, dec=crval2, unit='deg')
            log_message(f"Querying Simbad near RA={coord.ra.deg:.6f}, DEC={coord.dec.deg:.6f}")
            result_table = Simbad.query_region(coord, radius='0d1m0s')

            if result_table is not None and len(result_table) > 0:
                best_name = None
                constellation = "Unknown"
                for row in result_table:
                    main_id = str(row['MAIN_ID'])
                    ids = row['IDS'] if 'IDS' in result_table.colnames else ""
                    log_message(f"Checking candidate: MAIN_ID={main_id} | IDS={ids}")
                    if ids:
                        for alias in ids.split('|'):
                            alias = alias.strip().upper()
                            if alias.startswith("M") and alias[1:].strip().isdigit():
                                best_name = alias
                                break
                            elif alias.startswith("NGC"):
                                best_name = alias
                                break
                    if best_name:
                        constellation = row['COO_CONST'] if 'COO_CONST' in result_table.colnames else "Unknown"
                        break

                if not best_name:
                    # Attempt to identify from hardcoded Messier catalog even if Simbad returned rows
                    best_match = None
                    min_sep = float('inf')
                    target_coord = SkyCoord(ra=crval1, dec=crval2, unit='deg')
                    for name, (ra, dec) in MESSIER_CATALOG.items():
                        catalog_coord = SkyCoord(ra=ra, dec=dec, unit='deg')
                        sep = target_coord.separation(catalog_coord).degree
                        if sep < 1.0 and sep < min_sep:
                            best_match = name
                            min_sep = sep
                    if best_match:
                        object_name = best_match
                        full_title = f"Imaging Session: {object_name} (from catalog)"
                        root.after(0, session_title_var.set, full_title)
                        log_message(f"🗂️ Fallback to catalog: {full_title}")
                        return
                    best_name = str(result_table[0]['MAIN_ID'])
                    constellation = result_table[0]['COO_CONST'] if 'COO_CONST' in result_table.colnames else "Unknown"

                object_name = best_name
                all_ids = result_table[0]['IDS'] if 'IDS' in result_table.colnames else ""
                log_message(f"All aliases: {all_ids}")

                # Calculate total integration time from light frames
                total_integration = 0
                try:
                    for lf in light_files:
                        with fits.open(lf) as hdul:
                            hdr = hdul[0].header
                            total_integration += float(hdr.get('EXPTIME', 0))
                except Exception as e:
                    log_message(f"Warning: Failed to calculate total integration time: {e}")

                full_title = f"Imaging Session: {object_name} in {constellation} — {int(total_integration)}s"
                root.after(0, session_title_var.set, full_title)
                log_message(f"Updated session title to: {full_title}")
                log_message("About to set session title using root.after")
                root.after(0, session_title_var.set, f"Imaging Session: {object_name}")
                log_message("Called root.after to update session title")
                log_message(f"Identified object: {object_name}")
            else:
                # Attempt to identify from hardcoded Messier catalog
                best_match = None
                min_sep = float('inf')
                target_coord = SkyCoord(ra=crval1, dec=crval2, unit='deg')
                for name, (ra, dec) in MESSIER_CATALOG.items():
                    catalog_coord = SkyCoord(ra=ra, dec=dec, unit='deg')
                    sep = target_coord.separation(catalog_coord).degree
                    if sep < 1.0 and sep < min_sep:
                        best_match = name
                        min_sep = sep
                if best_match:
                    object_name = best_match
                    full_title = f"Imaging Session: {object_name} (from catalog)"
                    root.after(0, session_title_var.set, full_title)
                    log_message(f"🗂️ Fallback to catalog: {full_title}")
                else:
                    log_message(f"No object found near RA={coord.ra.deg:.6f}, DEC={coord.dec.deg:.6f}")

        except Exception as wcs_err:
            log_message(f"⚠️ Could not verify WCS or identify object for {os.path.basename(fits_path)}: {wcs_err}")

    except Exception as e:
        log_message(f"Plate solving failed for {os.path.basename(fits_path)}: {e}")
        

def plate_solve_all_light_frames():
    summary_window = tk.Toplevel(root)
    summary_window.title("Plate Solving Log")
    summary_text = tk.Text(summary_window, wrap='word', height=20, width=80, bg='black', fg='#D16FFF')
    summary_text.pack(expand=True, fill='both')
    cancel_flag = threading.Event()

    def log_to_summary(message):
        print(message)
        summary_text.insert(tk.END, message + '\n')
        summary_text.see(tk.END)
        log_message(message)

    def run_plate_solve():
        import time
        total = len(light_files)
        log_to_summary("Starting plate solving of light frames...")
        log_to_summary(f"Total frames to solve: {total}")
        success, failed = [], []
        start = time.time()

        for idx, path in enumerate(light_files, 1):
            if cancel_flag.is_set():
                log_to_summary("❌ Plate solving cancelled by user.")
                break
            log_to_summary(f"[{idx}/{total}] Solving: {os.path.basename(path)}")
            try:
                plate_solve_and_update_header(path)
                # Imaging session title now updated in plate_solve_and_update_header only
                log_to_summary(f"✅ Success: {os.path.basename(path)}")
                success.append(os.path.basename(path))
            except Exception as e:
                error_msg = f"❌ Failed: {os.path.basename(path)} — {str(e)}"
                log_to_summary(error_msg)
                failed.append(error_msg)
            elapsed = time.time() - start
            eta = (elapsed / idx) * (total - idx) if idx > 0 else 0
            log_to_summary(f"⏳ ETA: {int(eta)}s")

        if not cancel_flag.is_set():
            log_to_summary("Plate solving complete.")
            log_message("🟣 All selected light frames have been plate solved.")
            # Wiggle Run Calibration button
            def wiggle(button, count=6):
                if count <= 0:
                    button.place_forget()
                    button.pack(anchor='e', pady=5)
                    return
                offset = -5 if count % 2 == 0 else 5
                x = button.winfo_x() + offset
                button.place(x=x, y=button.winfo_y())
                root.after(50, lambda: wiggle(button, count - 1))

            for child in right_step.winfo_children():
                if isinstance(child, tk.Button) and child['text'] == 'Run Calibration':
                    child.place(x=child.winfo_x(), y=child.winfo_y())  # Temporarily place to allow wiggle
                    wiggle(child)
                    break
            log_to_summary(f"✅ Successful solves: {len(success)}")
            log_to_summary(f"❌ Failed solves: {len(failed)}")
            if failed:
                log_to_summary("Failures:\n" + "\n".join(failed))

    cancel_button = tk.Button(summary_window, text="Cancel", command=cancel_flag.set)
    cancel_button.pack(pady=5)
    threading.Thread(target=run_plate_solve).start()

import json

SETTINGS_FILE = os.path.join(os.path.expanduser('~'), '.astro_calibrator_settings.json')

def load_settings():
    try:
        with open(SETTINGS_FILE, 'r') as f:
            return json.load(f)
    except:
        return {}

def save_settings():
    settings = {
        'sigma_threshold': sigma_threshold.get(),
        'sigma_iterations': sigma_iterations.get(),
        'stack_lights': stack_lights_var.get(),
        'sigma_clip': sigma_clip_var.get()
    }
    try:
        with open(SETTINGS_FILE, 'w') as f:
            json.dump(settings, f)
    except Exception as e:
        log_message(f"⚠️ Failed to save settings: {e}")

root = tk.Tk()
session_title_var = tk.StringVar(value="Imaging Session: Unknown")
root.title("Astro FITS Calibrator")

session_title_entry = tk.Entry(root, textvariable=session_title_var, font=("Arial", 14, "bold"), justify='center')
session_title_entry.pack(pady=5, fill='x', padx=10)

output_folder_var = tk.StringVar()
max_threads_var = tk.IntVar(value=cpu_count())
progress_var = tk.DoubleVar()

# preview_label = tk.Label(root)  # Disabled image preview
log_text = tk.Text(root)

frame = tk.LabelFrame(root, text="Select Calibration Frames", padx=10, pady=10, font=("Arial", 10, "bold"))
frame.pack(padx=10, pady=10)

entry_labels = []
btns = [
    ("Select Light Frames", light_files),
    ("Select Dark Frames", dark_files),
    ("Select Flat Frames", flat_files),
    ("Select Dark Flats", dark_flat_files),
]
for i, (label, store) in enumerate(btns):
    var = tk.Label(frame, text="No files selected")
    btn = tk.Button(frame, text=label, command=lambda s=store, l=var: select_files(s, l))
    btn.tooltip_text = f"Select your {label.lower().replace('select ', '')}. Must be FITS files."
    btn.grid(row=i, column=0, padx=5, pady=5, sticky='w')
    var.grid(row=i, column=1, padx=5, pady=5, sticky='w')
    entry_labels.append(var)

output_display_label = tk.Label(root, text="Output Directory: Will be set after selecting light frames")
output_display_label.pack(pady=5)
tk.Button(root, text="Change Output Folder", command=select_output_folder).pack()

save_masters_var = tk.BooleanVar()
stack_lights_var = tk.BooleanVar()
sigma_clip_var = tk.BooleanVar()
save_cb = tk.Checkbutton(root, text="Save Master Calibration Files", variable=save_masters_var)
save_cb.pack()
save_cb.tooltip_text = "Save the master dark, flat, and dark-flat calibration frames separately."
stack_cb = tk.Checkbutton(root, text="Stack Light Frames into One", variable=stack_lights_var)
stack_cb.pack()
stack_cb.tooltip_text = "Combine all calibrated light frames into one stacked image per filter."
clip_cb = tk.Checkbutton(root, text="Use Sigma Clipping When Stacking", variable=sigma_clip_var)
clip_cb.pack()
clip_cb.tooltip_text = "Apply sigma clipping to remove outliers during stacking."

sigma_threshold = tk.DoubleVar(value=3.0)
settings = load_settings()
sigma_threshold = tk.DoubleVar(value=settings.get('sigma_threshold', 3.0))
sigma_iterations = tk.IntVar(value=settings.get('sigma_iterations', 5))
stack_lights_var.set(settings.get('stack_lights', False))
sigma_clip_var.set(settings.get('sigma_clip', False))

sigma_frame = tk.Frame(root)
sigma_frame.pack(pady=2)
sigma_label = tk.Label(sigma_frame, text="Sigma Threshold:")
sigma_label.grid(row=0, column=0, padx=5)
sigma_label.tooltip_text = "Standard deviation multiplier used to clip outlier pixels."
sigma_entry = tk.Entry(sigma_frame, textvariable=sigma_threshold, width=5)
sigma_entry.grid(row=0, column=1, padx=5)
sigma_entry.tooltip_text = "Enter the sigma threshold (e.g. 3.0) for outlier rejection."
iter_label = tk.Label(sigma_frame, text="Iterations:")
iter_label.grid(row=0, column=2, padx=5)
iter_label.tooltip_text = "Number of times to run sigma clipping on the stack."
iter_entry = tk.Entry(sigma_frame, textvariable=sigma_iterations, width=5)
iter_entry.grid(row=0, column=3, padx=5)
iter_entry.tooltip_text = "Enter how many iterations to run sigma clipping (e.g. 5)."

stacking_method = tk.StringVar(value='median')
tk.Label(root, text="Stacking Method:").pack()
method_box = ttk.Combobox(root, textvariable=stacking_method, values=['median', 'average'], state='readonly')
method_box.pack()
method_box.tooltip_text = "Choose stacking method: median is more robust to outliers."

tk.Label(root, text="Max Threads:").pack()
threads_spin = tk.Spinbox(root, from_=1, to=cpu_count()*2, textvariable=max_threads_var, width=5)
threads_spin.pack()
threads_spin.tooltip_text = "Set number of CPU threads to use for calibration."

step_frame = tk.Frame(root)
step_frame.pack(fill='x', padx=10, pady=10)

left_step = tk.Frame(step_frame)
left_step.pack(side='left', expand=True, padx=20)
tk.Label(left_step, text="Step 1: Plate Solve Light Frames").pack(anchor='w')
plate_solve_btn = tk.Button(left_step, text="Plate Solve", command=plate_solve_all_light_frames, font=('Arial', 12, 'bold'), height=2, width=25, state='disabled')
plate_solve_btn.pack(anchor='w', pady=5)

right_step = tk.Frame(step_frame)
right_step.pack(side='right', expand=True, padx=20)
tk.Label(right_step, text="Step 2: Run Calibration").pack(anchor='e')
run_calibration_btn = tk.Button(right_step, text="Run Calibration", command=run_calibration, font=('Arial', 12, 'bold'), height=2, width=25, state='disabled')
run_calibration_btn.pack(anchor='e', pady=5)

progress_bar = ttk.Progressbar(root, variable=progress_var, maximum=100)
progress_bar.pack(fill=tk.X, padx=10, pady=10)

log_frame = tk.Frame(root)
log_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
log_text = tk.Text(log_frame, wrap='word', height=10)
log_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
scrollbar = ttk.Scrollbar(log_frame, command=log_text.yview)
scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
log_text.config(yscrollcommand=scrollbar.set)

def reset_to_defaults():
    sigma_threshold.set(3.0)
    sigma_iterations.set(5)
    stack_lights_var.set(False)
    sigma_clip_var.set(False)
    save_settings()
    log_message("🔄 Settings reset to defaults.")

menu_bar = tk.Menu(root)
file_menu = tk.Menu(menu_bar, tearoff=0)
file_menu.add_command(label="View FITS Header", command=lambda: view_header(filedialog.askopenfilename(filetypes=[("FITS files", "*.fits")])))
menu_bar.add_cascade(label="Tools", menu=file_menu)
about_menu = tk.Menu(menu_bar, tearoff=0)
about_menu.add_command(label="Visit geekastro.dev", command=lambda: os.system('start https://geekastro.dev'))
menu_bar.add_cascade(label="About", menu=about_menu)

reset_menu = tk.Menu(menu_bar, tearoff=0)
reset_menu.add_command(label="Reset Settings to Defaults", command=reset_to_defaults)
menu_bar.add_cascade(label="Reset", menu=reset_menu)
root.config(menu=menu_bar)

save_settings()

# Apply tooltips
try:
    import idlelib.tooltip as tooltip
    def attach_tooltip(widget):
        if hasattr(widget, 'tooltip_text'):
            tooltip.Hovertip(widget, widget.tooltip_text, hover_delay=500)
    for widget in root.winfo_children():
        for child in widget.winfo_children():
            attach_tooltip(child)
            for grandchild in child.winfo_children():
                attach_tooltip(grandchild)
            attach_tooltip(child)
        attach_tooltip(widget)
except Exception:
    log_message("Tooltips not available - idlelib module missing")

root.mainloop()
