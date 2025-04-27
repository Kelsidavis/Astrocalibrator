import os
import subprocess
import shutil
import concurrent.futures
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
from astropy import units as u
from calibration import load_fits_by_filter
from messier_catalog import MESSIER_CATALOG

# Adjustable Settings
PIXEL_SIZE_MICRONS = 3.00   # <-- Default pixel size of camera in microns
FOCAL_LENGTH_MM = 700       # <-- Default telescope focal length in mm

def calculate_pixel_scale():
    """Calculate pixel scale (arcseconds per pixel)."""
    return (206.265 * PIXEL_SIZE_MICRONS) / FOCAL_LENGTH_MM

def query_object_name(ra_deg, dec_deg, log_message):
    coord = SkyCoord(ra=ra_deg*u.degree, dec=dec_deg*u.degree, frame='icrs')
    try:
        result = Simbad.query_region(coord, radius='2m')  # 2 arcminutes search radius
        if result is not None and len(result) > 0:
            main_id = result[0]['MAIN_ID']
            if isinstance(main_id, bytes):
                return main_id.decode('utf-8')
            return main_id
        else:
            return "Unknown Object"
    except Exception as e:
        log_message(f"⚠️ Simbad query failed: {e}")
        return "Unknown Object"

def plate_solve_and_update_header(fits_path, log_message):
    try:
        print("👣 Entered plate_solve_and_update_header()")
        cmd = [
            r"C:\\Program Files\\astap\\astap.exe" if os.path.exists(r"C:\\Program Files\\astap\\astap.exe") else "astap.exe",
            "-f", fits_path,
            "-wcs", os.path.splitext(fits_path)[0] + ".wcs"
        ]

        print(f"🗂️ Input FITS path received: {fits_path}")
        print(f"🛤️ ASTAP executable path: {cmd[0]}")
        print(f"🔧 Running ASTAP with command: {cmd}")

        try:
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=90)
        except FileNotFoundError as e:
            raise FileNotFoundError("Plate solver executable not found") from e

        # Check if ASTAP failed
        if result.returncode != 0:
            print(f"❌ ASTAP solver failed with exit code {result.returncode}")
            print(f"❌ ASTAP stderr: {result.stderr.strip()}")
            return  # Exit early if solver failed

        print(f"[INFO] ASTAP stdout: {result.stdout.strip()}")


        with fits.open(fits_path) as hdul:
            hdr = hdul[0].header
            crval1 = hdr.get('CRVAL1')
            crval2 = hdr.get('CRVAL2')
            if not crval1 or not crval2:
                print("No WCS headers found after solve.")
                wcs_file = os.path.splitext(fits_path)[0] + '.wcs'
                if os.path.exists(wcs_file):
                    print(f"🗂️ Trying to read WCS from {wcs_file}")
                    try:
                        with open(wcs_file, 'r') as f:
                            for line in f:
                                if 'RA center' in line:
                                    crval1 = float(line.split(':')[1].strip())
                                if 'DEC center' in line:
                                    crval2 = float(line.split(':')[1].strip())
                        if crval1 and crval2:
                            print(f"🧭 Sidecar WCS center: RA={crval1:.4f}°, Dec={crval2:.4f}°")
                        else:
                            print("⚠️ Sidecar WCS file missing RA/DEC information.")
                            return
                    except Exception as e:
                        print(f"⚠️ Failed to parse sidecar WCS file: {e}")
                        return
                else:
                    print("⚠️ No sidecar WCS file found.")
                    return

        # Use the newly added query_object_name() here
        session_name = query_object_name(crval1, crval2, log_message)

        print(f"📅 Session: Imaging Session: {session_name}")
        return session_name

    except Exception as e:
        import traceback
        print(f"💥 Fatal crash in plate_solve_and_update_header: {e} {traceback.format_exc()}")
        raise

def open_pixel_settings_window():
    import tkinter as tk
    from tkinter import simpledialog, messagebox

    global PIXEL_SIZE_MICRONS, FOCAL_LENGTH_MM

    root = tk.Tk()
    root.withdraw()  # Hide the main root window

    try:
        pixel_size = simpledialog.askfloat("Pixel Size", "Enter pixel size (microns):", initialvalue=PIXEL_SIZE_MICRONS)
        focal_length = simpledialog.askfloat("Focal Length", "Enter focal length (mm):", initialvalue=FOCAL_LENGTH_MM)

        if pixel_size and focal_length:
            PIXEL_SIZE_MICRONS = pixel_size
            FOCAL_LENGTH_MM = focal_length
            messagebox.showinfo("Settings Updated", f"New pixel scale: {calculate_pixel_scale():.2f} arcsec/pixel")
        else:
            messagebox.showwarning("Settings Cancelled", "No changes made.")

    except Exception as e:
        messagebox.showerror("Error", f"Failed to update settings: {e}")

    root.destroy()
