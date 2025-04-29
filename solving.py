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

def cleanup_wcs_file(wcs_path):
    try:
        if os.path.exists(wcs_path):
            os.remove(wcs_path)
            print(f"🧹 Deleted temporary WCS file: {wcs_path}")
    except Exception as e:
        print(f"⚠️ Failed to delete WCS file: {e}")

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

def parse_wcs_sidecar(wcs_path):
    """Extract WCS information (CRVAL, CD matrix) from ASTAP .wcs file."""
    wcs_data = {}
    try:
        with open(wcs_path, 'r') as f:
            for line in f:
                if 'RA center' in line:
                    wcs_data['CRVAL1'] = float(line.split(':')[1].strip())
                elif 'DEC center' in line:
                    wcs_data['CRVAL2'] = float(line.split(':')[1].strip())
                elif 'Pixel scale' in line:
                    pixel_scale_arcsec = float(line.split(':')[1].split()[0])
                    wcs_data['CDELT1'] = -pixel_scale_arcsec / 3600.0  # Negative for image orientation
                    wcs_data['CDELT2'] = pixel_scale_arcsec / 3600.0
                elif 'Field rotation' in line:
                    rotation = float(line.split(':')[1].strip())
                    wcs_data['CROTA2'] = rotation
    except Exception as e:
        print(f"⚠️ Failed to parse .wcs file: {e}")
    return wcs_data

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

        if result.returncode != 0:
            print(f"❌ ASTAP solver failed with exit code {result.returncode}")
            print(f"❌ ASTAP stderr: {result.stderr.strip()}")
            return

        print(f"[INFO] ASTAP stdout: {result.stdout.strip()}")

        # Parse the WCS sidecar file
        wcs_file = os.path.splitext(fits_path)[0] + '.wcs'
        if not os.path.exists(wcs_file):
            print("⚠️ No sidecar WCS file found.")
            return

        wcs_info = parse_wcs_sidecar(wcs_file)
        if not wcs_info:
            print("⚠️ Failed to parse WCS info.")
            return

        with fits.open(fits_path, mode='update') as hdul:
            hdr = hdul[0].header
            hdr['CRVAL1'] = wcs_info.get('CRVAL1', 0)
            hdr['CRVAL2'] = wcs_info.get('CRVAL2', 0)
            hdr['CTYPE1'] = 'RA---TAN'
            hdr['CTYPE2'] = 'DEC--TAN'
            hdr['EQUINOX'] = 2000.0
            hdr['RADECSYS'] = 'ICRS'
            hdr['CDELT1'] = wcs_info.get('CDELT1', -0.000277778)
            hdr['CDELT2'] = wcs_info.get('CDELT2', 0.000277778)
            hdr['CROTA2'] = wcs_info.get('CROTA2', 0.0)
            print(f"✅ FITS header updated with WCS information.")

        # Query object name
        session_name = query_object_name(wcs_info['CRVAL1'], wcs_info['CRVAL2'], log_message)
        print(f"📅 Session: Imaging Session: {session_name}")

        cleanup_wcs_file(wcs_file)
        return session_name

    except Exception as e:
        import traceback
        print(f"💥 Fatal crash in plate_solve_and_update_header: {e} {traceback.format_exc()}")
        raise
