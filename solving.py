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

def plate_solve_and_update_header(fits_path, log_message):
    try:
        print("👣 Entered plate_solve_and_update_header()")

        cmd = [r"C:\Program Files\astap\astap.exe" if os.path.exists(r"C:\Program Files\astap\astap.exe") else "astap.exe", "-f", fits_path, "-wcs"]
        try:
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=90)
        except FileNotFoundError as e:
            raise FileNotFoundError("Plate solver executable not found") from e
        print(f"[INFO] ASTAP stdout: {result.stdout.strip()}")

        with fits.open(fits_path) as hdul:
            hdr = hdul[0].header
            crval1 = hdr.get('CRVAL1')
            crval2 = hdr.get('CRVAL2')
            if not crval1 or not crval2:
                print("No WCS headers found after solve.")
                return

        coord = SkyCoord(ra=crval1, dec=crval2, unit='deg')
        print(f"Querying Simbad: RA={coord.ra.deg:.6f}, DEC={coord.dec.deg:.6f}")

        def safe_simbad_query(coord):
            try:
                return Simbad.query_region(coord, radius='0d1m0s')
            except Exception as e:
                print(f"Simbad query failed: {e}")
                return None

        with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
            future = executor.submit(safe_simbad_query, coord)
            try:
                result_table = future.result(timeout=10)
            except concurrent.futures.TimeoutError:
                print("⏰ Simbad query timed out.")
                result_table = None

            simbad_id = None
        if result_table and len(result_table) > 0:
            simbad_id = str(result_table[0]['MAIN_ID'])
            print(f"Identified object (Simbad): {simbad_id}")

        print("Checking Messier catalog...")
        
        best_match, min_sep = None, float('inf')
        for name, (ra, dec) in MESSIER_CATALOG.items():
            sep = SkyCoord(ra=ra*u.deg, dec=dec*u.deg).separation(coord).degree
            if sep < 2.0 and sep < min_sep:
                print(f"🔍 Close match: {name} at {sep:.3f}°")
                best_match = name
                min_sep = sep

        if best_match:
            print(f"🗂️ Fallback match: {best_match}")
            session_name = best_match
        else:
            print("No Messier fallback match found.")
            session_name = simbad_id if simbad_id else "Unknown"

        print(f"📅 Session: Imaging Session: {session_name}")
        return session_name
    except Exception as e:
        import traceback
        print(f"💥 Fatal crash in plate_solve_and_update_header: {e} {traceback.format_exc()}")
        raise
