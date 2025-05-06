import os
import subprocess
import shutil
import concurrent.futures
from astropy.io import fits
from astroquery.simbad import Simbad
from astropy import units as u
from object_info import object_info
from astropy.coordinates import SkyCoord, Angle


def cleanup_wcs_file(wcs_path):
    try:
        if os.path.exists(wcs_path):
            os.remove(wcs_path)
            print(f"☑ Deleted temporary WCS file: {wcs_path}")
    except Exception as e:
        print(f"⚠️ Failed to delete WCS file: {e}")


def query_object_name(ra_deg, dec_deg, log_message):
    print(f"DEBUG RAW: ra={ra_deg} ({type(ra_deg)}), dec={dec_deg} ({type(dec_deg)})")
    try:
        ra = float(ra_deg)
        dec = float(dec_deg)
        print(f"🔍 Validating declination before SkyCoord: {dec}")
        if not (-90.0 <= dec <= 90.0):
            raise ValueError(f"Invalid declination: {dec}")
    except Exception as e:
        import traceback
        print("💥 Coordinate validation failed:\n" + traceback.format_exc())
        log_message(f"⚠️ WCS matching failed: {e}")
        return "Unknown Object"

    try:
        skycoord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')

        closest_name = None
        closest_distance = None
        for name, (_, _, obj_ra, obj_dec) in object_info.items():
            obj_coord = SkyCoord(ra=obj_ra * u.deg, dec=obj_dec * u.deg, frame='icrs')
            sep = skycoord.separation(obj_coord)
            print(f"📏 Checking against {name}: sep = {sep.to_string('arcmin')}")
            if sep < Angle(10, unit=u.arcminute):
                if closest_distance is None or sep < closest_distance:
                    closest_name = name
                    closest_distance = sep

        if closest_name:
            print(f"🔭 Matched local object: {closest_name} at {closest_distance.to_string('arcmin')}")
            log_message(f"🔭 Matched local object: {closest_name} ({closest_distance.to_string('arcmin')})")
            return closest_name

        result = Simbad.query_region(skycoord, radius='2m')
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
    """Try reading ASTAP .wcs sidecar as FITS first, fallback to text."""
    def safe_float(val):
        try:
            return float(str(val).strip())
        except Exception:
            return None

    try:
        with fits.open(wcs_path) as hdul:
            hdr = hdul[0].header
            wcs_data = {
                'CRVAL1': safe_float(hdr.get('CRVAL1')),
                'CRVAL2': safe_float(hdr.get('CRVAL2')),
                'CRPIX1': safe_float(hdr.get('CRPIX1')),
                'CRPIX2': safe_float(hdr.get('CRPIX2')),
                'CD1_1': safe_float(hdr.get('CD1_1')),
                'CD1_2': safe_float(hdr.get('CD1_2')),
                'CD2_1': safe_float(hdr.get('CD2_1')),
                'CD2_2': safe_float(hdr.get('CD2_2')),
                'CROTA2': safe_float(hdr.get('CROTA2')),
                'CTYPE1': hdr.get('CTYPE1'),
                'CTYPE2': hdr.get('CTYPE2'),
            }
            return wcs_data
    except Exception as e:
        print(f"⚠️ Failed to read FITS .wcs sidecar: {e}")
        return {}


def plate_solve_and_update_header(fits_path, log_message):
    try:
        print("👣 Entered plate_solve_and_update_header()")
        cmd = [
            r"C:\\Program Files\\astap\\astap.exe" if os.path.exists(r"C:\\Program Files\\astap\\astap.exe") else "astap.exe",
            "-f", fits_path,
            "-wcs", os.path.splitext(fits_path)[0] + ".wcs"
        ]

        print(f"🗂️ Input FITS path received: {fits_path}")
        print(f"🚤 ASTAP executable path: {cmd[0]}")
        print(f"🔧 Running ASTAP with command: {cmd}")

        try:
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=90)
        except FileNotFoundError as e:
            raise FileNotFoundError("Plate solver executable not found") from e

        if result.returncode != 0:
            print(f"❌ ASTAP solver failed with exit code {result.returncode}")
            print(f"❌ ASTAP stderr: {result.stderr.strip()}")
            return None

        print(f"[INFO] ASTAP stdout: {result.stdout.strip()}")

        wcs_file = os.path.splitext(fits_path)[0] + '.wcs'
        if not os.path.exists(wcs_file):
            print("⚠️ No sidecar WCS file found.")
            return None

        wcs_info = parse_wcs_sidecar(wcs_file)
        if not wcs_info:
            print("⚠️ Failed to parse WCS info.")
            return None

        ra = wcs_info.get('CRVAL1')
        dec = wcs_info.get('CRVAL2')

        try:
            ra = float(ra)
            dec = float(dec)
        except (TypeError, ValueError):
            print(f"⚠️ Failed to convert RA/Dec to float: RA={ra}, Dec={dec}")
            log_message(f"⚠️ Invalid WCS values. RA={ra}, Dec={dec}")
            cleanup_wcs_file(wcs_file)
            return "UnknownObject"

        print(f"🔬 Attempting SkyCoord with RA={ra}, Dec={dec}")
        print(f"DEBUG RAW: ra={ra} ({type(ra)}), dec={dec} ({type(dec)})")
        print(f"🔍 Validating declination before SkyCoord: {dec}")

        if not (-90.0 <= dec <= 90.0):
            print(f"⚠️ Invalid declination from WCS: {dec}")
            log_message(f"⚠️ Skipping WCS injection due to invalid declination: {dec}")
            cleanup_wcs_file(wcs_file)
            return "UnknownObject"

        imaging_date = None

        with fits.open(fits_path, mode='update') as hdul:
            hdr = hdul[0].header
            naxis1 = hdr.get('NAXIS1', 2048)
            naxis2 = hdr.get('NAXIS2', 2048)

            hdr['CRVAL1'] = ra
            hdr['CRVAL2'] = dec
            hdr['CTYPE1'] = 'RA---TAN'
            hdr['CTYPE2'] = 'DEC--TAN'
            hdr['CUNIT1'] = 'deg'
            hdr['CUNIT2'] = 'deg'
            hdr['EQUINOX'] = 2000.0
            hdr['RADECSYS'] = 'ICRS'
            hdr['CRPIX1'] = naxis1 / 2
            hdr['CRPIX2'] = naxis2 / 2
            hdr['CD1_1'] = wcs_info.get('CD1_1', -0.000277778)
            hdr['CD1_2'] = wcs_info.get('CD1_2', 0.0)
            hdr['CD2_1'] = wcs_info.get('CD2_1', 0.0)
            hdr['CD2_2'] = wcs_info.get('CD2_2', 0.000277778)

            if 'CROTA2' in wcs_info:
                hdr['CROTA2'] = wcs_info['CROTA2']

            imaging_date = hdr.get('DATE-OBS', None)
            if imaging_date:
                print(f"📅 Captured imaging date: {imaging_date}")
            else:
                print("⚠️ DATE-OBS missing in FITS header.")

            print(f"✅ FITS header updated with full WCS fields for Tycho compatibility.")

        session_name = query_object_name(ra, dec, log_message)
        print(f"📅 Session: Imaging Session: {session_name}")

        cleanup_wcs_file(wcs_file)

        return session_name

    except Exception as e:
        import traceback
        print(f"💥 Fatal crash in plate_solve_and_update_header: {e} {traceback.format_exc()}")
        raise
