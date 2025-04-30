from astropy.io import fits
import numpy as np
import os
import zipfile
import concurrent.futures
from datetime import datetime
import gc
from wcs_utils import inject_wcs_from_sidecar
import shutil
from drizzle import drizzle_stack  # New CFA drizzle integration
from scipy.stats import skew
from collections import defaultdict

AUTO_DARK_SCALE = True
ENABLE_CFA_DRIZZLE = True

MASTER_NAMES = {
    'bias': 'master_bias',
    'dark': 'master_dark',
    'flat': 'master_flat',
    'dark_flat': 'master_dark_flat'
}

def inject_minimal_sip(header):
    try:
        header['A_ORDER'] = 0
        header['B_ORDER'] = 0
        header['AP_ORDER'] = 0
        header['BP_ORDER'] = 0

        for ctype_key in ['CTYPE1', 'CTYPE2']:
            if ctype_key in header and '-SIP' not in header[ctype_key]:
                header[ctype_key] = header[ctype_key].strip() + '-SIP'

        print("✅ Minimal SIP keywords and CTYPE updated for compatibility.")
    except Exception as e:
        print(f"⚠️ Failed to inject SIP headers: {e}")

def load_fits_data(path):
    try:
        with fits.open(path, memmap=False) as hdul:
            return hdul[0].data.astype(np.float32)
    except Exception as e:
        print(f"⚠️ Failed to load FITS data from {path}: {e}")
        return None

def estimate_dark_scale(light_data, dark_data):
    if not AUTO_DARK_SCALE:
        return 1.015

    best_scale = 1.0
    best_skew = float('inf')

    mask = light_data < 60000
    if not np.any(mask):
        return 1.015

    for scale in np.linspace(0.9, 1.2, 61):
        residual = (light_data - dark_data * scale)[mask]
        s = abs(skew(residual.flatten()))
        if s < best_skew:
            best_skew = s
            best_scale = scale
    return best_scale

def load_fits_grouped_by_filter(file_list, filter_sensitive=True):
    grouped = defaultdict(list)
    for path in file_list:
        try:
            hdr = fits.getheader(path)
            key = hdr.get('FILTER', 'UNKNOWN').strip().upper() if filter_sensitive else 'ALL'
            grouped[key].append(path)
        except Exception:
            grouped['UNKNOWN'].append(path)
    return grouped

def save_master_per_filter(groups, output_folder, base_name, dark_flat_paths=None, filter_sensitive=True):
    path_map = {}
    written = set()
    for filt, paths in groups.items():
        if not filter_sensitive and 'ALL' in written:
            continue

        stack = []
        for p in paths:
            data = fits.getdata(p).astype(np.float32)
            if dark_flat_paths and filt in dark_flat_paths:
                df = fits.getdata(dark_flat_paths[filt]).astype(np.float32)
                data -= df
            stack.append(data)

        if not stack:
            continue

        combined = np.median(stack, axis=0).astype(np.float32)
        header = fits.getheader(paths[0])
        suffix = f"_{filt}" if filter_sensitive else ""
        filename = os.path.join(output_folder, f"{base_name}{suffix}_master.fits")
        fits.writeto(filename, combined, header, overwrite=True)
        path_map[filt] = filename if filter_sensitive else filename
        written.add('ALL' if not filter_sensitive else filt)
    return path_map

def load_filter_master(path_map, filter_name):
    return path_map.get(filter_name) or path_map.get('UNKNOWN') or path_map.get('ALL')

def calibrate_image(light_path, use_master=False, master_dark_paths=None, master_flat_paths=None, master_bias_paths=None):
    with fits.open(light_path, memmap=False) as hdul:
        light_data = hdul[0].data.astype(float)
        light_header = hdul[0].header.copy()

    exposure_light = light_header.get('EXPTIME', None)
    filter_name = light_header.get('FILTER', 'UNKNOWN').strip().upper()

    master_dark = None

    if use_master and master_bias_paths:
        bias_path = load_filter_master(master_bias_paths, filter_name)
        if bias_path and os.path.exists(bias_path):
            with fits.open(bias_path, memmap=False) as bias_hdul:
                master_bias = bias_hdul[0].data.astype(float)
            light_data -= master_bias

    if use_master and master_dark_paths:
        dark_path = load_filter_master(master_dark_paths, filter_name)
        if dark_path and os.path.exists(dark_path):
            with fits.open(dark_path, memmap=False) as dark_hdul:
                master_dark = dark_hdul[0].data.astype(float)
                exposure_dark = dark_hdul[0].header.get('EXPTIME', None)

            scale_factor = estimate_dark_scale(light_data, master_dark)
            if exposure_light and exposure_dark and exposure_dark != 0:
                scale_factor *= exposure_light / exposure_dark

            light_data -= master_dark * scale_factor
            light_header['DARKSCL'] = (round(scale_factor, 4), 'Dark scale factor used')

    if use_master and master_flat_paths:
        flat_path = load_filter_master(master_flat_paths, filter_name)
        if flat_path and os.path.exists(flat_path):
            with fits.open(flat_path, memmap=False) as flat_hdul:
                master_flat = flat_hdul[0].data.astype(float)
            normalized_flat = master_flat / np.median(master_flat)
            light_data /= normalized_flat
            light_header['CALFLAT'] = (os.path.basename(flat_path), 'Flat field used')

    return light_data, light_header

def create_master_frame(image_list, method='median', dark_flat_path=None):
    if not image_list:
        return None

    stack = []
    for path in image_list:
        data = load_fits_data(path)
        if data is not None:
            if dark_flat_path and os.path.exists(dark_flat_path):
                dark_flat = load_fits_data(dark_flat_path)
                if dark_flat is not None:
                    data -= dark_flat
            stack.append(data)

    if not stack:
        return None

    if method == 'median':
        return np.median(stack, axis=0).astype(np.float32)
    else:
        return np.mean(stack, axis=0).astype(np.float32)

def save_master_frame(data, header, output_folder, name):
    output_path = os.path.join(output_folder, f"{name}_master.fits")
    fits.writeto(output_path, data, header=header, overwrite=True)
    return output_path

def load_fits_by_filter(file_list):
    from collections import defaultdict
    filtered = defaultdict(list)
    for path in file_list:
        try:
            hdr = fits.getheader(path)
            filter_name = hdr.get('FILTER', 'UNKNOWN')
            filtered[filter_name].append(path)
        except Exception:
            filtered['UNKNOWN'].append(path)
    return filtered

def calibrate_and_save(light_path, master_dark_paths, master_flat_paths, master_bias_paths, calibrated_folder):
    calibrated_data, header = calibrate_image(
        light_path,
        use_master=True,
        master_dark_paths=master_dark_paths,
        master_flat_paths=master_flat_paths,
        master_bias_paths=master_bias_paths
    )

    header['CALIB'] = (True, 'Frame has been calibrated')
    header.add_history('Calibrated using Astrocalibrator')

    important_fields = ['OBJECT', 'INSTRUME', 'FILTER', 'DATE-OBS', 'EXPTIME', 'TELESCOP']
    preserved = []
    missing = []

    for key in important_fields:
        if key in header:
            preserved.append(f"{key}='{header[key]}'")
        else:
            missing.append(key)

    likely_osc = ('BAYERPAT' in header) or ('BAYER' in str(header)) or ('FILTER' not in header)

    log_lines = [f"🧾 Preserved header fields for {os.path.basename(light_path)}:"]
    for item in preserved:
        log_lines.append(f"   - {item}")

    if missing:
        suppressed = False
        if likely_osc and 'FILTER' in missing:
            missing = [field for field in missing if field != 'FILTER']
            suppressed = True
        if missing:
            log_lines.append(f"⚠️ Missing expected fields: {', '.join(missing)}")
        elif suppressed:
            log_lines.append(f"ℹ️ Missing 'FILTER' suppressed for likely OSC image.")

    full_log_message = "\n".join(log_lines)

    base_name = os.path.basename(light_path)
    output_path = os.path.join(
        calibrated_folder,
        f"{os.path.splitext(base_name)[0]}_cal.fits"
    )

    fits.writeto(output_path, calibrated_data, header=header, overwrite=True)

    # 🛰️ Copy .wcs sidecar if available
    original_sidecar = os.path.splitext(light_path)[0] + '.wcs'
    calibrated_sidecar = os.path.splitext(output_path)[0] + '.wcs'

    if os.path.exists(original_sidecar):
        shutil.copyfile(original_sidecar, calibrated_sidecar)
        print(f"📄 Copied WCS sidecar: {os.path.basename(calibrated_sidecar)}")
    else:
        print(f"⚠️ No sidecar found for {os.path.basename(light_path)}")

    # 🛰️ Inject WCS from copied sidecar into the calibrated FITS
    from main import inject_wcs_from_sidecar
    success = inject_wcs_from_sidecar(output_path)

    # ✨ If WCS injection succeeded, also inject minimal fake SIP
    if success:
        try:
            with fits.open(output_path, mode='update') as hdul:
                header = hdul[0].header
                inject_minimal_sip(header)
                hdul.flush()
        except Exception as e:
            print(f"⚠️ Failed to reopen calibrated FITS for SIP injection: {e}")

    # 🧹 Free memory
    del calibrated_data
    del header
    gc.collect()

def build_and_save_master(image_list, output_folder, name, dark_flat_path=None):
    if not image_list:
        return None

    header = fits.getheader(image_list[0])
    master = create_master_frame(image_list, dark_flat_path=dark_flat_path)
    if master is not None:
        save_master_frame(master, header, output_folder, name)

    # 🧹 Free memory
    del master
    del header
    gc.collect()
    return None

def run_parallel_calibration(
    light_by_filter,
    dark_by_filter,
    flat_by_filter,
    bias_by_filter,
    output_folder,
    session_title="UnknownObject",
    log_callback=None,
    save_masters=False
):
    if log_callback is None:
        log_callback = print

    log_callback("🛠️ Starting master calibration frame creation...")

    calibrated_folder = os.path.join(output_folder, "calibrated")
    os.makedirs(calibrated_folder, exist_ok=True)

    master_dark_paths = {}
    master_flat_paths = {}
    master_bias_paths = {}

    for filter_name in light_by_filter:
        dark_paths = dark_by_filter.get(filter_name, [])
        flat_paths = flat_by_filter.get(filter_name, [])
        bias_paths = bias_by_filter.get(filter_name, [])

        log_callback(f"🔧 Calibrating filter: {filter_name} ({len(light_by_filter[filter_name])} frames)")

        if dark_paths:
            dark = create_master_frame(dark_paths)
            if dark is not None:
                path = save_master_frame(dark, fits.getheader(dark_paths[0]), output_folder, f"{filter_name}_master_dark")
                master_dark_paths[filter_name] = path

        if flat_paths:
            flat = create_master_frame(flat_paths)
            if flat is not None:
                path = save_master_frame(flat, fits.getheader(flat_paths[0]), output_folder, f"{filter_name}_master_flat")
                master_flat_paths[filter_name] = path

        if bias_paths:
            bias = create_master_frame(bias_paths)
            if bias is not None:
                path = save_master_frame(bias, fits.getheader(bias_paths[0]), output_folder, f"{filter_name}_master_bias")
                master_bias_paths[filter_name] = path

    # Now dispatch calibration threads using file paths
    futures = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for filter_name, light_paths in light_by_filter.items():
            for light_path in light_paths:
                futures.append(
                    executor.submit(
                        calibrate_and_save,
                        light_path,
                        master_dark_paths,
                        master_flat_paths,
                        master_bias_paths,
                        calibrated_folder
                    )
                )

        for i, future in enumerate(concurrent.futures.as_completed(futures)):
            try:
                result = future.result()
                if result:
                    log_callback(result)
                if (i + 1) % 5 == 0:
                    log_callback(f"🖼️ Calibrated {i + 1}/{len(futures)} frames...")
            except Exception as e:
                log_callback(f"💥 Calibration failed: {e}")

    log_callback("✅ Per-filter calibration complete.")

    # 🧹 ZIP master calibration frames if Save Masters is enabled
    if save_masters:
        log_callback("📦 Creating ZIP archive of master calibration frames...")

        session_date = datetime.now().strftime("%Y-%m-%d")
        object_safe = session_title.replace(' ', '_').replace(':', '_').replace('/', '_') or 'UnknownObject'
        zip_name = f"{object_safe}_{session_date}_masters.zip"
        zip_path = os.path.join(output_folder, zip_name)

        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for key, filename in MASTER_NAMES.items():
                path = os.path.join(output_folder, f"{filename}_master.fits")
                if os.path.exists(path):
                    zipf.write(path, arcname=f"{filename}_master.fits")

        log_callback(f"📦 Master calibration frames zipped successfully: {zip_name}")
    else:
        log_callback("ℹ️ Save Masters not enabled. Skipping master frames ZIP creation.")

    return master_dark_paths, master_flat_paths, master_bias_paths