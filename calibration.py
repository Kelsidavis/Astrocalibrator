from astropy.io import fits
import numpy as np
import os
import zipfile
import concurrent.futures
from datetime import datetime
import gc
from wcs_utils import inject_wcs_from_sidecar
import shutil
from scipy.stats import skew
from collections import defaultdict
import unicodedata
import re

AUTO_DARK_SCALE = True

MASTER_NAMES = {
    'bias': 'master_bias',
    'dark': 'master_dark',
    'flat': 'master_flat',
    'dark_flat': 'master_dark_flat'
}

def group_darks_by_exposure(dark_files):
    grouped = defaultdict(list)
    for path in dark_files:
        try:
            hdr = fits.getheader(path)
            exptime = float(hdr.get('EXPTIME', -1))
            if exptime > 0:
                grouped[exptime].append(path)
            else:
                log_message(f"⚠️ Skipping dark {path} — invalid EXPTIME")
        except Exception as e:
            log_message(f"⚠️ Failed to group dark {path}: {e}")
    return dict(grouped)

def group_flats_by_filter_and_exposure(flat_files):
    grouped = defaultdict(list)
    for path in flat_files:
        try:
            hdr = fits.getheader(path)
            filt = hdr.get('FILTER', 'UNKNOWN').strip().upper()
            exptime = float(hdr.get('EXPTIME', -1))
            if exptime > 0:
                grouped[(filt, exptime)].append(path)
        except Exception as e:
            print(f"⚠️ Failed to group flat: {path} – {e}")
    return dict(grouped)

def group_dark_flats_by_filter_and_exptime(dark_flat_files):
    grouped = defaultdict(list)
    for path in dark_flat_files:
        try:
            hdr = fits.getheader(path)
            filt = normalize_filter_name(hdr.get("FILTER", "UNKNOWN"))
            exptime_raw = hdr.get("EXPTIME", -1)
            try:
                exptime = round(float(exptime_raw), 1)
            except (TypeError, ValueError):
                print(f"⚠️ Invalid EXPTIME value: {exptime_raw} in {path}")
                continue
            if exptime > 0:
                grouped[(filt, exptime)].append(path)
        except Exception as e:
            print(f"⚠️ Failed to group dark flat: {path} – {e}")
    return dict(grouped)

def normalize_filter_name(raw):
    if not raw:
        return 'UNKNOWN'

    raw = str(raw)  # 🔒 Ensure input is a string

    # Normalize Unicode: strip accents, convert full-width & compatibility characters
    raw = unicodedata.normalize('NFKD', raw)
    raw = ''.join(c for c in raw if not unicodedata.combining(c))

    # Lowercase and strip non-alphanumerics
    raw = raw.lower().strip()
    raw = re.sub(r'[^a-z0-9]', '', raw)

    aliases = {
        'ha': ['ha', 'halpha', 'hα', 'hαlpha', 'h-alpha'],
        'oiii': ['oiii', 'o3'],
        'sii': ['sii', 's2'],
        'l': ['l', 'lum', 'luminance'],
        'r': ['r', 'red'],
        'g': ['g', 'green'],
        'b': ['b', 'blue']
    }

    for norm, variants in aliases.items():
        if raw in variants:
            return norm.upper()

    return raw.upper()

def create_master_flat_scaled(flat_paths, dark_flat_path=None):
    """Create a master flat by subtracting matching dark flat and normalizing."""
    stack = []

    if dark_flat_path and os.path.exists(dark_flat_path):
        dark_flat = load_fits_data(dark_flat_path)
    else:
        dark_flat = None

    for path in flat_paths:
        try:
            data = fits.getdata(path).astype(np.float32)
            if dark_flat is not None:
                data -= dark_flat
            norm = np.median(data)
            if norm == 0:
                continue
            stack.append(data / norm)
        except Exception as e:
            print(f"⚠️ Failed to normalize flat: {path} – {e}")

    if not stack:
        return None

    combined = np.mean(stack, axis=0).astype(np.float32)
    return combined

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

def load_filter_master(path_map, filter_name, log_callback=print):
    filter_name = normalize_filter_name(filter_name)
    path = path_map.get(filter_name)
    if not path:
        log_callback(f"⚠️ No master for filter '{filter_name}', trying fallback...")
        path = path_map.get('UNKNOWN') or path_map.get('ALL')
    if not path:
        log_callback(f"❌ No suitable master file found for filter '{filter_name}'")
    return path


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
            print(f"🧪 Applied master bias for filter '{filter_name}'")

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
            print(f"🧪 Applied master dark for filter '{filter_name}' (scale={round(scale_factor, 4)})")

    if use_master and master_flat_paths:
        flat_path = load_filter_master(master_flat_paths, filter_name)
        if flat_path and os.path.exists(flat_path):
            with fits.open(flat_path, memmap=False) as flat_hdul:
                master_flat = flat_hdul[0].data.astype(float)
                #Flats need a little bump at least on my test bed perhaps this needs to be dynamically calculated and applied
                median_flat = np.median(master_flat)
                if median_flat > 0:
                    normalized_flat = (master_flat / median_flat) ** 1.09  # Slight exaggeration
                    light_data /= normalized_flat
                    light_header['CALFLAT'] = (os.path.basename(flat_path), 'Flat field used (exp=1.1)')
                    print(f"🧪 Applied master flat for filter '{filter_name}'")
                else:
                    print(f"⚠️ Flat median is zero in {flat_path}, skipping flat normalization.")
    return light_data, light_header

def create_master_frame(image_list, method='median', dark_flat_path=None, master_bias=None):
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
            if master_bias is not None:
                data -= master_bias
            stack.append(data)

    if not stack:
        return None

    if method == 'median':
        return np.median(stack, axis=0).astype(np.float32)
    else:
        return np.mean(stack, axis=0).astype(np.float32)

def save_master_frame(data, header, output_folder, name):
    output_path = os.path.join(output_folder, f"{name}.fits")
    fits.writeto(output_path, data, header=header, overwrite=True)
    return output_path

def load_fits_by_filter(file_list):
    from collections import defaultdict
    filtered = defaultdict(list)
    for path in file_list:
        try:
            hdr = fits.getheader(path)
            filter_name = normalize_filter_name(hdr.get('FILTER'))
            filtered[filter_name].append(path)
        except Exception:
            filtered['UNKNOWN'].append(path)
    return filtered

def safe_load_by_exptime(fits_paths):
    from collections import defaultdict
    result = defaultdict(list)
    for path in fits_paths:
        try:
            hdr = fits.getheader(path)
            exptime = float(hdr.get('EXPTIME', -1))
            if exptime > 0:
                result[exptime].append(path)
            else:
                print(f"⚠️ Skipping dark: {path} – invalid or missing EXPTIME")
        except Exception as e:
            print(f"⚠️ Failed to read EXPTIME for {path}: {e}")
    return dict(result)


def calibrate_and_save(light_path, master_dark_paths, master_flat_paths, master_bias_paths, calibrated_folder):
    try:
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

        print("\n".join(log_lines))

        base_name = os.path.basename(light_path)
        output_path = os.path.join(
            calibrated_folder,
            f"{os.path.splitext(base_name)[0]}_cal.fits"
        )

        fits.writeto(output_path, calibrated_data, header=header, overwrite=True)

        # 🛰️ Inject WCS from sidecar only if not already present (skip object detection)
        try:
            with fits.open(output_path, mode='update') as hdul:
                hdr = hdul[0].header
                if not any(k in hdr for k in ('CD1_1', 'PC1_1', 'WCSAXES')):
                    if inject_wcs_from_sidecar(output_path):
                        inject_minimal_sip(hdr)
                        hdul.flush()
                    else:
                        print(f"⚠️ WCS injection from sidecar failed for {os.path.basename(output_path)}")
                else:
                    print(f"ℹ️ WCS already present in {os.path.basename(output_path)}, skipping injection.")
        except Exception as e:
            print(f"⚠️ Exception during WCS injection: {e}")

        return f"✅ Calibrated {os.path.basename(light_path)}"

    except Exception as e:
        return f"💥 Calibration failed for {os.path.basename(light_path)}: {e}"

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

def run_parallel_calibration(light_by_filter, dark_by_exptime, flat_by_filter,
                              bias_by_filter, output_folder, session_title,
                              log_callback, save_masters=False,
                              dark_flat_file_list=None):

    if log_callback is None:
        log_callback = print

    log_callback("🛠️ Starting master calibration frame creation...")

    def normalize_dict_keys(d):
        return {normalize_filter_name(k): v for k, v in d.items()}

    light_by_filter = normalize_dict_keys(light_by_filter)
    dark_by_exptime = normalize_dict_keys(dark_by_exptime)
    flat_by_filter = normalize_dict_keys(flat_by_filter)
    bias_by_filter = normalize_dict_keys(bias_by_filter)

    log_callback(f"📁 Filters in lights: {list(light_by_filter.keys())}")
    log_callback(f"📁 Exposures in darks: {[round(float(k), 2) if isinstance(k, (int, float, str)) and str(k).replace('.', '', 1).isdigit() else k for k in dark_by_exptime.keys()]}")
    log_callback(f"📁 Filters in flats: {list(flat_by_filter.keys())}")
    log_callback(f"📁 Filters in biases: {list(bias_by_filter.keys())}")

    calibrated_folder = os.path.join(output_folder, "calibrated")
    os.makedirs(calibrated_folder, exist_ok=True)

    used_filters = set(light_by_filter.keys()) if not save_masters else set(
        list(flat_by_filter.keys()) | list(bias_by_filter.keys())
    )

    master_bias_paths = {}
    all_biases = [f for k, v in bias_by_filter.items() if k in used_filters for f in v]
    if all_biases:
        log_callback(f"🛠️ Creating global master bias from {len(all_biases)} files...")
        bias = create_master_frame(all_biases)
        if bias is not None:
            bias_path = save_master_frame(bias, fits.getheader(all_biases[0]), output_folder, "master_bias")
            for f in light_by_filter.keys():
                master_bias_paths[f] = bias_path
            log_callback(f"✅ Master bias saved to {bias_path}")
        else:
            log_callback("⚠️ Failed to create global master bias frame.")
    else:
        log_callback("⚠️ No bias frames provided.")

    master_dark_paths = {}
    all_darks = [f for v in dark_by_exptime.values() for f in v]
    if all_darks:
        log_callback(f"🛠️ Creating global master dark from {len(all_darks)} files...")
        bias_path = next(iter(master_bias_paths.values()), None)
        master_bias = load_fits_data(bias_path) if bias_path else None
        dark = create_master_frame(all_darks, master_bias=master_bias)
        if dark is not None:
            dark_path = save_master_frame(dark, fits.getheader(all_darks[0]), output_folder, "master_dark")
            for f in light_by_filter.keys():
                master_dark_paths[f] = dark_path
            log_callback(f"✅ Master dark saved to {dark_path}")
        else:
            log_callback("⚠️ Failed to create global master dark frame.")
    else:
        log_callback("⚠️ No dark frames provided.")

    grouped_dark_flats = group_dark_flats_by_filter_and_exptime(dark_flat_file_list or [])
    log_callback(f"🔍 Dark flat keys: {list(grouped_dark_flats.keys())}")

    def get_matching_dark_flat(flat_path):
        try:
            hdr = fits.getheader(flat_path)
            filt = normalize_filter_name(hdr.get("FILTER", "UNKNOWN"))
            exptime = round(float(hdr.get("EXPTIME", -1)), 1)
            return grouped_dark_flats.get((filt, exptime), [None])[0]
        except Exception as e:
            print(f"⚠️ Failed to match dark flat: {e}")
            return None

    master_flat_paths = {}
    for filter_name, flat_paths in flat_by_filter.items():
        if not save_masters and filter_name not in used_filters:
            continue
        log_callback(f"🧪 Calibrating flats for {filter_name} ({len(flat_paths)} frames)")
        matched_dark = get_matching_dark_flat(flat_paths[0])
        if matched_dark:
            log_callback(f"🌑 Using dark flat: {os.path.basename(matched_dark)}")
        else:
            log_callback(f"⚠️ No matching dark flat for {filter_name} – bias-only flat calibration")
        flat = create_master_flat_scaled(flat_paths, dark_flat_path=matched_dark)
        if flat is not None:
            flat_median = np.median(flat)
            log_callback(f"📊 Master flat median for {filter_name}: {flat_median:.2f}")
            path = save_master_frame(flat, fits.getheader(flat_paths[0]), output_folder, f"{filter_name}_master_flat")
            master_flat_paths[filter_name] = path
        else:
            log_callback(f"⚠️ Failed to create master flat for {filter_name}")

    futures = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for filter_name, light_paths in light_by_filter.items():
            for light_path in light_paths:
                futures.append(executor.submit(
                    calibrate_and_save,
                    light_path,
                    master_dark_paths,
                    master_flat_paths,
                    master_bias_paths,
                    calibrated_folder
                ))

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

    if save_masters:
        log_callback("📦 Creating ZIP archive of master calibration frames...")
        session_date = datetime.now().strftime("%Y-%m-%d")
        object_safe = session_title.replace(' ', '_').replace(':', '_').replace('/', '_') or 'UnknownObject'
        zip_name = f"{object_safe}_{session_date}_masters.zip"
        zip_path = os.path.join(output_folder, zip_name)
        master_files = list(master_bias_paths.values()) + list(master_dark_paths.values()) + list(master_flat_paths.values())
        with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_STORED) as zipf:
            for path in master_files:
                if path and os.path.exists(path):
                    zipf.write(path, arcname=os.path.basename(path))
        log_callback(f"📦 Master calibration frames zipped successfully: {zip_name}")

    return master_dark_paths, master_flat_paths, master_bias_paths