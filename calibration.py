from astropy.io import fits
import numpy as np
import os
import zipfile
import concurrent.futures
from datetime import datetime

def load_fits_data(path):
    if os.path.exists(path):
        with fits.open(path) as hdul:
            return hdul[0].data.astype(np.float32)
    return None

def calibrate_image(light_path, use_master=False, master_dark_path=None, master_flat_path=None, master_bias_path=None):
    with fits.open(light_path) as hdul:
        light_data = hdul[0].data.astype(float)
        light_header = hdul[0].header

    exposure_light = light_header.get('EXPTIME', None)

    if use_master and master_bias_path and os.path.exists(master_bias_path):
        with fits.open(master_bias_path) as bias_hdul:
            master_bias = bias_hdul[0].data.astype(float)
        light_data -= master_bias

    if use_master and master_dark_path and os.path.exists(master_dark_path):
        with fits.open(master_dark_path) as dark_hdul:
            master_dark = dark_hdul[0].data.astype(float)
            dark_header = dark_hdul[0].header
            exposure_dark = dark_header.get('EXPTIME', None)

        if exposure_light and exposure_dark and exposure_dark != 0:
            scale_factor = exposure_light / exposure_dark
            light_data -= master_dark * scale_factor
        else:
            light_data -= master_dark

    if use_master and master_flat_path and os.path.exists(master_flat_path):
        with fits.open(master_flat_path) as flat_hdul:
            master_flat = flat_hdul[0].data.astype(float)
        normalized_flat = master_flat / np.median(master_flat)
        light_data /= normalized_flat

    return light_data

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

def calibrate_and_save(light_path, master_dark_path, master_flat_path, master_bias_path, calibrated_folder, log_callback=None):
    calibrated_data = calibrate_image(
        light_path,
        use_master=True,
        master_dark_path=master_dark_path,
        master_flat_path=master_flat_path,
        master_bias_path=master_bias_path
    )

    header = fits.getheader(light_path)

    # Add calibration information
    header['CALIB'] = (True, 'Frame has been calibrated')
    header.add_history('Calibrated using Astrocalibrator')

    # Log preserved important fields
    important_fields = ['OBJECT', 'OBSERVER', 'INSTRUME', 'FILTER', 'DATE-OBS', 'EXPTIME', 'TELESCOP']
    preserved = []
    missing = []

    for key in important_fields:
        if key in header:
            preserved.append(f"{key}='{header[key]}'")
        else:
            missing.append(key)

    # 🛠️ Suppress 'FILTER' warning if image is likely OSC
    likely_osc = ('BAYERPAT' in header) or ('BAYER' in str(header)) or ('FILTER' not in header)

    if log_callback:
        log_callback(f"🧾 Preserved header fields for {os.path.basename(light_path)}:")
        for item in preserved:
            log_callback(f"   - {item}")

        if missing:
            if not (likely_osc and missing == ['FILTER']):
                log_callback(f"⚠️ Missing expected fields: {', '.join(missing)}")

    base_name = os.path.basename(light_path)
    output_path = os.path.join(
        calibrated_folder,
        f"{os.path.splitext(base_name)[0]}_cal.fits"
    )

    fits.writeto(output_path, calibrated_data, header=header, overwrite=True)
    return output_path

def run_parallel_calibration(light_images, dark_images, flat_images, bias_images, output_folder, session_title="UnknownObject", log_callback=None):
    if log_callback is None:
        log_callback = print  # fallback if no logger passed

    master_dark = None
    master_flat = None
    master_bias = None
    master_dark_flat = None

    log_callback("🛠️ Starting master calibration frame creation...")

    if bias_images:
        master_bias = create_master_frame(bias_images)
        if master_bias is not None:
            save_master_frame(master_bias, fits.getheader(bias_images[0]), output_folder, "master_bias")
            log_callback(f"💾 Master Bias created from {len(bias_images)} files.")

    if dark_images:
        master_dark = create_master_frame(dark_images)
        if master_dark is not None:
            save_master_frame(master_dark, fits.getheader(dark_images[0]), output_folder, "master_dark")
            log_callback(f"💾 Master Dark created from {len(dark_images)} files.")

    dark_flat_images = [img for img in dark_images if 'flat' in img.lower()]
    if dark_flat_images:
        master_dark_flat = create_master_frame(dark_flat_images)
        if master_dark_flat is not None:
            save_master_frame(master_dark_flat, fits.getheader(dark_flat_images[0]), output_folder, "master_dark_flat")
            log_callback(f"💾 Master Dark-Flat created from {len(dark_flat_images)} files.")

    if flat_images:
        dark_flat_path = os.path.join(output_folder, "master_dark_flat_master.fits") if master_dark_flat is not None else None
        master_flat = create_master_frame(flat_images, method='median', dark_flat_path=dark_flat_path)
        if master_flat is not None:
            save_master_frame(master_flat, fits.getheader(flat_images[0]), output_folder, "master_flat")
            log_callback(f"💾 Master Flat created from {len(flat_images)} files.")

    master_dark_path = os.path.join(output_folder, "master_dark_master.fits") if master_dark is not None else None
    master_flat_path = os.path.join(output_folder, "master_flat_master.fits") if master_flat is not None else None
    master_bias_path = os.path.join(output_folder, "master_bias_master.fits") if master_bias is not None else None

    calibrated_folder = os.path.join(output_folder, "calibrated")
    os.makedirs(calibrated_folder, exist_ok=True)

    log_callback(f"🔧 Starting calibration of {len(light_images)} light frames...")

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(
                calibrate_and_save,
                light_path,
                master_dark_path,
                master_flat_path,
                master_bias_path,
                calibrated_folder,
                log_callback
            )
            for light_path in light_images   # <<< 🛠️ This was missing too!
        ]

        completed = 0
        total = len(futures)
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()
                completed += 1
                if completed % 5 == 0 or completed == total:
                    log_callback(f"🖼️ Calibrated {completed}/{total} frames...")
            except Exception as e:
                log_callback(f"💥 Error in calibration: {e}")

    log_callback("📦 Creating ZIP archive of calibrated frames...")

    session_date = datetime.now().strftime("%Y%m%d")
    object_safe = session_title.replace(' ', '_').replace(':', '_') or 'UnknownObject'
    zip_name = f"{object_safe}_{session_date}.zip"
    zip_path = os.path.join(output_folder, zip_name)

    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for file_name in os.listdir(calibrated_folder):
            full_path = os.path.join(calibrated_folder, file_name)
            zipf.write(full_path, arcname=file_name)

    log_callback(f"📦 Calibrated frames zipped successfully: {zip_name}")

    return master_dark, master_flat, master_bias
