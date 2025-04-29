from astropy.io import fits
import numpy as np
import os
import zipfile
import concurrent.futures
from datetime import datetime
import gc
from wcs_utils import inject_wcs_from_sidecar
import shutil

def inject_minimal_sip(header):
    """Inject minimal fake SIP headers for compatibility with picky software like Tycho."""
    try:
        header['A_ORDER'] = 0
        header['B_ORDER'] = 0
        header['AP_ORDER'] = 0
        header['BP_ORDER'] = 0
        print("✅ Minimal SIP keywords injected for compatibility.")
    except Exception as e:
        print(f"⚠️ Failed to inject minimal SIP keywords: {e}")


def load_fits_data(path):
    if os.path.exists(path):
        with fits.open(path, memmap=False) as hdul:
            return hdul[0].data.astype(np.float32)
    return None

def calibrate_image(light_path, use_master=False, master_dark_path=None, master_flat_path=None, master_bias_path=None):
    with fits.open(light_path, memmap=False) as hdul:
        light_data = hdul[0].data.astype(float)
        light_header = hdul[0].header.copy()

    exposure_light = light_header.get('EXPTIME', None)

    if use_master and master_bias_path and os.path.exists(master_bias_path):
        with fits.open(master_bias_path, memmap=False) as bias_hdul:
            master_bias = bias_hdul[0].data.astype(float)
        light_data -= master_bias

    if use_master and master_dark_path and os.path.exists(master_dark_path):
        with fits.open(master_dark_path, memmap=False) as dark_hdul:
            master_dark = dark_hdul[0].data.astype(float)
            exposure_dark = dark_hdul[0].header.get('EXPTIME', None)

        if exposure_light and exposure_dark and exposure_dark != 0:
            scale_factor = exposure_light / exposure_dark
            light_data -= master_dark * scale_factor
        else:
            light_data -= master_dark

    if use_master and master_flat_path and os.path.exists(master_flat_path):
        with fits.open(master_flat_path, memmap=False) as flat_hdul:
            master_flat = flat_hdul[0].data.astype(float)
        normalized_flat = master_flat / np.median(master_flat)
        light_data /= normalized_flat

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

def calibrate_and_save(light_path, master_dark_path, master_flat_path, master_bias_path, calibrated_folder):
    calibrated_data, header = calibrate_image(
        light_path,
        use_master=True,
        master_dark_path=master_dark_path,
        master_flat_path=master_flat_path,
        master_bias_path=master_bias_path
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

def run_parallel_calibration(light_images, dark_images, flat_images, bias_images, output_folder, session_title="UnknownObject", log_callback=None, save_masters=False):
    if log_callback is None:
        log_callback = print

    log_callback("🛠️ Starting master calibration frame creation...")

    master_dark = None
    master_flat = None
    master_bias = None
    master_dark_flat = None

    dark_flat_images = [img for img in dark_images if 'flat' in img.lower()]

    # Build masters in parallel
    tasks = {}
    with concurrent.futures.ThreadPoolExecutor() as executor:
        if bias_images:
            tasks['bias'] = executor.submit(build_and_save_master, bias_images, output_folder, "master_bias")
        if dark_images:
            tasks['dark'] = executor.submit(build_and_save_master, dark_images, output_folder, "master_dark")
        if dark_flat_images:
            tasks['dark_flat'] = executor.submit(build_and_save_master, dark_flat_images, output_folder, "master_dark_flat")
        if flat_images:
            tasks['flat'] = executor.submit(build_and_save_master, flat_images, output_folder, "master_flat",
                                            dark_flat_path=os.path.join(output_folder, "master_dark_flat_master.fits") if dark_flat_images else None)

        for name, task in tasks.items():
            result = task.result()
            if result is not None and log_callback:
                log_callback(f"💾 Master {name.replace('_', ' ').title()} created.")

            if name == 'bias':
                master_bias = result
            elif name == 'dark':
                master_dark = result
            elif name == 'flat':
                master_flat = result
            elif name == 'dark_flat':
                master_dark_flat = result

    master_dark_path = os.path.join(output_folder, "master_dark_master.fits") if master_dark is not None else None
    master_flat_path = os.path.join(output_folder, "master_flat_master.fits") if master_flat is not None else None
    master_bias_path = os.path.join(output_folder, "master_bias_master.fits") if master_bias is not None else None

    calibrated_folder = os.path.join(output_folder, "calibrated")
    os.makedirs(calibrated_folder, exist_ok=True)

    log_callback(f"🔧 Starting calibration of {len(light_images)} light frames...")

    with concurrent.futures.ProcessPoolExecutor(max_workers=os.cpu_count() - 1) as executor:
        futures = [
            executor.submit(
                calibrate_and_save,
                light_path,
                master_dark_path,
                master_flat_path,
                master_bias_path,
                calibrated_folder,
            )
            for light_path in light_images
        ]

        completed = 0
        total = len(futures)
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()
                if result is not None:
                    output_path, full_log_message = result
                    log_callback(full_log_message)
                else:
                    log_callback("⚠️ Calibration returned None for a frame.")

                completed += 1
                if completed % 5 == 0 or completed == total:
                    log_callback(f"🖼️ Calibrated {completed}/{total} frames...")
            except Exception as e:
                log_callback(f"💥 Error in calibration: {e}")

    # 🧹 ZIP master calibration frames if Save Masters is enabled
    if save_masters:
        log_callback("📦 Creating ZIP archive of master calibration frames...")

        session_date = datetime.now().strftime("%Y-%m-%d")
        object_safe = session_title.replace(' ', '_').replace(':', '_').replace('/', '_') or 'UnknownObject'
        zip_name = f"{object_safe}_{session_date}_masters.zip"
        zip_path = os.path.join(output_folder, zip_name)

        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for file_name in os.listdir(output_folder):
                if file_name.endswith("_master.fits"):
                    full_path = os.path.join(output_folder, file_name)
                    zipf.write(full_path, arcname=file_name)

        log_callback(f"📦 Master calibration frames zipped successfully: {zip_name}")
    else:
        log_callback("ℹ️ Save Masters not enabled. Skipping master frames ZIP creation.")

    return master_dark, master_flat, master_bias
