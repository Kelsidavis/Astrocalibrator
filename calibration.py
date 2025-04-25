from astropy.io import fits
import numpy as np
import os

def load_fits_data(path):
    if os.path.exists(path):
        with fits.open(path) as hdul:
            return hdul[0].data.astype(np.float32)
    return None

def calibrate_image(image_path, use_master=False, master_dark_path=None, master_flat_path=None, master_bias_path=None):
    from gui import master_dark_enabled, master_flat_enabled, master_bias_enabled

    with fits.open(image_path) as hdul:
        image_data = hdul[0].data.astype(np.float32)

    # Load master calibration files only if enabled via GUI
    master_bias = load_fits_data(master_bias_path) if use_master and master_bias_enabled.get() and master_bias_path else None
    master_dark = load_fits_data(master_dark_path) if use_master and master_dark_enabled.get() and master_dark_path else None
    master_flat = load_fits_data(master_flat_path) if use_master and master_flat_enabled.get() and master_flat_path else None

    # Apply master bias
    if master_bias is not None:
        image_data -= master_bias

    # Apply master dark
    if master_dark is not None:
        image_data -= master_dark

    # Apply flat correction
    if master_flat is not None:
        flat_corrected = master_flat / np.median(master_flat)
        image_data /= flat_corrected

    return image_data

def create_master_frame(image_list, method='median'):
    if not image_list:
        return None

    stack = []
    for path in image_list:
        data = load_fits_data(path)
        if data is not None:
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

def run_parallel_calibration(light_images, dark_images, flat_images, bias_images, output_folder):
    # Placeholder: define if using multiprocessing in the future
    pass
