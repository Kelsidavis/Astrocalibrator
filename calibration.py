import os
import numpy as np
from astropy.io import fits
from collections import defaultdict
from datetime import datetime
from multiprocessing import Pool


def load_fits_by_filter(file_list):
    from gui import log_message
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
    from gui import log_message
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


def run_parallel_calibration(tasks, threads, progress_var=None):
    from gui import log_message, root
    from time import time
    log_message("⚙️ Beginning parallel calibration using multiprocessing...")
    start_time = time()
    results = []
    with Pool(processes=threads) as pool:
        for i, result in enumerate(pool.imap_unordered(calibrate_image, tasks), 1):
            results.append(result)
            if result.startswith("Error"):
                log_message(result)
            else:
                log_message(f"Saved calibrated: {result}")
            elapsed = time() - start_time
            eta = (elapsed / i) * (len(tasks) - i) if i > 0 else 0
            if progress_var:
                progress_var.set((i / len(tasks)) * 100)
            root.update_idletasks()
    return results
