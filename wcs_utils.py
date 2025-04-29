import os
from astropy.io import fits

def inject_wcs_from_sidecar(fits_path):
    """
    Attempts to inject WCS info into a calibrated FITS file
    from an ASTAP-created sidecar file (same name, .wcs extension).
    """
    sidecar_path = os.path.splitext(fits_path)[0] + '.wcs'

    if not os.path.exists(sidecar_path):
        print(f"ℹ️ No sidecar WCS file found for {os.path.basename(fits_path)}")
        return

    try:
        with fits.open(sidecar_path) as sidecar_hdul:
            wcs_header = sidecar_hdul[0].header

        with fits.open(fits_path, mode='update') as hdul:
            target_header = hdul[0].header

            # Inject essential WCS fields
            for key in ['CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2',
                        'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
                        'CTYPE1', 'CTYPE2', 'CUNIT1', 'CUNIT2']:
                if key in wcs_header:
                    target_header[key] = wcs_header[key]

            target_header['WCS_INJ'] = (True, 'WCS injected from sidecar file')
            hdul.flush()

        print(f"✅ WCS injected from sidecar into {os.path.basename(fits_path)}")

    except Exception as e:
        print(f"⚠️ Failed to inject WCS for {os.path.basename(fits_path)}: {e}")
