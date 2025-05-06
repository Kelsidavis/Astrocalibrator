import os
from astropy.io import fits

def inject_wcs_from_sidecar(fits_path):
    """
    Injects WCS info from a plain-text ASTAP-generated .wcs file into the FITS header.
    """
    sidecar_path = os.path.splitext(fits_path)[0] + '.wcs'

    if not os.path.exists(sidecar_path):
        print(f"ℹ️ No sidecar WCS file found for {os.path.basename(fits_path)}")
        return

    wcs_info = {}
    try:
        with open(sidecar_path, 'r') as f:
            for line in f:
                if 'RA center' in line:
                    wcs_info['CRVAL1'] = float(line.split(':')[1].strip())
                elif 'DEC center' in line:
                    wcs_info['CRVAL2'] = float(line.split(':')[1].strip())
                elif 'Field rotation' in line:
                    wcs_info['CROTA2'] = float(line.split(':')[1].strip())
                elif 'Pixel scale' in line:
                    pix_scale = float(line.split(':')[1].split()[0])
                    wcs_info['CDELT1'] = -pix_scale / 3600.0
                    wcs_info['CDELT2'] = pix_scale / 3600.0

        ra = wcs_info.get('CRVAL1')
        dec = wcs_info.get('CRVAL2')
        if ra is None or dec is None:
            print(f"⚠️ Missing RA/Dec in WCS sidecar for {os.path.basename(fits_path)}")
            return

        if not (-90.0 <= dec <= 90.0):
            print(f"⚠️ Skipping WCS injection due to invalid declination: {dec}")
            return

        with fits.open(fits_path, mode='update') as hdul:
            hdr = hdul[0].header
            naxis1 = hdr.get('NAXIS1', 2048)
            naxis2 = hdr.get('NAXIS2', 2048)

            hdr['CRVAL1'] = ra
            hdr['CRVAL2'] = dec
            hdr['CRPIX1'] = naxis1 / 2
            hdr['CRPIX2'] = naxis2 / 2
            hdr['CD1_1'] = wcs_info.get('CDELT1', -0.000277778)
            hdr['CD1_2'] = 0.0
            hdr['CD2_1'] = 0.0
            hdr['CD2_2'] = wcs_info.get('CDELT2', 0.000277778)
            #Append -SIP
            ctype1 = hdr.get('CTYPE1', 'RA---TAN')
            ctype2 = hdr.get('CTYPE2', 'DEC--TAN')
            if '-SIP' not in ctype1:
                ctype1 += '-SIP'
            if '-SIP' not in ctype2:
                ctype2 += '-SIP'
            hdr['CTYPE1'] = ctype1
            hdr['CTYPE2'] = ctype2
            hdr['CUNIT1'] = 'deg'
            hdr['CUNIT2'] = 'deg'
            hdr['RADECSYS'] = 'ICRS'
            hdr['EQUINOX'] = 2000.0
            if 'CROTA2' in wcs_info:
                hdr['CROTA2'] = wcs_info['CROTA2']
            hdr['WCS_INJ'] = (True, 'WCS injected from ASTAP sidecar')

            hdul.flush()

        print(f"✅ WCS injected from sidecar into {os.path.basename(fits_path)}")

    except Exception as e:
        print(f"⚠️ Failed to inject WCS for {os.path.basename(fits_path)}: {e}")
