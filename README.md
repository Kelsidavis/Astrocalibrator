# Astrocalibrator

Astrocalibrator is a FITS image calibration and plate solving tool built in Python with a GUI interface. It supports multi-threaded processing, filter-matched dark/flat calibration, and integration with ASTAP and Simbad for solving and object data.

## 🚀 Features

- GUI interface using `tkinter`
- Calibration of FITS images using dark and flat frames
- Multi-threaded processing
- Plate solving with ASTAP
- Object lookup with `astroquery` (Simbad)
- Log output and session tracking

## 🛠️ Requirements

Make sure you have Python 3.8+ installed. You can then install the required dependencies using pip.

### Install dependencies:

```bash
pip install pillow astroquery
```

You may also need:

- **ASTAP** installed and added to your system's PATH for plate solving.
- Internet connection for `astroquery`-based Simbad queries.

## 📦 Running the App

```bash
python main.py
```

