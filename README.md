# Astrocalibrator

**Astrocalibrator** is a Python-based GUI tool for astronomers that streamlines the calibration of FITS images. Designed with simplicity and automation in mind, it supports full calibration workflows with optional plate solving and archive packaging.

## 🚀 Features
- Graphical interface using Tkinter
- Batch support for Light, Dark, Flat, and Dark Flat frames
- Automatic calibration with multiprocessing
- Optional plate solving using ASTAP
- Saves final output as a ZIP archive with session naming
- Option to save master calibration frames
- Progress bar with visual feedback and log console
- Session naming based on plate-solved object

## 📦 Output Structure
- Calibrated images saved in a temporary folder and then archived
- Final ZIP placed in the parent directory of the light frames folder
- File name format: `[ObjectName]_[Date].zip`

## 🛠 Requirements
- Python 3.x
- `astropy`, `tkinter`, `astap` (installed separately)
- Windows (full support), others may require adaptation

## 🧰 Usage
1. Launch `main.py`
2. Select FITS files for Light, Dark, Flat, and Dark Flats
3. Click **"Solve & Calibrate"**
4. Optionally check **"Save Calibration Masters"** to keep master frames
5. Output will be zipped and placed in the appropriate folder

## 🧠 Notes
- ASTAP must be installed and accessible via system PATH
- Solving failures will be logged and skipped gracefully
- Session title is inferred from the first successful solve result

## 📜 License
MIT License

