# Astrocalibrator

**Astrocalibrator** is a Python-based GUI tool for astronomers designed primarily for the archival of FITS data for later scientific processing. It streamlines image calibration while providing organized, compressed outputs for long-term storage. Designed with simplicity and automation in mind, it supports full calibration workflows with optional plate solving and archive packaging.

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
- Output ZIP contains calibrated light frames and optionally master calibration frames
- Temporary files are cleaned after zipping
- Logs are shown live in the application
- Calibrated images saved in a temporary folder and then archived
- Final ZIP placed in the parent directory of the light frames folder
- File name format: `[ObjectName]_[Date].zip`

## 🛠 Installation

1. **Clone the repository**:
   ```bash
   git clone https://github.com/Kelsidavis/astrocalibrator.git
   cd astrocalibrator
   ```

2. **Install required Python packages**:
   ```bash
   pip install astropy
   ```
   Note: `tkinter` is typically included with Python installations. If not, refer to your OS-specific installation method.

3. **Install ASTAP**:
   - Download and install ASTAP from [https://www.hnsky.org/astap.htm](https://www.hnsky.org/astap.htm).
   - Ensure that the ASTAP executable is added to your system's PATH.

4. **Run the application**:
   - Double-click the `start.bat` file provided in the repo.
   - Or launch manually with:
     ```bash
     python main.py
     ```


- Python 3.x
- `astropy`, `tkinter`, `astap` (installed separately)
- Windows (full support), others may require adaptation

## 🧰 Usage
1. Launch `main.py`
2. Select FITS files for Light, Dark, Flat, and Dark Flats
3. Click **"Solve & Calibrate"**
4. Optionally check **"Save Calibration Masters"** to keep master frames
5. Output will be zipped and placed in the appropriate folder

## 🤝 Contributing
Contributions, feature requests, and issue reports are welcome! Feel free to fork the repo and submit a pull request.

## 🧠 Notes
- ASTAP must be installed. By default, the application attempts to locate it in Program Files; adding it to your system PATH is optional but recommended
- Solving failures will be logged and skipped gracefully
- Session title is inferred from the first successful solve result

## 📜 License
MIT License

