# 🌌 AstroCalibrator

**AstroCalibrator** is a standalone FITS calibration and plate-solving GUI for amateur astrophotographers, featuring dark/flat correction, per-filter stacking with sigma clipping, automatic object identification via SIMBAD, and convenient archive output.

> 🪐 Designed for automation, ease of use, and integration with ASTAP and Simbad.

---

## 🚀 Features

- 📁 Load and organize light, dark, flat, and dark-flat frames
- 🧮 Create master calibration frames (optional save)
- 📊 Stack light frames per filter (median or mean)
- 🔍 Sigma clipping for stacking with adjustable threshold and iteration
- 📌 Automatic plate solving with **ASTAP**
- 🛰️ Simbad object identification (fallback to Messier catalog)
- 🧠 Integrated metadata: session title, RA/DEC, total exposure, constellation
- 📦 Output ZIP archive containing:
  - Calibrated & stacked frames
  - Readme metadata
- 🎛️ Multiprocessing (thread count control)
- 💾 Settings persistence (sigma values, stacking method, etc.)
- 🎨 Tooltips, adjustable options, and cheerful UI extras

---

## 🛠 Requirements

- Python 3.8+
- Packages:
