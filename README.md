# LC-OD – Layer Charge Calculator (OD Method)

## Overview

**LC-OD** is a Shiny-based application for calculating the layer charge of smectites using the spectroscopic OD method from ATR-IR spectra (Kuligiewicz et al., 2015).

The tool is designed for researchers working with clay minerals and infrared spectroscopy, providing a reproducible and user-friendly workflow.

Input Data Format

The application accepts FTIR spectra stored in `.csv` files.

Each file must contain exactly **two columns**:

1) wavenumber (in cm⁻¹)
2) signal as ATR absorbance

Requirements

* Files must **not contain headers or column names**
* Supported column separators: `;`, `,`, or tab
* If `,` is used as the **column separator**, numeric values must use `.` as the **decimal separator**
* If `;` or tab is used as the **column separator**, both `.` and `,` are accepted as decimal separators
* Data must be numeric and contain no missing or invalid values

### Notes

* Spectra do not need to share the same wavenumber grid — the application will automatically align them
* Files that cannot be read or processed correctly will be skipped and reported in the application

---

## 🚀 For End Users (No Installation Required)

A fully portable version of the application is available.

### How to run:

1. Go to the **Releases** section of this repository
2. Download the latest version:

   ```
   LC-OD_x.x_Windows.zip
   ```
3. Extract the ZIP archive to a folder (e.g. `C:\LC-OD\`)
4. Open the folder and double-click:

   ```
   START_LCOD.bat
   ```
5. The application will automatically open in your web browser

---

### Important Notes

* Do **not** run the application directly from the ZIP file
* Do **not** move or delete the `R` or `Application` folders
* Do **not** close the console window while using the app

---

## 🧪 For Developers

This repository contains the source code and a reproducible R environment managed with renv.

### Requirements

* R (recommended ≥ 4.x)
* Internet connection (for restoring packages)

### Setup

```r
install.packages("renv")   # if needed
renv::restore()
```

### Run the application

```r
shiny::runApp("app.R")
```

or:

```r
source("launch_app.R")
```

---

## 📁 Repository Structure

```
lc-od/
  app.R
  launch_app.R
  renv.lock
  renv/
  README.md
  LICENSE
  .gitignore
```

---

## 📦 Portable Version

The portable distribution includes:

* embedded R environment
* all required packages
* preconfigured runtime

This version is intended for:

* users not familiar with the R environment
* teaching and workshops

---

## ⚠️ Notes

* The full R environment is **not included** in this repository
* Large files and binaries are distributed via Releases

---

## 📄 License

This project is licensed under the Mozilla Public License 2.0 (MPL-2.0).

---

## 👤 Author

Artur Kuligiewicz, Ph.D.
© 2026 Institute of Geological Sciences, Polish Academy of Sciences

---

