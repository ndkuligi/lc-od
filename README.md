# LC-OD – Layer Charge Calculator (OD Method)

## Overview

**LC-OD** is a Shiny-based application for calculating the layer charge of smectites using the spectroscopic OD method from ATR-IR spectra (Kuligiewicz et al., 2015).

Kuligiewicz, A., Derkowski, A., Emmerich, K., Christidis, G. E., Tsiantos, C., Gionis, V., & Chryssikos, G. D. (2015). Measuring the Layer Charge of Dioctahedral Smectite by O–D Vibrational Spectroscopy. Clays and Clay Minerals, 63(6), 443-456. https://doi.org/10.1346/ccmn.2015.0630603 

The tool is designed for researchers working with clay minerals and infrared spectroscopy, providing a reproducible and user-friendly workflow.

Input Data Format

The application accepts FTIR spectra stored in `.csv`, '.dpt', or '.xy' files.

Each file must contain wavenumbers (in cm⁻¹) and signal (ATR absorbance) in the first two columns.

Requirements

* Files must **not contain headers or column names**
* Supported column separators: `;`, `,`, or tab
* If `,` is used as the **column separator**, numeric values must use `.` as the **decimal separator**
* If `;` or tab is used as the **column separator**, both `.` and `,` are accepted as decimal separators
* Data must be numeric and contain no missing or invalid values

### Notes

* Spectra do not need to share the same wavenumber grid — the application will automatically align them
* Files that cannot be read or processed correctly will be skipped and reported in the application

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

A portable Windows version is available from the latest release:

https://github.com/ndkuligi/lc-od/releases/latest/download/LC-OD-portable-Windows.zip

The portable distribution includes:

* embedded R environment
* all required packages
* preconfigured runtime

This version is intended for:

* users not familiar with the R environment
* teaching and workshops

---
### Important Notes

* Do **not** run the application directly from the ZIP file
* Do **not** move or delete the `R` or `Application` folders
* Do **not** close the console window while using the app

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

## Version 1.1.0 — main updates

### New features

* Added support for `.xy` and `.dpt` spectral files
* Improved automatic separator detection
* Added uncertainty estimates for layer-charge calculations
* Improved robustness of spectral import

### Methodological changes

* Expanded O–D minimum search range from `2650–2720 cm⁻¹` to `2600–2720 cm⁻¹`

### Export improvements

* Exported statistics now include:

  * `LC_SFM_ERROR`
  * `LC_AAM_ERROR`


