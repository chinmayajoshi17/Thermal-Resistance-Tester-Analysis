**README: Modified ASTM D5470 Thermal Resistance Tester**

**Paper Name: Quantification and Mitigation of Uncertainties in Thermal Conductivity Measurements Using a Modified ASTM D5470 Thermal Resistance Tester**

- **Sampling Rate:** 1000 Hz (1 ms interval), starting from 0 s  
- **Thermocouples:** OMEGA TJ36-CPSS-032G-3 T-type, 0.032” probe diameter  
- **DAQ:** NI-9210 with spring terminals + NI-cDAQ-9174 chassis  

Each material (PG, TiG2) contains:
- `WithTIM/` – Measurements with thermal interface material applied  
- `NoTIM/` – Measurements without TIM  
- Each has `Run1/`, `Run2/`, `Run3/` folders with individual sample files  
- PG: 9 samples -> 9 *.lvm files for each run
- TiG2: 4 samples -> 4 *.lvm files for each run

---

## Python Notebooks

> All code is written in Python 3 using Jupyter notebooks.

###  `SSCalc.ipynb`
Checks if the system has reached steady-state by analyzing ΔT (temperature change) over 30s, 60s, 120s, and 180s intervals. This helps verify stability and considers noise or small drifts.

###  `AnalysisCodeTRT.ipynb`
- Inputs: `materialName`, `TIMstate`, `RunNum` (defined in Cell 2)
- Outputs:  
- Thermal conductivity (W/m·K) using data from a run
- Uncertainty in thermal resistance (m²·K/W) for each test in a run

###  `ScedasticityCheck.ipynb`
Takes uncertainty values and sample thicknesses to numerically check for **Homoscedasticity** or **Heteroscedasticity** and visualizes the results.

###  `YorkRegression.ipynb`
Uses `YorkRegression.py`, a modified version of Mikko Pitkänen’s `fit_bivariate.py`, to compare:
- Ordinary Least Squares (OLS)
- **York’s Regression** (York et al., 2004)

- `fit_bivariate.py` from https://gist.github.com/mikkopitkanen/ce9cd22645a9e93b6ca48ba32a3c85d0

Outputs updated thermal conductivity and uncertainty values using York’s method to account for uncertainties in thickness measurement (x-axis) and thermal resistance measurement (y-axis).

---

To run the notebooks, install and import the following Python libraries:

```bash
pip install numpy
pip install pandas
pip install matplotlib
pip install scipy
pip install os
pip install glob
```


---
