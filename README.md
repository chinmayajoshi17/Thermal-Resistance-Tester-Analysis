#README: Quantification and Mitigation of Uncertainties in Thermal Conductivity Measurements Using a Modified ASTM D5470 Thermal Resistance Tester

Description of the data and file structure
This dataset includes temperature measurements over time for 6 thermocouples for steady state analysis of thermal resistance and thermal conductivity of Pyrolytic Graphite and Titanium Grade 2. The temperature data is measured using OMEGA TJ36-CPSS-032G-3 T-type thermocouples with a 0.032” probe diameter connected to a NI-9210 Thermocouple DAQ with Spring Terminals, attached to a NI-cDAQ-9174 DAQ Chassis which is connected to a computer running Windows 10. Time data is sampled at 1000 Hz (every 1 ms), starting at 0 s and continuing until manually stopped. The data acquisition setup is used in the facility described in [INSERT PAPER NAME/LINK]. 

Files and variables:
File: DatasetsAndCode.zip
Description: The zip file includes 2 folders, i.e. “DataForPaper” and “PythonCode” 
DatasetsAndCode.zip  DataForPaper/PythonCode

The “DataForPaper” folder further contains two folders, i.e., PG for the Pyrolytic Graphite data, and TiG2 for the Titanium Grade 2 data. 
DataForPaper  PG/TiG2  WithTIM/NoTIM  Run X (X ∈ [1,2,3])  *.lvm

Each folder is then divided into sub-folders for “WithTIM” for datasets with TIMs applied, and “NoTIM” for datasets without any TIMs applied.
Each “WithTIM” and “NoTIM” folder is then sub-divided into “Run1”, “Run2”, and “Run3” which contain the data for each material, i.e. 9 files for PG corresponding to the 9 PG samples, and 4 files for TiG2 corresponding to the 4 TiG2 samples.
All the Data files are .lvm files with 23 lines of headers before the actual data starts. The data is divided into 7 columns delimited using “tabs” as follows: Time, Temperature at TC1, Temperature at TC2, Temperature at TC3, Temperature at TC4, Temperature at TC5, Temperature at TC6. Time is measured in seconds (s) and the temperature is measured in Celsius (oC).

The “PythonCode” folder contains multiple .ipynb files that are used for the steady state analysis, uncertainty analysis, thermal conductivity, and thermal resistance calculations.
PythonCode  *.ipynb/*.py

“SSCalc.ipynb” checks if the data collected has reached steady-state and compares the ΔT (change in temperature) with respect to time over a period of 120s or 2 minutes. It also checks the ΔT with respect to time over a period of 30s, 60s, and 180s to ensure that minor changes in T or sensor noise does not affect long term steady-state temperatures.
“AnalysisCodeTRT.ipynb” requires a defined material, TIM state (With or No TIM), and the repeated “Run Number” for the tests in Cell 2 of the code. It outputs the measured thermal conductivity of the materials and the uncertainty in the thermal resistance measurements.
“ScedasticityCheck.ipynb” uses the Uncertainty in each Thermal Resistance measurement from “AnalysisCodeTRT.ipynb” and the sample thicknesses to numerically check for homoscedasticity or heteroscedasticity in the data and visualize it using plots.
“YorkRegression.ipynb” uses “YorkRegression.py”, a modified version of “fit_bivariate.py” by Mikko Pitkanen, to compare regular Least Squares Regression methods to York’s Regression mentioned in York et al. 2004 and provides thermal conductivity and uncertainty values based on York’s Regression.
