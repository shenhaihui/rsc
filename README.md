# Source Codes of Ranking and Selection with Covariates

This repository contains the source codes used in the following paper:
* Shen, Haihui, L. Jeff Hong, and Xiaowei Zhang (2020). Ranking and selection with covariates for personalized decision making. Submitted to *INFORMS Journal on Computing*.

### Disclaimer
LICENSE: Redistribution and use in source and binary forms, with or without modification, are permitted, under the terms of the [BSD license](./BSD_License.txt).

WARNING: These codes are written only for the purpose of demonstration and verification. While the correctness has been carefully checked, the quality such as standardability,
clarity, generality, and efficiency has not been well considered.

## 1 Introduction
The `\R&S-C2020` folder contains the entire MATLAB implementations of all the numerical experiments and case study in Shen et al. (2020). The codes are split into the following
three folders.
* `\NumericalExperiments`, contains codes for numerical experiments in Shen et al. (2020, §6, EC.6).
* `\RobustnessTest`, contains codes of the robustness test to linearity assumptions in Shen et al. (2020, §7.2).
* `\EsophagealCancer`, contains codes of the case study of personalized treatment for cancer prevention in Shen et al. (2020, §8).

## 2 Installation
The codes were written and run in MATLAB R2018b, on Windows 7 Enterprise 64-bit OS,
with Intel i9-9900K CPU @ 3.60 GHz, 16 GB RAM.

To install the MATLAB codes, just copy the entire folder `\R&S-C2020` into your MATLAB directory, or change the path of MATLAB to the folder `\R&S-C2020`

## 3 Details on Numerical Experiments and Case Study
### 3.1 Numerical Experiments in Shen et al. (2020, §6, EC.6)
Get into folder `\NumericalExperiments`.
* To run Procedure TS and TS+ for all the 13 problems with pre-specified target
<img src="https://latex.codecogs.com/svg.latex?{\text{PCS}_{\text{E}}\geq{1-\alpha}}">, run script `PCSE.m`, which is self-explanatory.
 Scripts `find_h_hom_PCSE.m` and `find_h_hom_PCSE_SA.m` are functions to calculates the constant <img src="https://latex.codecogs.com/svg.latex?{h}"> of Procedure TS for homoscedastic errors, where the the former uses numerical integration (`integral` or `trapz`) and MATLAB built-in root finding function `fzero`, while the latter uses stochastic approximation method.
 Similarly, scripts `find_h_het_PCSE.m` and `find_h_het_PCSE_SA.m` are functions to calculates the constant <img src="https://latex.codecogs.com/svg.latex?{h_{\text{Het}}}"> of Procedure TS+ for heteroscedastic errors.
 
* To run Procedure TS and TS+ for all the 13 problems with pre-specified target
<img src="https://latex.codecogs.com/svg.latex?{\text{PCS}_{\text{min}}\geq{1-\alpha}}">, run script `PCSmin.m`, which is self-explanatory.
 Scripts `find_h_hom_PCSmin.m` and `find_h_het_PCSmin.m` are functions to calculates the constant <img src="https://latex.codecogs.com/svg.latex?{h}"> of Procedure TS and constant <img src="https://latex.codecogs.com/svg.latex?{h_{\text{Het}}}"> of Procedure TS+ , respectively.

### 3.2 Robustness Test in Shen et al. (2020, §7.2)
Get into folder `\RobustnessTest`, and run script `RobustnessTest.m`, which is self-explanatory.

### 3.3 Case Study in Shen et al. (2020, §8)
Get into folder `\EsophagealCancer`.
* The script `EsophagealCancerSim.m` is the implementation of a discrete-time Markov chain model that simulates the progression of Barrett’s esophagus (BE) to esophageal
adenocarcinoma (EAC), which was developed by Hur et al. (2004) and Choi et al. (2014). See more details [here](https://simopt.github.io/ECSim).

* The script `BruteForceSim.m` runs the simulation model for 10^6 replications on a grid of points to approximate the true response surfaces, whose results are saved in `QALY_true.mat`.

* To select personalized treatment regimen with Procedure TS+ and compare the performance with traditional approach, run script `PersonalizedTreatment.m`, which
is self-explanatory.

## References
* Choi, Sung Eun, Katherine E. Perzan, Angela C. Tramontano, Chung Yin Kong, and Chin Hur (2014). Statins and aspirin for chemoprevention in Barrett’s esophagus: Results
of a cost-effectiveness analysis. *Cancer Prevention Research* **7**(3), 341–350.

* Hur, Chin, Norman S. Nishioka, and G. Scott Gazelle (2004). Cost-effectiveness of aspirin chemoprevention for Barrett’s esophagus. *Journal of the National Cancer Institute* **96**(4), 316–325.
