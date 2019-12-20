# Systematic differences in bucket sea surface temperatures caused by misclassification of engine room intake measurements

<br>

![GitHub last commit](https://img.shields.io/github/last-commit/duochanatharvard/LME_Offsets_vs_Diurnal_Amplitudes)
![GitHub repo size](https://img.shields.io/github/repo-size/duochanatharvard/LME_Offsets_vs_Diurnal_Amplitudes)


Matlab scripts associated with the paper "Systematic differences in bucket sea surface temperatures caused by misclassification of engine room intake measurements" by [Duo Chan](https://github.com/duochanatharvard) and Peter Huybers.

All codes are [Matlab](https://www.mathworks.com/products/matlab.html) .m files.  We provide a [main](DA_LME_main.m) script for fast reproduction of Figures and bucket model simulations in the main text.  If you are reproducing the [full analysis](#full-analysis), we also provide codes and step-by-step instructions for running these codes.  "Parallel" computations using multiple CPUs are encouraged are some steps, though fast-reproduction is runnable on a laptop.  

**Dependencies**
<!--
* [**Matlab m_map toolbox**](https://www.eoas.ubc.ca/~rich/map.html).
-->

* [**distinguishable_colors**](https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors) for generating figures.

* [**sun_position**](https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/5430/versions/1/previews/sun_position.m/index.html) for running the extended wooden bucket model.

If you have issues implementing these scripts, or identify any deficiencies, please contact Duo Chan (duochan@g.harvard.edu).

<br>

## Table of Contents
 * [Fast reproduction of Figures](#fast-reproduction-of-figures)
 * [Full Analysis](#full-analysis)
  * [A. Getting Started](#a-getting-started)
  * [B. Diurnal Cycles](#b-diurnal-cycles)
  * [C. LME Intercomparison](#c-lme-intercomparison)
  * [D. Merging DA and LME](#d-merging-da-and-lme)

<br>

## Fast reproduction of Figures

[<span style="color:gray">Back to Table of contents</span>](#table-of-contents)

After cloning this repository, open a Matlab interface and simply run the following lines of codes to reproduce figures in the main text.

```
addpath(genpath(pwd));
DA_LME_main
```

## Full Analysis
[<span style="color:gray">Back to Table of contents</span>](#table-of-contents)

Files are organized following the figure below.
![image](Graphics/Setup.png)

<br>

### A. Getting started
[<span style="color:gray">Back to Table of contents</span>](#table-of-contents)

To get started, run [DA_LME_init.m](DA_LME_init.m), which sets up directories structured following the figure above.  

```
dir_data = 'the directory for outputs';
DA_LME_init(dir_data);
```

After creating directories, you need to add pre-processed ICOAS3.0 data to `$DATA/ICOADS3/ICOADS_QCed/`, where data (30G) can be downloaded from [here](https://dataverse.harvard.edu/file.xhtml?persistentId=doi:10.7910/DVN/DXJIGA/KWDPTS&version=2.0).  Codes for downloading and preprocessing ICOADS3.0 data are published [here](https://github.com/duochanatharvard/Homogeneous_early_20th_century_warming#a-preprocess).

Moreover, the analysis also requires ships tracks from [Carella et al., (2015)](https://rmets.onlinelibrary.wiley.com/doi/abs/10.1002/joc.4492).  Because these data are not yet publicly available, please obtain a permission from Dr. [Elizabeth Kent](mailto:eck@noc.ac.uk) and forward her message to [Duo Chan](mailto:duochan@g.harvard.edu).  Duo will provide `.mat` files upon seeing Kent's permission.

<br>

### B. Diurnal Cycles

[<span style="color:gray">Back to Table of contents</span>](#table-of-contents)

In this step, diurnal anomalies of individual ships are extracted, and then binned and averaged for individual SST groups.  Outputs are diurnal amplitudes of individual bucket SST groups as a function of region-season combinations and 20-year periods.  Output data are placed in `$DATA/DIURNAL/DATA_for_figures/`.

To get started, first run the following lines to extract diurnal anomalies from individual ships and days,
```
P.relative = 'mean_SST';
for yr = 1880:2009
    DIURNAL_Step_01_ship_diurnal_signal(yr,P);
end
```
Computed diurnal anomalies have names `IMMA1_R3.0.0_YYYY-MM_Ship_Diurnal_Signal_relative_to_mean_SST.mat` and are saved in `$DATA/DIURNAL/Step_01_Ship_Signal/`.  Looping over individual years can take quite a while, thus it would be helpful to use multiple CPUs and "parallel" the computation.

Next, run the following lines to subset diurnal anomalies by SST methods and combine them across years.  Note that unlike HadSSTs, we do not treat hull sensor SSTs as ERI measurements.
```
DA_POST_sum_and_fitting(1,[],[],'bucket');
DA_POST_sum_and_fitting(1,[],[],'ERI');
```
The combined diurnal anomalies are saved in `SUM_METHOD_DA_signals_YYYY_YYYY_Annual_relative_to_mean_SST.mat` and placed in `$DATA/DIURNAL/DATA_for_figures/`.

Finally, run the following lines to compute the averaged diurnal cycles for each group over individual 20-year bins and region-season combinations.  The code also fits a sinusoidal function for the amplitude of diurnal cycles.  

<!--Please remember to replace method to `ERI` after running the analysis for bucket SSTs.-->

```
for yr_start = 1880:1990
    yr_end   = yr_start + 19;
    method   = 'bucket';
    DA_POST_sum_and_fitting(0,yr_start,yr_end,method);
end

for yr_start = 1920:1990
    yr_end   = yr_start + 19;
    method   = 'ERI';
    DA_POST_sum_and_fitting(0,yr_start,yr_end,method);
end
```
Output files will have names following `STATS_METHOD_DA_signals_YYYY_YYYY_REGION_relative_to_mean_SST.mat` and are placed in `$DATA/DIURNAL/DATA_for_figures/`.


<br>

###  C. LME Intercomparison

[<span style="color:gray">Back to Table of contents</span>](#table-of-contents)

Intercomparison of collocated SSTs using linear-mixed-effect (LME) model is the same as [here](https://github.com/duochanatharvard/Homogeneous_early_20th_century_warming#b-main-code), but we have improved the algorithm such that the code avoids gigantic matrix inversion and runs faster than the previous version.  The analysis consists of three steps: **1.** pairing SSTs, **2.** screening pairs, and **3.** run the LME model.

To get started, first run the following code to obtain all available pairs from different groups that are within 300km and 2days in time.
```
LME_Step_01_Run_Pairs
```
Computed initial pairs have names `IMMA1_R3.0.0_YYYY-MM_All_pairs.mat` and are saved in `$DATA/LME_intercomparison/Step_01_All_Pairs/`.  This code picks out not only pairs between buckets groups, but also bucket-ERI and ship-buoy pairs.  The code loops over individual years, which takes more than a week to run.  Thus it would be wise to use multiple CPUs and "parallel" the computation.

The second step is to screen pairs such that each measurement is used at most once.  The current code treats all ERI measurements as a single group and uses measurements that have valid diurnal anomaly estimates.  In other words, LME analysis should be performed after the diurnal cycle analysis.
```
LME_Step_02_Screen_Pairs
```
Output files are have names `IMMA1_R3.0.0_YYYY-MM_Bucket_vs_ERI_in_one_group_diurnal_points_mean_SST.mat` and are saved in `$DATA/LME_intercomparison/Step_02_Screen_Pairs/`.

[**IMPORTANT**] Before running the final step, please download `OI_SST_inter_annual_decorrelation_20180316.mat` from [here](https://dataverse.harvard.edu/file.xhtml?persistentId=doi:10.7910/DVN/DXJIGA/TQ8THW&version=2.0) and place it in `$DATA/Miscellaneous/`.  This file is used to compute the uncertainty of physical SSTs between pairs of measurements.  

<!--It is provided separately because of its 284MB size.-->

The final step is to run the LME model and estimate relative offsets between groups of SSTs by running,
```
LME_Step_03_LME_bucket_regional_ERI
```
This step outputs `LME_Bucket_vs_ERI_in_one_group_diurnal_points_mean_SST_YYYY_YYYY_Full_SST_Global.mat` files and saves them in `$DATA/LME_intercomparison/Step_04_LME_output/`.  Intermediate outputs, including the concatenated and aggregated pairs, are output to `$DATA/LME_intercomparison/Step_03_Binned_Pairs/`.


<br>

### D. Merging DA and LME
[<span style="color:gray">Back to Table of contents</span>](#table-of-contents)

Once step A-C are completed, simply run the following script to match outputs from Diurnal and LME analyses and prepare for generating figures and statistics.
```
DA_LME_Prepare_Data
```
The script will output a file named `All_lme_offsets_and_diurnal_amplitudes.mat` and place it in `$DATA/`.  When running `DA_LME_main.m` to generate figures, please copy `All_lme_offsets_and_diurnal_amplitudes.mat` to `$CODE/` such that Matlab can find it.


<br>

Maintained by __Duo Chan__ (duochan@g.harvard.edu)
