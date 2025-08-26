
# EnKF_CH4
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16947849.svg)](https://doi.org/10.5281/zenodo.16947849)


## Description

This repository is a global methane emission assimilation system based on an Ensemble Kalman Filter (EnKF) framework and the GEOS-Chem Global Chemical Transport Model (v12.5.0). It is an updated version of the ESA PyOSSE: Package for Observation System Simulation Experiments developed by [Liang Feng](https://www.geos.ed.ac.uk/~lfeng/). We have converted the original FORTRAN modules to Python scripts for better use and easier update.

Included in this repository are:
  - Tag run program to create the Jacobian matrix needed for the EnKF;
  - Zip file for GEOS-Chem v12.5.0 and user-updated global_ch4_mod.F for an easier way to define the tag run;
  - Scripts to create inversion run directory, GEOS-Chem rerun directory and diagnosis directories;
  - Configuration files that specify inversion options and diagnosis options;
  - Scripts to run GEOS-Chem tests;
  - Scripts for CH4 observations, including GOSAT, TROPOMI, NOAA, TCCON;
  - Scripts to correct GEOS-Chem CH4 latitudinal biases during inversion; and
  - Main Driver routines.


## Acknowlegdement
We acknowledge the GEOS-Chem community, in particular the Harvard University team that helps maintain the GEOS-Chem model, and the NASA Global Modelling and Assimilation Office (GMAO) for providing the MERRA2 data product. 

### Links to GEOS-Chem
  - [http://wiki.seas.harvard.edu/geos-chem/index.php/Main_Page/](http://wiki.seas.harvard.edu/geos-chem/index.php/Main_Page)
  - [https://github.com/geoschem/](https://geoschem.github.io)
  - [https://github.com/geoschem/GEOSChem-python-tutorial](https://github.com/geoschem/GEOSChem-python-tutorial)


## Table of Contents

- [Configuration](#Configuration)
- [Framework](#Framework)
- [InputData](#InputData)
- [Examples](#Examples)
- [References](#References)
- [License](#license)

## Configuration

#### Python setup
```bash
# Python environment installation 
$ conda env create -vv -n idp -f enkf_environment.yml
$ source activate idp
```

#### Tag run
```bash
# configuration
$ cd  .../global_tagrun/code/
$ vim geos_chem_def.py
# tag run 
$ ./enkf_drive_sh.py
```

#### Inversion run

```bash
# configuration
$ vim /obs/operator_config.py ## Observation options
$ cd  .../enkf/
$ vim geos_chem_def.py  ## inversion options and diagnosis path
$ vim restart_config.py 
# inversion run 
$ ./etkf_main.py
```

## Framework
 
<img src="https://github.com/Rainbow1994/EnKF_CH4/blob/master/images/enkf-ch4-framework.png" width="500">

## InputData

Data | Usage  | Official Access
| :--- | :--- |:---
NOAA ObsPack products | In-situ CH4 recordings for inversion| [Global Monitoring Laboratory (noaa.gov)](https://gml.noaa.gov/ccgg/obspack/release_notes.html#obspack_ch4_1_GLOBALVIEWplus_v5.0_2022-10-17)
GOSAT Proxy XCH4 data(v9.0)| GOSAT XCH4 retrievals for inversion | [University of Leicester GOSAT Proxy XCH4 v9.0](https://catalogue.ceda.ac.uk/uuid/18ef8247f52a4cb6a14013f8235cc1eb)
ACE-FTS CH4 profiles(v4.1)| CH4 profiles for GEOS-chem latitudinal bias correction| [Data access](https://databace.scisat.ca/) and [Data quality flags](https://borealisdata.ca/dataset.xhtml?persistentId=doi:10.5683/SP2/BC4ATC)
TCCON data| In-situ XCH4 measurements for validation | [TCCON Data Archive](https://tccondata.org/)

_We appreciate all of the scientists and professionals who contributed to the datasets listed above._

## Examples

#### CH4 concentration
The decadal mean difference between GOSAT-retrieved methane column concentrations (XCH4) and those simulated using GEOS-Chem with a priori emissions at 4x5 (R4, a) and 2x2.5 (R2, b) scales and using a posteriori emissions after inversion at grid scales of R4 (c) and R2 (d).

<img src="https://github.com/Rainbow1994/EnKF_CH4/blob/master/images/GOSAT_comparison.jpg" width="500">

#### CH4 flux
Decadal mean distribution of a priori and a posteriori methane emissions in using R4 (a, d) and the R2 (b, e) inversions and their differences (a posteriori minus a priori) (c, f).

<img src="https://github.com/Rainbow1994/EnKF_CH4/blob/master/images/flux.jpg" width="500">

#### Global total emissions
Annual mean variations of global total methane emissions (a) in R4 (blue) and R2 (orange) versions of the GEOS-Chem model and their monthly variations (b) from 2010 to 2019.

<img src="https://github.com/Rainbow1994/EnKF_CH4/blob/master/images/total_emission.jpg" width="500">

Global annual total emissions during the 2010s (Tg/yr).

<img src="https://github.com/Rainbow1994/EnKF_CH4/blob/master/images/Table1.png" width="500">

#### Validation
Taylor diagrams of statistical results (correlation coefficient, standard deviation, and root-mean-square deviation (RMSD)) between surface-measured methane column concentrations (XCH4) from the TCCON network and those simulated using GEOS-Chem with a priori emissions at R4 (a) and R2 (b), and using a posteriori emissions after inversion at R4 (c) and R2 (d) (Deep blue: latitudes of TCCON sites are larger than 60°N; Green: the latitudes of sites are within 45° – 60°N; Yellow: the latitudes of sites are within 30° – 45°N; Red: the latitudes of sites are within –15°S – 45°N; Blue: The sites located in the mid-latitudes of SH).

<img src="https://github.com/Rainbow1994/EnKF_CH4/blob/master/images/tccon.jpg" width="500">

## References
[1] Feng, L., Palmer, P.I., Bösch, H. and Dance, S., 2009. Estimating surface CO 2 fluxes from space-borne CO 2 dry air mole fraction observations using an ensemble Kalman Filter. Atmospheric chemistry and physics, 9(8), 2619−2633, https://doi.org/10.5194/acp-9-2619-2009.

[2] Feng, L., Palmer, P.I., Bösch, H., Parker, R.J., Webb, A.J., Correia, C.S., Deutscher, N.M., Domingues, L.G., Feist, D.G., Gatti, L.V. and Gloor, E., 2017. Consistent regional fluxes of CH 4 and CO 2 inferred from GOSAT proxy XCH 4: XCO 2 retrievals, 2010–2014. Atmospheric chemistry and physics, 17(7), 4781−4797, https://doi.org/10.5194/acp-17-4781-2017.

[3] Zhu, S., Feng, L., Liu, Y., Wang, J. and Yang, D., 2022. Decadal Methane Emission Trend Inferred from Proxy GOSAT XCH4 Retrievals: Impacts of Transport Model Spatial Resolution. Advances in Atmospheric Sciences, 39(8), 1343-1359, https://doi.org/10.1007/s00376-022-1434-6.

## License
This project is licensed under the [MIT License](https://github.com/Rainbow1994/EnKF_CH4/blob/master/MIT%20License).






