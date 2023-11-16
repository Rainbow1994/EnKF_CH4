
# EnKF_CH4

## Description

This repository is a global methane emission assimilation system based on an Ensemble Kalman Filter (EnKF) framework and the GEOS-Chem transport model (v12.5.0). It is an updated version of the ESA PyOSSE: Package for Observation System Simulation Experiments developed by [Liang Feng](https://www.geos.ed.ac.uk/~lfeng/) We have converted the original FORTRAN modules to Python scripts for better use and easier update.

Included in this repository are:
  - Tag run program to create Jacobi matrix;
  - Zip file for GEOS-Chem v12.5.0 and user-updated global_CH4.f for easier Tag run;
  - Scripts to create inversion run directory, GEOS-Chem rerun directory and diagnosis directories;
  - Configuration files that specify inversion options and diagnosis options;
  - Scripts to run GEOS-Chem tests;
  - Scripts for CH4 observations, including GOSAT, TROPOMI, NOAA, TCCON;
  - Scripts to correct GEOS-Chem CH4 latitudinal biases during inversion;
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

Python setup
```bash
# Python environment installation 
$ conda env create -vv -n enkf -f enkf_environment.yml
$ source activate enkf
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


## References
[1] Feng, L., Palmer, P.I., Bösch, H. and Dance, S., 2009. Estimating surface CO 2 fluxes from space-borne CO 2 dry air mole fraction observations using an ensemble Kalman Filter. Atmospheric chemistry and physics, 9(8), 2619−2633, https://doi.org/10.5194/acp-9-2619-2009.

[2] Feng, L., Palmer, P.I., Bösch, H., Parker, R.J., Webb, A.J., Correia, C.S., Deutscher, N.M., Domingues, L.G., Feist, D.G., Gatti, L.V. and Gloor, E., 2017. Consistent regional fluxes of CH 4 and CO 2 inferred from GOSAT proxy XCH 4: XCO 2 retrievals, 2010–2014. Atmospheric chemistry and physics, 17(7), 4781−4797, https://doi.org/10.5194/acp-17-4781-2017.

[3] Zhu, S., Feng, L., Liu, Y., Wang, J. and Yang, D., 2022. Decadal Methane Emission Trend Inferred from Proxy GOSAT XCH4 Retrievals: Impacts of Transport Model Spatial Resolution. Advances in Atmospheric Sciences, 39(8), 1343-1359, https://doi.org/10.1007/s00376-022-1434-6.

## License







