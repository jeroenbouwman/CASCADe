
[![pipeline status](https://gitlab.com/jbouwman/CASCADe/badges/master/pipeline.svg)](https://gitlab.com/jbouwman/CASCADe/commits/master)

#  <span style="color:#1F618D">CASCADe </span>: <span style="color:#1F618D">C</span>alibration of tr<span style="color:#1F618D">A</span>nsit <span style="color:#1F618D">S</span>pectroscopy using <span style="color:#1F618D">CA</span>usal <span style="color:#1F618D">D</span>ata

At present several thousand transiting exoplanet systems have been discovered.
For relatively few systems, however, a spectro-photometric characterization of
the planetary atmospheres could be performed due to the tiny photometric
signatures of the atmospheres and the large systematic noise introduced by the
used instruments or the earth atmosphere. Several methods have been developed to deal with instrument and atmospheric noise. These methods include high precision calibration and modeling of the instruments, modeling of the noise using methods like principle component analysis or Gaussian processes and the simultaneous observations of many reference stars. Though significant progress has been made, most of these methods have drawbacks as they either have to make too many assumptions or do not fully utilize all information available in the data to
negate the noise terms.

The <span style="color:#1F618D">CASCADe </span> package, developed within the EC
Horizons 2020 project <span style="color:#FF0000">Exoplanets A </span>, implements
a full spectroscopic pipeline for HST/WFC3 and Spitzer/IRS spectroscopic
timeseries observations as well as lightcurve calibration and fitting functionality.
The <span style="color:#1F618D">CASCADe </span> project implements a novel
“data driven” method, pioneered by Schoelkopf et al (2016) utilizing the causal connections within a data set, and uses this to calibrate the spectral timeseries data of a transiting systems, observed in a single object mode. The current code
has been tested successfully on spectroscopic data obtained with the Spitzer and
HST observatories, as well as JWST MIRI simulations.

## Installing <span style="color:#1F618D">CASCADe</span>

The easiest way to install the <span style="color:#1F618D">CASCADe </span>
package is to create an Anaconda environment, download the distribution from PyPi,
and install the package in the designated Anaconda environment with the following
commands:

```bash

conda create --name cascade python=3.8 ipython
conda activate cascade
pip install CASCADe-spectroscopy

```

This will install all code and scripts you need for the package to work.
> **_NOTE:_**  <span style="color:#1F618D">CASCADe </span>  The batman package is only guaranteed to work when using numpy version 1.22.1, and with this numpy version one should install version 0.53.1 of the numba package.

## Installing the required <span style="color:#1F618D">CASCADe</span> data and examples

All necessary data needed by <span style="color:#1F618D">CASCADe </span> for it to
work properly, such as calibration files for the different spectroscopic instruments
of HST and Spitzer and configuration templates, need be downloaded from the GitLab
repository before running the code. To initialize the data download and setup of
the <span style="color:#1F618D">CASCADe </span> data storage one can use the
following bash command:
```bash

setup_cascade.py

```

or alternatively from within the python interpreter:

```python

from cascade.initialize import initialize_cascade
initialize_cascade()

```

The additional downloaded data also includes examples and observational data to
try out the <span style="color:#1F618D">CASCADe </span> package, which are explained
below.

> **_NOTE:_**  The data files will be downloaded by default to a `CASCADeSTORAGE/` directory in the users home directory. If a different location is preferred, please read the section on how to set the <span style="color:#1F618D">CASCADe </span>
environment variables first.

## Installing alternatives for the <span style="color:#1F618D">CASCADe</span> package

The <span style="color:#1F618D">CASCADe </span> code can also be downloaded from
GitLab directly by either using git or pip. To download and install with a
single command using pip, type in a terminal the following command

```bash

pip install git+git://gitlab.com/jbouwman/CASCADe.git@master

```

which will download the latest version. For other releases replace the `master`
branch with one of the available releases on GitLab. Alternatively, one can first
clone the repository and then install, either using the HTTPS protocal:

```bash

git clone https://gitlab.com/jbouwman/CASCADe.git

```

or clone using SSH:

```bash

git clone git@gitlab.com:jbouwman/CASCADe.git

```

Both commands will download a copy of the files in a folder named after the
project's name. You can then navigate to the directory and start working on it
locally. After accessing the root folder from terminal, type

```bash

pip install .

```

to install the package.

In case one is installing <span style="color:#1F618D">CASCADe </span> directly from
GitLab, and one is using Anaconda,  make sure a cascade environment is created and
activated before using our package. For convenience, in the
<span style="color:#1F618D">CASCADe </span> main package directory an
environment.yml can be found. You can use this yml file to create or update the
cascade Anaconda environment. If you not already had created an cascade environment
execute the following command:

```bash

conda env create -f environment.yml

```

In case you already have an cascade environment, you can update the necessary
packages with the following command (also use this after updating
<span style="color:#1F618D">CASCADe </span> itself):

```bash

conda env update -f environment.yml

```

Make sure the <span style="color:#1F618D">CASCADe </span> package is in your path.
You can either set a `PYTHONPATH` environment variable pointing to the location of the
<span style="color:#1F618D">CASCADe </span> package on your system, or when using
anaconda with the following command:

```bash

conda develop <path_to_the_CASCADe_package_on_your_system>/CASCADe

```

## Using  <span style="color:#1F618D">CASCADe </span>

The <span style="color:#1F618D">CASCADe </span> distribution comes with a few
working examples and data sets which can be found in the examples directory of
the <span style="color:#1F618D">CASCADe </span> distribution on GitLab, and which
should have been copied to the storage directory specified by the
`CASCADE_STORAGE_PATH` environment variable. All needed parameters for running
the code are defined by  initialization files and a few environment variables to
be set by the user.

### <span style="color:#1F618D">CASCADe </span> Environment Variables

<span style="color:#1F618D">CASCADe </span> uses the following environment variables:

```

    CASCADE_STORAGE_PATH
        Default path to all functional data needed by the CASCADe package.
    CASCADE_DATA_PATH
        Default path to all observational data to be analyzed with CASCADe.
    CASCADE_SAVE_PATH
        Default path to where CASCADe saves all pipeline output.
    CASCADE_INITIALIZATION_FILE_PATH:
        Default path to CASCADe pipeline initialization files.
    CASCADE_SCRIPTS_PATH
        Default path to the CASCADe data reduction pipeline scripts for each
        observation.
    CASCADE_LOG_PATH:
        Default path to the saved CASCADe log files.
    CASCADE_WARNINGS
        Switch to show or not show warnings. Can either be 'on' or 'off'

```

These environment variables control where <span style="color:#1F618D">CASCADe </span>
searches and stores all required data. In case the environment variables are not
set by the user, <span style="color:#1F618D">CASCADe </span> uses default values
defined in the `initialize` module. The main environment variable is
`CASCADE_STORAGE_PATH`. This environment variable, if not set by the user has a
default value of `<user_home_directory>/CASCADeSTORAGE/`. All other path settings
stored in the other environment variables are set relative to this, so in principle
the user only has to set the `CASCADE_STORAGE_PATH` variable to the preferred location,
such as a directory on a larger file system rather then in the user home directory.
The `CASCADE_DATA_PATH`, `CASCADE_SAVE_PATH`, `CASCADE_INITIALIZATION_FILE_PATH`,
`CASCADE_SCRIPTS_PATH` and `CASCADE_LOG_PATH` are set by default to
`data/`, `results/`, `init_files/`, `scripts/` and `logs/`, respectively, all
relative to the path defined by `CASCADE_STORAGE_PATH`.
The `CASCADE_WARNINGS` variable is by default set to `"on"` to show possible
warnings generated by the <span style="color:#1F618D">CASCADe </span> code.


### <span style="color:#1F618D">CASCADe </span> Script Examples

The distribution comes with three examples demonstrating the use of the
<span style="color:#1F618D">CASCADe </span> package. These use cases, together with
the observational data can be found in the 'examples/' sub-directory of the GitLab
repository, and should be installed (see above) in the directory defined by the
`CASCADE_STORAGE_PATH` and other environment variables. The examples cover both HST and
Spitzer data as well as an example for a Generic spectroscopic dataset, in this
case an observation with the GMOS instrument on Gimini.

#### Example 1a: Extracting a Spectroscopic Timeseries from HST WFC3 spectral images.

In this example we demonstrate how extract a spectral timeseries from spectral images.
For this we use observations of WASP-19b with the Wide Field Camera 3 (WFC3) on
board the HST observatory. As these observations were made in `"Staring Mode"`,
we use as a start the `flt` data product. The spectral images can be found in the
`HST/WFC3/WASP-19b_ibh715/SPECTRAL_IMAGES` sub-directory of the main data directory
specified by the `CASCADE_DATA_PATH` environment variable.

The pipeline script for this example can be found in the `HST/WFC3/` sub-directory
of the scripts directory as specified by the `CASCADE_SCRIPTS_PATH` environment
variable. To run this example, execute the following script:

```bash

python3 run_CASCADe_WASP19b_extract_timeseries.py

```

The individual steps in this script are commented and a detailed explanation on the
reduction steps can be found in the <span style="color:#1F618D">CASCADe </span>
documentation (see below). Note that while running the script, interactive plots are
opened which the user needs to close before the script continues. One can prevent plots
from opening by un-commenting the line `matplotlib.use('AGG')` in the python script.

The initialization files for this example can be found in the
`HST/WFC3/WASP-19b_ibh715_example` sub-directory of the directory specified by the
`CASCADE_INITIALIZATION_FILE_PATH` environment variable.
The `cascade_WASP19b_object.ini` contains all parameters specifying the target star
and planet. The ephemeris and orbital period given in the `.ini` file are used
to calculated an orbital phase which is attached to the individual spectra in the
final timeseries. In case of HST WFC3 spectra, as is the case here, the stellar
parameters are used to calculated an expected spectrum, which is then used to
determine a global shift of the spectrum in the wavelength direction. It is,
therefore, important that the `effective temperature`, `logg` and `metallicity`
parameters are close as possible to the correct values of the system which the
observations are being analyzed.  All other parameters controlling the behavior of
<span style="color:#1F618D">CASCADe </span> are given in the
`cascade_WASP19b_extract_timeseries.ini` initialization file. The parameters are
grouped different sections in a logical way, such as parameters controlling the
processing steps, describing the observations and data and overall behavior of
the code. Detailed information on each of these parameters can be found in the
documentation (see below.)

Executing the script results in the extraction of the individual (1d) spectra.
These are stored as fits files in the sub-directory `SPECTRA` at the same level
as the `SPECTRAL_IMAGES` directory containing the spectral images. Two types of
files are written: `COE`, which are the <span style="color:#1F618D">CASCADe </span>
`Optimal Extraction` spectra and `CAE`, which are the
<span style="color:#1F618D">CASCADe </span> `Aperture Extraction` spectra, produced
using respectively, optimal extraction or an extraction aperture. For further details
we refer to the documentation.  Next to the spectral fits files, also several
diagnostic plots are produced which can be found in the
`HST/WFC3/WASP-19b_ibh715_transit_output_from_extract_timeseries` sub-directory
of the directory specified by the `CASCADE_SAVE_PATH` environment variable.


#### Example 1b: Calibrating a HST WFC3 Spectroscopic Timeseries and extracting a transit spectrum.

After creating a spectral timeseries from the spectral images in Example 1a,
we can proceed with the calibration of the spectral lightcurves and deriving the
planetary spectrum. To demonstrate how to use <span style="color:#1F618D">CASCADe </span>
for this, we will take the extracted `COE` spectra of the WASP 19b observation
from Example 1a and proceed to characterize the systematics and extract the transit
spectrum. To run Example 1b, execute the following script:

```bash

python3 run_CASCADe_WASP19b_calibrate_planet_spectrum.py

```

As with the first example, the individual steps in this script are commented and
a detailed explanation on the reduction steps can be found in the
<span style="color:#1F618D">CASCADe </span> documentation (see below).
The initialization files for this example can be found in the same directory as the
initialization files for example 1a. We use again the `cascade_WASP19b_object.ini` file
to specify the target star and orbital parameters. In contrast to the previous example,
all system parameters specified in this initialization file are relevant.
> **_NOTE:_** The current <span style="color:#1F618D">CASCADe </span> version only fits
for the transit depth. All other system parameters such as the ephemeris, period,
inclination and semi-major axis are fixed to the values specified in the initialization
file.

All other parameters controlling the behavior of
the <span style="color:#1F618D">CASCADe </span> pipline  are given in the
`cascade_WASP19b_calibrate_planet_spectrum.ini` initialization file. A detailed
explanation of the control parameters is given in the documentation.

The resulting transit spectrum and diagnostic plots are stored in the
`HST/WFC3/WASP-19b_ibh715_transit_from_hst_wfc3_spectra` sub-directory of the
directory specified by `CASCADE_SAVE_PATH` environment variable. The calibrated spectrum
of WASP 19b is stored in the `WASP-19b_ibh715_bootstrapped_exoplanet_spectrum.fits` and
the derived systematics for this dataset in `WASP-19b_ibh715_bootstrapped_systematics_model.fits`

#### Example 2: Calibrating a Spitzer/IRS Spectroscopic Timeseries and extracting a transit spectrum.

The <span style="color:#1F618D">CASCADe </span> package can not only calibrate
observations with the WFC3 instrument onboard HST, but can also handle transit
spectroscopy observations with the IRS instrument onboard the Spitzer Space
Observatory. As an example, we analyze Spitzer/IRS observations of an eclipse
of HD189733b, using the with the <span style="color:#1F618D">CASCADe </span>
package pre-extracted `COE` spectral data product. The data can be found
in the `SPITZER/IRS/HD189733b_AOR23439616/SPECTRA/` sub-directory of the main
data directory specified by the `CASCADE_DATA_PATH` environment variable.
The pipeline script for this example can be found in the `SPITZER/IRS/` sub-directory
of the scripts directory as specified by the `CASCADE_SCRIPTS_PATH` environment
variable. To run Example 2, execute the following script:

```bash

python3 run_CASCADe_HD189733b_calibrate_planet_spectrum.py

```

The pipeline steps used in this example are identical to the ones of Example 1b.
The initialization files for example 2 can be found in the
`SPITZER/IRS/HD189733b_AOR23439616_example` sub-directory of the directory
specified by the `CASCADE_INITIALIZATION_FILE_PATH` environment variable.
Similar to the first example, the `cascade_HD189733b_object.ini` file contains
all parameters specifying the target star and orbital parameters, while the
`cascade_HD189733b_calibrate_planet_spectrum.ini` initialization file specifies
all other parameters controlling the behavior of the
<span style="color:#1F618D">CASCADe </span> pipeline. The HD189733b eclipse
spectrum and diagnostic plots are stored in the
`SPITZER/IRS/HD189733b_AOR23439616_eclipse_from_spitzer_irs_spectra` sub-directory
of the directory specified by the `CASCADE_SAVE_PATH` environment variable.

#### Example 3: Calibrating a GIMINI/GMOS Spectroscopic Timeseries and extracting a transit spectrum.

As a final example we show how to use <span style="color:#1F618D">CASCADe </span>
for spectral timeseries extracted with another software package for a generic
instrument. Though spectral extraction from spectral images or cubes is currently only
implemented for the WFC3 instrument of HST and the IRS instrument of Spitzer,
the calibration of spectral lightcurves and derivation of the planetary spectrum
can be performed for any generic spectroscopic timeseries. The previous examples
showed how to use <span style="color:#1F618D">CASCADe </span> with HST and Spitzer
observations. In this example we use an observation of WASP-103b with the GMOS
instrument installed at the Gemini telescope (See Lendl et al 2017, A&A 606).

The spectral timeseries data for this example is located in the
`Generic/Gemini/GMOS/WASP103b/SPECTRA/` sub-directory of the main data directory.
To be able to run this example we stored the GMOS spectra as fits files with an identical format as the spectral fits files created by <span style="color:#1F618D">CASCADe </span>.
The pipeline script is located in the `Generic/Gemini/GMOS/` sub-directory in
the scripts directory.

To run this example, execute the following script:

```bash

python3 run_CASCADe_WASP103b_calibrate_planet_spectrum.py

```

The initialization files for example 2 can be found in the
`Generic/Gemini/GMOS/WASP-103b_example/` sub-directory of the directory specified
by the `CASCADE_INITIALIZATION_FILE_PATH` environment variable. Similar to
the other examples, the `cascade_WASP103b_object.ini` initialization file contains
all parameters defining the system, and the `cascade_WASP103b_calibrate_planet_spectrum.ini` file contains all other
parameters needed by the <span style="color:#1F618D">CASCADe </span> pipeline.
The WASP-103 b transit spectrum and diagnostic plots are stored in the
`Generic/Gemini/GMOS/WASP103b_transit_from_generic_instrument/`
sub-directory of the directory specified by the `CASCADE_SAVE_PATH` environment
variable.

## Documentation

The full documentation can be found online at:

```

https://jbouwman.gitlab.io/CASCADe/

```

The full documentation includes further descriptions of the pipeline,
initialization files and the <span style="color:#1F618D">CASCADe </span> modules.

Alternatively, the documentation can be found in the `docs`  directory of the
<span style="color:#1F618D">CASCADe </span> GitLab repository.
After cloning the git repository, the full documentation can be generated
by executing in the in the `docs` directory the following commands:

```bash

make html
make latexpdf

```

The generated HTML and PDF files will be located in the `build/html` and
`build/latex` sub-directories of the main documentation directory, respectively.

## Acknowledgments

The <span style="color:#1F618D">CASCADe </span> code was developed by
Jeroen Bouwman, with contributions from the following collaborators:

Fred Lahuis (SRON)\
Rene Gastaus (CEA)\
Raphael Peralta (CEA)\
Matthias Samland (MPIA)

This work was supported by the European Unions Horizon 2020 Research and
Innovation Programme, under Grant Agreement N 776403.

## Publications

https://ui.adsabs.harvard.edu/abs/2021AJ....161..284M/abstract

https://ui.adsabs.harvard.edu/abs/2021A%26A...646A.168C/abstract

https://ui.adsabs.harvard.edu/abs/2020ASPC..527..179L/abstract

https://exoplanet-talks.org/talk/271
