
[![pipeline status](https://gitlab.com/jbouwman/CASCADe/badges/master/pipeline.svg)](https://gitlab.com/jbouwman/CASCADe/commits/master)

#  <span style="color:#1F618D">CASCADe </span>: Calibration of trAnsit Spectroscopy using CAusal Data

At present several thousand transiting exoplanet systems have been discovered.
For relatively few systems, however, a spectro-photometric characterization of
the planetary atmospheres could be performed due to the tiny photometric signatures
of the atmospheres and the large systematic noise introduced by the used instruments
or the earth atmosphere. Several methods have been developed to deal with instrument
and atmospheric noise. These methods include high precision calibration and modeling
of the instruments, modeling of the noise using methods like principle component
analysis or Gaussian processes and the simultaneous observations of many reference
stars. Though significant progress has been made, most of these methods have drawbacks
as they either have to make too many assumptions or do not fully utilize all
information available in the data to negate the noise terms.

The <span style="color:#1F618D">CASCADe </span> package, developed within the EC
Horizons 2020 project <span style="color:#FF0000">Exoplanets A </span>, implements
a full spectroscopic pipeline for HST/WFC3 and Spitzer/IRS spectroscopic
timeseries observations as well as lightcurve calibration and fitting functionality.
The <span style="color:#1F618D">CASCADe </span> project implements a novel “data driven”
method, pioneered by Schoelkopf et al (2016) utilizing the causal connections within a
data set, and uses this to calibrate the spectral timeseries data of single transiting
systems. The current code has been tested successfully on spectroscopic data
obtained with the Spitzer and HST observatories.

## Installing Cascade

The <span style="color:#1F618D">CASCADe </span> code can be downloaded from gitlab either
using git or pip. To download and install with a single command using pip, type in a
terminal the following command

```bash

pip install git+git://gitlab.com/jbouwman/CASCADe.git@master

```

Alternatively, one can first clone the repository and then install, either using the
HTTPS protocal:

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

In case one is using Anaconda make sure a cascade environment is created and
activated before using our package In the <span style="color:#1F618D">CASCADe </span>
main package directory an environment.yml can be found. You can
use this yml file to create or update the cascade Anaconda environment. If you
not already had created an cascade environment execute the following command:

```bash

conda env create -f environment.yml

```

Incase you already have an cascade environment, you can update the necessary
packages with the following command (also use this after updating
<span style="color:#1F618D">CASCADe </span> itself):

```bash

conda env update -f environment.yml

```

Make sure the <span style="color:#1F618D">CASCADe </span> package is in your path.
You can either set a PYTHONPATH environment variable pointing to the location of the
<span style="color:#1F618D">CASCADe </span> package on your system, or when using
anaconda with the following command:

```bash

conda develop <path_to_the_CASCADe_package_on_your_system>/CASCADe

```

## Using  <span style="color:#1F618D">CASCADe </span>

The CASCADe distribution comes with a few working examples and data sets which can
be found in the examples and data sub-directories, respectively. All needed parameters
for running the code are defined in the provided initialization files and a few environment
variables. The user only has to change those, particularly the path settings to the data
and .ini files to run the examples.

### <span style="color:#1F618D">CASCADe </span> Environment Variables

<span style="color:#1F618D">CASCADe </span> uses the following environment variables:

```

    CASCADE_WARNINGS
        Switch to show or not show warnings. Can either be ‘on’ or ‘off’
    CASCADE_DEFAULT_PATH
        Default directory for CASCADe.
    CASCADE_DATA_PATH
        Default path to the input observational data and calibration files.
    CASCADE_SAVE_PATH
        Default path to where CASCADe saves output.
    CASCADE_INITIALIZATION_FILE_PATH:
        Default directory for CASCADe initialization files.
    CASCADE_LOG_PATH
        Default path to were log files are stored.    

```

in case the environment variables are not set by the user,
<span style="color:#1F618D">CASCADe </span> uses default values defined in the
initialize module. The `CASCADE_WARNINGS` variable is by default set to `"on"`.
The `CASCADE_DEFAULT_PATH` should be set to the installation directory of the
<span style="color:#1F618D">CASCADe </span> code. Per default, all other environment
variables, defining the location of the data, were results need to be stored,
were the initialization files are located and were log files are written, are set
relative to the `CASCADE_DEFAULT_PATH` environment variable by default to
`"data/"`,  `"examples/results/"`, `"examples/init_files/"` and  `"examples/logs/"`,
respectively. With these standard values the example scripts provided in the
<span style="color:#1F618D">CASCADe </span> distribution can be executed without problems.
However, the user is advised to set these variables to values other than the standard
values as not to clutter the <span style="color:#1F618D">CASCADe </span> source
directory. Especially for large data sets, changing the default to a directory on a
large storage disk is recommended.

### <span style="color:#1F618D">CASCADe </span> Script Examples

The distribution comes with three examples demonstrating the use of the
<span style="color:#1F618D">CASCADe </span> package. These use cases can be found
in the 'examples/' sub-directory. They cover both HST and Spitzer data as well as
an example for a Generic spectroscopic dataset.

#### Example 1a: Extracting a Spectroscopic Timeseries from HST WFC3 spectral images.

In this example we demonstrate how extract a spectral timeseries from spectral images.
For this we use observations of WASP-19b with the Wide Field Camera 3 (WFC3) on board the
HST observatory. As these observations were made in `"Staring Mode"`, we use as a start the `flt`
data product. The spectral images can be found in the data sub-directory under
`data/HST/WFC3/SPECTRAL_IMAGES`. If the user has changed the default settings of the
environment variables, these data will need to be copied to the corresponding user directory.

To run this example, execute the following script:

```bash

python3 run_CASCADe_WASP19b_extract_timeseries.py

```

The individual steps in this script are commented and a detailed explanation on the
reduction steps can be found in the <span style="color:#1F618D">CASCADe </span> documentation
(see below). The initialization files for this example can be found in the `examples/init_files`
sub-directory The `cascade_WASP19b_object.ini` contains all parameters specifying the target star
and planet. For this example, the ephemeris and orbital period given in the `.ini` file are used
to calculated an orbital phase which is attached to the individual spectra in the final timeseries,
though this has no impact on the spectra themselfs and is only used as a additional
service to the user. In case of HST WFC3 spectra, as is the case here, the stellar parameters
are used to calculated an expected spectrum, which is then used to determine a global shift of the
spectrum in the wavelength direction. It is,  therefore, important that these parameters
(`effective temperature`, `logg` and `metallicity`) are close as possible to the correct
values of the system which the observations are being analyzed.  All other parameters controlling
the behavior of <span style="color:#1F618D">CASCADe </span> are given in the
`cascade_WASP19b_extract_timeseries.ini` initialization file.  The parameters are grouped
different sections in a logical way, such as parameters controlling the processing steps, describing the
observations and data and overall behavior of the code. Detailed information on each of these parameters
can be found in the documentation (see below.)

Executing the script results in the extraction of the individual (1d) spectra. These are stored
in the sub-directory 'SPECTRA' as fits files at the same level as the SPECTRAL_IMAGES directory containing the
spectral images. Two types of files are written: `COE`, which are the <span style="color:#1F618D">CASCADe </span>
optimal extraction results and `CAE`, which are the <span style="color:#1F618D">CASCADe </span>
aperture extraction results, produced using respectively an extraction profile or an extraction aperture.
For further details we refer to the documentation.  Next to the spectral fits files,
also several diagnostic plots are produced which can be found in the
`WASP-19b_ibh715_transit_output_from_extract_timeseries` sub-directory of the `CASCADE_SAVE_PATH`
directory which, if not set by the user, is the `CASCADe\examples\results` directory of the <span style="color:#1F618D">CASCADe </span> distribution.

#### Example 1b: Calibrating a HST WFC3 Spectroscopic Timeseries and extracting a transit spectrum.

After extracting the spectral timeseries from spectral data cubes or images as shown in Example 1a,
we can proceed with the extraction of the planetary spectrum and the calibration of the spectral lightcurves.
To demonstrate how to use <span style="color:#1F618D">CASCADe </span> to do this, we will take the extracted spectra of the WASP 19b observation from Example 1a and proceed to characterize the systematics and extract the transit spectrum.
To run Example 1b, execute the following script:

```bash

python3 run_CASCADe_WASP19b_calibrate_planet_spectrum.py

```

The individual steps in this script are commented and a detailed explanation on the
reduction steps can be found in the <span style="color:#1F618D">CASCADe </span> documentation
(see below). The initialization files for this example can be found in the `examples/init_files`
sub-directory. The `cascade_WASP19b_object.ini` file contains all parameters specifying the target star
and planet. In contrast to the previous example, all system parameters specified in
this initialization file are relevant. Note, that the current <span style="color:#1F618D">CASCADe </span>
version only fits for the transit depth. All other system parameters such as the ephemeris, period, inclination and
semi-major axis are fixed to the values specified in the `.ini` file. All other parameters controlling
the behavior of <span style="color:#1F618D">CASCADe </span> are given in the
`cascade_WASP19b_calibrate_planet_spectrum.ini` initialization file.  A detailed explanation of
the control parameters is given in the documentation.

The resulting transit spectrum and diagnostic plots are stored in the  `WASP-19b_ibh715_transit_from_hst_wfc3_spectra`
sub-directory in the `CASCADE_SAVE_PATH` directory which, if not set by the user, is the `CASCADe\examples\results`
directory of the <span style="color:#1F618D">CASCADe </span> distribution.

#### Example 2: Calibrating a Spitzer/IRS Spectroscopic Timeseries and extracting a transit spectrum.

Apart from the observations with the WFC3 instrument onboard HST, <span style="color:#1F618D">CASCADe </span>
can also handle transit spectroscopy observations with the IRS instrument onboard the Spitzer Space Observatory.
As an example, we will use a Spitzer/IRS observation  of an eclipse of HD189733b.

To run Example 2, execute the following script:

```bash

python3 run_CASCADe_HD189733b_calibrate_planet_spectrum.py

```

The pipeline steps used in this example  are identical to the ones of Example 1b. Again, the the  
initialization files for this example can be found in the `examples/init_files`
sub-directory. The `cascade_HD189733b_object.ini` file contains all parameters specifying the target star
and planet, while the `cascade_HD189733b_calibrate_planet_spectrum.ini` initialization file specifies
all other parameters controlling the behavior of <span style="color:#1F618D">CASCADe </span>.

The HD189733b eclipse spectrum and diagnostic plots are stored in the `HD189733b_AOR23439616_eclipse_from_spitzer_irs_spectra`
sub-directory in the `CASCADE_SAVE_PATH` directory which, if not set by the user, is the `CASCADe\examples\results`
directory of the <span style="color:#1F618D">CASCADe </span> distribution.

#### Example 3: Calibrating a GIMINI/GMOS Spectroscopic Timeseries and extracting a transit spectrum.

As a final example we show how to use <span style="color:#1F618D">CASCADe </span> for spectral timeseries
which have been extracted with another software package for a generic instrument. Though
spectral extraction from spectral images or cubes is currently only implemented for HST/WFC3 and SpitzerIRS data,
the extraction of the planetary signal and the calibration of the spectral lightcurves can be performed for any generic spectroscopic timeseries.
The previous examples showed how to use <span style="color:#1F618D">CASCADe </span> with HST and Spitzer observations. In this example we use an observation with the GIMINI/GMOS instrument of WASP-103b (See Lendl et al 2017, A&A 606).

To run this example, execute the following script:

```bash

python3 run_CASCADe_WASP103b_calibrate_planet_spectrum.py

```
To be able to run this example we stored the GMOS spectra as fits files with an identical format as the fits files
created by <span style="color:#1F618D">CASCADe </span> to store the extracted spectra.  

The initialization files for this example are `cascade_WASP103b_object.ini`, which contains the system parameters and
`cascade_WASP103b_calibrate_planet_spectrum.ini`, which contains all other necessary parameters.
The WASP-103 b transit spectrum and diagnostic plots are stored in the `WASP103b_transit_from_generic_instrument`
sub-directory in the `CASCADE_SAVE_PATH` directory which, if not set by the user, is the `CASCADe\examples\results`
directory of the <span style="color:#1F618D">CASCADe </span> distribution.

## Documentation

Extended documentation can be found in the `docs`  directory.
To generate the documentation, make sure you have 'Sphinx' installed and than execute
in the in the `docs` directory the following command:

```bash

make html
make latexpdf

```

The generated HTML and PDF will be located in the location specified by `conf.py`,
in this case `build/html` and `build/latex`. Alternatively one can get the latest
documentation from:

```
https://jbouwman.gitlab.io/CASCADe/
```


## Acknowledgments

The <span style="color:#1F618D">CASCADe </span> code was developed by Jeroen Bouwman, with contributions
from the following collaborators:

Fred Lahuis (SRON)\
Rene Gastaus (CEA)\
Raphael Peralta (CEA)\
Matthias Samland (MPIA)

This work was supported by the European Unions Horizon 2020 Research and Innovation Programme,
under Grant Agreement N 776403.

## Publications

https://ui.adsabs.harvard.edu/abs/2021AJ....161..284M/abstract

https://ui.adsabs.harvard.edu/abs/2021A%26A...646A.168C/abstract

https://ui.adsabs.harvard.edu/abs/2020ASPC..527..179L/abstract

https://exoplanet-talks.org/talk/271
