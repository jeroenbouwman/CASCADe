
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

The CASCADe distribution comes with a few working examples and datasets which can
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
However, the user is adviced to set these variables to values other than the standard 
values as not to clutter the <span style="color:#1F618D">CASCADe </span> source
directory. Especially for large data sets, changing the default to a directory on a
large storage disk is recommended. 

### <span style="color:#1F618D">CASCADe </span> Script Examples

The distribution comes with three examples demonstrating the use of the
<span style="color:#1F618D">CASCADe </span> package. These use cases can be found
in the 'examples/' subdirectory. They cover both HST and Spitzer data as well as
an example for a Generic spectroscopic dataset.

#### Example1a: Extracting a Spectroscopic Timeseries from HST WFC3 spectral images.

In this example we demonstrate how extract a spectral timeseries from spectral images.
For this we use observations of WASP-19b with the Wide Field Camera 3 (WFC3) on board the
HST observatory. As these observations were made in `"Staring Mode"`, we use as a start the `flt`
data product. The spectral images can be found in the data subdirectory under
`data/HST/WFC3/SPECTRAL_IMAGES`. If the user has changed the default settings of the
environment variables, these data will need to be copied to the corresponding user directory.

To run this example, execute the following script:

```bash

python3 run_CASCADe_WASP19b_extract_timeseries.py

```

The individual steps in this script are commentend and a detailed explanation on the 
reduction steps can be found in the <span style="color:#1F618D">CASCADe </span> documentation
(see below). The initialization files for this example can be found in the `examples/init_files`
subdirectory The `cascade_WASP19b_object.ini` contains all parameters specifying the target star
and planet. For this example, the ephemeris and orbital period given in the .ini file are used
to calculated an orbital phase which is attached to the individual spactra in the final timeseries,
though this has no impact on the spectra themselfs and is only used as a additional
servvice to the user. In case of HST WFC3 spectra, as is the case here, the stellar parameters 
are used to calculated an expected spectrum, which is then used to determine a global shift of the
spectrum in the wavelength direction. It is,  therefore, important that these parameters
(`effective temperature`, `logg` and `metallicity`) are close as posible to the correct
values of the system which the observations are beeing analysed.  All other parameters controling
the behaviour of <span style="color:#1F618D">CASCADe </span> are given in the
`cascade_WASP19b_extract_timeseries.ini` initialization file.  The parameters are grouped
different sections in a logical way, such as parameters controling the processing steps, describing the
observations and data and overal behaviour of the code. Detailed information on each of these parameters
can be found in the documentation (see below.)

Executing the script results in the extraction of the individual (1d) spectra. These are stored
in the subdirectory 'SPECTRA' as fits files at the same level as the SPECTRAL_IMAGES directory containing the
spectral images. Two types of files are written: `COE`, which are the <span style="color:#1F618D">CASCADe </span>
optimal extraction results and 'CAE', which are the <span style="color:#1F618D">CASCADe </span>
aperture extraction results, produced using respectively an exration profile or an extraction aperture.
For furthre details we refer to the documentation.  Next to the spectral fits files,
also several diagnostic plots are produced which can be found in the 
'examples/results/WASP-19b_ibh715_transit_output_from_extract_timeseries' subdirectory.

#### Example1b: Calibrating a HST WFC3 Spectroscopic Timeseries and extracting a transit spectrum.

```bash

python3 run_CASCADe_WASP19b_calibrate_planet_spectrum.py

```

#### Example2: Calibrating a Spitzer/IRS Spectroscopic Timeseries and extracting a transit spectrum.

```bash

python3 run_CASCADe_HD189733b_calibrate_planet_spectrum.py

```

#### Example3: Calibrating a GIMINI/GMOS Spectroscopic Timeseries and extracting a transit spectrum.

```bash

python3 run_CASCADe_WASP103b_calibrate_planet_spectrum.py

```

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
