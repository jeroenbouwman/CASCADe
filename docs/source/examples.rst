
.. role:: blue

.. raw:: html

    <style> .blue {color:#1f618d} </style>
    <style> .red {color:red} </style>

:blue:`CASCADe` examples
========================

The distribution comes with three examples demonstrating the use of the
:blue:`CASCADe` package. These use cases, together with
the observational data can be found in the **examples/** sub-directory of the
GitLab repository, and should be installed (see above) in the directory defined
by the ``CASCADE_STORAGE_PATH`` and other environment variables. The examples
cover both HST and Spitzer data as well as an example for a Generic spectroscopic
dataset, in this case an observation with the GMOS instrument on Gimini.

Example 1a: Extracting a Spectroscopic Timeseries from HST WFC3 spectral images
-------------------------------------------------------------------------------

In this example we demonstrate how extract a spectral timeseries from spectral
images. For this we use observations of WASP-19b with the Wide Field Camera 3
(WFC3) on board the HST observatory. As these observations were made in
**Staring Mode**, we use as a start the **flt** data product. The spectral
images can be found in the **HST/WFC3/WASP-19b_ibh715/SPECTRAL_IMAGES/**
sub-directory of the main data directory specified by the ``CASCADE_DATA_PATH``
environment variable.

The pipeline script for this example can be found in the **HST/WFC3/**
sub-directory of the scripts directory as specified by the
``CASCADE_SCRIPTS_PATH`` environment variable. To run this example, execute the
following script:

.. code-block:: bash

  python3 run_CASCADe_WASP19b_extract_timeseries.py

The individual steps in this script are commented and a detailed explanation on
the reduction steps can be found in the :blue:`CASCADe` documentation (see below).
Note that while running the script, interactive plots are opened which the user
needs to close before the script continues. One can prevent plots from opening
by un-commenting the line ``matplotlib.use('AGG')`` in the python script.

The initialization files for this example can be found in the
**HST/WFC3/WASP-19b_ibh715_example/** sub-directory of the directory specified
by the ``CASCADE_INITIALIZATION_FILE_PATH`` environment variable.
The **cascade_WASP19b_object.ini** contains all parameters specifying the target
star and planet. The ephemeris and orbital period given in the **.ini** file are
used to calculated an orbital phase which is attached to the individual spectra
in the final timeseries. In case of HST WFC3 spectra, as is the case here, the
stellar parameters are used to calculated an expected spectrum, which is then
used to determine a global shift of the spectrum in the wavelength direction.
It is, therefore, important that the **effective temperature, logg and
metallicity** parameters are close as possible to the correct values of the
system which the observations are being analyzed.  All other parameters
controlling the behavior of :blue:`CASCADe` are given in the
**cascade_WASP19b_extract_timeseries.ini** initialization file. The parameters
are grouped different sections in a logical way, such as parameters controlling
the processing steps, describing the observations and data and overall behavior
of the code. Detailed information on each of these parameters can be found in the
documentation.

Executing the script results in the extraction of the individual (1d) spectra.
These are stored as fits files in the sub-directory **SPECTRA/** at the same
level as the **SPECTRAL_IMAGES/** directory containing the spectral images.
Two types of files are written: **COE**, which are the :blue:`CASCADe`
**Optimal Extraction** spectra and **CAE**, which are the :blue:`CASCADe`
**Aperture Extraction** spectra, produced using respectively, optimal
extraction or an extraction aperture. For further details we refer to the
documentation. Next to the spectral fits files, also several diagnostic plots
are produced which can be found in the
**HST/WFC3/WASP-19b_ibh715_transit_output_from_extract_timeseries/** sub-directory
of the directory specified by the ``CASCADE_SAVE_PATH`` environment variable.


Example 1b: Calibrating a HST WFC3 Spectroscopic Timeseries and extracting a transit spectrum
---------------------------------------------------------------------------------------------

After creating a spectral timeseries from the spectral images in Example 1a,
we can proceed with the calibration of the spectral lightcurves and deriving the
planetary spectrum. To demonstrate how to use :blue:`CASCADe` for this, we will
take the extracted **COE** spectra of the WASP 19b observation from Example 1a
and proceed to characterize the systematics and extract the transit spectrum. To
run Example 1b, execute the following script:

.. code-block:: bash

  python3 run_CASCADe_WASP19b_calibrate_planet_spectrum.py

As with the first example, the individual steps in this script are commented and
a detailed explanation on the reduction steps can be found in the :blue:`CASCADe`
documentation. The initialization files for this example can be found in the
same directory as the initialization files for example 1a. We use again the
**cascade_WASP19b_object.ini** file to specify the target star and orbital
parameters. In contrast to the previous example, all system parameters specified
in this initialization file are relevant.

.. note::
  The current :blue:`CASCADe` version only fits
  for the transit depth. All other system parameters such as the ephemeris,
  period, inclination and semi-major axis are fixed to the values specified in
  the initialization file.

All other parameters controlling the behavior of the :blue:`CASCADe` pipeline
are given in the **cascade_WASP19b_calibrate_planet_spectrum.ini**
initialization file. A detailed explanation of the control parameters is given
in the initialization section of the documentation.

The resulting transit spectrum and diagnostic plots are stored in the
**HST/WFC3/WASP-19b_ibh715_transit_from_hst_wfc3_spectra/** sub-directory of the
directory specified by ``CASCADE_SAVE_PATH`` environment variable. The calibrated
spectrum of WASP 19b is stored in the
**WASP-19b_ibh715_bootstrapped_exoplanet_spectrum.fits** and the derived
systematics for this dataset in the
**WASP-19b_ibh715_bootstrapped_systematics_model.fits** file.

Example 2: Calibrating a Spitzer/IRS Spectroscopic Timeseries and extracting a transit spectrum
-----------------------------------------------------------------------------------------------

The :blue:`CASCADe` package can not only calibrate observations with the WFC3
instrument onboard HST, but can also handle transit spectroscopy observations
with the IRS instrument onboard the Spitzer Space Observatory. As an example,
we analyze Spitzer/IRS observations of an eclipse of HD189733b, using the with
the :blue:`CASCADe` package pre-extracted **COE** spectral data product.
The data can be found in the **SPITZER/IRS/HD189733b_AOR23439616/SPECTRA/**
sub-directory of the main data directory specified by the ``CASCADE_DATA_PATH``
environment variable. The pipeline script for this example can be found in the
**SPITZER/IRS/** sub-directory of the scripts directory as specified by the
``CASCADE_SCRIPTS_PATH`` environment variable. To run Example 2, execute the
following script:

.. code-block:: bash

  python3 run_CASCADe_HD189733b_calibrate_planet_spectrum.py

The pipeline steps used in this example are identical to the ones of Example 1b.
The initialization files for example 2 can be found in the
**SPITZER/IRS/HD189733b_AOR23439616_example/** sub-directory of the directory
specified by the ``CASCADE_INITIALIZATION_FILE_PATH`` environment variable.
Similar to the first example, the **cascade_HD189733b_object.ini** file contains
all parameters specifying the target star and orbital parameters, while the
**cascade_HD189733b_calibrate_planet_spectrum.ini** initialization file
specifies all other parameters controlling the behavior of the :blue:`CASCADe`
pipeline. The HD189733b eclipse spectrum and diagnostic plots are stored in the
**SPITZER/IRS/HD189733b_AOR23439616_eclipse_from_spitzer_irs_spectra/**
sub-directory of the directory specified by the ``CASCADE_SAVE_PATH`` environment
variable.

Example 3: Calibrating a GIMINI/GMOS Spectroscopic Timeseries and extracting a transit spectrum
-----------------------------------------------------------------------------------------------

As a final example we show how to use :blue:`CASCADe` for spectral timeseries
extracted with another software package for a generic instrument. Though spectral
extraction from spectral images or cubes is currently only implemented for the
WFC3 instrument of HST and the IRS instrument of Spitzer, the calibration of
spectral lightcurves and derivation of the planetary spectrum can be performed
for any generic spectroscopic timeseries. The previous examples showed how to
use :blue:`CASCADe` with HST and Spitzer observations. In this example we use an
observation of WASP-103b with the GMOS instrument installed at the Gemini
telescope (See Lendl et al 2017, A&A 606).

The spectral timeseries data for this example is located in the
**Generic/Gemini/GMOS/WASP103b/SPECTRA/** sub-directory of the main data
directory. To be able to run this example we stored the GMOS spectra as fits
files with an identical format as the spectral fits files created by
:blue:`CASCADe`. The pipeline script is located in the **Generic/Gemini/GMOS/**
sub-directory in the scripts directory.

To run this example, execute the following script:

.. code-block:: bash

  python3 run_CASCADe_WASP103b_calibrate_planet_spectrum.py

The initialization files for example 2 can be found in the
**Generic/Gemini/GMOS/WASP-103b_example/** sub-directory of the directory
specified by the ``CASCADE_INITIALIZATION_FILE_PATH`` environment variable.
Similar to the other examples, the **cascade_WASP103b_object.ini**
initialization file contains all parameters defining the system, and the
**cascade_WASP103b_calibrate_planet_spectrum.ini** file contains all other
parameters needed by the :blue:`CASCADe` pipeline. The WASP-103 b transit
spectrum and diagnostic plots are stored in the
**Generic/Gemini/GMOS/WASP103b_transit_from_generic_instrument/** sub-directory
of the directory specified by the ``CASCADE_SAVE_PATH`` environment variable.
