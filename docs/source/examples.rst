
.. role:: blue

.. raw:: html

    <style> .blue {color:#1f618d} </style>
    <style> .red {color:red} </style>

:blue:`CASCADe` examples
========================

The distribution comes with three examples demonstrating the use of the
:blue:`CASCADEe` package. These use cases can be found
in the 'examples/' subdirectory. They cover both HST and Spitzer data as well as
an example for a Generic spectroscopic dataset.

Example1a: Extracting a Spectroscopic Timeseries from HST WFC3 spectral images.
-------------------------------------------------------------------------------

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
reduction steps can be found in the :blue:`CASCADe` documentation
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
the behavior of :blue:`CASCADe` are given in the
`cascade_WASP19b_extract_timeseries.ini` initialization file.  The parameters are grouped
different sections in a logical way, such as parameters controlling the processing steps, describing the
observations and data and overall behavior of the code. Detailed information on each of these parameters
can be found in the documentation (see below.)

Executing the script results in the extraction of the individual (1d) spectra. These are stored
in the sub-directory 'SPECTRA' as fits files at the same level as the SPECTRAL_IMAGES directory containing the
spectral images. Two types of files are written: `COE`, which are the :blue:`CASCADe`
optimal extraction results and `CAE`, which are the :blue:`CASCADe`
aperture extraction results, produced using respectively an extraction profile or an extraction aperture.
For further details we refer to the documentation.  Next to the spectral fits files,
also several diagnostic plots are produced which can be found in the
`WASP-19b_ibh715_transit_output_from_extract_timeseries` sub-directory of the `CASCADE_SAVE_PATH`
directory which, if not set by the user, is the `CASCADe\examples\results` directory of the :blue:`CASCADe` distribution.

Example1b: Calibrating a HST WFC3 Spectroscopic Timeseries and extracting a transit spectrum.
---------------------------------------------------------------------------------------------

After extracting the spectral timeseries from spectral data cubes or images as shown in Example 1a,
we can proceed with the extraction of the planetary spectrum and the calibration of the spectral lightcurves.
To demonstrate how to use :blue:`CASCADe` to do this, we will take the extracted spectra of the WASP 19b observation from Example 1a and proceed to characterize the systematics and extract the transit spectrum.
To run Example 1b, execute the following script:

```bash
python3 run_CASCADe_WASP19b_calibrate_planet_spectrum.py
```

The individual steps in this script are commented and a detailed explanation on the
reduction steps can be found in the :blue:`CASCADe` documentation
(see below). The initialization files for this example can be found in the `examples/init_files`
sub-directory. The `cascade_WASP19b_object.ini` file contains all parameters specifying the target star
and planet. In contrast to the previous example, all system parameters specified in
this initialization file are relevant. Note, that the current :blue:`CASCADe`
version only fits for the transit depth. All other system parameters such as the ephemeris, period, inclination and
semi-major axis are fixed to the values specified in the `.ini` file. All other parameters controlling
the behavior of :blue:`CASCADe` are given in the
`cascade_WASP19b_calibrate_planet_spectrum.ini` initialization file.  A detailed explanation of
the control parameters is given in the documentation.

The resulting transit spectrum and diagnostic plots are stored in the `WASP-19b_ibh715_transit_from_hst_wfc3_spectra`
sub-directory in the `CASCADE_SAVE_PATH` directory which, if not set by the user, is the `CASCADe\examples\results`
directory of the :blue:`CASCADe` distribution.

Example2: Calibrating a Spitzer/IRS Spectroscopic Timeseries and extracting a transit spectrum.
-----------------------------------------------------------------------------------------------

Apart from the observations with the WFC3 instrument onboard HST, :blue:`CASCADe`
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
all other parameters controlling the behavior of :blue:`CASCADe`.

The HD189733b eclipse spectrum and diagnostic plots are stored in the `HD189733b_AOR23439616_eclipse_from_spitzer_irs_spectra`
sub-directory in the `CASCADE_SAVE_PATH` directory which, if not set by the user, is the `CASCADe\examples\results`
directory of the :blue:`CASCADe` distribution.

Example3: Calibrating a GIMINI/GMOS Spectroscopic Timeseries and extracting a transit spectrum.
-----------------------------------------------------------------------------------------------

As a final example we show how to use :blue:`CASCADe` for spectral timeseries
which have been extracted with another software package for a generic instrument. Though
spectral extraction from spectral images or cubes is currently only implemented for HST/WFC3 and SpitzerIRS data,
the extraction of the planetary signal and the calibration of the spectral lightcurves can be performed for any generic spectroscopic timeseries.
The previous examples showed how to use :blue:`CASCADe` with HST and Spitzer observations. In this example we use an observation with the GIMINI/GMOS instrument of WASP-103b (See Lendl et al 2017, A&A 606).

To run this example, execute the following script:

```bash
python3 run_CASCADe_WASP103b_calibrate_planet_spectrum.py
```

To be able to run this example we stored the GMOS spectra as fits files with an identical format as the fits files
created by :blue:`CASCADe` to store the extracted spectra.

The initialization files for this example are `cascade_WASP103b_object.ini`, which contains the system parameters and
`cascade_WASP103b_calibrate_planet_spectrum.ini`, which contains all other necessary parameters.
The WASP-103 b transit spectrum and diagnostic plots are stored in the `WASP103b_transit_from_generic_instrument`
sub-directory in the `CASCADE_SAVE_PATH` directory which, if not set by the user, is the `CASCADe\examples\results`
directory of the :blue:`CASCADe` distribution.
