
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
data product. The spectral images can be found in the data subdirectory under
`data/HST/WFC3/SPECTRAL_IMAGES`. If the user has changed the default settings of the
environment variables, these data will need to be copied to the corresponding user directory.

To run this example, execute the following script:

```bash
python3 run_CASCADe_WASP19b_extract_timeseries.py
```

The individual steps in this script are commentend and a detailed explanation on the 
reduction steps can be found in the :blue:`CASCADEe` documentation
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
the behaviour of :blue:`CASCADEe` are given in the
`cascade_WASP19b_extract_timeseries.ini` initialization file.  The parameters are grouped
different sections in a logical way, such as parameters controling the processing steps, describing the
observations and data and overal behaviour of the code. Detailed information on each of these parameters
can be found in the documentation (see below.)

Executing the script results in the extraction of the individual (1d) spectra. These are stored
in the subdirectory 'SPECTRA' as fits files at the same level as the SPECTRAL_IMAGES directory containing the
spectral images. Two types of files are written: `COE`, which are the :blue:`CASCADEe`
optimal extraction results and `CAE`, which are the :blue:`CASCADEe`
aperture extraction results, produced using respectively an exration profile or an extraction aperture.
For furthre details we refer to the documentation.  Next to the spectral fits files,
also several diagnostic plots are produced which can be found in the 
'examples/results/WASP-19b_ibh715_transit_output_from_extract_timeseries' subdirectory.

Example1b: Calibrating a HST WFC3 Spectroscopic Timeseries and extracting a transit spectrum.
---------------------------------------------------------------------------------------------

```bash
python3 run_CASCADe_WASP19b_calibrate_planet_spectrum.py
```

Example2: Calibrating a Spitzer/IRS Spectroscopic Timeseries and extracting a transit spectrum.
-----------------------------------------------------------------------------------------------

```bash
python3 run_CASCADe_HD189733b_calibrate_planet_spectrum.py
```

Example3: Calibrating a GIMINI/GMOS Spectroscopic Timeseries and extracting a transit spectrum.
-----------------------------------------------------------------------------------------------

```bash
python3 run_CASCADe_WASP103b_calibrate_planet_spectrum.py
```
