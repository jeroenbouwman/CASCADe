
[![pipeline status](https://gitlab.com/jbouwman/CASCADe/badges/master/pipeline.svg)](https://gitlab.com/jbouwman/CASCADe/commits/master)

# CASCADe: Calibration of trAnsit Spectroscopy using CAusal Data

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

The <span style="color:#1F618D">CASCADe </span> project implements a novel “data driven” method, pioneered by
Schoelkopf et al (2016) utilizing the causal connections within a data set,
and uses this to calibrate the spectral timeseries data of single transiting
systems. The current code has been tested successfully to spectroscopic data
obtained with the Spitzer and HST observatories.


## Using Cascade

The CASCADe distribution comes with a few working examples and datasets which can be found in the examples and
data sub-directories, respectively. All needed parameters for running the code are defined in the provided 
initialization files. The user only has to change those, particularly the path settings to the data and .ini files
to run the examples. 

The basic structure of the CASCADe pipeline is as follows:
To run the code, first import the cascade module
```python
import cascade
```

Then, create transit spectroscopy object
```python
tso = cascade.TSO.TSOSuite()
```

To reset all previous divined or initialized parameters
```python
tso.execute("reset")
```

Initialize the TSO object using .ini files which define the data, system and causal noise model parameters 
needed to load all data, create a calibration model and to extract the planetary signal
```python
path = cascade.initialize.default_initialization_path
tso = cascade.TSO.TSOSuite("initialize", "cascade_cpm.ini",
                           "cascade_object.ini",
                           "cascade_data.ini", path=path_to_ini_files)
```

Load the observational data. In case of spectral images, any background data or background model
will also be loaded at this step.
```python
tso.execute("load_data")
```

Subtract the background from the science data. In case no background needs to be subtracted, 
a flag can be set in the configuration files and this step will be ignored.
```python
tso.execute("subtract_background")
```

Filter the input data to flag bad pixels and create a cleaned dataset.
```python
   tso.execute("filter_dataset")
```

Determine the relative position, rotation and scale change of the source spectrum from the input spectroscopic data set.
```python
   tso.execute("determine_source_movement")
```

Correct the wavelength assaciated with the detector pixels for telescope movements.
```python
   tso.execute("correct_wavelengths")
```

Set the extraction area within which the signal of the exoplanet will be determined
```python
   tso.execute("set_extraction_mask")
```

Extract the spectrum of the Star + planet using both optimal as well as aperture extraction.
```python
   tso.execute("extract_1d_spectra")
```

Setup the matrix of regressors used to model the noise
```python
tso.execute("select_regressors")
```

Define the eclipse model
```python
tso.execute("define_eclipse_model")
```

Derive the calibrated time series and fit for the planetary signal
```python
tso.execute("calibrate_timeseries")
```

Extract the planetary signal
```python
tso.execute("extract_spectrum")
```

Correct the extracted planetary signal for non uniform subtraction of average eclipse/transit signal
```python
tso.execute("correct_extracted_spectrum")
```

Save the planetary signal
```python
tso.execute("save_results")
```

Plot results (planetary spectrum, residual etc.)
```python
tso.execute("plot_results")
```

## Documentation

Extended documentation can be found in the `docs`  directory. 
To generate the documentation, make sure you have 'Sphinx' installed and than execute in the in the
`docs` directory the following command:
```bash
make html
make latexpdf
```

The generated HTML and PDF will be located in the location specified by `conf.py`,
in this case `build/html` and 'build\latex'.
