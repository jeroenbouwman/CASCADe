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

The <span style="color:#1f618d">CASCADe </span> project implements a novel “data driven” method, pioneered by
Schoelkopf et al (2016) utilizing the causal connections within a data set,
and uses this to calibrate the spectral timeseries data of single transiting
systems. The current code has been tested successfully to spectroscopic data
obtained with the Spitzer and HST observatories.



## Using Cascade

to run the code, first load all needed modules:
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

initialize the TSO object using ini files which define the data, model parameters and behavior of the causal pixel model implemented in CASCADe.
```python
path = cascade.initialize.default_initialization_path
tso = cascade.TSO.TSOSuite("initialize", "cascade_cpm.ini",
                           "cascade_object.ini",
                           "cascade_data_spectral_images.ini", path=path)
```

Load the observational data
```python
tso.execute("load_data")
```

Subtract the background
```python
tso.execute("subtract_background")
```

Sigma clip data
```python
tso.execute("sigma_clip_data")
```

Determine the position of source from the spectroscopic data set
```python
tso.execute("determine_source_position")
```

Set the extraction area within which the signal of the exoplanet will be determined
```python
tso.execute("set_extraction_mask")
```

Extract the spectrum of the Star + planet in an optimal way
```python
tso.execute("optimal_extraction")
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
