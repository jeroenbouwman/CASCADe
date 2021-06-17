
.. role:: blue

.. raw:: html

    <style> .blue {color:#1f618d} </style>
    <style> .red {color:red} </style>

The :blue:`CASCADe` pipelines
=============================
The :blue:`CASCADe` package combines basically two distinct pipelines: The first can be considered
a "classical" spectroscopic pipeline to extract spectra from spectral images or cubes. The focus of the second pipeline
is on the calibration of spectral timeseries and lightcurve fitting to extract a transit or eclipse spectrum. In the following
we will explain both pipelines in more detail.

Using :blue:`CASCADe` to extract a spectral timeseries
-------------------------------------------------------

In this section of the manual we will discuss all the pipeline steps needed to extract a
spectral timeseries. As a first step, to be able to run the :blue:`CASCADe` code,
first import the cascade package:

.. code-block:: python

   import cascade

This module will provide all needed functionality. To start with the spectral extraction,
one first creates a transit spectroscopy object (TSO) which will hold all data and intermediate pipeline
data products:

.. code-block:: python

   tso = cascade.TSO.TSOSuite()

Next we need to initialize the TSO object using the initialization files which define the data,
model parameters and behavior of the causal pixel model implemented in :blue:`CASCADe`. For details
on these files, please consult the :doc:`initialization` section of the documentation.
To initialize the TSO instance, execute the following command:


.. code-block:: python

   path = cascade.initialize.default_initialization_path
   tso.execute("initialize", "cascade_WASP19b_extract_timeseries.ini",
               "cascade_WASP19b_object.ini", path=path)

The files used in this example are for spectral data of the transit of WASP 19 b observed
with the WFC3 instrument of HST and are part of the example use cases part of the
:blue:`CASCADe` distribution. For further details see the :doc:`examples` sections of the documentation.
Before running the code, make sure that the paths specified in the initialization
files are correctly set in addition to the path settings of the :blue:`CASCADe` environment
variables (see the :doc:`environment` section of the documentation for further details).
If the initialization files are not in the standard directory where CASCADe expect these file to be,
an additional path keyword can be set. In the example above we set the path the the standard initialization
path. If the path is the standard initialization path or that specified by the `CASCADE_INITIALIZATION_FILE_PATH`
environment variable, the path keyword can be omitted. Also note that a relative path can be specified.
In this case the path will be relative to the the standard initialization path or that specified by the
`CASCADE_INITIALIZATION_FILE_PATH` environment variable.
The loaded paramters can be found under `tso.cascade_parameters`

In case the user whiches to re-initialize the TSO object, it is highly advised to first reset the object
with the following command:

.. code-block:: python

   tso.execute("reset")

After which all previously defined control and data parameters will be reset and a cleaned
re-initialization can be made.

After initialization of the TSO object we can load the observational dataset by executing the
following pipeline step:

.. code-block:: python

   tso.execute("load_data")

This will load all spectroscopic data specified in the initialization files, including any background data or,
like with the WFC3 data, will fit a background model to the data. The loaded data can be found `tso.observations`

The next pipeline step is the subtraction the background from the spectral images:

.. code-block:: python

   tso.execute("subtract_background")

This will subtract the observed or fitted background from the data under `tso.observations.dataset` and
set the `self.observation.dataset.isBackgroundSubtracted` flag to `True`.

The next step is to identify and flag bad pixels, and to create a cleaned dataset and a spectral
extraction profile. For this we use a directional filter. To execute this step use the following command:

.. code-block:: python

   tso.execute("filter_dataset")

This will update the mask of the background subtracted dataset found in `tso.observations.dataset`
and create two new data product: the cleaded dataset `tso.cpm.cleaned_dataset` and a smoothed,
filtered dataproduct `tso.cpm.filtered_dataset` on which the extraction profile for the optimal spectral
extraction will be based.

Continuing, next we need to determine the pointing movement of the telescope as this has a direct impact on the
wavelength registration. To determine the relative position, rotation and scale change of the source spectrum
from the input spectroscopic data set, run the following pipeline step:

.. code-block:: python

   tso.execute("determine_source_movement")

We use a cross-correlation in phase space in combination with a polar transform to determine the
translational and rotational movements of the telescope. The movements are relative to a reference
observation. To minimize the risk of using a reference spectral image which a problem, we compare
to a number of reference images and than take the median values. The results of this step are stored
in `tso.cpm.spectral_movement`.

The movement of the telescope can now be used to correct the initial wavelength associated with the detector
pixels for telescope movements:

.. code-block:: python

   tso.execute("correct_wavelengths")

This will correct the wavelengths for each time step of the data and the cleaned and filtered data product.

The last step before the spectral extraction is setting the extraction area on the detector within which the spectrum
will be extracted. This is done with the following command:

.. code-block:: python

   tso.execute("set_extraction_mask")

The extraction aperture is centered around the spectral trace with a certain width which is defined in the initialization
files, and takes into account the movements of the telescope determined in the `determine_source_movement` step.
This ensures that also for larger movements, like nodding of spatial scanning like with the HST/WFC3 observations,
the extraction aperture is always properly centered. The extraction mask in time is stored under `tso.cpm.extraction_mask`.

The final step is then the spectral extraction of the spectral timeseries of the star & planet:

.. code-block:: python

   tso.execute("extract_1d_spectra")

This step will extract the spectral timeseries using both optimal as well as aperture extraction,
and rebin the resulting spectra to a uniform wavelength grid, for which the rebin factor is specified
in the initialization files. The results are stored in `tso.observation.dataset_optimal_extracted`,
respectively `tso.observation.dataset_aperture_extracted`. The extraction profile and its associated
mask are stored under `tso.cpm.extraction_profile` and `tso.cpm.extraction_profile_mask`. Apart from
storing the results in the TSO object, the spectra are also saved as fits spectral tables at the
location specified by the `CASCADE_DATA_PATH` environment variable and the 'observations_path'
variable defined in the initialization files. The fits files of the optimal extracted spectra are
labeled with 'COE' (Cascade Optimal Extracted) and the aperture extracted spectra are labeled
with CAE (Cascade Aperture Extracted).



Using :blue:`CASCADe` to calibrate the spectral timeseries and determine the planetary spectrum
-----------------------------------------------------------------------------------------------

After extracting the spectral timeseries, we can now proceed to determine the systematics on the timeseries
and extract the transit or eclipse spectrum of the planet. The first few pipeline steps are identical to those of
the previous section. First we import the :blue:`CASCADe` package:

.. code-block:: python

   import cascade

than we create an instance of the timeseries object:

.. code-block:: python

   tso = cascade.TSO.TSOSuite()

next we initialize the TSO object:

.. code-block:: python

  path = cascade.initialize.default_initialization_path
  tso.execute("initialize", "cascade_WASP19b_calibrate_planet_spectrum.ini",
                  "cascade_WASP19b_object.ini")

note that again we use the WASP-19 b example provided in the examples coming with the :blue:`CASCADe` package.

We then load the spectra into the tso object:

.. code-block:: python

  tso.execute("load_data")

For completeness we also execute the  `subtract_background` step:

.. code-block:: python

  tso.execute("subtract_background")

as there is the possibility that the user provides not background subtracted spectra. If the
spectra are background subtracted, this can be indicated in the initialization files by switching
the `observations_has_background` variable to `False`. In that case the `subtract_background`
step will be skipped.

As with the spectral extraction pipeline, we also execute the `filter_dataset` pipeline step:

.. code-block:: python

  tso.execute("filter_dataset")

to flag any spurious spectral data points and create a cleaned dataset. Note that in this case
we use a simple median filtering in contrast to the directional filtering used for the spectral
data cubes and images.

In case of HST/WFC3 spectra, a check of the wavelength solution for an overall wavelength shift
is made with the following pipeline step:

.. code-block:: python

  tso.execute("check_wavelength_solution")

For this, a sample model of the observed spectrum is created using the stellar parameters
defined in the initialization files together with the sensitivity curve of the spectrograph.
Using a cross correlation between the model and the observed time averaged spectrum,
an overall wavelength shift is determined and corrected for. Note that this step is already
done during spectral extraction with the :blue:`CASCADe` pipeline. For other instruments this step
can be ignored by switching the `processing_determine_initial_wavelength_shift` parameter to `False`

After all the previous steps we now can run the main pipeline task for calibrating the spectral
timeseries and fitting the transit or eclipse:

.. code-block:: python

  tso.execute("calibrate_timeseries")

The derived transit spectra from the regression fit and bootstrap analysis are stored under
`tso.exoplanet_spectrum` The used transit model under `tso.model` and the results from the
regression analysis under `tso`.calibration_results'.

The final pipeline step will plot the resulting spectra and save the transit or eclipse spectrum
as a fits file:

.. code-block:: python

   tso.execute("save_results")

All results from this pipeline will be stored in the location defined by the `CASCADE_SAVE_PATH` environment variable
and the `cascade_save_path` parameter from the initialization files.
