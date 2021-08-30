
.. role:: blue

.. raw:: html

    <style> .blue {color:#1f618d} </style>
    <style> .red {color:red} </style>

:blue:`CASCADe` Initialization
==============================
:blue:`CASCADe` uses initialization files to control its behavior. This includes the definition of the
input data and observational setup, the regression modeling and lightcurve fitting and general settings for
the :blue:`CASCADe` package. In this section of the documentation we will explain the control parameters defined
in the initialization files in detail. For this we will use the WASP-19b example provided with this package
(see also the :doc:examples section of this manual).

Extracting a spectral timeseries
--------------------------------
The following initialization parameters are used with the `run_CASCADe_WASP19b_extract_timeseries.py` example script
and can be found in the `cascade_WASP19b_extract_timeseries.ini` file, both coming with the :blue:`CASCADe` package.
The configuration parameters are grouped in sections, depending on what particular aspect of they
control.

[CASCADE]
^^^^^^^^^

This initialization file section controls some of the general behavior of the package, the following parameters
can be specified here:

.. code-block:: python

 cascade_save_path = WASP-19b_ibh715_transit_output_from_extract_timeseries/
  cascade_use_multi_processes = True
  cascade_max_number_of_cpus = 22
  cascade_verbose = True
  cascade_save_verbose = True

The ``cascade_save_path``  specifies the path were verbose output will be saved. This path can either be
relative to the path specified by the `CASCADE_SAVE_PATH` environment variable, or absolute, in the latter case
overwriting the environment settings. The ``cascade_use_multi_processes`` flag specifies if the main calculations should be
done in parallel using the ray package. The ``cascade_max_number_of_cpus`` parameter controls the maximum number of CPUs
the parallel calculations can use. The ``cascade_verbose`` flag controls if verbose plots will be generated and
the ``cascade_save_verbose`` flag if these plots will be saved or not.


[PROCESSING]
^^^^^^^^^^^^

The parameters in this initialization file section control the processing and spectral extraction
of the :blue:`CASCADe` package, and as such will have to be chosen with care.

.. code-block:: python

 processing_source_selection_method = nearest
  processing_bits_not_to_flag = [0, 10, 12, 14]
  processing_compress_data = True
  processing_sigma_filtering = 4.0
  processing_max_number_of_iterations_filtering = 20
  processing_fractional_acceptance_limit_filtering = 0.005
  processing_quantile_cut_movement = 0.1
  processing_order_trace_movement = 1
  processing_nreferences_movement = 7
  processing_main_reference_movement = 3
  processing_upsample_factor_movement = 111
  processing_angle_oversampling_movement = 2
  processing_nextraction = 7
  processing_extend_roi = [1.0, 1.0, 1.0, 1.0]
  processing_drop_frames = {'up': [-1], 'down': [-1]}
  processing_rebin_factor_extract1d = 1.05
  processing_auto_adjust_rebin_factor_extract1d = True
  processing_renorm_spatial_scans = True
  processing_determine_initial_wavelength_shift = True

The ``processing_source_selection_method`` parameter controls how the target is found
in the target acquisition images taken with the HST/WFC3 spectroscopic observations.
Note that finding the target in these images is important as it determines the initial
wavelength solution and the placement of the region of interest on the spectroscopic
images.  If set to `nearest`, the target found nearest to the expected position is
used. If set to 'brightest' the brightest target in the FOV is used. Other possible
value are 'second nearest' and 'second brightest'. The ``processing_bits_not_to_flag``
is specific to HST observations and indicate which bit set in the data quality mask
should be flagged. If one of those bits are set, the corresponding pixel data is
flagged as 'bad' The ``processing_compress_data`` flag indicates if data rows flagged
for all time steps as 'bad' should be removed from the dataset. The default value is
'True' The ``processing_sigma_filtering`` parameter indicates the sigma level above
which to flag pixels as bad during the spatial filtering. The ``processing_max_number_of_iterations_filtering``
parameter indicates the maximum number of iterations used in the spatial filtering.
The standard number is 20. This parameter is advised to be used only to ensure that
the spatial filtering process terminates after some time. The filtering normally
terminates upon convergence which is reached if the fraction of pixels still considered
bad is lower than the relative fraction specified by the
``processing_fractional_acceptance_limit_filtering`` parameter, which is set to a
standard value of half a percent. To determine the source position in the cross-dispersion
direction and spectral trace trace the ''center-of-light'' (COL) of the dispersed
light is determined. For comparison to the dispersion profile and calibration files
a polynomial trace is fitted to the COL. The parameter ``processing_quantile_cut_movement``
indicates the lower quantile which is cut before determining the cross-dispersion
position. The ``processing_order_trace_movement`` parameter determines the polynomial
order of the fitted trace.
To determine the relative telescope movement, a cross correlation method is applied to
determine the relative movement (including rotation and scaling) of the spectral images to
a reference image. To ensure that a chosen reference image did not have by chance some
artifacts which influences the results a number, geven by the parameter ``processing_nreferences_movement``
of reference images are picked, after which a median movement is calculated. As the
main reference of the reference images the image with the the index number given by the
``processing_main_reference_movement`` is used. To ensure sub-pixel movement accuracy
the images are over-sampled. The parameter ``processing_upsample_factor_movement`` indicates the
oversampling of the spectral pixels, while the ``processing_angle_oversampling_movement`` indicates the
oversampling of the angle after transformation to polar coordinates to determine the
relative rotation and scaling.
For the spectral extraction both a region-of-interest (ROI) as well as a extraction aperture are defined.
The width of the latter is defined by the ``processing_nextraction`` parameter, which
gives the extraction aperture in number of pixels. The ROI is set based on the instrument and
expected position of the target. The ROI can be extended by changing the values of the
``processing_extend_roi`` parameter. The values are multiplicative and extend the lower wavelength bound,
the upper wavelength bound, the left aperture boundary and the the right boundary, respectively.
Values smaller then 1 indicate a shrinkage, larger then one an expansion.
In case of the use of spectral cubes (line with HST spatial scanning mode observations),
the user has the option to drop samples from the detector reads. The ``processing_drop_frames`` parameter
lists which samples on a detector ramp should be dropped. In case of HST spatial scans for both up and down
scans this can be indicated separately. A -1 value indicates that no samples are dropped.
After spectral extraction, the extracted 1d spectra are re-binned to a uniform
wavelength grid, with a down sampling factor given by the parameter ``processing_rebin_factor_extract1d``.
The value of this parameter should be equal or larger then 1.0. If the
``processing_auto_adjust_rebin_factor_extract1d`` is set, the wavelength re-binning
is outomatically adjusted such that there are fewer wavelength bins than time samples.
The value given by the ``processing_rebin_factor_extract1d`` parameter is in this
case used as a lower limit. In case of HST spatial scanning observations with both up and down scans, if the
``processing_renorm_spatial_scans`` flag is set, the up and down scan data are re-normalized
to the overall mean flux value. If the ``processing_determine_initial_wavelength_shift`` flag is set,the time
average spectrum is compared to a model spectrum to determine a possible general
wavelength shift. Currently only implemented for HST observations


[MODEL]
^^^^^^^

The parameters in this section define the limbdarkening and lightcurve model. However, for the spectral
extraction pipeline only the stellar model grids, part of the used limbdarkening code, are used to make a
simple estimate of the expected spectrum.

.. code-block:: python

 model_type_limb_darkening = exotethys
  model_stellar_models_grid = Phoenix_2018

The ``model_type_limb_darkening`` parameter selects which limbdarkening code is used with :blue:`CASCADe`.
At present only `exotethys` can be selected. The ``model_stellar_models_grid`` indicated which stellar grid to use.
Apart from the 'Phoenix_2018' grid, other options are 'Phoenix_2012_13' and 'Atlas2000', the latter being the
standard option. These grids come with the used Exotethys limbdarkening code.

[INSTRUMENT]
^^^^^^^^^^^^

The parameters in this section describe the used instrument, in this example the WFC3 instrument onboard HST.
The only other instrument currently implemented in the spectral extraction pipeline is the IRS instrument of Spitzer.
JWST instruments will follow in the near future.

.. code-block:: python

 instrument_observatory = HST
  instrument = WFC3
  instrument_filter = G141
  instrument_aperture = IRSUB128
  instrument_cal_filter = F139M
  instrument_cal_aperture = IRSUB512
  instrument_beam = A


[OBSERVATIONS]
^^^^^^^^^^^^^^

The parameters in this section of the initialization files describe the observational data.

.. code-block:: python

 observations_type = TRANSIT
  observations_mode = STARING
  observations_data = SPECTRAL_IMAGE
  observations_path = data/
  observations_target_name = WASP-19b_ibh715
  observations_cal_path = calibration/
  observations_id = ibh715
  observations_cal_version = 4.32
  observations_data_product = flt
  observations_has_background = True
  observations_uses_background_model = True


Calibrating the spectral timeseries and extracting the transit spectrum
-----------------------------------------------------------------------

[CASCADE]
^^^^^^^^^^^^^^

.. code-block:: python

 cascade_save_path = WASP-19b_ibh715_transit_from_hst_wfc3_spectra/
  cascade_use_multi_processes = True
  cascade_max_number_of_cpus = 6
  cascade_verbose = True
  cascade_save_verbose = True

[PROCESSING]
^^^^^^^^^^^^^^

The parameters in this initialization file section control the processing the spectral timeseries observations
before the regression analysis. As far less pre-processing needs to be done in this pipeline compared to the
spectral extraction pipeline only a few parameters need to be set.

.. code-block:: python

 processing_compress_data = True
  processing_sigma_filtering = 4.0
  processing_nfilter = 5
  processing_stdv_kernel_time_axis_filter = 0.4
  processing_nextraction = 1
  processing_determine_initial_wavelength_shift = True

The ``processing_sigma_filtering``, ``processing_nfilter`` and ``processing_stdv_kernel_time_axis_filter`` control
the filtering for deviating signals in the spectral time series. In the above setting a sigma clip is made for signals
deviating by more than 4 sigma in a 5 pixel box in the wavelength direction. The small kernel size in the time direction
ensures very little filtering in time. This is important as we are constructing a regression model of a time series and which
to preserve the systematics in time for proper characterization. The user is advised to leave these settings as the are.
The ``processing_nextraction`` in this pipeline is only here for legacy reasons and has no effect. It will be removed
in future releases. The ``processing_determine_initial_wavelength_shift`` switch controls the use of the
`check_wavelength_solution` pipeline step. it is currently only implemented for the slitless WFC3 observations,
for any other instrument this should be set to `False`

[CPM]
^^^^^^^^^^^^^^

This section of the initialization parameters contains all parameters controlling the regression model (or Causal Pixel Model)
applied to calibrate the spectral lightcurves and to extract the transit or eclipse spectrum.

.. code-block:: python

 cpm_lam0 = 0.001
  cpm_lam1 = 10000.0
  cpm_nlam = 140
  cpm_deltapix = 7
  cpm_ncut_first_integrations = 10
  cpm_nbootstrap = 250
  cpm_regularization_method = value
  cpm_add_time = True
  cpm_add_time_model_order = 1
  cpm_add_position = True
  cpm_regularize_depth_correction = True
  cpm_sigma_mse_cut = 5.0


[MODEL]
^^^^^^^^^^^^^^

The parameters in this section define the limbdarkening and lightcurve model used to fit the transit signal.

.. code-block:: python

 model_type = batman
  model_type_limb_darkening = exotethys
  model_limb_darkening = nonlinear
  model_stellar_models_grid = Phoenix_2018
  model_calculate_limb_darkening_from_model = True
  model_limb_darkening_coeff = [0.0, 0.0, 0.0, 0.0]
  model_nphase_points = 10000
  model_phase_range = 0.5
  model_apply_dilution_correcton = False


[INSTRUMENT]
^^^^^^^^^^^^^^

The parameters in this section describe the used instrument and observatory. The parameters are the same as
described in the previous section describing the initialization of the spectral extract pipeline. Note that some
of the parameters are here of legacy reasons and have no influence on the lightcurve calibration pipeline, such as
the ``instrument_cal_filter``, ``instrument_cal_aperture`` and ``instrument_beam`` parameters.

.. code-block:: python

 instrument_observatory = HST
  instrument = WFC3
  instrument_filter = G141
  instrument_aperture = IRSUB128
  instrument_cal_filter = F139M
  instrument_cal_aperture = IRSUB512
  instrument_beam = A


[OBSERVATIONS]
^^^^^^^^^^^^^^

The parameters in this section of the initialization files describe the observational data. They are the same
as used in the spectral extraction pipeline, though with slightly different values as described below.

.. code-block:: python

 observations_type = TRANSIT
  observations_mode = STARING
  observations_data = SPECTRUM
  observations_path = data/
  observations_target_name = WASP-19b_ibh715
  observations_cal_path = calibration/
  observations_id = ibh715
  observations_cal_version = 4.32
  observations_data_product = COE
  observations_has_background = False

For the calibration pipeline, as it is only working with spectral timeseries, the ``observations_mode`` and
``observations_data`` parameters should always be set to, respectively, `STARING` and `SPECTRUM`.
The ``observations_data_product`` parameter is set to `COE` in this example to use the optimal extracted spectra.
Alternative values are 'CAE' for the aperture extracted spectra. Note that the user can define/use their own data product.
The name of the data product should appear in the fits file names just before the .fits extension.
As the used spectra are already background subtracted, the ``observations_has_background`` switch needs to be set
to 'False'


The stellar and planetary parameters
------------------------------------
Finally, in addition to all parameters controlling the behavior of the spectral extraction and
regression analysis, we have the parameters specifying the star and planetary system. These parameters are
used in the lightcurve model and the `check_wavelength_solution` pipeline step.

[OBJECT]
^^^^^^^^

In this section all relevant stellar and planetary parameters are specified. These should be self explanatory.
Note the units for the different parameters. When modifying the parameters, these units should be correctly
specified such that the astropy package can handle them.

.. code-block:: python

 object_name = WASP-19 b
  object_ephemeris = 2454775.3366 d
  object_period = 0.7888396 d
  object_radius = 1.31 jupiterRad
  object_semi_major_axis = 0.01649 AU
  object_inclination = 80.5 deg
  object_eccentricity = 0.0061
  object_omega = 259.0 deg
  object_radius_host_star = 0.94 solRad
  object_temperature_host_star = 5568.0 K
  object_kmag = 10.48 Kmag
  object_metallicity_host_star = 0.0 dex
  object_logg_host_star = 4.45 dex(cm / s2)
  object_distance = 268.325 pc

[CATALOG]
^^^^^^^^^

.. code-block:: python

 catalog_use_catalog = False
  catalog_name = NASAEXOPLANETARCHIVE
  catalog_update = True

[DILUTION]
^^^^^^^^^^

.. code-block:: python

 dilution_temperature_star = 3700.0 K
  dilution_metallicity_star = 0.0 dex
  dilution_logg_star = 5.0 dex(cm / s2)
  dilution_flux_ratio = 0.01
  dilution_band_wavelength = 1.26 micron
  dilution_band_width = 0.05 micron
  dilution_wavelength_shift = 0.0 micron
