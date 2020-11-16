
.. role:: blue

.. raw:: html

    <style> .blue {color:#1f618d} </style>
    <style> .red {color:red} </style>

:blue:`CASCADe` Initialization
==============================
:blue:`CASCADEe` uses initialization files to control its behaviour.

CASCADE
-------

.. code-block:: python

 [CASCADE]
 cascade_save_path = WASP19b_transit_using_spectral_images/
 cascade_useMultiProcesses = True
 cascade_verbose = True
 cascade_save_verbose = True

PROCESSING
----------

.. code-block:: python

 [PROCESSING]
 processing_sigma_filtering = 3.5
 processing_max_number_of_iterations_filtering = 15
 processing_fractional_acceptance_limit_filtering = 0.005
 processing_quantile_cut_movement = 0.1
 processing_order_trace_movement = 1
 processing_nreferences_movement = 6
 processing_main_reference_movement = 4
 processing_upsample_factor_movement = 111
 processing_angle_oversampling_movement = 2
 processing_nextraction = 7
 processing_rebin_factor_extract1d = 1.05


CPM
---

.. code-block:: python

 cpm_cv_method = gcv
 cpm_lam0 = 1.0e-9
 cpm_lam1 = 1.0e3
 cpm_nlam = 150
 cpm_deltapix = 1
 cpm_nrebin = 1
 cpm_use_pca = False
 cpm_use_pca_filter = False
 cpm_number_of_pca_components = 39
 cpm_add_time = False
 cpm_add_postition = True
 cpm_add_calibration_signal = False
 cpm_calibration_signal_position = after
 cpm_clip_percentile_time = 0.00
 cpm_clip_percentile_regressors = 0.00
 cpm_calibration_signal_depth = 0.08
 cpm_relative_sig_value_limit = 4.e-1

MODEL
-----

.. code-block:: python

 [MODEL]
 model_type = batman
 model_limb darkening = quadratic
 model_limb_darkening_coeff = [0.153, 0.293]
 model_nphase_points = 10000
 model_phase_range = 0.5

OBJECT
------

.. code-block:: python

 [OBJECT]
 object_name = WASP-19 b
 object_radius = 1.386 Rjup
 object_radius_host_star = 1.004 Rsun
 object_temperature_host_star = 5500.0 K
 object_semi_major_axis = 0.01634 AU
 object_inclination = 78.78 deg
 object_eccentricity = 0.0020
 object_omega = 259 deg
 object_period = 0.788838989 d
 object_ephemeris = 2455168.96801 d
 object_kmag = 10.48 Kmag
 object_metallicity_host_star = 0.14 dex
 object_logg_host_star = 4.3932 dex(cm/s2)

CATALOG
-------

.. code-block:: python

 [CATALOG]
 catalog_use_catalog = False
 catalog_name = EXOPLANETS.ORG
 catalog_update = True


INSTRUMENT
----------     

.. code-block:: python

 [INSTRUMENT]
 instrument_observatory = HST
 instrument = WFC3
 instrument_filter = G141
 instrument_aperture = IRSUB128
 instrument_cal_filter = F139M
 instrument_cal_aperture = IRSUB512
 instrument_beam = A


OBSERVATIONS
------------

.. code-block:: python

 [OBSERVATIONS]
 observations_type = TRANSIT
 observations_mode = STARING
 observations_data = SPECTRAL_IMAGE
 observations_path = data/HST/
 observations_target_name = WASP19b
 observations_cal_path = calibration/HST/
 observations_id = ibh715
 observations_cal_version = 4.32
 observations_data_product = flt
 observations_has_background = True
 observations_uses_background_model = True
 observations_median_signal = 0.02

