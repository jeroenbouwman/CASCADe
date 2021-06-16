
.. role:: blue

.. raw:: html

    <style> .blue {color:#1f618d} </style>
    <style> .red {color:red} </style>

:blue:`CASCADe` Initialization
==============================
:blue:`CASCADEe` uses initialization files to control its behaviour.

Extracting a spectral Timeseries
--------------------------------

[CASCADE]
^^^^^^^^^

.. code-block:: python

  cascade_save_path = WASP-19b_ibh715_transit_output_from_extract_timeseries/
  cascade_use_multi_processes = True
  cascade_max_number_of_cpus = 6
  cascade_verbose = True
  cascade_save_verbose = True

[PROCESSING]
^^^^^^^^^^^^

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
  processing_determine_initial_wavelength_shift = True


[MODEL]
^^^^^^^

.. code-block:: python

  model_type_limb_darkening = exotethys
  model_stellar_models_grid = Atlas_2000


[INSTRUMENT]
^^^^^^^^^^^^

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


Calibrating the specral timeseries and extracting the transit spectrum
----------------------------------------------------------------------

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

.. code-block:: python

  processing_compress_data = True
  processing_sigma_filtering = 4.0
  processing_nfilter = 5
  processing_stdv_kernel_time_axis_filter = 0.4
  processing_nextraction = 1
  processing_determine_initial_wavelength_shift = True


[CPM]
^^^^^^^^^^^^^^

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


[MODEL]
^^^^^^^^^^^^^^

.. code-block:: python

  model_type = batman
  model_type_limb_darkening = exotethys
  model_limb_darkening = nonlinear
  model_stellar_models_grid = Atlas_2000
  model_calculate_limb_darkening_from_model = True
  model_limb_darkening_coeff = [0.0, 0.0, 0.0, 0.0]
  model_nphase_points = 10000
  model_phase_range = 0.5
  model_apply_dilution_correcton = False


[INSTRUMENT]
^^^^^^^^^^^^^^

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
