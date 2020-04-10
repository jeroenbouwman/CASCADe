
.. role:: blue

.. raw:: html

    <style> .blue {color:#1f618d} </style>
    <style> .red {color:red} </style>

Using :blue:`CASCADe`
=====================

to run the code, first load all needed modules:

.. code-block:: python

   import cascade

Then, create transit spectroscopy object

.. code-block:: python

   tso = cascade.TSO.TSOSuite()

To reset all previous defined or initialized parameters

.. code-block:: python

   tso.execute("reset")

Initialize the TSO object using ini files which define the data, model parameters and behavior of the causal pixel model implemented in :blue:`CASCADe`.

.. code-block:: python

   path = cascade.initialize.default_initialization_path
   tso = cascade.TSO.TSOSuite("initialize", "cascade_cpm.ini",
                              "cascade_object.ini",
                              "cascade_data_spectral_images.ini", path=path)

Load the observational dataset.

.. code-block:: python

   tso.execute("load_data")

Subtract the background from the spectral images.

.. code-block:: python

   tso.execute("subtract_background")

Filter the input data to flag bad pixels and create a cleaned dataset.

.. code-block:: python

   tso.execute("filter_dataset")

Determine the relative position, rotation and scale change of the source spectrum from the input spectroscopic data set.

.. code-block:: python

   tso.execute("determine_source_movement")

Correct the wavelength assaciated with the detector pixels for telescope movements.

.. code-block:: python

   tso.execute("correct_wavelengths")

Set the extraction area within which the signal of the exoplanet will be determined

.. code-block:: python

   tso.execute("set_extraction_mask")

Extract the spectrum of the Star + planet using both optimal as well as aperture extraction.

.. code-block:: python

   tso.execute("extract_1d_spectra")

Setup the matrix of regressors used to model the noise

.. code-block:: python

   tso.execute("select_regressors")

Define the eclipse model

.. code-block:: python

   tso.execute("define_eclipse_model")

Derive the calibrated time series and fit for the planetary signal

.. code-block:: python

   tso.execute("calibrate_timeseries")

Extract the planetary signal

.. code-block:: python

   tso.execute("extract_spectrum")

Correct the extracted planetary signal for non uniform subtraction of average eclipse/transit signal

.. code-block:: python

   tso.execute("correct_extracted_spectrum")

Save the planetary signal

.. code-block:: python

   tso.execute("save_results")

Plot results (planetary spectrum, residual etc.)

.. code-block:: python

   tso.execute("plot_results")

