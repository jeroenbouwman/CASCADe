
.. role:: blue

.. raw:: html

    <style> .blue {color:#1f618d} </style>
    <style> .red {color:red} </style>

:blue:`CASCADe` and HST data
============================
The :blue:`CASCADe` package comes with the functionality to download data from the
HST data archive (MAST), and simultaneously generate the correct ``.ini`` files
and pipeline scripts. files. To make use of this functionality, the user can make
use of the ``build_local_hst_archive.py`` python script. This script makes use of the
``WFC_files.pickle`` database file which can be found in the ``archive_databases``
sub-directory of the main :blue:`CASCADe` storage directory as defined by the
``CASCADE_STORAGE_PATH`` environment variable. This data base file contains most
spectroscopic transit observations made with the HST WFC3 instrument. In this file,
for each visit, the data files belonging to that visit are defined, together with
the main instrument settings.  To get an overview of the main functionality one can
execute of the command line the following command:

.. code-block:: bash

  build_local_hst_archive.py --help

which will print all options available with this python script. To list all
planets for which an HST data set is defined in the database file, one can use
the following command:

.. code-block:: bash

  build_local_hst_archive.py -lap

which will print a list of all available planets. One can then pick a planet and
download all data from the MAST archive. As an example, for the data of WASP-117b,
as published in
`Carone et al 2021 <https://ui.adsabs.harvard.edu/abs/2021A%26A...646A.168C/abstract>`_
one can execute the folowing command:

.. code-block:: bash

  build_local_hst_archive.py -avp WASP-117b

which will download all the spectroscopic data for this planet. In this case the
`ima` data product is downloaded as it is a spatial scanning observation,
together with several target acquisition images, for which we use the `flt` data
product. In addition to downloading the data, also the required initialization
files are generated, which will be written to the directory specified by the
``CASCADE_INITIALIZATION_FILE_PATH`` environment variable. for generating the
initialization files, :blue:`CASCADe` makes use of the ``processing_exceptions.ini``
file which can be found in the ``archive_databases`` sub directory where also the
database file is located. This file contains several predefined values for most
observations in the database which we found to be optimal, or where derived in
published studies. The user has the possibility to overwrite these values by
creating a ``user_processing_exceptions.ini`` file and define any initialization
parameter value in there in a similar manner as is done in the
``processing_exceptions.ini`` file. If an observation has no entry in one of
these files, standard values are used from one of the exoplanet archives and from
the fits headers of the data files.

Finally, also an command line pipeline execution script is generated and stored
in the  directory defined by the ``CASCADE_SCRIPTS_PATH`` environment variable.
Executing that script will run the entire :blue:`CASCADe` pipline, extracting the
spectral time series from the spectral images or cubes and fit the transit signal
to derive the planetary spectrum. 
