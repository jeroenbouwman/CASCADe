
.. role:: blue

.. raw:: html

    <style> .blue {color:#1f618d} </style>
    <style> .red {color:red} </style>

:blue:`CASCADe` Environment
===========================
:blue:`CASCADEe` uses the following environment variables::

    CASCADE_WARNINGS
        Switch to show or not show warnings. Can either be ‘on’ or ‘off’
    CASCADE_DEFAULT_PATH
        Default directory for CASCADe.
    CASCADE_DATA_PATH
        Default path to the input observational data and calibration files.
    CASCADE_SAVE_PATH
        Default path to where CASCADe saves output.
    CASCADE_INITIALIZATION_FILE_PATH:
        Default directory for CASCADe initialization files.
    CASCADE_LOG_PATH
        Default path to were log files are stored.    


in case the environment variables are not set by the user, :blue:`CASCADEe` uses default values defined in the
initialize module. The CASCADE_WARNINGS variable is by default set to 'on'. The CASCADE_DEFAULT_PATH should be
set to the installation directory of the :blue:`CASCADEe` code. Per default, all other environment variables,
defining the location of the data, were results need to be stored, were the initialization files are located
and were log files are written, are set relative to the CASCADE_DEFAULT_PATH environment variable by default to 
`data/`,  `examples/results/`, `examples/init_files/` and  `examples/logs/`, respectively. With these standard values the
example scripts provided in the :blue:`CASCADEe` distribution can be executed without problems. However,
the user is adviced to set these variables to values other than the standard values as not to clutter the
:blue:`CASCADEe` source directory. Especially for large data sets, changing the default to a directory on a
large storage disk is recommended. 


To set an environment veriable you can use the following shell commands on Linux or MAC OS:

.. code-block:: console

   export CASCADE_WARNINGS=off

Or in a Windows environment:

.. code-block:: console

   set CASCADE_WARNINGS=off

Alternatively one could set the environment variables in a python script before importing :blue:`CASCADe`:

.. code-block:: python

   import os
   os.environ['CASCADE_WARNINGS']='off'

   import cascade

Note that in the latter case the environment varables MUST be defined before importing cascade.