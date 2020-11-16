
.. role:: blue

.. raw:: html

    <style> .blue {color:#1f618d} </style>
    <style> .red {color:red} </style>

:blue:`CASCADe` Environment
===========================
:blue:`CASCADEe` uses the following environment variables::

    CASCADE_WARNINGS
        Switch to show or not show warnings. Can either be ‘on’ or ‘off’
    CASCADE_PATH
        Default directory for CASCADe.
    CASCADE_OBSERVATIONS
        Default path to the input observational data and calibration files.
    CASCADE_SAVE_PATH
        Default path to where CASCADe saves output.
    CASCADE_INITIALIZATION_FILE_PATH:
        Default directory for CASCADe initialization files.


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
