.. CASCADe documentation master file, created by
   sphinx-quickstart on Thu Jan 31 17:03:33 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. role:: red

.. role:: blue

.. raw:: html

    <style> .blue {color:#1f618d} </style>
    <style> .red {color:red} </style>

:blue:`CASCADe`: :blue:`C`\alibration of tr\ :blue:`A`\nsit :blue:`S`\pectroscopy using :blue:`CA`\usal :blue:`D`\ata
=====================================================================================================================

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

The :blue:`CASCADe` project implements a novel “data driven” method, pioneered by 
Schoelkopf et al (2016) utilizing the causal connections within a data set,
and uses this to calibrate the spectral timeseries data of single transiting
systems. The current code has been tested successfully to spectroscopic data
obtained with the Spitzer and HST observatories.

Document version:
|version|

CASCADe documentation
=====================

.. toctree::
   :maxdepth: 1

   install
   environment
   initialization
   howto
   examples
   cascade

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
