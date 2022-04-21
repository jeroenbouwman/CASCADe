
.. role:: blue

.. raw:: html

    <style> .blue {color:#1f618d} </style>
    <style> .red {color:red} </style>

Installing :blue:`CASCADe`
==========================

Using Anaconda
--------------

The easiest way is to create an anaconda environment
and to install and run :blue:`CASCADe` wihtin this environment.
In case anaconda is not yet installed on the local system, start with
downloading the installer, which can be found at:

	https://www.anaconda.com/download/

Once installed update conda to the latest version:

.. code-block:: bash

  conda update conda

and then create and environment for :blue:`CASCADe`:

.. code-block:: bash

  conda create --name cascade python=3.9 ipython

Note that specifying the python version and adding ipython to be installed are
optional. The :blue: `CASCADe` code is guaranteed to work with python version 3.8
but should also work for other versions like 3.7 and 3.9.
The conda cascade environment can then be activated with the following command:

.. code-block:: bash

  source activate cascade

One can now install all necessary packages for :blue:`CASCADe` within this environment
either with conda install or pip install. The :blue:`CASCADe` package itself can be
downloaded from gitlab.


Installing :blue:`CASCADe` with pip
-----------------------------------

The easiest way to install the :blue:`CASCADe`
package is to download the distribution from PyPi,
and install the package a designated Anaconda environment with the following
commands:

.. code-block:: bash

   conda activate cascade
   pip install CASCADe-spectroscopy

.. note::
  NOTE: With :blue:`CASCADe` version 1.1.5, the batman package version 2.4.8 is
  only guaranteed to work when using numpy version 1.22.1, and with this numpy
  version one should install version 0.53.1 of the numba package to ensure the
  functionality of the latter.


Installing the required :blue:`CASCADe` data and examples
---------------------------------------------------------

All necessary data needed by :blue:`CASCADe` for it to work properly, such as
calibration files for the different spectroscopic instruments of HST and Spitzer
and configuration templates, need be downloaded from the GitLab repository before
running the code. To initialize the data download and setup of the
:blue:`CASCADe` data storage one can use the following bash command:

.. code-block:: bash

  setup_cascade.py

or alternatively from within the python interpreter:

.. code-block:: python

  from cascade.initialize import initialize_cascade
  initialize_cascade()


The additional downloaded data also includes examples and observational data to
try out the :blue:`CASCADe` package, which are explained
below.

.. note::
  NOTE: The data files will be downloaded by default to a **CASCADeSTORAGE/**
  directory in the users home directory. If a different location is preferred,
  please read the section on how to set the :blue:`CASCADe`
  environment variables first.


Installing alternatives for the :blue:`CASCADe` package
-------------------------------------------------------

The :blue:`CASCADe` code can also be downloaded from
GitLab directly by either using git or pip. To download and install with a
single command using pip, type in a terminal the following command

.. code-block:: bash

  pip install git+git://gitlab.com/jbouwman/CASCADe.git@master

which will download the latest version. For other releases replace the ``master``
branch with one of the available releases on GitLab. Alternatively, one can first
clone the repository and then install, either using the HTTPS protocal:

.. code-block:: bash

  git clone https://gitlab.com/jbouwman/CASCADe.git

or clone using SSH:

.. code-block:: bash

  git clone git@gitlab.com:jbouwman/CASCADe.git

Both commands will download a copy of the files in a folder named after the
project's name. You can then navigate to the directory and start working on it
locally. After accessing the root folder from terminal, type

.. code-block:: bash

  pip install .

to install the package.

In case one is installing :blue:`CASCADe` directly from
GitLab, and one is using Anaconda,  make sure a cascade environment is created
and activated before using our package. For convenience, in the :blue:`CASCADe`
main package directory an **environment.yml** can be found. You can use this yml
file to create or update the cascade Anaconda environment. If you not already
had created an cascade environment execute the following command:

.. code-block:: bash

  conda env create -f environment.yml

In case you already have an cascade environment, you can update the necessary
packages with the following command (also use this after updating :blue:`CASCADe`
itself):

.. code-block:: bash

  conda env update -f environment.yml

Make sure the :blue:`CASCADe` package is in your path. You can either set a
``PYTHONPATH`` environment variable pointing to the location of the
:blue:`CASCADe` package on your system, or when using anaconda with the
following command:

.. code-block:: bash

  conda develop <path_to_the_CASCADe_package_on_your_system>/CASCADe
