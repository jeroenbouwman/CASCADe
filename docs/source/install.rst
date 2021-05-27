
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

```bash
conda update conda
```

and then create and environment for :blue:`CASCADe`:

```bash
conda create --name cascade
```

The conda cascade environment can then be activated with the following command:

```bash
source activate cascade
```

One can now install all necessary packages for :blue:`CASCADe` within this environment
either with conda install or pip install. The :blue:`CASCADe` package itself can be
downloaded from gitlab.


Getting :blue:`CASCADe` from gitlab
-----------------------------------

The :blue:`CASCADe` code can be found at gitlab and can be downloaded with git. Note
that the gitlab repository is private and to download the code, one needs to
create a gitlab account after which the account can be added to the list of
users which have acces to the repository. 

To download and install with a single command, type on the terminal

```bash
pip install git+git://github.com/ucl-exoplanets/ExoTETHyS.git@master
```

Alternatively, one can first clone the repository and then install.
The :blue:`CASCADe` repository path is, depending if the HTTPS or SSH protocol is used:

- HTTPS: `https://gitlab.com/jbouwman/CASCADe`
- SSH: `` git@gitlab.com:jbouwman/CASCADe ``


To get started, open a terminal window in the directory
you wish to clone the repository files into, and run one
of the following commands:

Clone via HTTPS:

```bash
git clone https://gitlab.com/jbouwman/CASCADe.git
```

Clone via SSH:

```bash
git clone git@gitlab.com:jbouwman/CASCADe.git
```

Both commands will download a copy of the files in a folder named after the
project's name. You can then navigate to the directory and start working on it
locally. After accessing the root folder from terminal, type

```bash
pip install .
```

to install the package.
In case one is using Anaconda make sure the cascade environment is
activated before using the :blue:`CASCADe` package


Installing all necessary packages within an Anaconda environment
-----------------------------------------------------------------

In the :blue:`CASCADe` main package directory an environment.yml can be found. You can
use this yml file to create or update the cascade Anaconda environment. If you
not already had created an cascade environment execute the following command:
    
```bash
conda env create -f environment.yml
```

Incase you already have an cascade environment, you can update the necessary 
packages with the following command (also use this after updating :blue:`CASCADe`
itself):

```bash
conda env update -f environment.yml
```

Make sure the :blue:`CASCADe` package is in your path. You can either set a PYTHONPATH
environment variable pointing to the location of the :blue:`CASCADe` package on your
system, or when using anaconda with the following command:

```bash
conda develop <path_to_the_CASCADe_package_on_your_system>/CASCADe
```

Updating the :blue:`CASCADe` code
---------------------------------

The insalled :blue:`CASCADe` code can be kept up to date by regularly checking for
updates. The latest version can be pulled from gitlab with the command:

```bash
git pull master
```
The latest development version can be pulled from:

```bash
git pull development
```
