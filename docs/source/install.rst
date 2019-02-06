Installing CASCADe
==================

Using Anaconda
--------------

The easiest way is to create an anaconda environment
and to install and run CASCADe wihtin this environment.
In case anaconda is not yet installed on the local system, start with 
downloading the installer, which can be found at:

	https://www.anaconda.com/download/

Once installed update conda to the latest version:

```bash
conda update conda
```

and then create and environment for CASCADe:

```bash
conda create --name cascade
```

The conda cascade environment can then be activated with the following command:

```bash
source activate cascade
```

One can now install all necessary packages for CASCADe within this environment either with conda install or
pip install. The CASCADe package itself can be downloaded from gitlab.


Getting CASCADe from gitlab
---------------------------
The CASCADe code can be found at gitlab and can be downloaded with git. Note that the gitlab repository
is private and to download the code, one needs to create a gitlab account after which the account can be 
added to the list of users which have acces to the repository. 

The CASCADe repository path is, depending if the HTTPS or SSH protocol is used:

- HTTPS: `https://gitlab.com/jbouwman/CASCADe`
- SSH: `` git@gitlab.com:jbouwman/CASCADe ``

To get started, open a terminal window in the directory
you wish to clone the repository files into, and run one
of the following commands:

Clone via HTTPS:

```bash
git clone https://gitlab.com/jbouwman/CASCADe
```

Clone via SSH:

```bash
git clone git@gitlab.com:jbouwman/CASCADe 
```

Both commands will download a copy of the files in a
folder named after the project's name. You can then navigate to the directory and start working
on it locally. In case one is using Anaconda make sure the environment is activated before installing CASCADe.

Updating the CASCADe code
-------------------------

The insalled CASCADe code can be kept up to date by regularly checking for updates. The latest version can be
pulled from gitlab with the command:

```bash
git pull master
```

