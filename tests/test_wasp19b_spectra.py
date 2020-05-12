#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  Check the comput of Spectra from Spectral Image for wasp19b
#  BEWARE : it takes 700 seconds on my computer
#
"""
    This test uses py.test.
    
    1) move the reference spectra
    (the environnement variable  CASCADE_DATA_PATH is used to read the spectral image and write the spectra)

    2) define CASCADE_SAVE_PATH

    3) launch run_CASCADe_WASP19b_extract_timeseries.py

    4) read the reference *COE.fits and compare to tso.observation.dataset_optimal_extracted.data

    5) read the reference *CAE.fits and compare to tso.observation.dataset_aperture_extracted.data

    6) Erase the output FITS files, rename the reference directory for spectra

"""

import os
import cascade
import dill as pickle
import shutil
import pytest
#
from .readwrite_spectra import read_spectra_dir, compare_spectra
DEBUG = False
#
#  rtol=1.1e-3  2.62e-2
#  1e-13 failed
def test_wasp19b_spectra(verbose = False, rtol=1.e-12):
    
    ################  Prepare the Environnement ##########################
    in_dir = os.path.join(os.path.dirname(cascade.__path__[0]), 'examples')
    my_command = os.path.join(in_dir, "run_CASCADe_WASP19b_extract_timeseries.py")

    spectra_dir = os.path.join(os.path.dirname(cascade.__path__[0]), 'data/data/HST/WFC3/WASP-19b_ibh715/SPECTRA')
    ref_dir = os.path.join(os.path.dirname(cascade.__path__[0]), 'data/data/HST/WFC3/WASP-19b_ibh715/SPECTRA_reference')
    os.rename(spectra_dir, ref_dir)
    os.environ["CASCADE_SAVE_PATH"] = os.path.join(os.getcwd(), 'temporary')

    ####  this is for testing the test : read tso, 577 MO, instead of computing it
    if(DEBUG):
        my_command='read_wasp19b_tso.py'
        shutil.copytree('/Users/gastaud/test_dap/cascade/simulation_LRS/jeroen/SPECTRA/', spectra_dir)

    ################  Launch the script ##########################
    ###  wait 11 minutes on my computer

    exec(open(my_command).read())

    ################  Read and Compare Optimal ##########################

    flux1, ferror1, wavelength1, time1 = read_spectra_dir(ref_dir, pattern='*COE.fits', verbose=verbose)
    flux2, ferror2, wavelength2, time2 = read_spectra_dir(spectra_dir, pattern='*COE.fits', verbose=verbose)

    statusO, messageO = compare_spectra(flux1, ferror1, wavelength1, time1, flux2, ferror2, wavelength2, time2, rtol=rtol,atol=0, verbose=verbose)
    print("O", statusO, messageO)

    ################  Read and Compare Aperture ##########################

    flux1, ferror1, wavelength1, time1 = read_spectra_dir(ref_dir, pattern='*CAE.fits', verbose=verbose)
    flux2, ferror2, wavelength2, time2 = read_spectra_dir(spectra_dir, pattern='*CAE.fits', verbose=verbose)
    
    statusA, messageA = compare_spectra(flux1, ferror1, wavelength1, time1, flux2, ferror2, wavelength2, time2, rtol=rtol,atol=0, verbose=verbose)
    print("A", statusA, messageA)

    ################  Clean the Environnement ##########################
    shutil.rmtree(spectra_dir)
    os.rename(ref_dir, spectra_dir)
    ########################################
    assert(statusA)
    assert(statusO)

    return
