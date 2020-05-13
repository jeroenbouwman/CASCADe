# script add  noises
# RG Jeu 11 jul 2019

#
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import math
from numpy.random import poisson


def add_noise(data, ron=None):
    '''
    Add photon noise to data, and readout noise if given.

    :param data: data not noisy
    :param ron: Readout Noise
    :type data: ndarray
    :type ron: float

    :return: Noisy data
    :rtype: ndarray
    '''
    data_poisson = poisson(data)

    if type(ron) !=type(None):
        data_poisson = np.random.normal(0, ron, data.shape) + data_poisson

    return data_poisson


def add_sine_noise(lightcurves, mean_amplitude):
    '''
    Add sinus to light curves to simulate a causal noise. Sinus have the same phase and pulsation for each wavelengths.
    Amplitude is sightly different for each wavelengths.

    :param lightcurves: Matrix of light curves taken at different wavelengths
    :param mean_amplitude: Mean amplitude of sinus
    :type lightcurves: ndarray; dim: (nw,nt); value in electron
    :type mean_amplitude: float
    :return: Light curves + sinus
    '''
    ny, nt = lightcurves.shape
    nper = 1
    phase = math.pi

    sin_amplitude = mean_amplitude / 3 + (mean_amplitude / 5) * np.random.uniform(-1, 1, ny)
    sin_time_vector = (np.arange(nt) - nt / 2) * 4 * math.pi / nt  # Time vector in [-2pi, 2pi]

    sin_noise = 2 * mean_amplitude / 3 + sin_amplitude.reshape((ny, 1)) * np.sin(
        sin_time_vector.reshape((1, nt)) * (nper / 2) + phase)

    return lightcurves + poisson(sin_noise)
