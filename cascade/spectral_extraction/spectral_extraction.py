#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This file is part of the CASCADe package which has been
# developed within the ExoplANETS-A H2020 program.
#
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright (C) 2018, 2019, 2020, 2021  Jeroen Bouwman
"""Module defining the spectral extraction functionality used in cascade."""

import math
from functools import partial
import collections
import warnings
import copy
from itertools import zip_longest
import multiprocessing as mp
from asyncio import Event
from typing import Tuple
from psutil import virtual_memory, cpu_count
from tqdm import tqdm
import ray
from ray.actor import ActorHandle
import statsmodels.api as sm
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
from scipy import ndimage
from scipy.ndimage import binary_dilation
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import Kernel2D
from astropy.modeling.models import Gaussian2D
from astropy.modeling.parameters import Parameter
from astropy.convolution import interpolate_replace_nans
from astropy.convolution import Gaussian1DKernel
from astropy.stats import sigma_clip
from skimage.registration import phase_cross_correlation
from skimage.transform import warp
from skimage._shared.utils import safe_as_int
from skimage.transform import rotate
from skimage.transform import SimilarityTransform
from sklearn.preprocessing import RobustScaler

from ..data_model import SpectralDataTimeSeries
from ..data_model import MeasurementDesc
from ..data_model import AuxilaryInfoDesc
from ..exoplanet_tools import SpectralModel
from ..utilities import _define_band_limits
from ..utilities import _define_rebin_weights
from ..utilities import _rebin_spectra

__all__ = ['directional_filters', 'create_edge_mask',
           'determine_optimal_filter', 'define_image_regions_to_be_filtered',
           'filter_image_cube', 'iterative_bad_pixel_flagging',
           'extract_spectrum', 'create_extraction_profile',
           'determine_relative_source_position',
           'warp_polar', 'highpass', '_log_polar_mapping',
           '_determine_relative_source_shift', 'register_telescope_movement',
           '_determine_relative_rotation_and_scale', '_derotate_image',
           '_pad_to_size', '_pad_region_of_interest_to_square',
           'correct_wavelength_for_source_movent',
           'rebin_to_common_wavelength_grid',
           'determine_center_of_light_posision',
           'determine_absolute_cross_dispersion_position',
           'combine_scan_samples', 'sigma_clip_data',
           'sigma_clip_data_cosmic', 'create_cleaned_dataset',
           'compressROI', 'compressSpectralTrace',
           'compressDataset', 'correct_initial_wavelength_shift',
           'renormalize_spatial_scans']


def _round_up_to_odd_integer(value):
    i = math.ceil(value)
    if i % 2 == 0:
        return i + 1
    return i


class Banana(Gaussian2D):
    """
    Modification of astrpy gaussian2D to get banana distribution.

    Notes
    -----
    https://tiao.io/post/building-probability-distributions-
    with-tensorflow-probability-bijector-api/
    """

    amplitude = Parameter(default=1)
    x_mean = Parameter(default=0)
    y_mean = Parameter(default=0)
    x_stddev = Parameter(default=1)
    y_stddev = Parameter(default=1)
    theta = Parameter(default=0.0)
    power = Parameter(default=1.0)
    sign = Parameter(default=1)

    def __init__(self, amplitude=amplitude.default, x_mean=x_mean.default,
                 y_mean=y_mean.default, x_stddev=None, y_stddev=None,
                 theta=None, cov_matrix=None, power=power.default,
                 sign=sign.default, **kwargs):
        if power is None:
            power = power.default
        if sign is None:
            sign = sign.default
        sign = np.sign(sign)
        super().__init__(
            amplitude=amplitude, x_mean=x_mean, y_mean=y_mean,
            x_stddev=x_stddev, y_stddev=y_stddev, theta=theta,
            cov_matrix=cov_matrix, power=power, sign=sign, **kwargs)

    @staticmethod
    def evaluate(x_in, y_in, amplitude, x_mean, y_mean, x_stddev, y_stddev,
                 theta, power, sign):
        """Two dimensional Gaussian function."""
        x = x_in
        y = y_in - sign*(np.abs(x_in)**power + 0.0)
        cost2 = np.cos(theta) ** 2
        sint2 = np.sin(theta) ** 2
        sin2t = np.sin(2. * theta)
        xstd2 = x_stddev ** 2
        ystd2 = y_stddev ** 2
        xdiff = x - x_mean
        ydiff = y - y_mean
        a = 0.5 * ((cost2 / xstd2) + (sint2 / ystd2))
        b = 0.5 * ((sin2t / xstd2) - (sin2t / ystd2))
        c = 0.5 * ((sint2 / xstd2) + (cost2 / ystd2))
        return amplitude * np.exp(-((a * xdiff ** 2) + (b * xdiff * ydiff) +
                                    (c * ydiff ** 2)))


class Banana2DKernel(Kernel2D):
    """
    Modification of astropy Gaussian2DKernel to get a babana shaped kernel.

    This class defines a banana shaped convolution kernel mimicking the shape
    of the dispersion pattern on the detector near the short and long
    wavelength ends.
    """

    _separable = True
    _is_bool = False

    def __init__(self, sigma, power=None, sign=None, **kwargs):
        self._model = Banana(1. / (2 * np.pi * sigma[0, 0] * sigma[1, 1]),
                             0., 0., cov_matrix=sigma, power=power, sign=sign)
        self._default_size = _round_up_to_odd_integer(
            8 * np.max([np.sqrt(sigma[0, 0]), np.sqrt(sigma[1, 1])]))
        super().__init__(**kwargs)
        self._truncation = np.abs(1. - self._array.sum())


def _define_covariance_matrix(x_stddev=None, y_stddev=None, theta=None):
    """
    Define covariance matrix.

    Define 2D covariance matrix based on standard deviation in x and y and
    rotation angle
    """
    if x_stddev is None:
        x_stddev = 1.0
    if y_stddev is None:
        y_stddev = 1.0
    if theta is None:
        theta = 0.0
    R = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta), np.cos(theta)]])
    S = np.array([[x_stddev, 0.0], [0.0, y_stddev]])
    sigma = np.linalg.multi_dot([R, S, S, R.T])
    return sigma


def directional_filters(return_valid_angle_range=False):
    """
    Directional filters for smoothing and filtering.

    These filters can be used in a Nagao&Matsuyama like edge preserving
    smoothing approach and are apropriate for dispersed spectra with a
    vertical dispersion direction. If the angle from vertical of the
    spectral trace of the dispersed light exceeds +- max(angle) radians,
    additional larger values need to be added to the angles list.

    Parameters
    ----------
    return_valid_angle_range : 'bool'
        optional, if True it returns the maximum angle range of the directional
        filters

    Returns
    -------
    nm_mask : numpy.ndarray of 'bool' type
        Array containing all oriented masks used for edge preserving smooting.
    maximum_angle_range : 'tuple' of 'float'
        If return_valid_angle_range is True, the maximum range of angles
        from vertical is returned

    Notes
    -----
    When adding kernels, make sure the maximum is in the central pixel
    """
    # note that the angels are in radians
    angles = [np.radians(0.0), np.radians(-1.5), np.radians(1.5),
              np.radians(-3.0), np.radians(3.0), np.radians(-4.5),
              np.radians(4.5), np.radians(-6.0), np.radians(6.0),
              np.radians(-9.0), np.radians(9.0), np.radians(-12.0),
              np.radians(12.0), np.radians(0.0), np.radians(90)-np.radians(60),
              np.radians(90)+np.radians(60), np.radians(90)-np.radians(60),
              np.radians(90)+np.radians(60)]
    if return_valid_angle_range:
        return (np.min(angles), np.max(angles))
    x_stddev = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                0.1, 0.1, 0.1, 0.1, 2.0, 0.1, 0.1, 0.1, 0.1]
    y_stddev = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
                3.0, 3.0, 3.0, 3.0, 2.0, 3.0, 3.0, 3.0, 3.0]
    sign = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1]
    power = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]

    x_kernel_size = 9
    y_kernel_size = 9

    Filters = np.zeros((x_kernel_size, y_kernel_size, len(angles)))

    for ik, (omega, xstd, ystd, p, s) in enumerate(zip(angles, x_stddev,
                                                       y_stddev, power, sign)):
        sigma = _define_covariance_matrix(xstd, ystd, omega)
        kernel = Banana2DKernel(sigma, x_size=x_kernel_size,
                                y_size=y_kernel_size, power=p, sign=s,
                                mode='oversample')
        kernel.normalize()
        Filters[..., ik] = kernel.array

    return Filters


def create_edge_mask(kernel, roi_mask):
    """
    Create an edge mask.

    Helper function for the optimal extraction task. This function
    creates an edge mask to mask all pixels for which the convolution
    kernel extends beyond the region of interest.

    Parameters
    ----------
    kernel : `array_like`
        Convolution kernel specific for a given instrument and observing
        mode, used in tasks such as replacing bad pixels and
        spectral extraction.
    roi_mask : `ndarray`
        Mask defining the region of interest from which the speectra will
        be extracted.

    Returns
    -------
    edge_mask : 'array_like'
        The edge mask based on the input kernel and roi_mask
    """
    # dilation_mask = np.ones(kernel.shape)
    dilation_mask = kernel
    edge_mask = (binary_dilation(roi_mask, structure=dilation_mask,
                                 border_value=True) ^ roi_mask)
    return edge_mask


def define_image_regions_to_be_filtered(ROI, filterShape):
    """
    Define image regions to be filtered.

    This function defines all pixels and corresponding sub-regions in the
    data cube to be be filtered.

    Parameters
    ----------
    ROI : 'ndarray' of 'float'
        Region of interest on the detector images
    fiterShape : 'tuple' of 'int'
        y (dispersion direction) and x (spatial direstion) size of the
        directional filters.

    Returns
    -------
    enumerated_sub_regions : 'list'
        enumerated definition of all regions in the spectral image cube and
        corresponding regions of the direction filter.
    """
    # find all images indices of the pixels of interest which are not flagged
    # in the region of interest
    y, x = np.where(ROI == False)
    indices_poi = [(yidx, xidx, ...) for (yidx, xidx) in zip(y, x)]

    # support boundary values used to define the sub-array ranges
    xmax = ROI.shape[1]
    xl = (filterShape[1]-1)//2+1
    xs = xl-1
    ymax = ROI.shape[0]
    yl = (filterShape[0]-1)//2+1
    ys = yl-1

    # all regions around the pixels not flaged in the roi with a size equal
    # to the filter size
    image_sub_regions = [(slice(np.max([0, yidx-ys]),
                                np.min([ymax, yidx+yl]), None),
                          slice(np.max([0, xidx-xs]),
                                np.min([xmax, xidx+xl]), None),
                          tidx) for (yidx, xidx, tidx) in indices_poi]

    # region of the filter used in the region around the pixels defined in roi
    filter_sub_regions = [(slice(np.max([0, yidx-ys])-yidx+ys,
                                 np.min([ymax, yidx+yl])-yidx+yl-1, None),
                           slice(np.max([0, xidx-xs])-xidx+xs,
                                 np.min([xmax, xidx+xl])-xidx+xl-1, None),
                           tidx) for (yidx, xidx, tidx) in indices_poi]

    # defines all regions in image and corresponding region of filter
    sub_regions = list(zip(image_sub_regions, filter_sub_regions))
    # eneumerated definition of all regions in image and corresponding region
    # of filter
    enumerated_sub_regions = \
        [(i, sub, poi) for (i, (sub, poi)) in enumerate(zip(sub_regions,
                                                            indices_poi))]
    return enumerated_sub_regions


def determine_optimal_filter(ImageCube, Filters, ROIcube, selector):
    """
    Determine optimal fileter.

    Determine the optimal Filter for the image cube using a procedure similar
    to Nagao & Matsuyama edge preserving filtering.

    Parameters
    ----------
    ImageCube : 'ndarray'
        Cube of Spectral images.
    Filters : 'ndarray'
        Cube of directional filters
    ROIcube : 'ndarray' of  'bool'
        Region of Interests for the input ImageCube
    selector : 'list'
        list containing all relevant information for each pixel within the
        ROI to select the sub region in the image cube and filter cube on which
        the filtering will be applied.

    Returns
    -------
    selectorNumber : 'int'
       id number of the selector (pixel)
    optimum_filter_index : 'ndarray' of 'int'
        array containing the optimal filter number for each pixel within the
        ROI in the image cube.
    SubImageOptimalMean : 'ndarray'
        Optimal mean for each sub image in the sub image cube defined by the
        selector.
    SubImageOptimaVariance : 'ndarray'
        Variance for each sub image in the sub image cube defined by the
        selector using the optimal Filter.
    """
    selectorNumber = selector[0]
    SubFilters = Filters[selector[1][1]]
    SubImageCube = ImageCube[selector[1][0]].data
    SubImageMask = ImageCube[selector[1][0]].mask
    mask = ROIcube[selector[1][0]]
    mask = SubImageMask | mask
    maskedCube = np.ma.array(SubImageCube, mask=mask)

    SubImageVariance = np.ma.zeros((Filters.shape[-1], SubImageCube.shape[-1]))
    SubImageMean = np.ma.zeros(SubImageVariance.shape)
    SubImageSquaredMean = np.ma.zeros(SubImageVariance.shape)

    for j, filterKernel in enumerate(SubFilters.T):
        weights = np.tile(filterKernel,
                          (maskedCube.shape[-1], 1, 1)).T
        SubImageMean[j, :] = np.ma.average(maskedCube, weights=weights,
                                           axis=(0, 1))
        SubImageSquaredMean[j, :] = np.ma.average(maskedCube**2,
                                                  weights=weights, axis=(0, 1))
        SubImageVariance[j, :] = (SubImageSquaredMean[j, :] -
                                  SubImageMean[j, :]**2)

    optimum_filter_index = np.ma.argmin(SubImageVariance, axis=0)
    idx = [optimum_filter_index, np.arange(len(optimum_filter_index))]
    SubImageOptimalMean = SubImageMean[tuple(idx)]
    SubImageOptimaVariance = SubImageVariance[tuple(idx)]
    return selectorNumber, optimum_filter_index, SubImageOptimalMean, \
        SubImageOptimaVariance


def chunks(lst, n):
    """
    Yield successive n-sized chunks from lst.

    Parameters
    ----------
    lst : 'list'
        Input list
    n : 'integer'
        Chunck size
    """
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


@ray.remote
def split_work(ImageCube, Filters, ROIcube, selector):
    """
    Split work.

    Ray wrapper to be able to devide data in chunks which
    can be filetered in parallel.

    Parameters
    ----------
    ImageCube : TYPE
        DESCRIPTION.
    Filters : TYPE
        DESCRIPTION.
    ROIcube : TYPE
        DESCRIPTION.
    selector : TYPE
        DESCRIPTION.

    Returns
    -------
    list
        DESCRIPTION.

    """
    return [determine_optimal_filter(ImageCube, Filters, ROIcube, x)
            for x in selector]


def filter_image_cube(data_in, Filters, ROIcube, enumeratedSubRegions,
                      useMultiProcesses=True, ReturnCleanedData=True):
    """
    Filter image cube.

    This routine filters in input data clube using an optimal choise of a
    directional filter and returns the filtered (smoothed) data. Optionally it
    also returns a cleaned dataset, where the masked data is replaced by the
    fitered data.

    Parameters
    ----------
    data_in : 'np.ndarray' of 'float'
        Input spectral data cube.

    Filters : 'ndarray' of 'float'
        Directional filter kernels.
    ROIcube : 'ndarray' of 'bool'
        Region of interest for each spectral image.
    enumeratedSubRegions : 'list'
        List containing all information to define the regions around the pixels
        within the ROI used in the filtering. This list is created with the
        define_image_regions_to_be_filtered function.
    useMultiProcesses : 'bool' (optional|True)
        Flag to choose to use the multiprocessing module or not.
    ReturnCleanedData : 'bool' (optional|True)
        Flag to choose to return a cleaned data set or not.

    Returns
    -------
    optimalFilterIndex : 'MaskedArray' of 'int'
        Index number of the optimal filter for each pixel in the input data.
    filteredImage : 'MaskedArray' of 'float'
        The filtered data set.
    filteredImageVariance : 'MaskedArray' of 'float'
        The variance of the filtered data set.
    cleanedData : 'MaskedArray' of 'float'
        The cleaned data set.
    """
    optimalFilterIndex = np.ma.zeros(data_in.shape, dtype='int')
    optimalFilterIndex.mask = ROIcube
    filteredImage = np.ma.zeros(data_in.shape)
    filteredImage.mask = ROIcube
    filteredImageVariance = np.ma.zeros(data_in.shape)
    filteredImageVariance.mask = ROIcube

    indices_poi = [poi for (_, _, poi) in enumeratedSubRegions]

    if not useMultiProcesses:
        # create new function with all fixed inout variables fixed.
        func = partial(determine_optimal_filter, data_in, Filters, ROIcube)
        filter_find_iterator = map(func, enumeratedSubRegions)
        for j in tqdm(filter_find_iterator, total=len(enumeratedSubRegions),
                      dynamic_ncols=True):
            optimalFilterIndex[indices_poi[j[0]]] = j[1]
            filteredImage[indices_poi[j[0]]] = j[2]
            filteredImageVariance[indices_poi[j[0]]] = j[3]
    else:
        ncpus = int(ray.cluster_resources()['CPU'])
        chunksize = len(enumeratedSubRegions)//ncpus + 1
        while chunksize > 256:
            chunksize = chunksize//ncpus + 1

        data_id = ray.put(data_in)
        filters_id = ray.put(Filters)
        roi_id = ray.put(ROIcube)

        result_ids = \
            [split_work.remote(data_id, filters_id, roi_id, x) for x in
             list(chunks(enumeratedSubRegions, chunksize))]

        pbar = tqdm(total=len(result_ids), dynamic_ncols=True,
                    desc='Filtering image cube')
        while len(result_ids):
            done_id, result_ids = ray.wait(result_ids)
            k = ray.get(done_id[0])
            for j in k:
                optimalFilterIndex[indices_poi[j[0]]] = j[1]
                filteredImage[indices_poi[j[0]]] = j[2]
                filteredImageVariance[indices_poi[j[0]]] = j[3]
            pbar.update(1)
        pbar.close()

    if ReturnCleanedData:
        cleanedData = data_in.copy()
        cleanedData[data_in.mask] = filteredImage[data_in.mask]
        cleanedData.mask = cleanedData.mask | ROIcube
        return optimalFilterIndex, filteredImage, filteredImageVariance, \
            cleanedData
    else:
        return optimalFilterIndex, filteredImage, filteredImageVariance


def iterative_bad_pixel_flagging(dataset, ROIcube, Filters,
                                 enumeratedSubRegions, sigmaLimit=4.0,
                                 maxNumberOfIterations=12,
                                 fractionalAcceptanceLimit=0.005,
                                 useMultiProcesses=True,
                                 maxNumberOfCPUs=2):
    """
    Flag bad pixels.

    This routine flags the bad pixels found in the input dataset and creates
    a cleaned dataset.

    Parameters
    ----------
    dataset : 'SpectralDataTimeSeries'
    ROIcube : 'ndarray' of 'bool'
        Region of interest for each spectral image.
    Filters : 'ndarray' of 'float'
        Directional filter kernels.
    enumeratedSubRegions : 'list'
        List containing all information to define the regions around the pixels
        within the ROI used in the filtering. This list is created with the
        define_image_regions_to_be_filtered function.
    sigmaLimit : 'float' (optional|4.0)
        Standard diviation limit used to identify bad pixels
    maxNumberOfIterations : 'int' (optional|12)
        The maximum number of iterations used in bad pixel masking
    fractionalAcceptanceLimit : 'float' (optional|0.005)
        Fractional number of bad pixels still found in the dataset below
        which the iteration can be terminated.
    useMultiProcesses : 'bool' (optional|True)
        Use multiprocessing or not.
    max_number_of_cpus : 'int' (optional|12)
        Maximum number of CPU's which can be used.

    Returns
    -------
    dataset : 'SpectralDataTimeSeries'
        Input data set containing the observed spectral time series. After
        completinion of this function, the dataset is returned with an
        updated mask with all diviating pixels flaged. Tho indicate this the
        flag isSigmaCliped=True is set.
    filteredDataset : 'SpectralDataTimeSeries'
        The optimal filterd dataset. Includied in this dataset is
        the optimal Filter Index, i.e. the index indicating the used
        optimal filter for each pixel.  To indicate that this is a
        filtered dataset the flag isFilteredData=True is set.
    cleanedDataset
        The cleaned dataset. To indicte that this is a cleaned dataset,
        the isCleanedData=True is set.

    Notes
    -----
    In the current version no double progress bar is used as in some IDE's
    double progress bars do not work properly.
    """
    acceptanceLimit = int(len(enumeratedSubRegions)*fractionalAcceptanceLimit)
    initialData = dataset.return_masked_array('data').copy()

    if useMultiProcesses:
        mem = virtual_memory()
        ncpu = int(np.min([maxNumberOfCPUs,
                           np.max([1, cpu_count(logical=True)//2])
                           ])
                   )
        mem_store = np.max([int(initialData.nbytes*3.0), int(1.1*78643200)])
        mem_workers = np.max([int(initialData.nbytes*5.0), int(1.1*52428800)])
        if mem.available < (mem_store+mem_workers):
            warnings.warn("WARNING: Not enough memory for Ray to start. "
                          "Required free memory: {} bytes "
                          "Available: {} bytes".
                          format(mem_store+mem_workers, mem.available))
        # ray.disconnect()
#        ray.init(num_cpus=ncpu, object_store_memory=mem_store,
#                 memory=mem_workers)
# bug fix
        ray.init(num_cpus=ncpu, object_store_memory=mem_store)

    numberOfFlaggedPixels = np.sum(~ROIcube & dataset.mask)
    tqdm.write('Initial bad pixel flagging. '
               '# of pixels masked: {}'.format(numberOfFlaggedPixels))

    optimalFilterIndex, filteredImage, filteredImageVariance, cleanedData = \
        filter_image_cube(initialData, Filters, ROIcube, enumeratedSubRegions,
                          useMultiProcesses=useMultiProcesses,
                          ReturnCleanedData=True)

    iiteration = 1
    while (numberOfFlaggedPixels > acceptanceLimit) & \
            (iiteration <= maxNumberOfIterations):
        mask = (cleanedData-filteredImage)**2 > (sigmaLimit *
                                                 filteredImageVariance)
        numberOfFlaggedPixels = np.sum(~ROIcube & mask)
        tqdm.write('Iteration #{} for bad pixel flagging. '
                   '# of pixels masked: {}'.format(iiteration,
                                                   numberOfFlaggedPixels))
        cleanedData.mask = cleanedData.mask | mask
        dataset.mask = dataset.mask | (~ROIcube & mask)
        (optimalFilterIndex, filteredImage, filteredImageVariance,
         cleanedData) = filter_image_cube(cleanedData, Filters, ROIcube,
                                          enumeratedSubRegions,
                                          useMultiProcesses=useMultiProcesses,
                                          ReturnCleanedData=True)
        iiteration += 1
    if (iiteration > maxNumberOfIterations) & \
            (numberOfFlaggedPixels < acceptanceLimit):
        warnings.warn("Iteration not converged in "
                      "iterative_bad_pixel_flagging. {} mask values not "
                      "converged. An increase of the maximum number of "
                      "iteration steps might be advisable.".
                      format(numberOfFlaggedPixels-acceptanceLimit))

    #  ray.disconnect()
    ray.shutdown()

    cleanedUncertainty = np.ma.array(dataset.uncertainty.data.value.copy(),
                                     mask=dataset.mask.copy())
    cleanedUncertainty[cleanedUncertainty.mask] = \
        np.sqrt(filteredImageVariance[cleanedUncertainty.mask])
    cleanedUncertainty.mask = filteredImage.mask

    ndim = dataset.data.ndim
    selection = tuple((ndim-1)*[0]+[Ellipsis])

    filteredDataset = \
        SpectralDataTimeSeries(wavelength=dataset.wavelength,
                               wavelength_unit=dataset.wavelength_unit,
                               data=filteredImage,
                               data_unit=dataset.data_unit,
                               time=dataset.time.data.value[selection],
                               time_unit=dataset.time_unit,
                               time_bjd=dataset.time_bjd.data.value[selection],
                               time_bjd_unit=dataset.time_bjd_unit,
                               uncertainty=np.sqrt(filteredImageVariance),
                               target_name=dataset.target_name,
                               dataProduct=dataset.dataProduct,
                               dataFiles=dataset.dataFiles,
                               isFilteredData=True)
    filteredDataset.optimalFilterIndex = optimalFilterIndex

    cleanedDataset = \
        SpectralDataTimeSeries(wavelength=dataset.wavelength,
                               wavelength_unit=dataset.wavelength_unit,
                               data=cleanedData,
                               data_unit=dataset.data_unit,
                               time=dataset.time.data.value[selection],
                               time_unit=dataset.time_unit,
                               time_bjd=dataset.time_bjd.data.value[selection],
                               time_bjd_unit=dataset.time_bjd_unit,
                               uncertainty=cleanedUncertainty,
                               target_name=dataset.target_name,
                               dataProduct=dataset.dataProduct,
                               dataFiles=dataset.dataFiles,
                               isCleanedData=True)

    dataset.isSigmaCliped = True

    return (dataset, filteredDataset, cleanedDataset)


def create_extraction_profile(fiteredSpectralDataset, ROI=None):
    """
    Create an extraction profile.

    Parameters
    ----------
    fiteredSpectralDataset : 'SpectralDataTimeSeries'
        Input filtered dataset on which the extraction profile is based
    ROI : 'ndarray' of 'bool', optional

    Returns
    -------
    spectralExtractionProfile : 'MaskedArray'
    """
    fiteredSpectralData = fiteredSpectralDataset.return_masked_array('data')
    if ROI is None:
        dataUse = fiteredSpectralData
    else:
        # newMask = np.logical_or(fiteredSpectralData.mask, ROI)
        newMask = fiteredSpectralData.mask | ROI
        dataUse = np.ma.array(fiteredSpectralData.data, mask=newMask)

    spectralExtractionProfile = dataUse/np.sum(dataUse, axis=1, keepdims=True)

    return spectralExtractionProfile


def extract_spectrum(dataset, ROICube, extractionProfile=None, optimal=False,
                     verbose=False, verboseSaveFile=None):
    """
    Extract 1D spectra.

    This routine extracts the spectrum both optimally as well as using an
    aperture.

    Parameters
    ----------
    dataset : 'SpectralDataTimeSeries'
        Input spectral dataset
    ROIcube : 'ndarray' of 'bool'
        Region of interest for each spectral image.
    extractionProfile : 'MaskedArray', optinonal
        Normalized extraction profile. Has to be set if optimal=True
    verbose : 'bool', optional
        If true diagnostic plots will be generated.
    verboseSaveFile : 'str', optional
        If not None, verbose output will be saved to the specified file.

    Returns
    -------
    extracted1dDataset : SpectralDataTimeSeries'
        Extracted 1D spectral timeseries dataset

    Raises
    ------
    ValueError
    """
    data = dataset.return_masked_array('data').copy()
    variance = (dataset.return_masked_array('uncertainty').copy())**2
    wavelength = dataset.return_masked_array('wavelength').copy()
    if optimal:
        mask = (~data.mask).astype(int) * (~ROICube).astype(int)
        extractedSpectra = \
            np.ma.sum(mask*extractionProfile*data/variance, axis=1) / \
            np.ma.sum(mask*extractionProfile**2/variance, axis=1)
        varianceExtractedSpectra = \
            np.ma.sum(mask*extractionProfile, axis=1) / \
            np.ma.sum(mask*extractionProfile**2/variance, axis=1)
        wavelengthExtractedSpectrum = \
            np.ma.sum(mask*extractionProfile**2*wavelength/variance, axis=1) /\
            np.ma.sum(mask*extractionProfile**2/variance, axis=1)
    else:
        mask = (~ROICube).astype(int)
        extractedSpectra = \
            np.ma.sum(mask*data, axis=1)
        varianceExtractedSpectra = np.ma.sum(mask*variance, axis=1)
        wavelengthExtractedSpectrum = np.ma.sum(mask*wavelength, axis=1) /\
            np.ma.sum(mask, axis=1)

    uncertaintyExtractedSpectra = np.sqrt(varianceExtractedSpectra)
    # As the canculations are done independently, the masks might be different
    # which causes problems when compressing rows. Make sure the final
    # mask is identical. This fixes bug #82
    new_mask = np.ma.logical_or(uncertaintyExtractedSpectra.mask,
                                extractedSpectra.mask)
    new_mask = np.ma.logical_or(new_mask, wavelengthExtractedSpectrum.mask)
    extractedSpectra.mask = new_mask
    uncertaintyExtractedSpectra.mask = new_mask
    wavelengthExtractedSpectrum.mask = new_mask
    extractedSpectra = np.ma.asanyarray(np.ma.compress_rows(extractedSpectra))
    uncertaintyExtractedSpectra = \
        np.ma.asanyarray(np.ma.compress_rows(uncertaintyExtractedSpectra))
    wavelengthExtractedSpectrum = \
        np.ma.asanyarray(np.ma.compress_rows(wavelengthExtractedSpectrum))

    dataProductOld = dataset.dataProduct
    if optimal:
        dataProduct = 'COE'
    else:
        dataProduct = 'CAE'
    dataFilesOld = dataset.dataFiles
    dataFiles = [fname.split("/")[-1].replace(dataProductOld, dataProduct)
                 for fname in dataFilesOld]
    extracted1dDataset = \
        SpectralDataTimeSeries(wavelength=wavelengthExtractedSpectrum,
                               wavelength_unit=dataset.wavelength_unit,
                               data=extractedSpectra,
                               data_unit=dataset.data_unit,
                               time=dataset.time.data.value[0, 0, :],
                               time_unit=dataset.time_unit,
                               time_bjd=dataset.time_bjd.data.value[0, 0, :],
                               time_bjd_unit=dataset.time_bjd_unit,
                               uncertainty=uncertaintyExtractedSpectra,
                               target_name=dataset.target_name,
                               dataProduct=dataProduct,
                               dataFiles=dataFiles,
                               isExtractedSpectra=True)
    if verbose:
        sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
        sns.set_style("white", {"xtick.bottom": True, "ytick.left": True})
        fig, ax0 = plt.subplots(figsize=(6, 6), nrows=1, ncols=1)
        for iwave in range(1, 8):
            ax0.plot(extractedSpectra[iwave, :])
        ax0.set_title('Extracted 1D spectral timeseries')
        ax0.set_xlabel('Integration #')
        ax0.set_ylabel('Flux [{}]'.format(dataset.data_unit))
        fig.subplots_adjust(hspace=0.3)
        fig.subplots_adjust(wspace=0.45)
        plt.show()
        if verboseSaveFile is not None:
            fig.savefig(verboseSaveFile, bbox_inches='tight')

    return extracted1dDataset


def determine_relative_source_position(spectralImageCube, ROICube,
                                       refIntegration,
                                       upsampleFactor=111,
                                       AngleOversampling=2):
    """
    Determine the shift of the spectra relative to the first integration.

    Parameters
    ----------
    spectralImageCube : 'ndarray'
        Input spectral image data cube. Fist dimintion is dispersion direction,
        second dimintion is cross dispersion direction and the last dimension
        is time.
    ROICube : 'ndarray' of 'bool'
        Cube containing the Region of interest for each integration.
        If not given, it is assumed that the mask of the spectralImageCube
        contains the region of interest.
    refIntegration : 'int'
        Index number of of integration to be taken as reference.
    upsampleFactor : 'int', optional
        integer factor by which to upsample image for FFT analysis to get
        sub-pixel accuracy. Default value is 111
    AngleOversampling : 'int', optional
        Oversampling factor for angle determination, Default value 2

    Returns
    -------
    relativeSourcePosition : 'collections.OrderedDict'
        relative rotation angle, scaling and x and y position as a
        function of time.

    Raises
    ------
    ValueError
        In case refIntegration exceeds number of integrations

    Notes
    -----
    Note that for sub-pixel registration to work correctly it should
    be performed on cleaned data i.e. bad pixels have been identified and
    corrected using the tools in this module.
    """
    nintegrations = spectralImageCube.shape[-1]
    refIntegration = int(refIntegration)
    if (refIntegration > nintegrations-1) | (refIntegration < 0):
        raise ValueError("Index number of reference integration exceeds \
                         limits. Aborting position determination")
    upsampleFactor = np.max([1, int(upsampleFactor)])
    AngleOversampling = np.max([1, int(AngleOversampling)])

    ImageCube = spectralImageCube.copy()
    refImage = ImageCube[..., refIntegration]

    if ROICube is None:
        ROICube = np.zeros((ImageCube.shape), dtype=bool)
    ROIref = ROICube[..., refIntegration]

    # First determine the rotation and scaling
    relativeAngle = np.zeros((nintegrations))
    relativeScale = np.ones((nintegrations))

    for it in range(1, nintegrations):
        relativeAngle[it], relativeScale[it] = \
         _determine_relative_rotation_and_scale(
             refImage, ROIref,
             ImageCube[..., it],
             ROICube[..., it],
             upsampleFactor=upsampleFactor,
             AngleOversampling=AngleOversampling)

    # Second, determine the shift
    yshift = np.zeros((nintegrations))
    xshift = np.zeros((nintegrations))

    for it, image in enumerate(ImageCube.T):
        if not np.allclose(relativeAngle[it], 0.0):
            derotateRefImage = _derotate_image(refImage, 0.0, ROI=ROIref,
                                               order=3)
            derotatedImage = _derotate_image(image.T, relativeAngle[it],
                                             ROI=ROICube[..., it], order=3)
            derotatedROIref = np.zeros_like(derotateRefImage).astype(bool)
            derotatedROI = np.zeros_like(derotatedImage).astype(bool)
        else:
            derotateRefImage = refImage
            derotatedImage = image.T
            derotatedROIref = ROIref
            derotatedROI = ROICube[..., it]

        shift = _determine_relative_source_shift(derotateRefImage,
                                                 derotatedImage,
                                                 referenceROI=derotatedROIref,
                                                 ROI=derotatedROI,
                                                 upsampleFactor=upsampleFactor)
        yshift[it] = shift[0]
        xshift[it] = shift[1]
    relativeSourcePosition = \
        collections.OrderedDict(relativeAngle=relativeAngle,
                                relativeScale=relativeScale,
                                cross_disp_shift=xshift,
                                disp_shift=yshift)
    return relativeSourcePosition


@ray.remote
def ray_determine_relative_source_position(spectralImageCube, ROICube,
                                           refIntegration, pba,
                                           upsampleFactor=111,
                                           AngleOversampling=2):
    """
    Ray wrapper for determine_relative_source_position.

    Parameters
    ----------
    spectralImageCube : 'ndarray'
        Input spectral image data cube. Fist dimintion is dispersion direction,
        second dimintion is cross dispersion direction and the last dimension
        is time.
    ROICube : 'ndarray' of 'bool'
        Cube containing the Region of interest for each integration.
        If not given, it is assumed that the mask of the spectralImageCube
        contains the region of interest.
    refIntegration : 'int'
        Index number of of integration to be taken as reference.
    upsampleFactor : 'int', optional
        integer factor by which to upsample image for FFT analysis to get
        sub-pixel accuracy. Default value is 111
    AngleOversampling : 'int', optional
        Oversampling factor for angle determination, Default value 2

    Returns
    -------
    movement : 'collections.OrderedDict'
        relative rotation angle, scaling and x and y position as a
        function of time.

    """
    movement = determine_relative_source_position(
        spectralImageCube, ROICube, refIntegration,
        upsampleFactor=upsampleFactor,
        AngleOversampling=AngleOversampling)
    pba.update.remote(1)
    return movement


def _determine_relative_source_shift(reference_image, image,
                                     referenceROI=None, ROI=None,
                                     upsampleFactor=111, space='real'):
    """
    Determine the relative shift of the spectral images.

    This routine determine the relative shift between a reference spectral
    image and another spectral image.

    Parameters
    ----------
    reference_image : 'ndarray or np.ma.MaskedArray' of 'float'
        Reference spectral image
    image : 'ndarray or np.ma.MaskedArray' of 'float'
        Spectral image
    referenceROI : ndarray' of 'bool', optional
    ROI : 'ndarray' of 'bool', optional
    upsampleFactor : 'int', optional
        Default value is 111
    space : 'str', optional
        Default value is 'real'
    Returns
    -------
    relativeImageShiftY
        relative shift compared to the reference image in the dispersion
        direction of the light (from top to bottom, shortest wavelength should
        be at row 0. Note that this shift is defined such that shifting a
        spectral image by this amound will place the trace at the exact same
        position as that of the reference image
    relativeImageShiftX
        relative shift compared to the reference image in the cross-dispersion
        direction of the light (from top to bottom, shortest wavelength should
        be at row 0. Note that this shift is defined such that shifting a
        spectral image by this amound will place the trace at the exact same
        position as that of the reference image.
    """
    ref_im = _pad_region_of_interest_to_square(reference_image, referenceROI)
    im = _pad_region_of_interest_to_square(image, ROI)

    # convolve with gaussian with sigma of 1 pixel to esnure that undersampled
    # spectra are properly registered.
    kernel = Gaussian2DKernel(1.0)
    ref_im = convolve(ref_im, kernel, boundary='extend')
    im = convolve(im, kernel, boundary='extend')

    # subpixel precision by oversampling image by upsampleFactor
    # returns shift, error and phase difference
    shift, _, _ = \
        phase_cross_correlation(ref_im, im, upsample_factor=upsampleFactor,
                                space=space)
    relativeImageShiftY = -shift[0]
    relativeImageShiftX = -shift[1]
    return relativeImageShiftY, relativeImageShiftX


def _determine_relative_rotation_and_scale(reference_image, referenceROI,
                                           image, ROI,
                                           upsampleFactor=111,
                                           AngleOversampling=2):
    """
    Determine rotation and scalng changes.

    This routine determines the relative rotation and scale change between
    an reference spectral image and another spectral image.

    Parameters
    ----------
    reference_image : 'ndarray or np.ma.MaskedArray' of 'float'
        Reference image
    referenceROI : 'ndarray' of 'float'
    image : 'ndarray or np.ma.MaskedArray' of 'float'
        Image for which the rotation and translation relative to the reference
        image will be determined
    ROI : 'ndarray' of 'float'
        Region of interest.
    upsampleFactor : 'int', optional
        Upsampling factor of FFT image used to determine sub-pixel shift.
        By default set to 111.
    AngleOversampling : 'int', optional
        Upsampling factor of the FFT image in polar coordinates for the
        determination of sub-degree rotation. Set by default to 2.

    Returns
    -------
    relative_rotation
        Relative rotation angle in degrees. The angle is defined such that the
        image needs to be rotated by this angle to have the same orientation
        as the reference spectral image
    relative_scaling
        Relative image scaling
    """
    AngleOversampling = int(AngleOversampling)
    nAngles = 360
    NeededImageSize = 2*AngleOversampling*nAngles
    ref_im = _pad_region_of_interest_to_square(reference_image, referenceROI)
    ref_im = _pad_to_size(ref_im, NeededImageSize, NeededImageSize)
    im = _pad_region_of_interest_to_square(image, ROI)
    im = _pad_to_size(im, NeededImageSize, NeededImageSize)

    # convolve with gaussian with sigma of 1 pixel to esnure that undersampled
    # spectra are properly registered.
    kernel = Gaussian2DKernel(1.0)
    ref_im = convolve(ref_im, kernel, boundary='extend')
    im = convolve(im, kernel, boundary='extend')

    h = np.hanning(im.shape[0])
    han2d = np.outer(h, h)  # 2D Hanning window

    fft_ref_im = np.abs(np.fft.fftshift(np.fft.fftn(ref_im*han2d)))**2
    fft_im = np.abs(np.fft.fftshift(np.fft.fftn(im*han2d)))**2

    h, w = fft_ref_im.shape
    radius = 0.8*np.min([w/2, h/2])

    hpf = highpass((h, w))

    fft_ref_im_filtered = fft_ref_im * hpf
    warped_fft_ref_im = warp_polar(fft_ref_im_filtered, scaling='log',
                                   radius=radius, output_shape=None,
                                   multichannel=None,
                                   AngleOversampling=AngleOversampling)
    fft_im_filtered = fft_im * hpf
    warped_fft_im = warp_polar(fft_im_filtered, scaling='log', radius=radius,
                               output_shape=None, multichannel=None,
                               AngleOversampling=AngleOversampling)

    tparams = phase_cross_correlation(warped_fft_ref_im, warped_fft_im,
                                      upsample_factor=upsampleFactor,
                                      space='real')

    shifts = tparams[0]
    # calculate rotation
    # note, only look for angles between +- 90 degrees,
    # remove any flip of 180 degrees due to search
    shiftr, shiftc = shifts[:2]
    shiftr = shiftr/AngleOversampling
    if shiftr > 90.0:
        shiftr = shiftr-180.0
    if shiftr < -90.0:
        shiftr = shiftr+180.0
    relative_rotation = -shiftr

    # Calculate scale factor from translation
    klog = radius / np.log(radius)
    relative_scaling = 1 / (np.exp(shiftc / klog))

    return relative_rotation, relative_scaling


def _derotate_image(image, angle, ROI=None, order=3):
    """
    Derotate image.

    Parameters
    ----------
    image : '2-D ndarray' of 'float'
        Input image to be de-rotated by 'angle' degrees.
    ROI : '2-D ndarray' of 'bool'
        Region of interest (default None)
    angle : 'float'
        Rotaton angle in degrees
    order : 'int'
        Order of the used interpolation in the rotation function of the
        skimage package.

    Returns
    -------
    derotatedImage : '2-D ndarray' of 'float'
        The zero padded and derotated image.
    """
    h, w = image.shape
    NeededImageSize = np.int(np.sqrt(h**2 + w**2))

    im = _pad_region_of_interest_to_square(image, ROI)
    im = _pad_to_size(im, NeededImageSize, NeededImageSize)
    derotatedImage = rotate(im, angle, order=order)

    return derotatedImage


def _pad_region_of_interest_to_square(image, ROI=None):
    """
    Pad ROI to square.

    Zero pad the extracted Region Of Interest of a larger image such that the
    resulting image is squared.

    Parameters
    ----------
    image : '2-D ndarray' of 'float'
        Input image to be de-rotated by 'angle' degrees.
    ROI : '2-D ndarray' of 'bool'
        Region of interest (default None)

    Returns
    -------
    padded_image : '2-D ndarray' of 'float'
    """
    if ROI is not None:
        label_im, _ = ndimage.label(ROI)
    elif isinstance(image, np.ma.core.MaskedArray):
        label_im, _ = ndimage.label(image.mask)
    else:
        raise AttributeError("For image 0 padding either use MaskedArray as \
                              input or provide ROI. Aborting 0 padding")
    slice_y, slice_x = ndimage.find_objects(label_im != 1)[0]
    padded_image = image[slice_y, slice_x]

    if isinstance(image, np.ma.core.MaskedArray):
        padded_image.set_fill_value(0.0)
        padded_image = padded_image.filled()

    h, w = padded_image.shape
    if h == w:
        return padded_image

    im_size = np.max([h, w])
    delta_h = im_size - h
    delta_w = im_size - w
    padding = ((delta_h//2, delta_h-(delta_h//2)),
               (delta_w//2, delta_w-(delta_w//2)))
    padded_image = np.pad(padded_image,
                          padding, 'constant', constant_values=(0))

    return padded_image


def _pad_to_size(image, h, w):
    """
    Zero pad the input image to an image of hight h and width w.

    Parameters
    ----------
    image : '2-D ndarray' of 'float'
        Input image to be zero-padded to size (h, w).
    h : 'int'
        Hight (number of rows) of output image.
    w : 'int'
        Width (number of columns) of output image.

    Returns
    -------
    padded_image : '2-D ndarray' of 'float'
    """
    padded_image = image.copy()
    if isinstance(padded_image, np.ma.core.MaskedArray):
        padded_image.set_fill_value(0.0)
        padded_image = padded_image.filled()
    h_image, w_image = padded_image.shape
    npad_h = np.max([1, (h-h_image)//2])
    npad_w = np.max([1, (w-w_image)//2])
    padding = ((npad_h, npad_h), (npad_w, npad_w))
    padded_image = np.pad(padded_image,
                          padding, 'constant', constant_values=(0))
    return padded_image


def _log_polar_mapping(output_coords, k_angle, k_radius, center):
    """
    Inverse mapping function to convert from cartesion to polar coordinates.

    Parameters
    ----------
    output_coords : ndarray
        `(M, 2)` array of `(col, row)` coordinates in the output image
    k_angle : float
        Scaling factor that relates the intended number of rows in the output
        image to angle: ``k_angle = nrows / (2 * np.pi)``
    k_radius : float
        Scaling factor that relates the radius of the circle bounding the
        area to be transformed to the intended number of columns in the output
        image: ``k_radius = width / np.log(radius)``
    center : tuple (row, col)
        Coordinates that represent the center of the circle that bounds the
        area to be transformed in an input image.

    Returns
    -------
    coords : ndarray
        `(M, 2)` array of `(col, row)` coordinates in the input image that
        correspond to the `output_coords` given as input.
    """
    angle = output_coords[:, 1] / k_angle
    rr = ((np.exp(output_coords[:, 0] / k_radius)) * np.sin(angle)) + center[0]
    cc = ((np.exp(output_coords[:, 0] / k_radius)) * np.cos(angle)) + center[1]
    coords = np.column_stack((cc, rr))
    return coords


def _linear_polar_mapping(output_coords, k_angle, k_radius, center):
    """
    Inverse mapping function to convert from cartesion to polar coordinates.

    Parameters
    ----------
    output_coords : ndarray
        `(M, 2)` array of `(col, row)` coordinates in the output image
    k_angle : float
        Scaling factor that relates the intended number of rows in the output
        image to angle: ``k_angle = nrows / (2 * np.pi)``
    k_radius : float
        Scaling factor that relates the radius of the circle bounding the
        area to be transformed to the intended number of columns in the output
        image: ``k_radius = ncols / radius``
    center : tuple (row, col)
        Coordinates that represent the center of the circle that bounds the
        area to be transformed in an input image.

    Returns
    -------
    coords : ndarray
        `(M, 2)` array of `(col, row)` coordinates in the input image that
        correspond to the `output_coords` given as input.
    """
    angle = output_coords[:, 1] / k_angle
    rr = ((output_coords[:, 0] / k_radius) * np.sin(angle)) + center[0]
    cc = ((output_coords[:, 0] / k_radius) * np.cos(angle)) + center[1]
    coords = np.column_stack((cc, rr))
    return coords


def warp_polar(image, center=None, *, radius=None, AngleOversampling=1,
               output_shape=None, scaling='linear', multichannel=False,
               **kwargs):
    """
    Remap image to polor or log-polar coordinates space.

    Parameters
    ----------
    image : ndarray
        Input image. Only 2-D arrays are accepted by default. If
        `multichannel=True`, 3-D arrays are accepted and the last axis is
        interpreted as multiple channels.
    center : tuple (row, col), optional
        Point in image that represents the center of the transformation (i.e.,
        the origin in cartesian space). Values can be of type `float`.
        If no value is given, the center is assumed to be the center point
        of the image.
    radius : float, optional
        Radius of the circle that bounds the area to be transformed.
    AngleOversampling : int
        Oversample factor for number of angles
    output_shape : tuple (row, col), optional
    scaling : {'linear', 'log'}, optional
        Specify whether the image warp is polar or log-polar. Defaults to
        'linear'.
    multichannel : bool, optional
        Whether the image is a 3-D array in which the third axis is to be
        interpreted as multiple channels. If set to `False` (default), only 2-D
        arrays are accepted.
    **kwargs : keyword arguments
        Passed to `transform.warp`.

    Returns
    -------
    warped : ndarray
        The polar or log-polar warped image.

    Examples
    --------
    Perform a basic polar warp on a grayscale image:
    >>> from skimage import data
    >>> from skimage.transform import warp_polar
    >>> image = data.checkerboard()
    >>> warped = warp_polar(image)
    Perform a log-polar warp on a grayscale image:
    >>> warped = warp_polar(image, scaling='log')
    Perform a log-polar warp on a grayscale image while specifying center,
    radius, and output shape:
    >>> warped = warp_polar(image, (100,100), radius=100,
    ...                     output_shape=image.shape, scaling='log')
    Perform a log-polar warp on a color image:
    >>> image = data.astronaut()
    >>> warped = warp_polar(image, scaling='log', multichannel=True)
    """
    if image.ndim != 2 and not multichannel:
        raise ValueError("Input array must be 2 dimensions "
                         "when `multichannel=False`,"
                         " got {}".format(image.ndim))

    if image.ndim != 3 and multichannel:
        raise ValueError("Input array must be 3 dimensions "
                         "when `multichannel=True`,"
                         " got {}".format(image.ndim))

    if center is None:
        center = (np.array(image.shape)[:2] / 2) - 0.5

    if radius is None:
        w, h = np.array(image.shape)[:2] / 2
        radius = np.sqrt(w ** 2 + h ** 2)

    if output_shape is None:
        height = 360*AngleOversampling
        width = int(np.ceil(radius))
        output_shape = (height, width)
    else:
        output_shape = safe_as_int(output_shape)
        height = output_shape[0]
        width = output_shape[1]

    if scaling == 'linear':
        k_radius = width / radius
        map_func = _linear_polar_mapping
    elif scaling == 'log':
        k_radius = width / np.log(radius)
        map_func = _log_polar_mapping
    else:
        raise ValueError("Scaling value must be in {'linear', 'log'}")

    k_angle = height / (2 * np.pi)
    warp_args = {'k_angle': k_angle, 'k_radius': k_radius, 'center': center}

    warped = warp(image, map_func, map_args=warp_args,
                  output_shape=output_shape, **kwargs)

    return warped


def highpass(shape):
    """
    Return highpass filter to be multiplied with fourier transform.

    Parameters
    ----------
    shape : 'ndarray' of 'int'
        Input shape of 2d filter

    Returns
    -------
    filter
        high pass filter
    """
    x = np.outer(
        np.cos(np.linspace(-math.pi/2., math.pi/2., shape[0])),
        np.cos(np.linspace(-math.pi/2., math.pi/2., shape[1])))
    return (1.0 - x) * (2.0 - x)


def grouper(iterable, n, fillvalue=None):
    """
    Collect data into fixed-length chunks or blocks.

    Parameters
    ----------
    iterable : TYPE
        DESCRIPTION.
    n : TYPE
        DESCRIPTION.
    fillvalue : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


@ray.remote
class ProgressBarActor:
    counter: int
    delta: int
    event: Event

    def __init__(self) -> None:
        self.counter = 0
        self.delta = 0
        self.event = Event()

    def update(self, num_items_completed: int) -> None:
        """
        Updates the ProgressBar with the incremental
        number of items that were just completed.
        """
        self.counter += num_items_completed
        self.delta += num_items_completed
        self.event.set()

    async def wait_for_update(self) -> Tuple[int, int]:
        """
        Blocking call.

        Waits until somebody calls `update`, then returns a tuple of
        the number of updates since the last call to
        `wait_for_update`, and the total number of completed items.
        """
        await self.event.wait()
        self.event.clear()
        saved_delta = self.delta
        self.delta = 0
        return saved_delta, self.counter

    def get_counter(self) -> int:
        """
        Returns the total number of complete items.
        """
        return self.counter


class ProgressBar:
    progress_actor: ActorHandle
    total: int
    description: str
    pbar: tqdm

    def __init__(self, total: int, description: str = ""):
        # Ray actors don't seem to play nice with mypy, generating
        # a spurious warning for the following line,
        # which we need to suppress. The code is fine.
        self.progress_actor = ProgressBarActor.remote()  # type: ignore
        self.total = total
        self.description = description

    @property
    def actor(self) -> ActorHandle:
        """
        Returns a reference to the remote `ProgressBarActor`.

        When you complete tasks, call `update` on the actor.
        """
        return self.progress_actor

    def print_until_done(self) -> None:
        """Blocking call.

        Do this after starting a series of remote Ray tasks, to which you've
        passed the actor handle. Each of them calls `update` on the actor.
        When the progress meter reaches 100%, this method returns.
        """
        pbar = tqdm(desc=self.description, total=self.total)
        while True:
            delta, counter = ray.get(self.actor.wait_for_update.remote())
            pbar.update(delta)
            if counter >= self.total:
                pbar.close()
                return


def ray_loop(dataCube, ROICube=None, upsampleFactor=111,
             AngleOversampling=2, nreference=6, maxNumberOfCPUs=2,
             useMultiProcesses=True):
    """
    Ray wrapper around determine_relative_source_position function.

    Performs parallel loop for different reference integrations to determine
    the relative source movement on the detector.

    Parameters
    ----------
    dataCube : 'ndarray'
        Input spectral image data cube. Fist dimention is dispersion direction,
        second dimintion is cross dispersion direction and the last dimension
        is time. The shortest wavelengths are at the first row
    ROICube : 'ndarray' of 'bool', optional
        Region of Interest
    nreferences : 'int', optional
        Number of reference times used to determine the relative movement
    upsampleFactor : 'int, optional
        Upsample factor for translational movement
    AngleOversampling : 'int, optional
        Upsample factor for determination of rotational movement.
    max_number_of_cpus : 'int', optional
        Maximum number of CPU used when using parallel calculations.
    useMultiProcesses : 'bool', optional
        If True, calculations will be done in parallel.

    Returns
    -------
    relativeSourcePosition : 'collections.OrderedDict'
        Ordered dict containing the relative rotation angle,
        scaling and x and y position as a function of time.
    """
    ntime = dataCube.shape[-1]
    if not useMultiProcesses:
        # create new function with all fixed inout variables fixed.
        func = partial(determine_relative_source_position, dataCube,
                       ROICube, upsampleFactor=upsampleFactor,
                       AngleOversampling=AngleOversampling)
        ITR = list(np.linspace(0, ntime-1, nreference, dtype=int))
        movement_iterator = map(func, ITR)

        for j in tqdm(movement_iterator, total=len(ITR),
                      dynamic_ncols=True):
            yield j
    else:
        ncpu = int(np.min([maxNumberOfCPUs, np.max([1, mp.cpu_count()-3])]))
        ray.init(num_cpus=ncpu)

        dataCube_id = ray.put(dataCube)
        ROICube_id = ray.put(ROICube)
        upsampleFactor_id = ray.put(upsampleFactor)
        AngleOversampling_id = ray.put(AngleOversampling)
        ITR = iter(np.linspace(0, ntime-1, nreference, dtype=int))

        pb = ProgressBar(nreference,
                         'Determine Telescope movement for '
                         '{} reference times'.format(nreference))
        actor = pb.actor
        result_ids = \
            [ray_determine_relative_source_position.remote(
                dataCube_id,
                ROICube_id,
                x,
                actor,
                upsampleFactor=upsampleFactor_id,
                AngleOversampling=AngleOversampling_id) for x in ITR]
        pb.print_until_done()
        MPITR = ray.get(result_ids)
        for relativeSourcePosition in MPITR:
            yield relativeSourcePosition

        ray.shutdown()


def register_telescope_movement(cleanedDataset, ROICube=None,  nreferences=6,
                                mainReference=4, upsampleFactor=111,
                                AngleOversampling=2, verbose=False,
                                verboseSaveFile=None, maxNumberOfCPUs=2,
                                useMultiProcesses=True):
    """
    Register the telescope movement.

    Parameters
    ----------
    cleanedDataset : 'SpectralDataTimeSeries'
        Input dataset. Note that for image registration to work properly,
        bad pixels need ti be removed (cleaned) first. This routine checks if
        a cleaned dataset is used by checking for the isCleanedData flag.
    ROICube : 'ndarray' of 'bool', optional
        Cube containing the Region of interest for each integration.
        If not given, it is assumed that the mask of the cleanedDataset
        contains the region of interest.
    nreferences : 'int', optional
        Default is 6.
    mainReference : 'int', optional
        Default is 4.
    upsampleFactor : 'int, optional
        Upsample factor of FFT images to determine relative movement at
        sub-pixel level. Default is 111
    AngleOversampling : 'int, optional
        Upsampling factor of the angle in the to polar coordinates transformed
        FFT images to determine the relative rotation adn scale change.
        Default is 2.
    verbose : 'bool', optional
        If true diagnostic plots will be generated. Default is False
    verboseSaveFile : 'str', optional
        If not None, verbose output will be saved to the specified file.
    max_number_of_cpus : 'int', optional
        Maxiumum bumber of cpu's to be used.

    Returns
    -------
    spectralMovement : 'OrderedDict'
        Ordered dict containing the relative rotation, scaling,
        and movement in the dispersion and cross dispersion direction.

    Raises
    ------
    ValueError, TypeError
        Errors are raised if certain data is not present of from the wrong
        type.
    """
    try:
        if cleanedDataset.isCleanedData is False:
            raise ValueError
    except (ValueError, AttributeError):
        raise TypeError("Input dataset is not recognized as cleaned dataset")

    maskeddata = cleanedDataset.return_masked_array('data').copy()
    dataUse = maskeddata.data
    if ROICube is None:
        ROICube = maskeddata.mask
    else:
        ROICube = np.logical_or(maskeddata.mask, ROICube)

    ntime = dataUse.shape[-1]
    if (nreferences < 1) | (nreferences > ntime):
        raise ValueError("Wrong nreferences value")
    if (mainReference < 0) | (mainReference > nreferences):
        raise ValueError("Wrong mainReference value")

    determinePositionIterator = \
        ray_loop(dataUse, ROICube=ROICube,
                 upsampleFactor=upsampleFactor,
                 AngleOversampling=AngleOversampling,
                 nreference=nreferences,
                 maxNumberOfCPUs=maxNumberOfCPUs,
                 useMultiProcesses=useMultiProcesses)
    iteratorResults = list(determinePositionIterator)

    referenceIndex = np.linspace(0, ntime-1, nreferences, dtype=int)
    testAngle = np.zeros((nreferences, ntime))
    testScale = np.zeros((nreferences, ntime))
    testCrossDispShift = np.zeros((nreferences, ntime))
    testDispShift = np.zeros((nreferences, ntime))
    for i in range(nreferences):
        testAngle[i, :] = iteratorResults[i]['relativeAngle'] - \
            iteratorResults[i]['relativeAngle'][referenceIndex[mainReference]]
        testScale[i, :] = iteratorResults[i]['relativeScale'] / \
            iteratorResults[i]['relativeScale'][referenceIndex[mainReference]]
        testCrossDispShift[i, :] = iteratorResults[i]['cross_disp_shift'] - \
            iteratorResults[i]['cross_disp_shift'][referenceIndex[
                mainReference]]
        testDispShift[i, :] = iteratorResults[i]['disp_shift'] - \
            iteratorResults[i]['disp_shift'][referenceIndex[mainReference]]

    relativeAngle = np.median(testAngle, axis=0)
    relativeScale = np.median(testScale, axis=0)
    crossDispersionShift = np.median(testCrossDispShift, axis=0)
    dispersionShift = np.median(testDispShift, axis=0)
    # shift to first time index
    testAngle = testAngle - relativeAngle[0]
    testScale = testScale / relativeScale[0]
    testCrossDispShift = testCrossDispShift - crossDispersionShift[0]
    testDispShift = testDispShift - dispersionShift[0]
    relativeAngle = relativeAngle - relativeAngle[0]
    relativeScale = relativeScale / relativeScale[0]
    crossDispersionShift = crossDispersionShift - crossDispersionShift[0]
    dispersionShift = dispersionShift - dispersionShift[0]

    if verbose:
        sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
        sns.set_style("white", {"xtick.bottom": True, "ytick.left": True})
        fig, axes = plt.subplots(figsize=(14, 12), nrows=2, ncols=2)
        ax0, ax1, ax2, ax3 = axes.flatten()
        ax0.plot(testAngle.T)
        ax0.plot(relativeAngle, lw=5)
        ax0.set_title('Relative Angle')
        ax0.set_xlabel('Integration #')
        ax0.set_ylabel('Angle [degrees]')
        ax1.plot(testScale.T)
        ax1.plot(relativeScale, lw=5)
        ax1.set_title('Relative Scale')
        ax1.set_xlabel('Integration #')
        ax1.set_ylabel('Scaling Factor')
        ax2.plot(testCrossDispShift.T)
        ax2.plot(crossDispersionShift, lw=5)
        ax2.set_title('Relative Cross-dispersion shift')
        ax2.set_ylabel('Shift [pixles]')
        ax2.set_xlabel('Integration #')
        ax3.plot(testDispShift.T)
        ax3.plot(dispersionShift, lw=5)
        ax3.set_title('Relative Dispersion shift')
        ax3.set_xlabel('Integration #')
        ax3.set_ylabel('Shift [pixles]')
        fig.subplots_adjust(hspace=0.3)
        fig.subplots_adjust(wspace=0.45)
        plt.show()
        if verboseSaveFile is not None:
            fig.savefig(verboseSaveFile, bbox_inches='tight')

    spectralMovement = \
        collections.OrderedDict(relativeAngle=relativeAngle,
                                relativeScale=relativeScale,
                                crossDispersionShift=crossDispersionShift,
                                dispersionShift=dispersionShift,
                                referenceIndex=referenceIndex[mainReference])
    return spectralMovement


def determine_center_of_light_posision(cleanData, ROI=None, verbose=False,
                                       quantileCut=0.5, orderTrace=2):
    """
    Determine the center of light position.

    This routine determines the center of light position (cross-dispersion)
    of the dispersed light. The center of light  is defined in a similar
    way as the center of mass.  This routine also fits a polynomial to the
    spectral trace.

    Parameters
    ----------
    cleanData : 'maskedArray' or 'ndarray'
        Input data
    ROI : 'ndarray' of 'bool', optional
        Region of interest
    verbose : 'bool'
        Default is False
    quantileCut : 'float', optional
        Default is 0.5
    orderTrace : 'int'
        Default is 2

    Returns
    -------
    total_light : 'ndarray'
        Total summed signal on the detector as function of wavelength.
    idx : 'int'
    COL_pos : 'ndarray'
        Center of light position of the dispersed light.
    ytrace : 'ndarray'
        Spectral trace position in fraction of pixels in the dispersion
        direction
    xtrace : 'ndarray'
        Spectral trace position in fraction of pixels in the cross dispersion
        direction
    """
    if isinstance(cleanData, np.ma.core.MaskedArray):
        data_use = cleanData.data
        if ROI is not None:
            mask_use = cleanData.mask | ROI
        else:
            mask_use = cleanData.mask
    else:
        data_use = cleanData
        if ROI is not None:
            mask_use = ROI
        else:
            mask_use = np.zeros_like(cleanData, dtype='bool')

    data = np.ma.array(data_use, mask=mask_use)
    npix, mpix = data.shape

    position_grid = np.mgrid[0:npix, 0:mpix]
    total_light = np.ma.sum(data, axis=1)

    COL = np.ma.sum(data*position_grid[1, ...], axis=1) / \
        total_light

    treshhold = np.quantile(total_light[~total_light.mask], quantileCut)
    idx_use = np.ma.where(total_light > treshhold)[0]
    ytrace = np.arange(npix)
    idx = ytrace[idx_use]
    X = []
    for i in range(orderTrace+1):
        X.append(idx**i)
    X = np.array(X).T
    robust_fit = sm.RLM(COL[idx_use], X).fit()
    z = robust_fit.params[::-1]
    f = np.poly1d(z)
    xtrace = f(ytrace)

    if verbose:
        sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
        sns.set_style("white", {"xtick.bottom": True, "ytick.left": True})
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.plot(idx_use, total_light[idx_use])
        ax.set_title('Integrated Signal')
        ax.set_xlabel('Pixel Position Dispersion Direction')
        ax.set_ylabel('Integrated Signal')
        plt.show()

    return total_light[idx_use], idx, COL[idx_use], ytrace, xtrace


def correct_initial_wavelength_shift(referenceDataset, cascade_configuration,
                                     *otherDatasets):
    """
    Determine if there is an initial wavelength shift and correct.

    Parameters
    ----------
    referenceDataset : 'cascade.data_model.SpectralDataTimeSeries'
        Dataset who's wavelength is used as refernce of the wavelength
        correction.
    cascade_configuration : 'cascade.initialize.initialize.configurator'
        Singleton containing the confifuration parameters of cascade.
    **otherDatasets : 'cascade.data_model.SpectralDataTimeSeries'
        Optional.
        Other datasets assumed to have the same walengths as the reference
        dataset and which will be corrected simultaneously with the reference.

    Returns
    -------
    referenceDataset : 'list' of 'cascade.data_model.SpectralDataTimeSeries'
        Dataset with corrected wavelengths.
    otherDatasets_list : 'list' of 'cascade.data_model.SpectralDataTimeSeries'
        Optinal output.
    modeled_observations : 'list' of 'ndarray'
    stellar_model : 'list' of 'ndarray'
    corrected_observations : 'list' of 'ndarray'
    """
    model_spectra = SpectralModel(cascade_configuration)
    wave_shift, error_wave_shift = \
        model_spectra.determine_wavelength_shift(referenceDataset)
    referenceDataset.wavelength = referenceDataset.wavelength+wave_shift
    referenceDataset.add_auxilary(wave_shift=wave_shift.to_string())
    referenceDataset.add_auxilary(error_wave_shift=error_wave_shift.to_string())
    otherDatasets_list = list(otherDatasets)
    for i, dataset in enumerate(otherDatasets_list):
        dataset.wavelength = dataset.wavelength+wave_shift
        dataset.add_auxilary(wave_shift=wave_shift.to_string())
        dataset.add_auxilary(error_wave_shift=error_wave_shift.to_string())
        otherDatasets_list[i] = dataset
    modeled_observations = \
        [model_spectra.model_wavelength, model_spectra.model_observation,
         model_spectra.scaling, model_spectra.relative_distanc_sqr,
         model_spectra.sensitivity]
    stellar_model = \
        [model_spectra.model_wavelength, model_spectra.rebinned_stellar_model]
    input_stellar_model = [model_spectra.sm[2], model_spectra.sm[3]]
    corrected_observations = \
        [model_spectra.corrected_wavelength, model_spectra.observation,
         wave_shift, error_wave_shift]
    stellar_model_parameters = model_spectra.par
    if len(otherDatasets_list) > 0:
        return [referenceDataset] + otherDatasets_list, modeled_observations,\
            stellar_model, corrected_observations, input_stellar_model, \
            stellar_model_parameters
    return referenceDataset,  modeled_observations, stellar_model, \
        corrected_observations, input_stellar_model, \
        stellar_model_parameters


def renormalize_spatial_scans(referenceDataset, *otherDatasets):
    """
    bla.

    Parameters
    ----------
    referenceDataset : TYPE
        DESCRIPTION.
    *otherDatasets : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    otherDatasets_list = list(otherDatasets)

    try:
        scan_direction = np.array(referenceDataset.scan_direction)
    except AttributeError:
        if len(otherDatasets_list) > 0:
            return [referenceDataset] + otherDatasets_list
        return referenceDataset

    unique_scan_directions = np.unique(scan_direction)
    if len(unique_scan_directions) != 2:
        if len(otherDatasets_list) > 0:
            return [referenceDataset] + otherDatasets_list
        return referenceDataset

    idx = scan_direction == 0.0
    med0 = np.median(referenceDataset.data[...,idx]).value
    med1 = np.median(referenceDataset.data[...,~idx]).value
    med = np.median(referenceDataset.data).value
    scaling0 = med / med0
    scaling1 = med / med1

    reference_data = copy.deepcopy(referenceDataset.data)
    reference_data[..., idx] = reference_data[..., idx]*scaling0
    reference_data[..., ~idx] = reference_data[..., ~idx]*scaling1
    referenceDataset.data = reference_data
    reference_uncertainty = copy.deepcopy(referenceDataset.uncertainty)
    reference_uncertainty[..., idx] = reference_uncertainty[..., idx]*scaling0
    reference_uncertainty[...,~idx] = reference_uncertainty[...,~idx]*scaling1
    referenceDataset.uncertainty = reference_uncertainty

    for i, dataset in enumerate(otherDatasets_list):
        reference_data = copy.deepcopy(dataset.data)
        reference_data[..., idx] = reference_data[..., idx]*scaling0
        reference_data[..., ~idx] = reference_data[..., ~idx]*scaling1
        dataset.data = reference_data
        reference_uncertainty = copy.deepcopy(dataset.uncertainty)
        reference_uncertainty[..., idx] = reference_uncertainty[..., idx]*scaling0
        reference_uncertainty[...,~idx] = reference_uncertainty[...,~idx]*scaling1
        dataset.uncertainty = reference_uncertainty
        otherDatasets_list[i] = dataset
    if len(otherDatasets_list) > 0:
        return [referenceDataset] + otherDatasets_list
    return referenceDataset


def determine_absolute_cross_dispersion_position(cleanedDataset, initialTrace,
                                                 ROI=None,
                                                 verbose=False,
                                                 verboseSaveFile=None,
                                                 quantileCut=0.5,
                                                 orderTrace=2):
    """
    Determine the initial cross dispersion position.

    This routine updates the initial spectral trace for positional shifts in
    the cross dispersion direction for the first exposure of the the time
    series observation.

    Parameters
    ----------
    cleanedDataset : 'SpectralDataTimeSeries'
    initialTrace : 'OrderedDict'
        input spectral trace.
    ROI : 'ndarray' of 'bool'
        Region of interest.
    verbose : 'bool', optional
        If true diagnostic plots will be generated. Default is False
    verboseSaveFile : 'str', optional
        If not None, verbose output will be saved to the specified file.
    quantileCut : 'float', optional
        Default is 0.5
    orderTrace : 'int', optional
        Default is 2

    Returns
    -------
    newShiftedTrace : 'OrderedDict'
        To the observed source poisiton shifted spectral trace
    newFittedTrace : 'OrderedDict'
        Trace determined by fit to the center of light position.
    initialCrossDispersionShift : 'float'
        Shift between initial guess for spectral trace position and
        fitted trace position of the first spectral image.
    """
    cleanedData = cleanedDataset.return_masked_array('data')

    newShiftedTrace = copy.copy(initialTrace)
    newFittedTrace = copy.copy(initialTrace)

    if ROI is not None:
        roiUse = cleanedData[..., 0].mask | ROI
    else:
        roiUse = cleanedData[..., 0].mask

    kernel = Gaussian2DKernel(1.0)
    convolvedFirstImage = convolve(cleanedData[..., 0], kernel,
                                   boundary='extend')

    _, idx, col, yTrace, xTrace = \
        determine_center_of_light_posision(convolvedFirstImage, ROI=roiUse,
                                           quantileCut=quantileCut,
                                           orderTrace=orderTrace)

    medianCrossDispersionPosition = np.ma.median(xTrace[idx])
    medianCrossDispersionPositionInitialTrace = \
        np.ma.median(initialTrace['positional_pixel'].value[idx])

    initialCrossDispersionShift = \
        medianCrossDispersionPosition-medianCrossDispersionPositionInitialTrace

    newShiftedTrace['positional_pixel'] = \
        newShiftedTrace['positional_pixel'] + \
        initialCrossDispersionShift*newShiftedTrace['positional_pixel'].unit
    newFittedTrace['positional_pixel'] = \
        xTrace*newFittedTrace['positional_pixel'].unit

    if verbose:
        sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
        sns.set_style("white", {"xtick.bottom": True, "ytick.left": True})
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.plot(idx, col, label='COL')
        ax.plot(yTrace, xTrace, label='Fitted Trace')
        ax.plot(newShiftedTrace['positional_pixel'].value,
                label='Shifted Instrument Trace')
        ax.plot(initialTrace['positional_pixel'].value,
                label='Initial Instrument Trace')
        ax.legend(loc='best')
        ax.set_title('Initial Trace Position')
        ax.set_xlabel('Pixel Position Dispersion Direction')
        ax.set_ylabel('Pixel Position Cross-Dispersion Direction')
        plt.show()
        if verboseSaveFile is not None:
            fig.savefig(verboseSaveFile, bbox_inches='tight')

    return newShiftedTrace, newFittedTrace, medianCrossDispersionPosition, \
        initialCrossDispersionShift


def correct_wavelength_for_source_movent(datasetIn, spectral_movement,
                                         useScale=False,
                                         useCrossDispersion=False,
                                         verbose=False, verboseSaveFile=None):
    """
    Correct wavelengths for source movement.

    This routine corrects the wavelength cube attached to the spectral
    image data cube for source (telescope) movements

    Parameters
    ----------
    datasetIn : 'SpectralDataTimeSeries'
        Input dataset for which the waveength will be corrected for telescope
        movement
    spectral_movement : 'OrderedDict'
        Ordered dict containing the relative rotation, scaling,
        and movement in the dispersion and cross dispersion direction.
     useScale : 'bool', optional
         If set the scale parameter is used to correct the wavelength.
         Default is False.
    useCrossDispersion : 'bool', optional
        If set the coress dispersion movement is used to correct the
        wavelength. Default is False.
    verbose : 'bool', optional
        If true diagnostic plots will be generated. Default is False
    verboseSaveFile : 'str', optional
        If not None, verbose output will be saved to the specified file.

    Returns
    -------
    dataset_out : 'SpectralDataTimeSeries'
        The flag isMovementCorrected=True is set to indicate that this dataset
        is corrected

    Notes
    -----
    Scaling changes are not corrected at the moment. Note that the used
    rotation and translation to correct the wavelengths is the relative
    source movement defined such that shifting the observed spectral image by
    these angles and shifts the position would be identical to the reference
    image. The correction of the wavelength using the reference spectral image
    is hence in the oposite direction.
    """
    dataset_out = copy.deepcopy(datasetIn)

    correctedWavelength = dataset_out.return_masked_array('wavelength').copy()
    # no need for mask here as wavekength should be difined for all pixels
    correctedWavelength = correctedWavelength.data

    ntime = correctedWavelength.shape[-1]
    for it in range(ntime):
        rows, cols = (np.array(correctedWavelength.shape)[:2] / 2) - 0.5
        center = np.array((cols, rows)) / 2. - 0.5
        tform1 = SimilarityTransform(translation=center)
        angle_rad = np.deg2rad(-spectral_movement['relativeAngle'][it])
        scale = spectral_movement['relativeScale'][it]
        tform2 = SimilarityTransform(rotation=angle_rad,
                                     scale=(1.0/scale-1.0)*int(useScale)+1.0)
        tform3 = SimilarityTransform(translation=-center)
        tform_rotate = tform3 + tform2 + tform1
        translation = (-spectral_movement['crossDispersionShift'][it] *
                       int(useCrossDispersion),
                       -spectral_movement['dispersionShift'][it])
        tform_translate = SimilarityTransform(translation=translation)
        tform_combined = tform_translate + tform_rotate
        correctedWavelength[..., it] = warp(correctedWavelength[..., it],
                                            tform_combined, order=3,
                                            cval=np.nan)

    # mask those regions of the images which are on the edge and migth
    # not be present at all times.
    correctedWavelength = np.ma.masked_invalid(correctedWavelength)
    ncols = correctedWavelength.shape[1]
    for ic in range(ncols):
        correctedWavelength[:, ic, :] = \
            np.ma.mask_rows(correctedWavelength[:, ic, :])
    # replace old wavelengths and update mask.
    dataset_out._wavelength = correctedWavelength.data
    dataset_out.mask = np.logical_or(dataset_out.mask,
                                     correctedWavelength.mask)
    dataset_out.isMovementCorrected = True

    if verbose:
        wnew = dataset_out.return_masked_array('wavelength')
        index_valid = np.ma.all(wnew[..., 0].mask, axis=1)
        index_valid = ~index_valid.data
        wnew = np.ma.median(wnew[index_valid, ...][1:8, ...], axis=1)
        wold = np.ma.median(datasetIn.wavelength[index_valid, ...][1:8, ...],
                            axis=1)
        sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
        sns.set_style("white", {"xtick.bottom": True, "ytick.left": True})
        fig, ax0 = plt.subplots(figsize=(6, 5), nrows=1, ncols=1)
        ax0.plot(wnew.T, zorder=5, lw=3)
        ax0.plot(wold.T, color='gray', zorder=4)
        ax0.set_ylabel('Wavelength [{}]'.format(datasetIn.wavelength_unit))
        ax0.set_xlabel('Integration #')
        ax0.set_title('Wavelength shifts')
        plt.show()
        if verboseSaveFile is not None:
            fig.savefig(verboseSaveFile, bbox_inches='tight')
    return dataset_out


def rebin_to_common_wavelength_grid(dataset, referenceIndex, nrebin=None,
                                    verbose=False, verboseSaveFile=None):
    """
    Rebin the spectra to single wavelength per row.

    Parameters
    ----------
    dataset : 'SpectralDataTimeSeries'
        Input dataset
    referenceIndex : 'int'
        Exposure index number which will be used as reference defining the
        uniform wavelength grid.
    nrebin : 'float', optional
        rebinning factor for the new wavelength grid compare to the old.
    verbose : 'bool', optional
        If True, diagnostic plots will be created
    verboseSaveFile : 'str', optional
        If not None, verbose output will be saved to the specified file.

    Returns
    -------
    rebinnedDataset : 'SpectralDataTimeSeries'
        Output to common wavelength grid rebinned dataset
    """
    if not isinstance(dataset, SpectralDataTimeSeries):
        raise TypeError("the input data to rebin_to_common_wavelength_grid "
                        "function needs to be a SpectralDataTimeSeries. "
                        "Aborting rebin to a common wavelength grid.")
    # all data with wavelength dependency + time
    spectra = dataset.return_masked_array('data')
    uncertainty = dataset.return_masked_array('uncertainty')
    wavelength = dataset.return_masked_array('wavelength')
    time = dataset.return_masked_array('time')

    # A pixel row (time) does not have the same wavelength in time
    # Need to find the miximum-lowest or minimum-higest wavelength for a
    # proper rebinning.
    min_wavelength = np.ma.min(np.ma.max(wavelength, axis=-1))
    max_wavelength = np.ma.max(np.ma.min(wavelength, axis=-1))

    referenceWavelength = np.sort(np.array(wavelength[1:-1, referenceIndex]))
    idx_min_select = np.where(referenceWavelength >= min_wavelength)[0][0]
    idx_max_select = np.where(referenceWavelength <= max_wavelength)[0][-1]
    referenceWavelength = referenceWavelength[idx_min_select:idx_max_select]

    lr, ur = _define_band_limits(wavelength)

    if nrebin is not None:
        referenceWavelength = \
            np.linspace(referenceWavelength[0+int(nrebin/2)],
                        referenceWavelength[-1-int(nrebin/2)],
                        int(len(referenceWavelength)/nrebin))
    lr0, ur0 = _define_band_limits(referenceWavelength)
    weights = _define_rebin_weights(lr0, ur0, lr, ur)
    rebinnedSpectra, rebinnedUncertainty = \
        _rebin_spectra(spectra, uncertainty, weights)
    rebinnedWavelength = np.tile(referenceWavelength,
                                 (rebinnedSpectra.shape[-1], 1)).T

    ndim = dataset.data.ndim
    selection = tuple((ndim-1)*[0]+[Ellipsis])

    dictTimeSeries = {}

    dictTimeSeries['data'] = rebinnedSpectra
    dictTimeSeries['data_unit'] = dataset.data_unit
    dictTimeSeries['uncertainty'] = rebinnedUncertainty
    dictTimeSeries['wavelength'] = rebinnedWavelength
    dictTimeSeries['wavelength_unit'] = dataset.wavelength_unit
    dictTimeSeries['time'] = time[selection]
    dictTimeSeries['time_unit'] = dataset.time_unit
    dictTimeSeries['isRebinned'] = True

    # get everything else apart from data, wavelength, time and uncertainty
    for key in vars(dataset).keys():
        if key[0] != "_":
            if isinstance(vars(dataset)[key], MeasurementDesc):
                measurement = getattr(dataset, key)
                dictTimeSeries[key] = measurement[selection]
            else:
                # print('can be added withour rebin')
                dictTimeSeries[key] = getattr(dataset, key)

    rebinnedDataset = SpectralDataTimeSeries(**dictTimeSeries)

    if verbose:
        index_valid = np.ma.all(wavelength.mask, axis=1)
        index_valid = ~index_valid.data
        sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
        sns.set_style("white", {"xtick.bottom": True, "ytick.left": True})
        fig, axes = plt.subplots(figsize=(10, 5), nrows=1, ncols=2)
        ax0, ax1 = axes.flatten()
        ax0.plot(rebinnedWavelength[1:6].T, zorder=5, lw=3)
        ax0.plot(wavelength[index_valid, :][2:7].T, color='gray', zorder=4)
        ax0.set_ylabel('Wavelength [{}]'.format(dataset.wavelength_unit))
        ax0.set_xlabel('Integration #')
        ax0.set_title('Rebinned Wavelengths')
        ax1.plot(rebinnedSpectra[1:6].T, zorder=5, lw=3)
        ax1.plot(spectra[index_valid, :][2:7].T, color='gray', zorder=4)
        ax1.set_ylabel('Flux [{}]'.format(dataset.data_unit))
        ax1.set_xlabel('Integration #')
        ax1.set_title('Rebinned Signal')
        fig.subplots_adjust(hspace=0.3)
        fig.subplots_adjust(wspace=0.45)
        plt.show()
        if verboseSaveFile is not None:
            fig.savefig(verboseSaveFile, bbox_inches='tight')
    return rebinnedDataset


def combine_scan_samples(datasetIn, scanDictionary, verbose=False):
    """
    Combine all (scan) samples.

    This routine creates a new SpectralDataTimeSeries of integration averaged
    spectra from time series data per sample-up-the-ramp.

    Parameters
    ----------
    datasetIn : 'SpectralDataTimeSeries'
        Input dataset
    scanDictionary : 'dict'
        Dictionary containg relevant indormation about the scans
    verbose : 'bool', optional
        If True, diagnostic plots will be created (default False).

    Returns
    -------
    combinedDataset : 'SpectralDataTimeSeries'
        Output dataset with average signals per integration

    Raises
    ------
    AttributeError
    """
    if not isinstance(datasetIn, SpectralDataTimeSeries):
        raise TypeError("the input data to combine_scan_sample function needs "
                        "to be a SpectralDataTimeSeries. Aborting combining "
                        "scan samples.")

    dataIn = datasetIn.return_masked_array('data').copy()
    dataInShape = dataIn.shape
    errorIn = datasetIn.return_masked_array('uncertainty').copy()
    dataUnit = datasetIn.data_unit
    waveIn = datasetIn.return_masked_array('wavelength').copy()
    wavelengthUnit = datasetIn.wavelength_unit
    timeIn = datasetIn.return_masked_array('time').copy()
    timeUnit = datasetIn.time_unit

    dictTimeSeries = {}

    def reshape_integration(data, shape, nreads):
        reshapedData = \
            np.ma.mean(np.ma.reshape(data, (shape[0], shape[1]//nreads,
                                            nreads)),
                       axis=-1)
        return reshapedData

    def reshape_error(error, shape, nreads):
        reshapedError = \
            np.ma.sqrt(np.ma.sum(np.ma.reshape(error,
                                               (shape[0], shape[1]//nreads,
                                                nreads))**2,
                       axis=-1))/nreads
        return reshapedError

    def reshape_auxilary(data, shape, nreads):
        reshapedData = \
            np.mean(np.reshape(data, (len(data)//nreads, nreads)), axis=-1)
        return list(reshapedData)
    
    def combine_list_of_strings(data, scanDictionary, sort_index):
        reshapedData = []
        for scan_dir, scan_par in scanDictionary.items():
            data_scan = np.ma.compress(scan_par['index'], data, axis=-1)
            reshapedData.append(data_scan[::scan_par['nsamples']])
        reshapedData = np.hstack(reshapedData) 
        reshapedData = np.take_along_axis(reshapedData,
                                          sort_index.mean(axis=0, dtype=int),
                                          axis=-1)
        reshapedData = list(reshapedData)
        base = [j for i in reshapedData
                for j in i.split("_") if 'sample' in j]
        if len(base) != 0:
            reshapedData = \
                [i.replace(base[j], "RESAMPLED{:04d}".format(j))
                 for j, i in enumerate(reshapedData)]
        return reshapedData

    def combine_scans_auxilary(data, scanDictionary, sort_index):
        reshapedData = []
        for scan_dir, scan_par in scanDictionary.items():
            data_scan = np.ma.compress(scan_par['index'], data, axis=-1)
            reshapedData.append(
                reshape_auxilary(data_scan, data_scan.shape,
                                 scan_par['nsamples']))
        reshapedData = np.hstack(reshapedData) 
        reshapedData = np.take_along_axis(reshapedData,
                                          sort_index.mean(axis=0, dtype=int),
                                          axis=-1)
        if isinstance(data, list):
            reshapedData = list(reshapedData)
        return reshapedData
    
    def reshape_data(data, scanDictionary, sort_index):
        reshapedData = []
        for scan_dir, scan_par in scanDictionary.items():
            data_scan = np.ma.compress(scan_par['index'], data, axis=-1)
            reshapedData.append(
                reshape_integration(data_scan, data_scan.shape,
                                    scan_par['nsamples']))
        reshapedData = np.hstack(reshapedData)    
        reshapedData = np.take_along_axis(reshapedData, sort_index, axis=-1)
        return reshapedData
    
    def reshape_primary_data(data, wave, error,  time,  scanDictionary):
        reshapedData = []
        reshapedTime = []
        reshapedWave = []
        reshapedError = []
        for scan_dir, scan_par in scanDictionary.items():
            data_scan = np.ma.compress(scan_par['index'], data, axis=-1)
            wave_scan = np.ma.compress(scan_par['index'], wave, axis=-1)
            time_scan = np.ma.compress(scan_par['index'], time, axis=-1)
            error_scan = np.ma.compress(scan_par['index'], error, axis=-1)
            reshapedData.append(
                reshape_integration(data_scan, data_scan.shape,
                                    scan_par['nsamples']))
            reshapedTime.append(
               reshape_integration(time_scan, data_scan.shape,
                                   scan_par['nsamples']))  
            reshapedWave.append(
               reshape_integration(wave_scan, data_scan.shape,
                                   scan_par['nsamples']))    
            reshapedError.append(
                reshape_error(error_scan, data_scan.shape,
                                   scan_par['nsamples']))
        reshapedData = np.hstack(reshapedData)
        reshapedTime = np.hstack(reshapedTime)
        reshapedWave = np.hstack(reshapedWave)
        reshapedError = np.hstack(reshapedError)
        idx_time_sort = np.argsort(reshapedTime, axis=-1)
        reshapedData = np.take_along_axis(reshapedData, idx_time_sort, axis=-1)
        reshapedTime = np.take_along_axis(reshapedTime, idx_time_sort, axis=-1)
        reshapedWave = np.take_along_axis(reshapedWave, idx_time_sort, axis=-1)
        reshapedError = np.take_along_axis(reshapedError, idx_time_sort, axis=-1)
        print(reshapedData.shape)
        return (reshapedData, reshapedWave, reshapedError, reshapedTime,
                idx_time_sort)

    (combinedData, combinedWavelength, combinedError, combinedTime,
     idx_time_sort) = reshape_primary_data(dataIn, waveIn, errorIn, timeIn,
                                           scanDictionary)

    dictTimeSeries['data'] = combinedData
    dictTimeSeries['data_unit'] = dataUnit
    dictTimeSeries['uncertainty'] = combinedError
    dictTimeSeries['wavelength'] = combinedWavelength
    dictTimeSeries['wavelength_unit'] = wavelengthUnit
    dictTimeSeries['time'] = combinedTime
    dictTimeSeries['time_unit'] = timeUnit

    # get everything appart from data, wavelength, time and uncertainty
    for key in vars(datasetIn).keys():
        if key[0] != "_":
            if isinstance(vars(datasetIn)[key], MeasurementDesc):
                measurement = getattr(datasetIn, key)
                # will be rebinned
                dictTimeSeries[key] = \
                    reshape_data(measurement, scanDictionary, idx_time_sort)
            elif isinstance(vars(datasetIn)[key], AuxilaryInfoDesc):
                aux = getattr(datasetIn, key)
                if isinstance(aux, list):
                    # list
                    if (len(aux) == dataInShape[-1]) & \
                      (isinstance(aux[0], str)):
                        # list of str needs to be rebinned
                        dictTimeSeries[key] = \
                            combine_list_of_strings(aux, scanDictionary,
                                                    idx_time_sort)
                    elif len(aux) == dataInShape[-1]:
                        dictTimeSeries[key] = \
                            combine_scans_auxilary(aux, scanDictionary,
                                                   idx_time_sort)
                    else:
                        dictTimeSeries[key] = aux
                elif isinstance(aux, np.ndarray):
                    if len(aux) == dataInShape[-1]:
                        # no list, no number
                        dictTimeSeries[key] = \
                            combine_scans_auxilary(aux, scanDictionary,
                                                   idx_time_sort)
                    else:
                        # no list but not an array
                        dictTimeSeries[key] = aux
                else:
                    # other aux, no rebin
                    dictTimeSeries[key] = aux
            else:
                # print('can be added withour rebin')
                dictTimeSeries[key] = getattr(datasetIn, key)

    combinedDataset = \
        SpectralDataTimeSeries(**dictTimeSeries)

    if verbose:
        sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
        sns.set_style("white", {"xtick.bottom": True, "ytick.left": True})
        fig, ax0 = plt.subplots(figsize=(6, 5), nrows=1, ncols=1)
        ax0.plot(np.ma.mean(dataIn, axis=-1), label='Signal per sample',
                 zorder=5, lw=3)
        ax0.plot(np.ma.mean(combinedData, axis=-1),
                 label='Signal per integration', zorder=6, lw=3)
        ax0.legend(loc='best')
        ax0.set_ylabel('Flux [{}]'.format(dataUnit))
        ax0.set_xlabel('Wavelength [{}]'.format(wavelengthUnit))
        ax0.set_title('Ramp averaged mean signal spectrum')
        plt.show()

    return combinedDataset


def sigma_clip_data_cosmic(data, sigma):
    """
    Sigma clip of time series data cube allong the time axis.

    Parameters
    ----------
    data : `ndarray`
        Input data to be cliped, last axis of data to be assumed the time
    sigma : `float`
        Sigma value of sigmaclip

    Returns
    -------
    sigmaClipMask : `ndarray` of 'bool'
        Updated mask for input data with bad data points flagged `(=1)`
    """
    # time axis always the last axis in data,
    # or the first in the transposed array
    filtereData = sigma_clip(data.T, sigma=sigma, axis=0)
    sigmaClipMask = filtereData.mask.T
    return sigmaClipMask


def sigma_clip_data(datasetIn, sigma, nfilter):
    """
    Perform sigma clip on science data to flag bad data.

    Parameters
    ----------
    datasetIn : 'SpectralDataTimeSeries'
        Input dataset
    sigma : `float`
        Sigma value of sigmaclip.
    nfilter : 'int'
        Filter length for sigma clip.

    Returns
    -------
    datasetOut : 'SpectralDataTimeSeries'
        Input data set containing the observed spectral time series. After
        completinion of this function, the dataset is returned with an
        updated mask with all diviating pixels flaged. Tho indicate this the
        flag isSigmaCliped=True is set.

    """
    if nfilter % 2 == 0:  # even
        nfilter += 1

    # mask cosmic hits
    data = datasetIn.return_masked_array("data")
    sigmaClipedMask = sigma_clip_data_cosmic(data, sigma)
    # update mask
    sigmaClipedMask = np.ma.mask_or(datasetIn.mask, sigmaClipedMask,
                                    shrink=False)
    datasetIn.mask = sigmaClipedMask

    dim = datasetIn.data.shape
    ndim = datasetIn.data.ndim
    newMask = datasetIn.mask.copy()

    for il in range(0+(nfilter-1)//2, dim[0]-(nfilter-1)//2):
        filter_index = \
            [slice(il - (nfilter-1)//2, il+(nfilter-1)//2+1, None)] + \
            [slice(None)]*(ndim-1)
        filter_index = tuple(filter_index)
        # reformat to masked array without quantity
        data = datasetIn.return_masked_array("data")
        # median along time axis
        data = np.ma.median(data[filter_index].T, axis=0)
        # filter in box in the wavelength direction
        data = sigma_clip(data, sigma=sigma, axis=ndim-2)
        # specra:  tiling=(dim[1], 1)
        # spectral images:  tiling=(dim[2], 1, 1)
        # spectral data cubes: tiling=(dim[3], 1, 1, 1)
        tiling = dim[ndim-1:] + tuple(np.ones(ndim-1).astype(int))
        mask = np.tile(data.mask, tiling)
        # add to mask
        newMask[filter_index] = np.ma.mask_or(newMask[filter_index], mask.T,
                                              shrink=False)

    newMask = np.ma.mask_or(datasetIn.mask, newMask, shrink=False)

    # update mask and set flag
    datasetIn.mask = newMask
    datasetIn.isSigmaCliped = True
    return datasetIn


def create_cleaned_dataset(datasetIn, ROIcube, kernel, stdvKernelTime):
    """
    Create a cleaned dataset to be used in regresion analysis.

    Parameters
    ----------
    datasetIn : 'SpectralDataTimeSeries'
        Input dataset
    ROIcube : 'ndarray' of 'bool'
        Region of interest.
    kernel : 'ndarray'
        Instrument convolution kernel
    stdvKernelTime : 'float'
        Standeard devistion in time direction used in convolution.

    Returns
    -------
    cleanedDataset : `SpectralDataTimeSeries`
        A cleaned version of the spectral timeseries data of the transiting
        exoplanet system

    """
    dataToBeCleaned = datasetIn.return_masked_array('data')
    uncertaintyToBeCleaned = datasetIn.return_masked_array('uncertainty')
    dataToBeCleaned.set_fill_value(np.nan)
    uncertaintyToBeCleaned.set_fill_value(np.nan)

    ndim = dataToBeCleaned.ndim

    if ndim == 2:
        RS = RobustScaler(with_scaling=True)
        data_scaled = RS.fit_transform(dataToBeCleaned.filled().T)
        dataToBeCleaned = \
            np.ma.array(data_scaled.T, mask=dataToBeCleaned.mask)
        RS2 = RobustScaler(with_scaling=True)
        data_scaled2 = RS2.fit_transform(uncertaintyToBeCleaned.filled().T)
        uncertaintyToBeCleaned = \
            np.ma.array(data_scaled2.T, mask=uncertaintyToBeCleaned.mask)

    dataToBeCleaned[ROIcube] = 0.0
    dataToBeCleaned.mask[ROIcube] = False
    dataToBeCleaned.set_fill_value(np.nan)
    uncertaintyToBeCleaned[ROIcube] = 1.0
    uncertaintyToBeCleaned.mask[ROIcube] = False
    uncertaintyToBeCleaned.set_fill_value(np.nan)

    kernel_size = kernel.shape[0]
    kernel_1d = Gaussian1DKernel(stdvKernelTime, x_size=kernel_size)
    kernel = np.repeat(np.expand_dims(kernel, axis=ndim-1),
                       (kernel_size), axis=ndim-1)
    selection = tuple([slice(None)])+tuple([None])*(ndim-1)
    kernel = kernel*kernel_1d.array[selection].T
    kernel = kernel/np.sum(kernel)

    cleanedData = \
        interpolate_replace_nans(dataToBeCleaned.filled(),
                                 kernel, boundary='extend')
    cleanedUncertainty = \
        interpolate_replace_nans(uncertaintyToBeCleaned.filled(),
                                 kernel, boundary='extend')
    if ndim == 2:
        cleanedData = \
            RS.inverse_transform(cleanedData.T).T
        cleanedUncertainty = \
            RS2.inverse_transform(cleanedUncertainty.T).T

#    cleanedData.mask = cleanedData.mask | ROI

    selection = tuple((ndim-1)*[0]+[Ellipsis])

    cleanedDataset = SpectralDataTimeSeries(
        wavelength=datasetIn.wavelength,
        wavelength_unit=datasetIn.wavelength_unit,
        data=cleanedData,
        data_unit=datasetIn.data_unit,
        mask=ROIcube,
        time=datasetIn.time.data.value[selection],
        time_unit=datasetIn.time_unit,
        time_bjd=datasetIn.time_bjd.data.value[selection],
        time_bjd_unit=datasetIn.time_bjd_unit,
        uncertainty=cleanedUncertainty,
        target_name=datasetIn.target_name,
        dataProduct=datasetIn.dataProduct,
        dataFiles=datasetIn.dataFiles,
        isCleanedData=True)

    return cleanedDataset


def compressROI(ROI, compressMask):
    """
    Remove masked wavelengths from ROI.

    Parameters
    ----------
    ROI : 'ndarray'
        Region of interest on detector.
    compressMask : 'ndarray'
        Compression mask indicating all valid data.

    Returns
    -------
    compressedROI : 'ndarray'
        Row (wavelength) compressed region of interest.

    """
    compressedROI = ROI[compressMask]
    return compressedROI


def compressSpectralTrace(spectralTrace, compressMask):
    """
    Remove masked wavelengths from spectral trace.

    Parameters
    ----------
    spectralTrace : 'dict'
        Spectral trace of the dispersed light on the detector.
    compressMask : 'ndarray'
        Compression mask indicating all valid data.

    Returns
    -------
    compressedsSpectralTrace : 'dict'
        Row (wavelength) compressed spectral trace.

    """
    compressedsSpectralTrace = spectralTrace.copy()
    for key in compressedsSpectralTrace.keys():
        compressedsSpectralTrace[key] = \
           compressedsSpectralTrace[key][compressMask]
    compressedsSpectralTrace['wavelength_pixel'] = \
        compressedsSpectralTrace['wavelength_pixel'] - \
        compressedsSpectralTrace['wavelength_pixel'][0]
    return compressedsSpectralTrace


def compressDataset(datasetIn, ROI):
    """
    Remove all flaged wavelengths from data set.

    Parameters
    ----------
    datasetIn : 'SpectralDataset'
        Spectral dataset.
    ROI : 'ndarray'
        Region of interest.

    Returns
    -------
    compressedDataset : SpectralDataset'
        Row (wavelength) compressed dataset.

    """
    dataIn = datasetIn.return_masked_array('data').copy()
    dataInShape = dataIn.shape
    errorIn = datasetIn.return_masked_array('uncertainty').copy()
    dataUnit = datasetIn.data_unit
    waveIn = datasetIn.return_masked_array('wavelength').copy()
    wavelengthUnit = datasetIn.wavelength_unit
    timeIn = datasetIn.return_masked_array('time').copy()
    timeUnit = datasetIn.time_unit

    fullMask = \
        np.ma.mask_or(dataIn.mask,
                      np.repeat(ROI[..., np.newaxis],
                                dataInShape[-1],
                                axis=-1),
                      shrink=False)
    compressMask = ~fullMask.all(axis=tuple(np.arange(1, dataIn.ndim)))

    dictTimeSeries = {}
    for key in vars(datasetIn).keys():
        if key[0] != "_":
            if isinstance(vars(datasetIn)[key], MeasurementDesc):
                dictTimeSeries[key] = \
                    getattr(datasetIn, key)[compressMask, ...]
            else:
                # print('can be added withour rebin')
                dictTimeSeries[key] = getattr(datasetIn, key)

    dictTimeSeries['data'] = dataIn[compressMask, ...]
    dictTimeSeries['data_unit'] = dataUnit
    dictTimeSeries['uncertainty'] = errorIn[compressMask, ...]
    dictTimeSeries['wavelength'] = waveIn[compressMask, ...]
    dictTimeSeries['wavelength_unit'] = wavelengthUnit
    dictTimeSeries['time'] = timeIn[compressMask, ...]
    dictTimeSeries['time_unit'] = timeUnit

    compressedDataset = \
        SpectralDataTimeSeries(**dictTimeSeries)
    return compressedDataset, compressMask
