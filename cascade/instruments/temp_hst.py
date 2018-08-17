#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 17:28:01 2018

@author: bouwman
"""
import cascade
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
import os
import fnmatch
from skimage.feature import register_translation
from skimage.feature.register_translation import _upsampled_dft
from astropy.time import Time
from astropy import coordinates as coord
from astropy import units as u
from astropy.coordinates import solar_system_ephemeris
import pandas
from astropy.wcs import WCS
from photutils import DAOStarFinder
from photutils import IRAFStarFinder

from astropy.stats import sigma_clipped_stats
from astropy.visualization import (MinMaxInterval, SqrtStretch,
                                   ImageNormalize, LogStretch)

from photutils import CircularAperture

def wfc3Trace(xc, yc, xref=522, yref=522, xref_grism=522, yref_grism=522,
                   subarray=256, subarray_grism=256):
    """
    """
    # adjust position in case different subarrays are used.
    xc = xc - (xref - xref_grism)
    yc = yc - (yref - yref_grism)

    coord0 = (1014 - subarray) // 2
    xc = xc + coord0
    yc = yc + coord0

    DYDX0 = [2.03396481352,
             -0.00019752130624389416,
             -0.002202066565067532,
             3.143514082596283e-8,
             4.3212786880932414e-7,
             1.210435999122636e-7]

    DYDX1 = [0.010205281672977665,
             -6.06056923866002e-6,
             -3.2485600412356953e-6,
             4.2363866304617406e-10,
             1.230956851333159e-8,
             1.6123073931033502e-9]

    dx = np.arange(1014) - xc
    M = np.sqrt(1.0 + (DYDX1[0] + DYDX1[1] * xc + DYDX1[2] * yc +
                DYDX1[3] * xc**2 + DYDX1[4] * xc * yc + DYDX1[5] * yc**2)**2)
    dp = dx * M

    trace = (DYDX0[0] + DYDX0[1] * xc + DYDX0[2] * yc +
             DYDX0[3] * xc**2 + DYDX0[4] * xc * yc + DYDX0[5] * yc**2) + \
        dp * (DYDX1[0] + DYDX1[1] * xc + DYDX1[2] * yc +
              DYDX1[3] * xc**2 + DYDX1[4] * xc * yc + DYDX1[5] * yc**2) / M
    if subarray < 1014:
        i0 = (1014 - subarray) // 2
        trace = trace[i0: i0 + subarray]

    idx_min = (subarray-subarray_grism)//2
    idx_max = (subarray-subarray_grism)//2 + subarray_grism
    trace = trace[idx_min:idx_max]
    return trace + yc - (1014 - subarray_grism) // 2

def wfc3Dispersion(xc, yc, xref=522, yref=522, xref_grism=522, yref_grism=522,
                   subarray=256, subarray_grism=256):
    """
    Convert pixel coordinate to wavelength. Method and coefficient
    adopted from Kuntschner et al. (2009), Wilkins et al. (2014)

    In case the direct image and spectral image are not taken with the
    same aperture, the centroid measurement is adjusted according to the
    table in: http://www.stsci.edu/hst/observatory/apertures/wfc3.html

    Input:
    ------
        xc:
            X coordinate of direct image centroid
        yc:
            Y coordinate of direct image centroid
    xref
    yref
    xref_grism
    yref_grism
    subarray
    subarray_grism

    Output:
    -------
        wavelength:
            return wavelength mapping of x coordinate in micron
    """

    # adjust position in case different subarrays are used.
    xc = xc - (xref - xref_grism)
    yc = yc - (yref - yref_grism)

    coord0 = (1014 - subarray) // 2
    xc = xc + coord0
    yc = yc + coord0

    DYDX1 = [0.010205281672977665,
             -6.06056923866002e-6,
             -3.2485600412356953e-6,
             4.2363866304617406e-10,
             1.230956851333159e-8,
             1.6123073931033502e-9]

    DLDP0 = [8949.40742544, 0.08044032819916265,
             -0.009279698766495334, 0.000021856641668116504,
             -0.000011048008881387708, 0.00003352712538187608]
    DLDP1 = [44.97227893276267,
             0.0004927891511929662,
             0.0035782416625653765,
             -9.175233345083485e-7,
             2.2355060371418054e-7,  -9.258690000316504e-7]
    # calculate field dependent dispersion coefficient
    p0 = DLDP0[0] + DLDP0[1] * xc + DLDP0[2] * yc + \
         DLDP0[3] * xc**2 + DLDP0[4] * xc * yc + DLDP0[5] * yc**2
    p1 = DLDP1[0] + DLDP1[1] * xc + DLDP1[2] * yc + \
         DLDP1[3] * xc**2 + DLDP1[4] * xc * yc + DLDP1[5] * yc**2
    dx = np.arange(1014) - xc

    M = np.sqrt(1.0 + (DYDX1[0] + DYDX1[1] * xc + DYDX1[2] * yc +
                DYDX1[3] * xc**2 + DYDX1[4] * xc * yc + DYDX1[5] * yc**2)**2)
    dp = dx * M
    wavelength = (p0 + dp * p1)
    if subarray < 1014:
        i0 = (1014 - subarray) // 2
        wavelength = wavelength[i0: i0 + subarray]

    idx_min = (subarray-subarray_grism)//2
    idx_max = (subarray-subarray_grism)//2 + subarray_grism
    wavelength = wavelength[idx_min:idx_max] * u.Angstrom
    return wavelength.to(u.micron)


tso = cascade.TSO.TSOSuite()

# reset
tso.execute("reset")

# initialization with ini files
path = cascade.initialize.default_initialization_path
tso = cascade.TSO.TSOSuite("initialize", "cascade_test_HST_cpm.ini",
                           "cascade_test_HST_object.ini",
                           "cascade_test_HST_data_spectra.ini", path=path)

tso.cascade_parameters.instrument

tso.cascade_parameters.instrument_filter
tso.cascade_parameters.instrument_aperture
tso.cascade_parameters.observations_path
tso.cascade_parameters.observations_target_name
tso.cascade_parameters.observations_cal_path
tso.cascade_parameters.instrument_cal_aperture
tso.cascade_parameters.instrument_cal_filter
tso.cascade_parameters.observations_id

tso.cascade_parameters.observations_mode

path_cal = "/home/bouwman/HST_CALIBRATION/WFC3/"
cal_file_name='wavelength_ref_pixel_wfc3.txt'

ptable_cal = pandas.read_table(path_cal+cal_file_name,delim_whitespace=True,
                               low_memory=False, skiprows=1,
                               names=['APERTURE', 'FILTER', 'XREF', 'YREF'])


bla = ptable_cal[(ptable_cal.APERTURE == 'IRSUB128')]
bla.shape[0]
bla = ptable_cal[(ptable_cal.APERTURE == 'IRSUB128') & (ptable_cal.FILTER == 'F139M')]
bla.shape[0]
bla = ptable_cal[(ptable_cal.APERTURE == 'IRSUB128') & (ptable_cal.FILTER.isnull())]
bla.size
bla.XREF.values
bla.YREF.values


path = "/home/bouwman/HST_OBSERVATIONS/WFC3/WASP19b/SPECTRA/"
header_spec = fits.getheader(path+"ibh715dhq_flt_xp_2.SPC.fits",ext=0)

header_spec1 = fits.getheader(path+"ibh715dhq_flt_xp_2.SPC.fits",ext=1)
spectral_data = fits.getdata(path+"ibh715dhq_flt_xp_2.SPC.fits",ext=1)

hdul = fits.open(path+"ibh715dhq_flt_xp_2.SPC.fits")


path = "/home/bouwman/HST_OBSERVATIONS/WFC3/WASP19b/SPECTRAL_IMAGES/"
header_im = fits.getheader(path+"ibh715zlq_flt.fits",ext=0)

spectral_image = fits.getdata(path+"ibh715zlq_flt.fits",ext=1)

hdu = fits.open(path+"ibh715zlq_flt.fits")[1]
ra_target = fits.getval(path+"ibh715zlq_flt.fits", "RA_TARG")
dec_target = fits.getval(path+"ibh715zlq_flt.fits", "DEC_TARG")
target_coord = coord.SkyCoord(ra_target, dec_target,
                              unit=(u.deg, u.deg), frame='icrs')

w = WCS(hdu.header)
lon, lat = w.all_pix2world(30, 40, 0)
print(lon, lat)
target_expected_coordinates = w.all_world2pix(ra_target, dec_target,0)
print(target_expected_coordinates)

mean, median, std = sigma_clipped_stats(spectral_image, sigma=3.0, iters=5)

# Create an ImageNormalize object
norm = ImageNormalize(np.log10(spectral_image),
                      stretch=LogStretch())
plt.imshow(np.log10(spectral_image), cmap='Greys', origin='lower', norm=norm)
plt.show()

# daofind = DAOStarFinder(fwhm=4.0, threshold=10.*std)
iraffind = IRAFStarFinder(fwhm=2.0, threshold=30.*std)
# sources = daofind(spectral_image - median)
sources = iraffind(spectral_image - median )


idx_max = sources['flux'].argmax()
positions = (sources[idx_max]['xcentroid'], sources[idx_max]['ycentroid'])
# positions = (sources['xcentroid'], sources['ycentroid'])
apertures = CircularAperture(positions, r=6.)
plt.imshow(np.log10(spectral_image), cmap='Greys', origin='lower', norm=norm)
apertures.plot(color='red', lw=3.5, alpha=0.5)
plt.show()

wave_cal = wfc3Dispersion(positions[0], positions[1],
                          xref=522, yref=522, xref_grism=522, yref_grism=522,
                          subarray=512, subarray_grism=128)  # * u.Angstrom
#idx_min = (512-128)//2
#idx_max = (512-128)//2 + 128
#print(wave_cal[idx_min:idx_max].to(u.micrometer))
print(wave_cal)
print((hdul[1].data['LAMBDA']*u.Angstrom).to(u.micrometer))

plt.plot(wave_cal)
plt.plot((hdul[1].data['LAMBDA']*u.Angstrom).to(u.micrometer)[1:])
plt.show()
plt.plot(wave_cal - (hdul[1].data['LAMBDA']*u.Angstrom).to(u.micrometer)[1:])
plt.show()

trace = wfc3Trace(positions[0], positions[1],
                  xref=522, yref=522, xref_grism=522, yref_grism=522,
                  subarray=512, subarray_grism=128)

plt.plot(trace)
plt.show()

# make list of all spectral images
def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return sorted(result)


data_files = find('*.fits', path)

# get the data
im_couter_corr = 0
nintegrations = len(data_files)
image_cube = np.zeros((256,256, nintegrations))
#image_cube = np.zeros((128, 128, nintegrations))
time = np.zeros((nintegrations))
for im, image_file in enumerate(data_files):
    instrument_fiter = fits.getval(image_file, "FILTER", ext=0)
    instrument_aperture = fits.getval(image_file, "APERTURE", ext=0)
    print(instrument_fiter, instrument_aperture)
    if instrument_fiter != 'G141':
        im_couter_corr += 1
        print(instrument_fiter)
        break
        continue
    # WARNING fits data is single precision!!
    spectral_image = fits.getdata(image_file, ext=1)
    if spectral_image.shape[0] != 256:
        im_couter_corr += 1
        continue
    image_cube[:, :, im-im_couter_corr] = spectral_image
    exptime = fits.getval(image_file, "EXPTIME", ext=0)
    expstart = fits.getval(image_file, "EXPSTART", ext=0)
    time[im-im_couter_corr] = expstart + 0.5*(exptime/(24.0*3600.0))

nintegrations = nintegrations-im_couter_corr
image_cube = image_cube[:,:,:-im_couter_corr]
time = time[:-im_couter_corr]
idx = np.argsort(time)
time = time[idx]
image_cube = image_cube[:,:, idx]

header = fits.getheader(data_files[0],ext=0)

ra_target = fits.getval(data_files[0], "RA_TARG")
dec_target = fits.getval(data_files[0], "DEC_TARG")
target_coord = coord.SkyCoord(ra_target, dec_target,
                              unit=(u.deg, u.deg), frame='icrs')
time_obs = Time(time, format='mjd', scale='utc',
                location=('120d', '40d'))
ltt_bary_jpl = time_obs.light_travel_time(target_coord, ephemeris='jpl')
time_barycentre = time_obs.tdb + ltt_bary_jpl
time = time_barycentre.jd

image0 = image_cube[:, :, 0]
shift_store = np.zeros((2, nintegrations))
for it in range(nintegrations):
    # subpixel precision
    shift, error, diffphase = register_translation(image0, image_cube[:, :, it], 100)
    shift_store[:,it] = shift

fig = plt.figure(figsize=(8, 3))
ax1 = plt.subplot(1, 2, 1, adjustable='box-forced')
ax2 = plt.subplot(1, 2, 2, sharex=ax1, sharey=ax1, adjustable='box-forced')
ax1.plot(time, shift_store[0,:])
ax1.set_title('Y offset')
ax2.plot(time, shift_store[1,:])
ax2.set_title('X offset')
plt.show()


image1 = image_cube[:, :, 0]
image2 = image_cube[:, :, 100]

fig = plt.figure(figsize=(4, 4))
ax1 = plt.subplot(1, 1, 1, adjustable='box-forced')
ax1.imshow(image1, cmap='gray')
ax1.set_axis_off()
ax1.set_title('Reference image + Trace')
ax1.plot(trace, '-r')
plt.show()

plt.imshow(image1)
plt.show()
plt.imshow(image2)
plt.show()

# pixel precision first
shift, error, diffphase = register_translation(image1, image2)

fig = plt.figure(figsize=(8, 3))
ax1 = plt.subplot(1, 3, 1, adjustable='box-forced')
ax2 = plt.subplot(1, 3, 2, sharex=ax1, sharey=ax1, adjustable='box-forced')
ax3 = plt.subplot(1, 3, 3)

ax1.imshow(image1, cmap='gray')
ax1.set_axis_off()
ax1.set_title('Reference image')

ax2.imshow(image2, cmap='gray')
ax2.set_axis_off()
ax2.set_title('Offset image')

# Show the output of a cross-correlation to show what the algorithm is
# doing behind the scenes
image_product = np.fft.fft2(image1) * np.fft.fft2(image2).conj()
cc_image = np.fft.fftshift(np.fft.ifft2(image_product))
ax3.imshow(cc_image.real)
ax3.set_axis_off()
ax3.set_title("Cross-correlation")

plt.show()

print("Detected pixel offset (y, x): {}".format(shift))

# subpixel precision
shift, error, diffphase = register_translation(image1, image2, 100)

fig = plt.figure(figsize=(8, 3))
ax1 = plt.subplot(1, 3, 1, adjustable='box-forced')
ax2 = plt.subplot(1, 3, 2, sharex=ax1, sharey=ax1, adjustable='box-forced')
ax3 = plt.subplot(1, 3, 3)

ax1.imshow(image1, cmap='gray')
ax1.set_axis_off()
ax1.set_title('Reference image')

ax2.imshow(image2, cmap='gray')
ax2.set_axis_off()
ax2.set_title('Offset image')

# Calculate the upsampled DFT, again to show what the algorithm is doing
# behind the scenes.  Constants correspond to calculated values in routine.
# See source code for details.
cc_image = _upsampled_dft(image_product, 150, 100, (shift*100)+75).conj()
ax3.imshow(cc_image.real)
ax3.set_axis_off()
ax3.set_title("Supersampled XC sub-area")


plt.show()

print("Detected subpixel offset (y, x): {}".format(shift))

from skimage.transform import (hough_line, hough_line_peaks,
                               probabilistic_hough_line)
from matplotlib import cm
h, theta, d = hough_line(image1)
print(h,theta, d)
for _, angle, dist in zip(*hough_line_peaks(h, theta, d)):
    y0 = (dist - 0 * np.cos(angle)) / np.sin(angle)
    y1 = (dist - image1.shape[1] * np.cos(angle)) / np.sin(angle)
    print(y0,y1)

# Generating figure 1
fig, axes = plt.subplots(1, 3, figsize=(7.5, 3),
                         subplot_kw={'adjustable': 'box-forced'})
ax = axes.ravel()
ax[0].imshow(image1, cmap=cm.gray)
ax[0].set_title('Input image')
ax[0].set_axis_off()

ax[1].imshow(np.log(1 + h),
             extent=[np.rad2deg(theta[-1]), np.rad2deg(theta[0]), d[-1], d[0]],
             cmap=cm.gray, aspect=1/1.5)
ax[1].set_title('Hough transform')
ax[1].set_xlabel('Angles (degrees)')
ax[1].set_ylabel('Distance (pixels)')
ax[1].axis('image')

ax[2].imshow(image1, cmap=cm.gray)
for _, angle, dist in zip(*hough_line_peaks(h, theta, d)):
    y0 = (dist - 0 * np.cos(angle)) / np.sin(angle)
    y1 = (dist - image1.shape[1] * np.cos(angle)) / np.sin(angle)
    ax[2].plot((0, image1.shape[1]), (y0, y1), '-r')
ax[2].set_xlim((0, image1.shape[1]))
ax[2].set_ylim((image1.shape[0], 0))
ax[2].set_axis_off()
ax[2].set_title('Detected lines')

plt.tight_layout()
plt.show()

from skimage.feature import canny
# Line finding using the Probabilistic Hough Transform
image = image1.copy()
edges = canny(image, 2, 1, 25)
edges = canny(image, 2,0.2,0.95, use_quantiles=True)
lines = probabilistic_hough_line(edges, threshold=10, line_length=5,
                                 line_gap=0)

# Generating figure 2
fig, axes = plt.subplots(1, 3, figsize=(7.5, 2.5), sharex=True, sharey=True)
ax = axes.ravel()

ax[0].imshow(image, cmap=cm.gray)
ax[0].set_title('Input image')

ax[1].imshow(edges, cmap=cm.gray)
ax[1].set_title('Canny edges')

ax[2].imshow(edges * 0)
for line in lines:
    p0, p1 = line
    ax[2].plot((p0[0], p1[0]), (p0[1], p1[1]))
ax[2].set_xlim((0, image.shape[1]))
ax[2].set_ylim((image.shape[0], 0))
ax[2].set_title('Probabilistic Hough')

for a in ax:
    a.set_axis_off()
    a.set_adjustable('box-forced')

plt.tight_layout()
plt.show()

