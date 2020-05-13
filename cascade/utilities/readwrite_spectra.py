#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This is an utilty file for easy testing :
#     reading is already coded in get_spectra in HST.py
#           and in utilities get_data_from_fits
#     writing is done in utilities  write_timeseries_to_fits
# But they are not easy to use
#
import numpy as np
import os
import glob
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table, Column
import astropy.units as u
import datetime
import matplotlib.pyplot as plt


##########################  READ/WRITE/COMPARE GENERIC SPECTRA #########
def read_spectra_dir(in_dir, pattern='*.fits', verbose=False):
    '''
    Stacks cascade generic 1D spectra  into a single matrix (time, wavelength),
    concatenate bjd, check that the wavelength vector is the same for all spectra.
    
    :param in_dir: Directory where fits files  are stored.
    :param pattern: pattern , e.g.  '*.fits'
    :param verbose: If True: print what is doing
    :type in_dir: str
    :type pattern: str;
    :type verbose: bool
    
    :return: wavelength , micron
    :rtype: ndarray; dim: (nw); value in electron/s
    :return: Matrix of light curves
    :rtype: ndarray; dim: (nw, nt); value in electron/s
    :return: Matrix of error
    :rtype: ndarray; dim: (nw nt); value in electron/s
    :return: Matrix of mask
    :rtype: ndarray; dim: (nw nt); value zero if good
    :return: Matrix of time unit day
    :rtype: ndarray; dim: (nt);
    '''
    data_files = glob.glob(os.path.join(in_dir, pattern))
    n_file = data_files.__len__()
    if (verbose):print(n_file)
    # read the first to get the shape
    filename = data_files[0]
    hdulist = fits.open(filename)
    waves0 =  hdulist[1].data['lambda']
    nw = waves0.size
    data_names =  hdulist[1].data.names
    hdr = hdulist[1].header
    i = data_names.index('LAMBDA')
    unit_lambda = u.Unit(hdr['TUNIT{:01}'.format(i+1)])
    
    i = data_names.index('FLUX')
    unit_flux = u.Unit(hdr['TUNIT{:01}'.format(i+1)])

    i = data_names.index('FERROR')
    unit_ferror = u.Unit(hdr['TUNIT{:01}'.format(i+1)])

    time_unit_name = hdulist[0].header.get('TBJDUNIT')
    if ( time_unit_name is None): time_unit_name='day'
    unit_time = u.Unit(time_unit_name)
    hdulist.close()
    #
    flux = np.zeros([nw, n_file])
    ferror = np.zeros([nw, n_file])
    mask = np.full([nw, n_file], True, dtype=bool)
    times = np.zeros([n_file])
    i = 0
    for filename in data_files:
        if (verbose): print(i, filename)
        # see also https://docs.astropy.org/en/stable/generated/examples/io/fits-tables.html
        # table =  Table.read(filename, hdu=1)
        # flux = events['FLUX']  # it is a column with unit
        hdulist = fits.open(filename)
        times[i]=hdulist[0].header['TIME_BJD']
        waves = hdulist[1].data['lambda']
        if not(np.isclose(waves, waves0, atol=1e-15).all()):
            diff = np.abs(waves-waves0)
            raise Exception('waves ne waves0: {}'.format(diff.max()))
        flux[:,i] = hdulist[1].data['FLUX']
        ferror[:,i] = hdulist[1].data['FERROR']
        mask[:,i] = hdulist[1].data['MASK']
        hdulist.close()
        i = i+1

    # now sort in time
    idx = np.argsort(times)
    times = times[idx]*unit_time
    flux = flux[:, idx]*unit_flux
    mask = mask[:, idx]
    mx = np.ma.masked_array(flux, mask)
    waves = waves*unit_lambda
    ferror = ferror[:, idx]*unit_ferror

    return mx, ferror, waves, times
# waves*unit_lambda, flux*unit_flux, ferror*unit_ferror, mask, times*unit_time

def write_spectra_dir(in_dir, pattern, wavelength, flux, ferror, mask, times, verbose=False, overwrite=False):
    '''
    Write several cascade generic 1D spectra
    :param in_dir: str, name of the output directory
    :param pattern: str, pattern of the name of the output files
    :param wavelength: ndarray: dim: (nw)  quantity (with unit)
    :param flux: ndarray: dim: (nw, nt)  quantity (with unit)
    :param error: ndarray: dim: (nw, nt) quantity (with unit)
    :param mask: ndarray: dim: (nw, nt) bolean, False means good value 
    :param time: ndarray: dim: (nt)   numpy array,quantity (with unit), barycenter julian date
    :param verbose: If True: print what is doing
    '''
    #
    nw, n_file = flux.shape
    my_names = ['LAMBDA', 'FLUX', 'FERROR', 'MASK']
    # should be UPPER CASE, and be carreful lambda is a python language keyword

    if (wavelength.shape[0] != nw):
        print('error with wavelength shape', wavelength.shape, nw )
        return -1
    if not(wavelength.unit.is_equivalent(u.meter)):
        print('error with wavelength.unit', wavelength.unit)
        return -1
    
    if (times.shape[0] != n_file):
        print('error with times', times.shape, n_file)
        return -1
    if not(times.unit.is_equivalent(u.day)):
        print('error with times.unit', times.unit)
        return -1
    bjd = times.to(u.day)
    
    for i in np.arange(n_file):
        filename = os.path.join(in_dir, pattern+'{:03}.fits'.format(i))
        print(i, filename)
        #t = Table([wavelength, flux[:,i], ferror[:,i], mask[:,i]], names = my_names)
        col1 = fits.Column(name='FLUX'   , format='D', array=flux[:,i].value,  unit=flux.unit.to_string())
        col2 = fits.Column(name='FERROR' , format='D', array=ferror[:,i].value, unit=ferror.unit.to_string())
        col3 = fits.Column(name='LAMBDA' , format='D', array=wavelength.value, unit=wavelength.unit.to_string())
        col4 = fits.Column(name='MASK'   , format='L', array=mask[:,i])
        #t.meta['TIME_BJD'] = times[i]
        primary_hdu = fits.PrimaryHDU()
        primary_hdu.header['TIME_BJD'] = bjd.data[i]
        table_hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4])
        #t.write(filename, format='fits')
        hdul = fits.HDUList([primary_hdu, table_hdu])
        hdul.writeto(filename, overwrite=overwrite)
        i = i+1
    return

def compare_spectra(flux1, ferror1, wavelength1, times1, flux2, ferror2, wavelength2, times2, rtol=5e-8,atol=0, verbose=False):
    '''
    Compare several cascade generic 1D spectra
    :param flux1: ndarray: dim: (nw, nt)  masked array of quantity (with unit)
    :param error1: ndarray: dim: (nw, nt) quantity (with unit)
    :param wavelength1: ndarray: dim: (nw)  quantity (with unit)
    :param time1: ndarray: dim: (nt)   numpy array, unit day
    #
    :param flux2: ndarray: dim: (nw, nt)  masked array of quantity (with unit)
    :param error2: ndarray: dim: (nw, nt) quantity (with unit)
    :param wavelength2: ndarray: dim: (nw)  quantity (with unit)
    :param time2: ndarray: dim: (nt)   numpy array, unit day
    :param verbose: If True: print what is doing
    '''
    #  atol has to be set too !
    # absolute(a - b) <= (atol + rtol * absolute(b))
    result = True
    message = ''
    if not(np.allclose(wavelength1.value, wavelength2.value, rtol=rtol, atol=atol)):
        message += "error on wavelength "
        result = False
    if(wavelength1.unit != wavelength2.unit):
        result = False
        message += "error on wavelength.unit "
    
    if not(np.allclose(flux1.data.value, flux2.data.value, rtol=rtol, atol=atol)):
        message += "error on flux "
        result = False
    if(flux1.data.unit != flux2.data.unit):
        result = False
        message += "error on flux.unit "

    if not(np.allclose(ferror1.value, ferror2.value, rtol=rtol, atol=atol)):
        message += "error on ferror "
        result = False
    if(ferror1.unit != ferror2.unit):
        message += "error on ferror.unit "
        result = False

    if not(np.array_equal(flux1.mask, flux2.mask)):
        message += "error on mask "
        result = False
    
    
    if not(np.allclose(times1.value, times2.value, rtol=rtol, atol=atol)):
        message += "error on times "
        result = False
    if(times1.unit != times2.unit):
        result = False
        message += "error on times.unit "

    if(result): message='OK'
    return result, message

##########################  READ/WRITE DEPTH ##############
######  depth (absorption) can be exonoodle input or cascade output ####
def write_depth_ecsv(fileName, array_depth, array_wl):
    col1 = Column(array_wl, name='wavelength', unit=u.micron, dtype='float32')
    col2 = Column(array_depth, name='depth', dtype='float32')
    
    t = Table([col1, col2])
    t.write(fileName, format='ascii.ecsv')
    ff = open(fileName, 'a')
    ff.write('# DATE = ' + str(datetime.datetime.now()) + '\n')
    ff.close()
    return

def write_depth_fits(fileName, array_depth, array_wl, mask=None, meta=None):
    col1 = Column(array_wl, name='wavelength', unit=u.micron, dtype='float32')
    col2 = Column(array_depth, name='depth', dtype='float32')
    t = Table([col1, col2])
    if mask is not None:
        col3 = Column(mask, name='mask', dtype='unit32')
        t = Table([col1, col2, col3])
    if meta is not None:
        t.meta=meta
    t.meta['DATE'] = str(datetime.datetime.now())
    #
    t.write(fileName, format='fits')
    return


def read_depth_ecsv(file_name):
    '''
        Read depths file used by ExoNoodle.
        
        :param file_name: Name of the file. It has to be calibrated the same way as for ExoNoodle
        :return: array_depth = f(array_wl)
        '''
    array_wl, array_depth = np.asarray(ascii.read(file_name)['wavelength']), np.asarray(ascii.read(file_name)['depth'])
    ## depth and not absorption
    return array_wl, array_depth


####
#################  PLOT ###################################

def plot_spectrum_2d(spectrum, wavelengths, times,  title, kw=None, normalise=False):
    nw, nt = spectrum.shape
    if (normalise):
        image = spectrum/(spectrum.mean(1)).reshape([nw, 1])
        my_title = title+ ' normalised image'
    else:
        image = spectrum
        my_title = title
    if(kw is None): kw = nw//2
    #kw = 375 # MaxLine
    fig, (ax1, ax2, ax3 ) = plt.subplots(3) #plt.subplots(1,3) #p
    #plt.title(title)

    fig.suptitle(my_title)
    ax1.imshow(image, origin='lower')
    ax1.set(xlabel='time', ylabel='line <==> wavelength')
    #
    #plt.title( title +' time fixed')
    title2 = 'cut at time= {0:8.3f}'.format(times[nt//2] )
    ax2.set_title(title2)
    ax2.plot(wavelengths, spectrum[:, nt//2])
    ax2.set(xlabel='wavelength', ylabel='flux')
   #plt.title(title+'line fixed')
    print(kw, wavelengths[kw])
    title3 = 'cut at wavelength= {0:8.3f}'.format(wavelengths[kw] )
    ax3.set_title(title3)
    ax3.plot(times, spectrum[kw, :])
    ax3.set(xlabel='time', ylabel='flux')
    # plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    plt.tight_layout()
    plt.show()
    return
#



