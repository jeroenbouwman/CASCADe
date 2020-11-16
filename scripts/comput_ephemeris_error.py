# CHECK THE EPEMERID
#  R Gastaud   14 mai 2020

import os
import numpy as np
import requests
import shutil
import pandas as pd
import pickle
from astropy.io import fits
import astropy.units as u
#
import cascade # just to get the path !
#
from scripts.write_parameters_for_ini_files import get_catalogue
from scripts.write_parameters_for_ini_files import reformat_nasa_exoplanet_archive_data
from scripts.write_parameters_for_ini_files import get_hst_parameters
#
############################################################################################
#  This function is copied from the function cascade.build_archive import return_header_info
def get_header_info(data_file):
    
    fits_keywords = ['TARGNAME', 'RA_TARG', 'DEC_TARG', 'PROPOSID',
                     'TELESCOP', 'INSTRUME', 'FILTER', 'APERTURE',
                     'NSAMP', 'EXPTIME','EXPSTART']
        
    url_data_archive = \
            'https://mast.stsci.edu/portal/Download/file/HST/product/{0}'
    temp_download_dir = 'temporary/'
 
    os.makedirs(temp_download_dir, exist_ok=True)
    ###  download data
    df = requests.get(url_data_archive.format(data_file), stream=True)
    with open(os.path.join(temp_download_dir, data_file), 'wb') as file:
        for chunk in df.iter_content(chunk_size=1024):
            file.write(chunk)
    # read header info
    header_info = {}
    with fits.open(os.path.join(temp_download_dir, data_file)) as hdul:
        for keyword in fits_keywords:
            header_info[keyword] = hdul['PRIMARY'].header[keyword]

    shutil.rmtree(temp_download_dir)
    return header_info

############################################################################################
#  CHECK THE EPEMERID for one example
def comput_ephemeris_error(period, period_error, ephemeris, expstart, exptime, verbose=False):
#    compare
    if(exptime < period):
        print('Good :observation covers less than one period ')
    else:
        print('Bad : observation covers more than one period ', exptime/period)
    #
    #print('MJD = JD - 2400000.5 ')
    delta_t = expstart + 2400000.5*u.day - ephemeris
    k = (delta_t/(period)).decompose().value

    ### need the error on the period to compute the error on the final phase ...
    phase_error = (period_error*k/period).decompose().value
    if( verbose):
        print('delta_t', delta_t, 'number of periods', k,'phase error=', phase_error)
    return phase_error, k

#exp_start = hdul[0].header['EXPSTART']*u.d
#exp_time = hdul[0].header['EXPTIME']*u.s
#period_error = 3e-07*u.d

############################################################################################
#  CHECK THE EPEMERID for a list of planets


def main(short_list, out_filename, verbose=False):
    root_dir = os.path.dirname(cascade.__path__[0])
    ###  first read the HST catalog
    ## see write_parameter_for_ini_files.get_hst_parameters(path_catalogue, name_catalogue)
    list_planets_hst_wfc3_path = os.path.join(root_dir, 'data/archive_databases/HST/WFC3')
    list_planets_hst_wfc3_filename = 'WFC3_files.pickle'
    hst_catalog_file =  os.path.join(list_planets_hst_wfc3_path, list_planets_hst_wfc3_filename )
    with open(hst_catalog_file, 'rb') as f:
        hst_data_catalog = pickle.load(f)

    ## now read the NASA catalog
    path_catalogue_nasa_exoplanet_archive = os.path.join(root_dir, 'data/exoplanet_data/NASAEXOPLANETARCHIVE')
    name_catalogue_nasa_exoplanet_archive = 'nasaexoplanetarchive_full.csv'
    url_catalogue_nasa_exoplanet_archive = ''
    nea_catalogue = get_catalogue(path_catalogue_nasa_exoplanet_archive,
                                  name_catalogue_nasa_exoplanet_archive,
                                  url_catalogue=url_catalogue_nasa_exoplanet_archive, update=False)
    catalogue_data = reformat_nasa_exoplanet_archive_data(nea_catalogue)

    planet_name = catalogue_data['PLANET'].values
    #stripped_planet_name = np.char.replace(planet_name.tolist(), " ", "")
    period = catalogue_data['PER'].values
    period_error = catalogue_data['PERUPPER'].values
    ephemeris =  catalogue_data['TT'].values
    ephemeris_error =  catalogue_data['TTUPPER'].values

    if (len(short_list) == 0):
        print(' no short list, take all planets in hst')
        hst_data = get_hst_parameters(list_planets_hst_wfc3_path, list_planets_hst_wfc3_filename)
        short_list = hst_data['PLANET'].values
    #
    short_catalog = {}
    for name in short_list:
        ii = np.where(planet_name == name)
        if(len(ii[0]) >0):
            i1 = ii[0][0]
            #print(name, planet_name[i1], period[i1], period_error[i1], ephemeris[i1])
            short_catalog[name]= (period[i1], period_error[i1], ephemeris[i1], ephemeris_error[i1])
        else:
            print('not found1', name)

    counter = 0
    for key in short_catalog:
        if(verbose): print(key, short_catalog[key])
        counter = counter+1

    f = open(out_filename,"w+")
    f.write("## planet_name, visit, observation, phase_error, number of periods \n")

    for name in short_catalog:
        mes_visites = []
        for visit in hst_data_catalog:
            if(hst_data_catalog[visit]['planet'] == name):
                mes_visites.append(visit)
                obs = hst_data_catalog[visit]
                filename = hst_data_catalog[mes_visites[0]]['observations_id_raw'].split(',')[0]
                # obs['observations_id_ima'].split(',')[0], obs['calibrations_id'].split(',')[0]
                if (verbose): print(visit, obs['planet'], obs['observation'],  filename)
                info = get_header_info(filename)
                a_period = short_catalog[name][0]*u.day
                a_period_error = short_catalog[name][1]*u.day
                a_ephemeris = short_catalog[name][2]*u.day
                expstart = info['EXPSTART']*u.day
                exptime = info['EXPTIME']*u.second
                phase_error, k = comput_ephemeris_error(a_period, a_period_error, a_ephemeris, expstart, exptime, verbose=verbose)
                #logging.info(visit, obs['planet'], obs['observation'])
                #f.write(visit+", "+obs['planet']+", "+ obs['observation']+ ", phase_error={:.2e}, k={:.2e} \n".format(phase_error, k))
                chaine = obs['planet']+", "+visit+", "+ obs['observation']+ ", {:.2e}, {:.0f} \n".format(phase_error, np.round(k))
                if(verbose): print(chaine)
                f.write(chaine)
                print(" ")
            else:
                print("not found2", name, visit)
    f.close()
    return

if __name__ == "__main__":
    short_list = ['K2-3d', 'HD97658b', 'HAT-P-26b', 'HAT-P-38b', 'HD3167c', 'HIP41378b', 'Kepler-16b', 'Kepler-51d',
                  'Kepler-9b', 'Kepler-9c', 'WASP-17b', '55Cnce']
    short_list = []
    out_filename = "guru99.txt"
    main(short_list, out_filename, verbose=False)






