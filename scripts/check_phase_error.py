#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy as np
import astropy.units as u

from scripts.write_parameters_for_ini_files import get_catalogue
from scripts.write_parameters_for_ini_files import reformat_nasa_exoplanet_archive_data
from scripts.write_parameters_for_ini_files import reformat_exoplanets_a_data
from scripts.write_parameters_for_ini_files import reformat_exoplanets_org_data
from scripts.write_parameters_for_ini_files import reformat_tepcat_data
from scripts.write_parameters_for_ini_files import merge_catalogue
from scripts.write_parameters_for_ini_files import remove_binary_name
from scripts.write_parameters_for_ini_files import remove_space
from scripts.write_parameters_for_ini_files import write_diff


def get_jeroen_parameters(path_catalogue, name_catalogue):
    """
    Read the file containing the HST observations data
    and extract the list of planet, type of observations and the filter used

    Parameters
    ----------
    path_catalogue : string
    name_catalogue : string

    Returns
    -------
    jeroen_cat : pandas object
         planet name, type of observation and filter used

    """
    catalog = pd.read_csv(os.path.join(path_catalogue, name_catalogue), sep=';')

    data = catalog.rename({'planet_name': 'PLANET'}, axis='columns')

    data['PLANET'] = remove_binary_name(data['PLANET'])
    data['PLANET'] = remove_space(data['PLANET'])

    jeroen_cat = data.sort_values(by=['PLANET'])

    return jeroen_cat


def compute_number_period(period, ephemeris, expstart, exptime, verbose=True):

    if verbose:
        if (exptime < period):
            print('Good :observation covers less than one period')
        else:
            print('Bad : observation covers more than one period', exptime/period)

    delta_t = expstart + 2400000.5*u.day - ephemeris    # 'MJD = JD - 2400000.5'
    nbr_period = (delta_t / period).decompose().value

    if verbose:
        print('delta_t = ', delta_t, 'number of periods = ', np.round(nbr_period))

    return np.round(nbr_period)


def compute_phase_error(period, period_error, nbr_period, verbose=False):

    phase_error = (period_error * nbr_period / period).decompose().value

    if verbose:
        print('phase error = ', phase_error)

    return phase_error


def get_phase_error(visit_cat, param_cat, best_value_only=True):

    visit_planet = visit_cat['PLANET'].values
    visit = visit_cat['visit'].values
    observation = visit_cat['observation'].values
    expstart = visit_cat['expstart'].values * u.day
    exptime = visit_cat['exptime'].values * u.s

    param_planet = param_cat['PLANET'].values
    ephemeris = param_cat['TT'].values * u.day  # primary transit [JD]
    ephemeris_err = param_cat['TTUPPER'].values * u.day  # primary transit [JD]
    period = param_cat['PER'].values * u.day
    period_err = param_cat['PERUPPER'].values * u.day

    result = []

    for i in range(visit.size):
        j = np.where(visit_planet[i] == param_planet)[0]

        nbr_period = compute_number_period(period[j], ephemeris[j], expstart[i], exptime[i], verbose=False)
        phase_error = compute_phase_error(period[j], period_err[j], nbr_period, verbose=False)

        if best_value_only:
            try:
                k = np.nanargmin(abs(phase_error))  # ignore Nan
            except:
                k = 0   # if phase_error contains only Nan values, returns the first item (a NaN value)
        else:
            k = 0

        if 'TTV' in param_cat:
            result.append({'planet_name': visit_planet[i], 'visit': visit[i], 'observation': observation[i],
                           'expstart': expstart[i].value, 'exptime': exptime[i].value,
                           'ephemeris': ephemeris[j].value[k], 'ephemeris_error': ephemeris_err[j].value[k],
                           'period': period[j].value[k], 'period_error': period_err[j].value[k],
                           'nbr_period': nbr_period[k], 'phase_error': phase_error[k],
                           'ttv': param_cat['TTV'][j].values[k]
                           })

        else:
            result.append({'planet_name': visit_planet[i], 'visit': visit[i], 'observation': observation[i],
                           'expstart': expstart[i].value, 'exptime': exptime[i].value,
                           'ephemeris': ephemeris[j].value[k], 'ephemeris_error': ephemeris_err[j].value[k],
                           'period': period[j].value[k], 'period_error': period_err[j].value[k],
                           'nbr_period': nbr_period[k], 'phase_error': phase_error[k]
                           })

    result = pd.DataFrame(result)
    return result


def main(root_dir, output_filename, verbose=False):

    # Jeroen catalogues compilation
    path_jeroen_compilation = os.path.join(root_dir, 'scripts')
    name_jeroen_compilation = 'jeroen_cat_compilation.csv'

    # Exoplanet.a catalogue
    path_catalogue_exoplanets_a = os.path.join(root_dir, 'data/exoplanet_data/EXOPLANETS_A')
    name_catalogue_exoplanets_a = 'exoplanets_a_full.csv'
    url_catalogue_exoplanets_a = 'http://svo2.cab.inta-csic.es/vocats/v2/exostars/cs.php?RA=180.000000&DEC=0.000000&S' \
                                 'R=180.000000&VERB=1&objlist=-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,' \
                                 '22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,' \
                                 '51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,' \
                                 '80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,' \
                                 '106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,' \
                                 '127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,' \
                                 '148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,' \
                                 '169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,' \
                                 '190,191,192,193,194,195,196,197,198,199,200,201&fldlist=-1,3,4,5,6,7,8,9,12,13,14,15,' \
                                 '16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,' \
                                 '45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,' \
                                 '74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,' \
                                 '102,103,104,105,106,107,108,109,110,111,112,113&nocoor=1&format=ascii'

    # Exoplanets.org catalogue
    path_catalogue_exoplanets_org = os.path.join(root_dir, 'data/exoplanet_data/exoplanets.org/')
    name_catalogue_exoplanets_org = 'exoplanets.csv'
    url_catalogue_exoplanets_org = 'http://www.exoplanets.org/csv-files/exoplanets.csv'

    # Nasa Exoplanet Archive catalogue
    path_catalogue_nasa_exoplanet_archive = os.path.join(root_dir, 'data/exoplanet_data/NASAEXOPLANETARCHIVE')
    name_catalogue_nasa_exoplanet_archive = 'PS_2020.06.09_02.17.21.csv'  #'nasaexoplanetarchive_full.csv'
    url_catalogue_nasa_exoplanet_archive = ''
    #https://exoplanetarchive.ipac.caltech.edu/workspace/TMP_XTDT5D_5033/TblSearch/2020.06.02_03.25.58_005033/index.html

    # Tepcat allplanet catalogue
    path_allplanet_tepcat = os.path.join(root_dir, 'data/exoplanet_data/tepcat/')
    name_allplanet_tepcat = 'allplanets.csv'
    url_allplanet_tepcat = 'http://www.astro.keele.ac.uk/jkt/tepcat/allplanets-csv.csv'

    # Tepcat observables catalogue
    path_observable_tepcat = os.path.join(root_dir, 'data/exoplanet_data/tepcat/')
    name_observable_tepcat = 'observables.csv'
    url_observable_tepcat = 'http://www.astro.keele.ac.uk/jkt/tepcat/observables.csv'

    # reading catalogues
    jeroen_compilation = get_jeroen_parameters(path_jeroen_compilation, name_jeroen_compilation)

    exoplanet_a_catalogue = get_catalogue(path_catalogue_exoplanets_a,
                                          name_catalogue_exoplanets_a,
                                          url_catalogue=url_catalogue_exoplanets_a, update=False)

    exoplanet_org_catalogue = get_catalogue(path_catalogue_exoplanets_org,
                                            name_catalogue_exoplanets_org,
                                            url_catalogue=url_catalogue_exoplanets_org, update=False)

    nea_catalogue = get_catalogue(path_catalogue_nasa_exoplanet_archive,
                                  name_catalogue_nasa_exoplanet_archive,
                                  url_catalogue=url_catalogue_nasa_exoplanet_archive,
                                  skiprows=303, update=False)

    tepcat_allplanets_catalogue = get_catalogue(path_allplanet_tepcat,
                                                name_allplanet_tepcat,
                                                url_catalogue=url_allplanet_tepcat, update=False)

    tepcat_observable_catalogue = get_catalogue(path_observable_tepcat,
                                                name_observable_tepcat,
                                                url_catalogue=url_observable_tepcat, update=False)

    # getting the list of planets
    list_planets_hst = jeroen_compilation.drop_duplicates('PLANET')['PLANET']

    # reformating catalogues
    exo_a_data = reformat_exoplanets_a_data(exoplanet_a_catalogue)
    exo_org_data = reformat_exoplanets_org_data(exoplanet_org_catalogue)
    nea_data = reformat_nasa_exoplanet_archive_data(nea_catalogue)
    tepcat_data = reformat_tepcat_data(tepcat_allplanets_catalogue, tepcat_observable_catalogue)

    exo_a_merged_hst = merge_catalogue(list_planets_hst, exo_a_data, verbose=False)
    exo_org_merged_hst = merge_catalogue(list_planets_hst, exo_org_data, verbose=False)
    nea_merged_hst = merge_catalogue(list_planets_hst, nea_data, verbose=False)
    tepcat_merged_hst = merge_catalogue(list_planets_hst, tepcat_data, verbose=False)

    # Write merged catalogue with parameters comparison
    #write_diff(exo_a_merged_hst, exo_org_merged_hst, nea_merged_hst, tepcat_merged_hst, output_filename='diff')

    # computing phase error
    exo_a = get_phase_error(jeroen_compilation, exo_a_merged_hst, best_value_only=True)
    exo_org = get_phase_error(jeroen_compilation, exo_org_merged_hst, best_value_only=False)
    nea = get_phase_error(jeroen_compilation, nea_merged_hst, best_value_only=True)
    tepcat = get_phase_error(jeroen_compilation, tepcat_merged_hst, best_value_only=False)

    # Writing results
    writer = pd.ExcelWriter(output_filename + '.xlsx', engine='xlsxwriter')

    exo_a.to_excel(writer, sheet_name='exo_a', index=False)
    exo_org.to_excel(writer, sheet_name='exo_org', index=False)
    nea.to_excel(writer, sheet_name='nea', index=False)
    tepcat.to_excel(writer, sheet_name='tepcat', index=False)

    writer.save()


if __name__ == "__main__":
    root_dir = '../'
    # root_dir = os.path.dirname(cascade.__path__[0])   #need to import cascade module
    main(root_dir=root_dir, output_filename='check_cat_phase_error', verbose=False)
