# from urllib.request import urlretrieve
import pickle
import os
import astropy.units as u
from astropy import constants as const
import numpy as np
from astropy.table import MaskedColumn

os.environ["CASCADE_WARNINGS"] = 'off'

from cascade.exoplanet_tools import Kmag
from cascade.initialize import cascade_default_data_path
from cascade.exoplanet_tools import parse_database
from cascade.exoplanet_tools import extract_exoplanet_data

# catalog_name = "EXOPLANETS_A"
catalog_name = "NASAEXOPLANETARCHIVE"
catalog1 = parse_database(catalog_name, update=True)
catalog_name = "TEPCAT"
catalog2 = parse_database(catalog_name, update=True)

cat_dict = {"NASAEXOPLANETARCHIVE": catalog1,
            "TEPCAT":  catalog2}
primary_cat = "NASAEXOPLANETARCHIVE"


path = os.path.join(cascade_default_data_path,
                    "data/HST/archive_data/")
hst_data = pickle.load(open(path+'WFC3_files.pickle', 'rb'))

all_observed_planets = [i['planet'] for i in hst_data.values()]

observables = ['RSTAR', 'TEFF', 'R', 'TT', 'PER', 'ECC', 'A', 'OM',
               'I', 'KMAG', 'LOGG', 'FE']
observables_units = [const.R_sun, u.K, const.R_jup, u.day, u.day,
                     u.dimensionless_unscaled, u.AU, u.deg, u.deg,
                     Kmag, u.dex(u.cm/u.s**2), u.dex]


def return_system_parameters(name, catalogs, observables, observables_units,
                             primary_cat):
    """
    Return observed system parameters.

    Parameters
    ----------
    name : TYPE
        DESCRIPTION.
    catalogs : TYPE
        DESCRIPTION.
    observables : TYPE
        DESCRIPTION.
    observables_units : TYPE
        DESCRIPTION.
    primary_cat : TYPE
        DESCRIPTION.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    values : TYPE
        DESCRIPTION.

    """
    search_results = {}
    for cat in catalogs:
        try:
            dr = extract_exoplanet_data(catalogs[cat],
                                        name, search_radius=2*u.arcsec)
            no_match = False
        except ValueError as e:
            print("No match in {}".format(cat))
            print(e)
            no_match = True
            dr = []
        search_results[cat] = {'no_match': no_match, 'record': dr}

    if np.all([i['no_match'] for i in search_results.values()]):
        raise ValueError("Planet not found")

    if not search_results[primary_cat]['no_match']:
        dr = search_results[primary_cat]['record'][0].copy()
    else:
        for cat in search_results:
            if not cat == primary_cat:
                if not search_results[cat]['no_match']:
                    dr = search_results[cat]['record'][0].copy()
                    break
    for observable, unit in zip(observables, observables_units):
        if observable not in dr.keys():
            dr[observable] = MaskedColumn(data=[0.0], mask=True, unit=unit)
            for cat in search_results:
                if observable in search_results[cat]['record'][0].keys():
                    dr[observable] = \
                        search_results[cat]['record'][0][observable]
                    break
    values = tuple([dr[key].filled(0.0)[0] * dr[key].unit
                    for key in observables])
    return values


for visit in hst_data:
    name = hst_data[visit]['planet']
    print(name)
    try:
        results = return_system_parameters(name, cat_dict, observables,
                                           observables_units, primary_cat)
    except ValueError as e:
        print(e)
        continue
    print(" Rstar: {} \n Teff: {} \n Rplanet: {} \n "
          "Ephemeris: {} \n Period: {} \n "
          "Eccentricity: {} \n Semi Major Axis: {} \n "
          "Omega: {} \n Inclination: {} \n "
          "Kmag: {} \n log g: {} \n FE: {} \n"
          " ".format(*results))


# for visit in hst_data:
#     name = hst_data[visit]['planet']
#     print(name)
#     try:
#         dr1 = extract_exoplanet_data(catalog1, name, search_radius=2*u.arcsec)
#         no_match1 = False
#     except ValueError as e:
#         print("No match in NASAEXOPLANETARCHIVE")
#         print(e)
#         no_match1 = True
#     try:
#         dr2 = extract_exoplanet_data(catalog2, name, search_radius=2*u.arcsec)
#         no_match2 = False
#     except ValueError as e:
#         print("No match in TEPCAT")
#         print(e)
#         no_match2 = True
#     if not no_match1:
#         dr = dr1
#     elif not no_match2:
#         dr = dr2
#     else:
#         continue
#     if not no_match2:
#         dr[0]['KMAG'] = dr2[0]['KMAG']
#     else:
#         dr[0]['KMAG'] = MaskedColumn(data=[0.0], mask=True, unit=Kmag)
#     print(" Rstar: {} \n Teff: {} \n Rplanet: {} \n "
#           "Ephemeris: {} \n Period: {} \n "
#           "Eccentricity: {} \n Semi Major Axis: {} \n "
#           "Omega: {} \n Inclination: {} \n "
#           "Kmag: {} \n log g: {} \n FE: {} \n"
#           " ".format(*tuple([dr[0][key].filled(0.0)[0] * dr[0][key].unit
#                             for key in observables])))

        
    # with fits.open(image_file) as hdul:
    #     instrument_fiter = hdul['PRIMARY'].header['FILTER']
    #     nsamp = hdul['PRIMARY'].header['NSAMP']
    #     scanAng = hdul['PRIMARY'].header['SCAN_ANG']
    #     paV3 = hdul['PRIMARY'].header['PA_V3']
    #     isSparse = 'SPARS' in hdul['PRIMARY'].header['SAMP_SEQ']
    #     instrument_beam
    #     inst_aperture
    #     scanspeed
    #     instrument
    #     obs_cal_version

# for visit in data:
#     print('\n', visit)
#     if not os.path.isdir(visit):
#         os.mkdir(visit)

#     for file in data[visit]['calibrations_id'].split(','):
#         print(file)
#         urlretrieve('https://mast.stsci.edu/portal/Download/file/HST/product/{0}'.format(file),
#                     os.path.join(visit, file))

#     for file in data[visit]['observations_id_raw'].split(','):
#         print(file)
#         urlretrieve('https://mast.stsci.edu/portal/Download/file/HST/product/{0}'.format(file),
#                     os.path.join(visit, file))

#     for file in data[visit]['observations_id_ima'].split(','):
#         print(file)
#         urlretrieve('https://mast.stsci.edu/portal/Download/file/HST/product/{0}'.format(file),
#                     os.path.join(visit, file))
