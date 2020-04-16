# from urllib.request import urlretrieve
import pickle
import os
import astropy.units as u
import numpy as np
from astropy.table import MaskedColumn

os.environ["CASCADE_WARNINGS"] = 'off'

from cascade.exoplanet_tools import Kmag
from cascade.initialize import cascade_default_data_path
from cascade.exoplanet_tools import parse_database
from cascade.exoplanet_tools import extract_exoplanet_data

# catalog_name = "EXOPLANETS_A"
catalog_name1 = "NASAEXOPLANETARCHIVE"
catalog1 = parse_database(catalog_name1, update=True)
catalog_name2 = "TEPCAT"
catalog2 = parse_database(catalog_name2, update=True)
catalog_name3 = "EXOPLANETS_A"
catalog3 = parse_database(catalog_name3, update=True)
catalog_name4 = "EXOPLANETS.ORG"
catalog4 = parse_database(catalog_name3, update=True)

cat_dict = {catalog_name1: catalog1,
            catalog_name2: catalog2,
            catalog_name3: catalog3,
            catalog_name4: catalog4}
primary_cat = catalog_name3


path = os.path.join(cascade_default_data_path,
                    "data/HST/archive_data/")
hst_data = pickle.load(open(path+'WFC3_files.pickle', 'rb'))


def remove_space(object_names):
    """
    Remove spaces.

    Parameters
    ----------
    object_names : TYPE
        DESCRIPTION.

    Returns
    -------
    new_name : TYPE
        DESCRIPTION.

    """
    new_name = np.asarray(object_names)
    for i in range(new_name.size):
        new_name[i] = (''.join(object_names[i].split(' ')))
    return new_name


def remove_duplicate(object_names):
    """
    Remove duplicates.

    Parameters
    ----------
    object_names : TYPE
        DESCRIPTION.

    Returns
    -------
    new_name : TYPE
        DESCRIPTION.

    """
    new_name = np.asarray(object_names)
    new_name = np.unique(new_name)
    return new_name


all_observed_planets = [i['planet'] for i in hst_data.values()]
all_observed_planets = remove_duplicate(all_observed_planets)

observables = ['RSTAR', 'TEFF', 'R', 'TT', 'PER', 'ECC', 'A', 'OM',
               'I', 'KMAG', 'LOGG', 'FE']
observables_units = [u.Rsun, u.K, u.Rjup, u.day, u.day,
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
        except (ValueError, KeyError) as e:
            print("No match in {}".format(cat))
            print(e)
            no_match = True
            dr = [{}]
        search_results[cat] = {'no_match': no_match, 'record': dr}
    # check if system is found in any catalog
    if np.all([i['no_match'] for i in search_results.values()]):
        raise ValueError("Planet not found")
    # does the primary found the target?  if not pick another
    if not search_results[primary_cat]['no_match']:
        dr = search_results[primary_cat]['record'][0].copy()
    else:
        for cat in search_results:
            if not cat == primary_cat:
                if not search_results[cat]['no_match']:
                    dr = search_results[cat]['record'][0].copy()
                    break
    # check for missing observable
    for observable, unit in zip(observables, observables_units):
        if observable not in dr.keys():
            dr[observable] = MaskedColumn(data=[0.0], mask=True, unit=unit)
            for cat in search_results:
                if observable in search_results[cat]['record'][0].keys():
                    dr[observable] = \
                        search_results[cat]['record'][0][observable]
                    break
    # check for missing values
    for observable in observables:
        if dr[observable].mask[0]:
            for cat in search_results:
                if observable in search_results[cat]['record'][0].keys():
                    if not search_results[cat]['record'][0][observable].mask[0]:
                        dr[observable] = \
                             search_results[cat]['record'][0][observable]
                        break
    # Make sure the units are those we use as standard units.
    values = [dr[key].filled(0.0)[0] * dr[key].unit for key in observables]
    for i, (value, unit) in enumerate(zip(values, observables_units)):
        values[i] = value.to(unit)
    values = tuple(values)
    return values


# loop over all entries in HST observations catalog
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
