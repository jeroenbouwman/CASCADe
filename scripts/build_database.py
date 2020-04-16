from urllib.request import urlretrieve
import pickle
import os
import astropy.units as u
import numpy as np
from astropy.table import MaskedColumn
import shutil
from astropy.io import fits

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
with open(path+'WFC3_files.pickle', 'rb') as f:
    hst_data = pickle.load(f)


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

temp_dir_path = os.path.join(cascade_default_data_path,
                             "data/HST/mastDownload/")

fits_keywords = ['TARGNAME', 'RA_TARG', 'DEC_TARG', 'PROPOSID',
                 'TELESCOP', 'INSTRUME', 'FILTER', 'APERTURE',
                 'EXPTIME', 
                 'SCAN_TYP', 'SCAN_RAT', 'SCAN_LEN']


def long_substr(data):
    """
    Find longest common substring.

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.

    Returns
    -------
    substr : TYPE
        DESCRIPTION.

    """
    substr = ''
    if len(data) > 1 and len(data[0]) > 0:
        for i in range(len(data[0])):
            for j in range(len(data[0])-i+1):
                if j > len(substr) and all(data[0][i:i+j] in x for x in data):
                    substr = data[0][i:i+j]
    return substr


# loop over all entries in HST observations catalog
for visit in hst_data:
    name = hst_data[visit]['planet']
    print("Target: {}, Visit: {}".format(name, visit))
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

    if hst_data[visit]['observation'] == 'transit':
        observation_type = 'TRANSIT'
    elif hst_data[visit]['observation'] == 'eclipse':
        observation_type = 'ECLIPSE'
    else:
        observation_type = 'PROBLEM'

    print("Observation type: {} \n".format(observation_type))
    if observation_type == 'PROBLEM':
        continue

    data_files = hst_data[visit]['observations_id_ima'].split(',')
    cal_data_files = hst_data[visit]['calibrations_id'].split(',')
    obsid = [file.split('_')[0] for file in data_files]
    cal_obsid = [file.split('_')[0] for file in cal_data_files]

    common_id = long_substr(obsid)
    print("Common ID: {} \n".format(common_id))

    os.makedirs(temp_dir_path, exist_ok=True)
    urlretrieve('https://mast.stsci.edu/portal/Download/file/HST/product/{0}'.
                format(data_files[0]), os.path.join(temp_dir_path,
                                                    data_files[0]))
    header_info = {}
    with fits.open(os.path.join(temp_dir_path, data_files[0])) as hdul:
        for keyword in fits_keywords:
            header_info[keyword] = hdul['PRIMARY'].header[keyword]
    obs_info = [header_info[key] for key in fits_keywords]

    print("FITS Header Info data: \n "
          "Target: {} \n "
          "RA {} \n "
          "DEC {} \n "
          "Proposal ID {} \n "
          "Telescope: {} \n "
          "Instrument: {} \n "
          "Filter: {} \n "
          "Aperture {} \n "
          "Exposure time: {} \n "
          "Scan type: {} \n "
          "Scan rate: {} \n "
          "Scan length: {} \n"
          " ".format(*tuple(obs_info)))

    if header_info['SCAN_TYP'] == 'N':
        observations_mode = 'STARING'
        observations_data = 'SPECTRAL_IMAGE'
        observations_data_product = 'flt'
    else:
        observations_mode = 'SCANNING'
        observations_data = 'SPECTRAL_CUBE'
        observations_data_product = 'ima'

    print("Mode: {} \n"
          "Data: {} \n"
          "Data product: {} \n"
          " ".format(observations_mode, observations_data,
                     observations_data_product))

    urlretrieve('https://mast.stsci.edu/portal/Download/file/HST/product/{0}'.
                format(cal_data_files[0]), os.path.join(temp_dir_path,
                                                        cal_data_files[0]))
    header_info = {}
    with fits.open(os.path.join(temp_dir_path, cal_data_files[0])) as hdul:
        for keyword in fits_keywords:
            header_info[keyword] = hdul['PRIMARY'].header[keyword]
    cal_info = [header_info[key] for key in fits_keywords]
    print("FITS Header Info aquisition image data: \n "
          "Target: {} \n "
          "RA {} \n "
          "DEC {} \n "
          "Proposal ID {} \n "
          "Telescope: {} \n "
          "Instrument: {} \n "
          "Filter: {} \n "
          "Aperture {} \n "
          "Exposure time: {} \n "
          "Scan type: {} \n "
          "Scan rate: {} \n "
          "Scan length: {} \n"
          " ".format(*tuple(cal_info)))

    shutil.rmtree(os.path.join(temp_dir_path))

                  
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
