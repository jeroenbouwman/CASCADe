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

#catalog_name = "EXOPLANETS_A"
catalog_name = "NASAEXOPLANETARCHIVE"
catalog1 = parse_database(catalog_name, update=True)
catalog_name = "TEPCAT"
catalog2 = parse_database(catalog_name, update=True)

path = os.path.join(cascade_default_data_path,
                    "data/HST/archive_data/")
hst_data = pickle.load(open(path+'WFC3_files.pickle', 'rb'))

observables = ['RSTAR', 'TEFF', 'R', 'TT', 'PER', 'ECC', 'A', 'OM',
               'I', 'KMAG', 'LOGG', 'FE']

for visit in hst_data:
    name = hst_data[visit]['planet']
    print(name)
    try:
        dr1 = extract_exoplanet_data(catalog1, name, search_radius=2*u.arcsec)
        no_match1 = False
    except ValueError as e:
        print("No match in NASAEXOPLANETARCHIVE")
        print(e)
        no_match1 = True
    try:
        dr2 = extract_exoplanet_data(catalog2, name, search_radius=2*u.arcsec)
        no_match2 = False
    except ValueError as e:
        print("No match in TEPCAT")
        print(e)
        no_match2 = True
    if not no_match1:
        dr = dr1
    elif not no_match2:
        dr = dr2
    else:
        continue
    if not no_match2:
        dr[0]['KMAG'] = dr2[0]['KMAG']
    else:
        dr[0]['KMAG'] = MaskedColumn(data=[0.0], mask=True, unit=Kmag)
    print(" Rstar: {} \n Teff: {} \n Rplanet: {} \n "
          "Ephemeris: {} \n Period: {} \n "
          "Eccentricity: {} \n Semi Major Axis: {} \n "
          "Omega: {} \n Inclination: {} \n "
          "Kmag: {} \n log g: {} \n FE: {} \n"
          " ".format(*tuple([dr1[0][key].filled(0.0)[0] * dr1[0][key].unit
                            for key in observables])))

        
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
