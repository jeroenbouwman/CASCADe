# from urllib.request import urlretrieve
import pickle
import os
import astropy.units as u

os.environ["CASCADE_WARNINGS"] = 'off'

from cascade.initialize import cascade_default_data_path
from cascade.exoplanet_tools import parse_database
from cascade.exoplanet_tools import extract_exoplanet_data

catalog_name = "EXOPLANETS_A"
# catalog_name = "NASAEXOPLANETARCHIVE"

catalog = parse_database(catalog_name, update=True)

path = os.path.join(cascade_default_data_path,
                    "data/HST/archive_data/")
hst_data = pickle.load(open(path+'WFC3_files.pickle', 'rb'))

for visit in hst_data:
    name = hst_data[visit]['planet']
    print(name)
    try:
        dr = extract_exoplanet_data(catalog, name, search_radius=2*u.arcsec)
        print("Ephemeris: {}, Period: {}".format(dr[0]['TT'].data,
                                                 dr[0]['PER'].data))
    except ValueError as e:
        print(e)


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
