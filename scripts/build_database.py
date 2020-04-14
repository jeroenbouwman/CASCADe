from urllib.request import urlretrieve
import pickle
import os

data = pickle.load(open('wfc3_files.pickle', 'rb'))

for visit in data:
    print('\n', visit)
    if not os.path.isdir(visit):
        os.mkdir(visit)

    for file in data[visit]['calibrations_id'].split(','):
        print(file)
        urlretrieve('https://mast.stsci.edu/portal/Download/file/HST/product/{0}'.format(file),
                    os.path.join(visit, file))

    for file in data[visit]['observations_id_raw'].split(','):
        print(file)
        urlretrieve('https://mast.stsci.edu/portal/Download/file/HST/product/{0}'.format(file),
                    os.path.join(visit, file))

    for file in data[visit]['observations_id_ima'].split(','):
        print(file)
        urlretrieve('https://mast.stsci.edu/portal/Download/file/HST/product/{0}'.format(file),
                    os.path.join(visit, file))
