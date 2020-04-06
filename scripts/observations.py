import os
import numpy as np
from astroquery.mast import Observations
import shutil

# this is the stadard format of the data structure in the data base file
# provided by Angelos defining the observational data for a single transit
db={'13467.01':
{'instrument_beam': 'A',
 'planet': 'WASP-12 b',
 'instrument_cal_filter': 'F126N',
 'observation': 'transit',
 'instrument_target': 'WASP-12',
 'instrument_filter': 'G141',
 'instrument_cal_root': 'iccz16010',
 'instrument_observatory': 'HST',
 'instrument': 'WFC3',
 'observations_id':
'iccz16010,iccz16pxq,iccz16pyq,iccz16pzq,iccz16q0q,iccz16q2q,iccz16q3q,\
iccz16q4q,iccz16q6q,iccz16q7q,iccz16q8q,iccz16qaq,iccz16qbq,iccz16qcq,iccz16qeq,\
iccz16qfq,iccz16qgq,iccz16qiq,iccz16qjq,iccz16qlq,iccz16qnq,iccz16qoq,\
iccz16qpq,iccz16qrq,iccz16qsq,iccz16qtq,iccz16qvq,iccz16qwq,iccz16qxq,\
iccz16qzq,iccz16r0q,iccz16r1q,iccz16r3q,iccz16r4q,iccz16r5q,iccz16r7q,\
iccz16r8q,iccz16r9q,iccz16rcq,iccz16rdq,iccz16req,iccz16rgq,iccz16rhq,\
iccz16riq,iccz16rkq,iccz16rlq,iccz16rmq,iccz16roq,iccz16rpq,iccz16rqq,\
iccz16rsq,iccz16rtq,iccz16ruq,iccz16rwq,iccz16rxq,iccz16ryq,iccz16s0q,\
iccz16s2q,iccz16s3q,iccz16s5q,iccz16s6q,iccz16s7q,iccz16s9q,iccz16saq,\
iccz16sbq,iccz16sdq,iccz16seq,iccz16sfq,iccz16shq,iccz16siq,iccz16sjq,\
iccz16slq,iccz16smq,iccz16snq,iccz16spq,iccz16sqq,iccz16ssq,iccz16suq,\
iccz16svq,iccz16swq,iccz16sxq,iccz16syq,iccz16szq,iccz16t0q,iccz16t1q,\
iccz16t2q,iccz16t3q,iccz16t4q,iccz16t5q,iccz16t6q,iccz16t7q,iccz16t8q,\
iccz16t9q,iccz16taq,iccz16tbq',
 'instrument_cal_aperture': 'GRISM256',
 'proposal': '13467',
 'instrument_aperture': 'GRISM256'}}
    
# WASP12b transit
db_id = '13467.01'
print(db[db_id])

# we need to agree on how to define the directory structure of the downloaded
# data. The save path is easy, just define the main path to the directory 
# where all HST data is stored.
save_path = '/home/bouwman/HST_OBSERVATIONS'
# I devided the sub directory to instrument, target name , spectral_images
# (in condtrast to spectra). The target name should be taken from the 
# source name and perhaps the instument cal route? or perhaps the longest sub
# string in tcommon n the observations id's (in this example iccz16) of the
# db_id. Misshien is dat laatste het makkelijkst.
sub_save_path = 'WFC3/WASP12b_SCAN1/SPECTRAL_IMAGES/'

target_path = os.path.join(save_path, sub_save_path)
os.makedirs(target_path, exist_ok=True)

#token='71647c4c0e634852bb3b829d7641a732'
#my_session = Observations.login(token=token)
#sessioninfo = Observations.session_info()

# Need to download both flt and ima data products. The ima will be used for
# scanning observations, the flt data will be used for staring and target
# aquisition images.
product = 'IMA'
calProduct = 'ima'
#product = 'FLT'
#calProduct = 'flt'
cal_id = db[db_id]['instrument_cal_root']
observations_id = (db[db_id]['observations_id']).split(',')
observations_id = observations_id + [cal_id]
observations_id = [i.strip() for i in observations_id]
for obs_id in observations_id:
    obsTable = Observations.query_criteria(obs_id=[obs_id])
    obsid_column = obsTable['obsid'][:1]
    dataProductsOverview = Observations.get_product_list(obsid_column)
    flt_index = np.where(dataProductsOverview['productSubGroupDescription']
                         == product)
    dataProductsByID = dataProductsOverview['obsID'][flt_index]
    dataProductsByObsID = dataProductsOverview['obs_id'][flt_index]
    manifest = \
        Observations.download_products(dataProductsByID,
                                       download_dir=target_path,
                                       productSubGroupDescription=[product],
                                       extension="fits")
    try:
        for id_in_product in dataProductsByObsID:
            shutil.move(os.path.join(target_path, 'mastDownload/HST',
                                     id_in_product,
                                     id_in_product+'_'+calProduct+'.fits'),
                        os.path.join(target_path,
                                     id_in_product+'_'+calProduct+'.fits'))
        shutil.rmtree(os.path.join(target_path, 'mastDownload/'))
    except FileNotFoundError:
        pass
