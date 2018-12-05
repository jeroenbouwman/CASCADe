#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 10:46:08 2018

@author: bouwman
"""

from astroquery import sha
# from astropy import coordinates as coord
# from astropy import units as u
from matplotlib import pyplot as plt
from astroquery.mast import Observations

####################################################
# Spitzer

rqk_t = sha.query(reqkey=23439616, dataset=1)
url = rqk_t['accessUrl'][0].strip()
img = sha.get_file(url)
print(img[0].header['HMJD_OBS'])

plt.imshow(img[0].data)
url = rqk_t['accessWithAnc1Url'][0].strip()
sha.save_file(url, out_dir='/home/bouwman/', out_name='test')


url = 'http://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/ProductDownload?DATASET=level1&ID=329292604&OPTIONS=anc2'
sha.save_file(url, out_dir='/home/bouwman/', out_name='test')
# changing anc1 to anc2 get the other data products.

###########################
# HST

obsTable = Observations.query_object("WASP12b", radius=".02 deg")
print(obsTable[:10])
print(obsTable['obs_id'][:10])

obs_by_name = obsTable.group_by(['obs_collection', 'instrument_name',
                                 'dataproduct_type', 'filters', 'obs_id',
                                 'obsid', 'proposal_id'])
obs_by_name.groups.keys
obs_by_name.groups.keys['dataproduct_type'] == 'spectrum'

bla = obsTable['instrument_name'].data.data
np.unique(bla)

obsTable = Observations.query_criteria(obs_id=['IBH715*'])
print(obsTable['obsid'][0])
obsids = obsTable['obsid'][:1]
dataProductsByID = Observations.get_product_list(obsids)
print(dataProductsByID)
print(dataProductsByID['productFilename'])
print(dataProductsByID['productSubGroupDescription'])
manifest = Observations.download_products(dataProductsByID,
                                          download_dir='/home/bouwman/',
                                          productSubGroupDescription=['FLT'],
                                          extension="fits")


obsids = obsTable['obsid']
for iobs in obsids:
    dataProductsByID = Observations.get_product_list(iobs)
    bla = ([i for i in dataProductsByID['productFilename'].data.data if 'jif' in i][0])
    bla2 = ([i for i in dataProductsByID['productFilename'].data.data if 'flt.fits' in i][0])
    print(bla, bla2)
