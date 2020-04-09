#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 15:30:09 2018

@author: bouwman
"""
# import xlsx
from astroquery.simbad import Simbad
import os.path
import numpy as np
import cascade
from cascade.exoplanet_tools import extract_exoplanet_data
import xlrd
import xlwt
from xlutils.copy import copy
from contextlib import closing
from urllib.request import Request
from urllib.request import urlopen
import warnings
import pandas as pd
import openpyxl
from astropy.table import Table
from astropy.table import MaskedColumn
from astropy import units as u
from astroquery import sha
from astropy.coordinates import SkyCoord
from astroquery import open_exoplanet_catalogue as oec
# from astroquery.open_exoplanet_catalogue import findvalue
from requests.exceptions import MissingSchema
import dill as pickle
import exoplanetsa
from exoplanetsa.tools import fovid_to_band
from exoplanetsa.tools import RetrieveAORs
from exoplanetsa.tools import write_orbit_parameters_to_file
from exoplanetsa.tools import read_orbit_parameters_from_file

###################
# input parameters
####################
path_flle = '/home/bouwman/'
file_name = 'EXOPLANETSA_test_HST_list.xlsx'

search_radius = 11.0*u.arcsec

file_name_orbit = 'EXOPLANETSA_test_Orbit_Parameters.xls'

####################################
# Read HST target list
####################################
wb_hst, targets, targets_stripped = \
    exoplanetsa.tools.read_target_list(file_name, path=path_flle)

#########################################################
# create table with orbital parameters and write to file
###########################################################
ParTable = \
    exoplanetsa.tools.get_orbital_parameters('NASAEXOPLANETARCHIVE',
                                             targets, targets_stripped,
                                             search_radius=search_radius)

write_orbit_parameters_to_file(ParTable, file_name_orbit, path=path_flle)

###########################################
# Read orbital parameter file
###########################################
ParTable = read_orbit_parameters_from_file(file_name_orbit, path=path_flle)


##############################
# Get all observations
##############################
searchRegion = 180.0*u.arcsec
searchOffsetWarning = 2.0*u.arcminute
minumum_aor_duration = 60.0*u.minute
spitzer_dict = {}
for itarget, system in enumerate(ParTable):
    print("searching spitzer archive for target: ", targets[itarget])
    c = SkyCoord(system["RA"], system["DEC"], frame='icrs',
                 unit=(u.hourangle, u.deg))
    try:
        query_error = False
        if (np.ma.is_masked(system["INCL"]) |
           np.ma.is_masked(system["OMEGA"])):
            query_error_orbit = True
        else:
            query_error_orbit = False
        ObsLog = RetrieveAORs(system["RA"]+','+system["DEC"],
                              searchRegion.value)
        if type(ObsLog) == bool:
            raise KeyError("Search Failed")
        idx_timeseries_aors = (ObsLog['status'].values == 'observed') & \
            (ObsLog['duration'].values.astype(np.float) >
             minumum_aor_duration.value)
        # aor numbers for target and duriation and offset
        aor_keys_in = ObsLog['aor_key'].values[idx_timeseries_aors]
        aor_offsets_in = np.asarray([np.float(i.split('(')[0])
                                     for i in ObsLog['_search_offset'].
                                     values])[idx_timeseries_aors]
        aor_offsets_flag_in = aor_offsets_in > searchOffsetWarning.value
        aor_duration_in = \
            ObsLog['duration'].values.astype(np.float)[idx_timeseries_aors]

        # initialize lists holding basic properties of observations
        aor_keys_out = []
        aor_offsets_out = []
        aor_offsets_flag_out = []
        aor_duration_out = []
        aor_phase = []
        aor_error_lv1 = []
        aor_tansit_flag = []
        aor_eclipse_flag = []
        aor_phase_curve_flag = []
        aor_band = []
        for iobs, obsid in enumerate(aor_keys_in):
            try:
                query_error_lv1 = [False]
                print('checking: ', obsid)
                rqk_t = sha.query(reqkey=obsid, dataset=1)
                idx_time_sorted = \
                    np.argsort((rqk_t['scet'].data).astype(np.datetime64))
                idx_first = idx_time_sorted[0]
                # wavelength ranges of observations
                bands = np.array(rqk_t['wavelength'].data)[idx_time_sorted]
                bands_uniq = np.unique(bands)
                # Field of view IDs of observations
                fovid = np.array(rqk_t['fovname'].data)[idx_time_sorted]
                fovid_uniq = np.unique(fovid)
                # flag if observations was intented for target
                is_primary_field = \
                    np.array(rqk_t['primaryfield'].data)[idx_time_sorted]
                is_primary_field_uniq = np.unique(is_primary_field)
                print('bands: {}'.format(bands_uniq))
                print('FOVIDs: {}'.format(fovid_uniq))
                print('primary: {}'.format(is_primary_field_uniq))
                # check bands agains fovids
                if (len(bands_uniq) == 1):
                    # only one band and bands give correct info
                    bands_in_obs = [bands_uniq[0].decode().strip()]
                elif (len(bands_uniq) > 1) & (len(fovid_uniq) == 1):
                    # only 1 fovid, use that to select band
                    bands_in_obs = fovid_to_band(fovid_uniq)
                elif (len(bands_uniq) > 1) & (len(fovid_uniq) > 1):
                    if len(is_primary_field_uniq) == 1:
                        bands_in_obs = [i.decode().strip() for i in bands_uniq]
                    else:
                        # strip PU
                        fovid_uniq = \
                            np.asarray([i for i in fovid_uniq
                                        if 'Peak-Up' not in i.decode()])
                        bands_in_obs = fovid_to_band(fovid_uniq)
                # get spitzer data
                url = rqk_t['accessUrl'][idx_first].strip()
                data = sha.get_file(url)
                # get time and calculate phase
                obs_time = data[0].header['HMJD_OBS'] + 2400000.5
                phase = ((obs_time-system["TT"])/system["PERIOD"]) % 1
                duration_in_phase = \
                    (aor_duration_in[iobs] /
                     (system["PERIOD"]*u.day).to(u.minute).value)
                # calculate transit and eclipse time in phase
                T1 = 0
                if np.ma.is_masked(system["ECC"]):
                    T2 = 0.5
                    query_error_orbit = True
                elif system["ECC"] < 0.01:
                    T2 = 0.5
                elif (np.ma.is_masked(system["INCL"]) |
                      np.ma.is_masked(system["OMEGA"])):
                    T2 = 0.5
                    query_error_orbit = True
                else:
                    T2 = (0.5*np.pi +
                          (1.0 + 1.0/np.sin(np.radians(system["INCL"]))**2) *
                          system["ECC"] *
                          np.cos(np.radians(system["OMEGA"]))) / \
                          np.pi
                print(phase, phase+duration_in_phase, T1, T2)
                # check if transit, eclispe or phase curve
                if (phase+duration_in_phase > 1):
                    isTransit = True
                else:
                    isTransit = False
                if ((phase+duration_in_phase > T2) &
                        (phase < T2)) | (phase+duration_in_phase > 1.0+T2):
                    isEclipse = True
                else:
                    isEclipse = False
                if isEclipse and isTransit:
                    isPhaseCurve = True
                else:
                    isPhaseCurve = False
                # store results
                isEclipse_in_obs = []
                isTransit_in_obs = []
                isPhaseCurve_in_obs = []
                phase_in_obs = []
                for iband, band in enumerate(bands_in_obs):
                    aor_keys_out.append(obsid)
                    aor_offsets_out.append(aor_offsets_in[iobs])
                    aor_offsets_flag_out.append(aor_offsets_flag_in[iobs])
                    aor_duration_out.append(aor_duration_in[iobs])
                    phase_in_obs.append(phase)
                    isEclipse_in_obs.append(isEclipse)
                    isTransit_in_obs.append(isTransit)
                    isPhaseCurve_in_obs.append(isPhaseCurve)
            except ValueError:
                query_error_lv1 = [True]
                aor_keys_out.append(obsid)
                aor_offsets_out.append(aor_offsets_in[iobs])
                aor_offsets_flag_out.append(aor_offsets_flag_in[iobs])
                aor_duration_out.append(aor_duration_in[iobs])
                bands_in_obs = [None]
                phase_in_obs = [None]
                isEclipse_in_obs = [False]
                isTransit_in_obs = [False]
                isPhaseCurve_in_obs = [False]
            aor_band += bands_in_obs
            aor_phase += phase_in_obs
            aor_error_lv1 += query_error_lv1
            aor_tansit_flag += isTransit_in_obs
            aor_eclipse_flag += isEclipse_in_obs
            aor_phase_curve_flag += isPhaseCurve_in_obs
        # convert to numpy array and store in dict
        aor_band = np.asarray(aor_band)
        aor_phase = np.asarray(aor_phase)
        aor_error_lv1 = np.asarray(aor_error_lv1)
        aor_tansit_flag = np.asarray(aor_tansit_flag)
        aor_eclipse_flag = np.asarray(aor_eclipse_flag)
        aor_phase_curve_flag = np.asarray(aor_phase_curve_flag)
        aor_keys_out = np.asarray(aor_keys_out)
        aor_offsets_flag_out = np.asarray(aor_offsets_flag_out)
        aor_duration_out = np.asarray(aor_duration_out)
        aor_offsets_out = np.asarray(aor_offsets_out)
        results_dict = {"TT": system["TT"], "PERIOD": system["PERIOD"],
                        "ECC": system["ECC"], "INCL": system["INCL"],
                        "OMEGA": system["OMEGA"],
                        "aor": aor_keys_out,
                        "band": aor_band,
                        "phase": aor_phase,
                        "tflag": aor_tansit_flag,
                        "eflag": aor_eclipse_flag,
                        "pflag": aor_phase_curve_flag,
                        "error": query_error,
                        "errorLV1": aor_error_lv1,
                        "errorOrbit": query_error_orbit,
                        "errorOffset": aor_offsets_flag_out,
                        "duration": aor_duration_out,
                        "searchOffset": aor_offsets_out}
    except (MissingSchema, KeyError):
        query_error = True
        results_dict = {"TT": system["TT"], "PERIOD": system["PERIOD"],
                        "ECC": system["ECC"], "INCL": system["INCL"],
                        "OMEGA": system["OMEGA"],
                        "aor": np.array([], dtype=np.int64),
                        "band": np.array([], dtype=np.byte),
                        "phase": np.array([], dtype=np.float64),
                        "tflag": np.array([], dtype=np.bool),
                        "eflag": np.array([], dtype=np.bool),
                        "pflag": np.array([], dtype=np.bool),
                        "error": query_error,
                        "errorLV1": np.array([], dtype=np.bool),
                        "errorOrbit": query_error_orbit,
                        "errorOffset": np.array([], dtype=np.bool),
                        "duration": np.array([], dtype=np.float64),
                        "searchOffset": np.array([], dtype=np.float64)}
    spitzer_dict[targets[itarget]] = results_dict

# save data
f = open('/home/bouwman/Spitzer_Search_Results.pickle', 'wb')
pickle.dump([targets, ParTable, spitzer_dict], f, -1)
f.close()

##################################
# create xls file to store resutls
###################################
alignWrap = xlwt.Alignment()
alignWrap.horz = xlwt.Alignment.HORZ_CENTER
alignWrap.vert = xlwt.Alignment.VERT_CENTER
alignWrap.wrap = xlwt.Alignment.WRAP_AT_RIGHT
style1 = xlwt.easyxf('font: bold 1, color black;')
style1.alignment = alignWrap
style1b = xlwt.easyxf('font: bold 1, color black; border: left thick')
style1b.alignment = alignWrap
style1c = xlwt.easyxf('font: bold 1, color black; border: right thick')
style1c.alignment = alignWrap
style2 = xlwt.easyxf('font: bold 1, color blue;')
style2.alignment = alignWrap
style2b = xlwt.easyxf('font: bold 1, color blue; border: left thick')
style2b.alignment = alignWrap
style2c = xlwt.easyxf('font: bold 1, color blue; border: right thick')
style2c.alignment = alignWrap
style2d = xlwt.easyxf('font: bold 1, color blue; border: right thick, left thick')
style2d.alignment = alignWrap
style3 = xlwt.easyxf('font: bold 1, color black; pattern: pattern solid, fore_colour red;')
style3.alignment = alignWrap
style3c = xlwt.easyxf('font: bold 1, color black; border: right thick; pattern: pattern solid, fore_colour red;')
style3c.alignment = alignWrap
style4 = xlwt.easyxf('font: bold 1, color red;')
style4.alignment = alignWrap
style4b = xlwt.easyxf('font: bold 1, color red; border: left thick')
style4b.alignment = alignWrap
style4c = xlwt.easyxf('font: bold 1, color red; border: right thick')
style4c.alignment = alignWrap


file_name_spitzer_query_results = "EXOPLANETSA_match_Spitzer_HST.xls"
# copy the xls file with the target list and etc parameters
wb_hst = xlrd.open_workbook(os.path.join(path_flle, file_name))
wb_spitzer = copy(wb_hst)
sheet_hst_targets = wb_spitzer.get_sheet(0)

# reformat target lists and create and format sheets with results
sheet_hst_targets.set_name('HST Archive Data')
sheet_hst_targets.col(0).width = 256 * 30
styleindex1 = wb_spitzer.add_style(style1)
styleindex1c = wb_spitzer.add_style(style1c)
styleindex2 = wb_spitzer.add_style(style2)
styleindex2c = wb_spitzer.add_style(style2c)
styleindex3c = wb_spitzer.add_style(style3c)
rows = sheet_hst_targets.get_rows()
for icol in [0, 1, 2, 4, 5]:
    rows[0]._Row__cells[icol].xf_idx = styleindex2c
for icol in [0, 3, 6]:
    rows[1]._Row__cells[icol].xf_idx = styleindex2c
for icol in [1, 2, 4, 5]:
    rows[1]._Row__cells[icol].xf_idx = styleindex2
for irow in list(rows.keys())[2:]:
    for icol in [0, 3, 6]:
        rows[irow]._Row__cells[icol].xf_idx = styleindex1c
    for icol in [1, 2, 4, 5]:
        rows[irow]._Row__cells[icol].xf_idx = styleindex1
for irow in list(rows.keys())[2:]:
    total = 0
    for icol in range(1, 7):
        total += rows[irow]._Row__cells[icol].number
    if total == 0:
        rows[irow]._Row__cells[0].xf_idx = styleindex3c

sheet_orbit_parameters = wb_spitzer.add_sheet("Orbit Parameters")
sheet_orbit_parameters.row(0).set_style(style2)
sheet_orbit_parameters.row(1).set_style(style2)
sheet_orbit_parameters.write_merge(0, 1, 0, 0, 'TARGET', style2)
sheet_orbit_parameters.col(0).width = 256 * 30
sheet_orbit_parameters.write_merge(0, 1, 1, 1, 'SIMBADNAME', style2)
sheet_orbit_parameters.col(1).width = 256 * 30
sheet_orbit_parameters.write(0, 2, 'RA', style2)
sheet_orbit_parameters.write(0, 3, 'DEC', style2)
sheet_orbit_parameters.col(2).width = 256 * 20
sheet_orbit_parameters.col(3).width = 256 * 20
sheet_orbit_parameters.write(0, 4, 'TT', style2)
sheet_orbit_parameters.col(4).width = 256 * 20
sheet_orbit_parameters.write(0, 5, 'PERIOD', style2)
sheet_orbit_parameters.write(0, 6, 'ECC', style2)
sheet_orbit_parameters.write(0, 7, 'INCL', style2)
sheet_orbit_parameters.write(0, 8, 'OMEGA', style2)
sheet_orbit_parameters.write(1, 2, '[hourangle]', style2)
sheet_orbit_parameters.write(1, 3, '[degrees]', style2)
sheet_orbit_parameters.write(1, 4, '[BJD]', style2)
sheet_orbit_parameters.write(1, 5, '[days]', style2)
sheet_orbit_parameters.write(1, 6, '', style2)
sheet_orbit_parameters.write(1, 7, '[degree]', style2)
sheet_orbit_parameters.write(1, 8, '[degree]', style2)
for ic, cn in enumerate(ParTable.keys()):
    for ir, rval in enumerate(ParTable[cn]):
        if isinstance(rval, bytes):
            sheet_orbit_parameters.write(ir+2, ic, str(rval.decode()), style1)
        else:
            sheet_orbit_parameters.write(ir+2, ic, str(rval), style1)

sheet_spitzer_targets = wb_spitzer.add_sheet("Spitzer Archive Data")
sheet_spitzer_targets.row(0).set_style(style2)
sheet_spitzer_targets.row(1).set_style(style2)

sheet_spitzer_targets.write_merge(0, 1, 0, 0, 'TARGET', style2)
sheet_spitzer_targets.col(0).width = 256 * 30
sheet_spitzer_targets.write_merge(0, 0, 1, 7, 'TRANSIT', style2b)
sheet_spitzer_targets.write_merge(0, 0, 8, 14, 'ECLIPSE', style2b)
sheet_spitzer_targets.write_merge(0, 0, 15, 21, 'PHASE CURVE', style2b)
sheet_spitzer_targets.write(1, 1, 'IRAC 3.6um', style2b)
sheet_spitzer_targets.write(1, 2, 'IRAC 4.5um', style2)
sheet_spitzer_targets.write(1, 3, 'IRAC 5.8um', style2)
sheet_spitzer_targets.write(1, 4, 'IRAC 8.0um', style2)
sheet_spitzer_targets.write(1, 5, 'IRS SL', style2)
sheet_spitzer_targets.write(1, 6, 'IRS PU Blue', style2)
sheet_spitzer_targets.write(1, 7, 'MIPS 24um', style2)
sheet_spitzer_targets.write(1, 8, 'IRAC 3.6um', style2b)
sheet_spitzer_targets.write(1, 9, 'IRAC 4.5um', style2)
sheet_spitzer_targets.write(1, 10, 'IRAC 5.8um', style2)
sheet_spitzer_targets.write(1, 11, 'IRAC 8.0um', style2)
sheet_spitzer_targets.write(1, 12, 'IRS SL', style2)
sheet_spitzer_targets.write(1, 13, 'IRS PU Blue', style2)
sheet_spitzer_targets.write(1, 14, 'MIPS 24um', style2)
sheet_spitzer_targets.write(1, 15, 'IRAC 3.6um', style2b)
sheet_spitzer_targets.write(1, 16, 'IRAC 4.5um', style2)
sheet_spitzer_targets.write(1, 17, 'IRAC 5.8um', style2)
sheet_spitzer_targets.write(1, 18, 'IRAC 8.0um', style2)
sheet_spitzer_targets.write(1, 19, 'IRS SL', style2)
sheet_spitzer_targets.write(1, 20, 'IRS PU Blue', style2)
sheet_spitzer_targets.write(1, 21, 'MIPS 24um', style2)

sheet_spitzer_targets.write_merge(0, 0, 22, 26, 'WARNINGS', style2d)
sheet_spitzer_targets.write(1, 22, 'ERROR', style2b)
sheet_spitzer_targets.write(1, 23, 'LVL1 ERROR', style2)
sheet_spitzer_targets.write(1, 24, 'ORBIT PAR ERROR', style2)
sheet_spitzer_targets.write(1, 25, 'CLASSIFICATION ERROR', style2)
sheet_spitzer_targets.write(1, 26, 'OFFSET ERROR', style2c)

sheet_spitzer_targets.col(22).width = 256 * 14
sheet_spitzer_targets.col(23).width = 256 * 14
sheet_spitzer_targets.col(24).width = 256 * 14
sheet_spitzer_targets.col(25).width = 256 * 18
sheet_spitzer_targets.col(26).width = 256 * 16
sheet_spitzer_targets.row(1).height = 255*2

band_search_list = ['IRAC 3', 'IRAC 4', 'IRAC 5', 'IRAC 8',  'IRS SL',
                    'IRS PU Blue', 'MIPS 24']

for itarget, target in enumerate(spitzer_dict.keys()):
    idx_good = ~(spitzer_dict[target]['band'] == None)
    aor_array = spitzer_dict[target]['aor'][idx_good]
    tflag_array = spitzer_dict[target]['tflag'][idx_good]
    eflag_array = spitzer_dict[target]['eflag'][idx_good]
    pflag_array = spitzer_dict[target]['pflag'][idx_good]
    total_aors = len(tflag_array)
    total_classified_aors = np.sum(tflag_array+eflag_array+pflag_array)
    error_not_classified_aors = total_aors - total_classified_aors
    error_flag = np.sum(spitzer_dict[target]['error'])
    error_flag_lvl1 = np.sum(spitzer_dict[target]['errorLV1'])
    error_flag_orbit = np.sum(spitzer_dict[target]['errorOrbit'])
    error_flag_offsets = np.sum(spitzer_dict[target]['errorOffset'])
    band_list = [band.strip() for band in spitzer_dict[target]['band'][idx_good].tolist()]
    for iss, str_search in enumerate(band_search_list):
        idx_match = [i for i, s in enumerate(band_list) if str_search in s]
        if len(idx_match) > 0:
            number_of_transit_obs = np.sum(tflag_array[idx_match])
            number_of_eclipse_obs = np.sum(eflag_array[idx_match])
            number_of_phase_curve_obs = np.sum(pflag_array[idx_match])
        else:
            number_of_transit_obs = 0
            number_of_eclipse_obs = 0
            number_of_phase_curve_obs = 0
        if iss == 0:
            style_use = style1b
        else:
            style_use = style1
        sheet_spitzer_targets.write(2+itarget, 1+iss,
                                    str(number_of_transit_obs), style_use)
        sheet_spitzer_targets.write(2+itarget, 1+iss+7,
                                    str(number_of_eclipse_obs), style_use)
        sheet_spitzer_targets.write(2+itarget, 1+iss+14,
                                    str(number_of_phase_curve_obs), style_use)
    if error_flag == 0:
        style_use = style1b
    else:
        style_use = style4b
    sheet_spitzer_targets.write(2+itarget, 22, str(error_flag), style_use)
    if error_flag_lvl1 == 0:
        style_use = style1
    else:
        style_use = style4
    sheet_spitzer_targets.write(2+itarget, 23, str(error_flag_lvl1), style_use)
    if error_flag_orbit == 0:
        style_use = style1
    else:
        style_use = style4
    sheet_spitzer_targets.write(2+itarget, 24, str(error_flag_orbit), style_use)
    if error_not_classified_aors == 0:
        style_use = style1
    else:
        style_use = style4
    sheet_spitzer_targets.write(2+itarget, 25, str(error_not_classified_aors),
                                style_use)
    if error_flag_offsets == 0:
        style_use = style1c
    else:
        style_use = style4c
    sheet_spitzer_targets.write(2+itarget, 26, str(error_flag_offsets),
                                style_use)
    if total_classified_aors > 0:
        sheet_spitzer_targets.write(2+itarget, 0, target, style1)
    else:
        sheet_spitzer_targets.write(2+itarget, 0, target, style3)

wb_spitzer.save(os.path.join(path_flle, file_name_spitzer_query_results))
