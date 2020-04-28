#!/bin/bash

conda activate cascade

export CASCADE_INITIALIZATION_FILE_PATH=/home/bouwman/CASCADeTest/init_files/
export CASCADE_DATA_PATH=/home/bouwman/CASCADeTest/
export CASCADE_SAVE_PATH=/home/bouwman/CASCADeTest/results/HST/WFC3/

python build_local_hst_archive.py -v 12181.14 -cio -nw

python run_cascade.py cascade_WASP-19b_ibh715_object.ini cascade_WASP-19b_ibh715_extract_timeseries.ini -ip HST/WFC3/WASP-19b_ibh715/ -nw

python run_cascade.py cascade_WASP-19b_ibh715_object.ini cascade_WASP-19b_ibh715_calibrate_planet_spectrum.ini -ip HST/WFC3/WASP-19b_ibh715/ -nw

unset CASCADE_INITIALIZATION_FILE_PATH
unset CASCADE_DATA_PATH
unset CASCADE_SAVE_PATH

conda deactivate cascade
