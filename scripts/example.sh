#!/bin/bash

conda activate cascade
#  to be customised
#export CASCADE_ROOT=/home/bouwman/CASCADeTest
export CASCADE_ROOT=/Users/gastaud/tmp/jeroen
export CASCADE_SCRIPTS=/Users/gastaud/test_dap/cascade/CASCADe/examples/init_files/
#
export CASCADE_INITIALIZATION_FILE_PATH=$CASCADE_ROOT/init_files/
export CASCADE_DATA_PATH=$CASCADE_ROOT
export CASCADE_SAVE_PATH=$CASCADE_ROOT/results/HST/WFC3/

echo $CASCADE_INITIALIZATION_FILE_PATH
echo $CASCADE_DATA_PATH
echo $CASCADE_SAVE_PATH

python $CASCADE_SCRIPTS/build_local_hst_archive.py -v 12181.14 -cio -nw

python $CASCADE_SCRIPTS/run_cascade.py cascade_WASP-19b_ibh715_object.ini cascade_WASP-19b_ibh715_extract_timeseries.ini -ip HST/WFC3/WASP-19b_ibh715/ -nw

python $CASCADE_SCRIPTS/run_cascade.py cascade_WASP-19b_ibh715_object.ini cascade_WASP-19b_ibh715_calibrate_planet_spectrum.ini -ip HST/WFC3/WASP-19b_ibh715/ -nw

unset CASCADE_INITIALIZATION_FILE_PATH
unset CASCADE_DATA_PATH
unset CASCADE_SAVE_PATH
unset CASCADE_ROOT
unset CASCADE_SCRIPTS

conda deactivate cascade
