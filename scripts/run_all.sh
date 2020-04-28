#!/bin/bash

#
# script to run run_cascade for all observations (if defined in the init directory)
# logs are created for each step. Using env var. CASCADE_LOG for this.
#
# can be run offline, e.g. on a server machine
#    screen -L ./run_all 
# then with Ctrl-A D you detatch the screen and log off.
# -L creates an overall log screenlog.0
#
#
# TBD remove WFC3 move after update from Jeroen
#

d=`date "+%Y-%m-%d"`
echo $d

mkdir -p $CASCADE_LOGS/run_cascade/

start=`date +%s`
planets=`ls $CASCADE_INITIALIZATION_FILE_PATH/HST/WFC3/`
for planet in $planets; do
  echo $planet

  start_planet=`date +%s`

  iof=cascade_${planet}_object.ini
  ief=cascade_${planet}_extract_timeseries.ini
  itf=cascade_${planet}_calibrate_planet_spectrum.ini
  ip=$CASCADE_INITIALIZATION_FILE_PATH/HST/WFC3/${planet}

  echo "logfile=$CASCADE_LOGS/run_cascade/extract_${planet}_$d.log"
  logfile=$CASCADE_LOGS/run_cascade/extract_${planet}_$d.log
  rm -f $logfile
  echo "$CASCADE_PATH/scripts/run_cascade.py $iof $ief -ip $ip &> $logfile"
  SECONDS=0
  $CASCADE_PATH/scripts/run_cascade.py $iof $ief -ip $ip &> $logfile
  echo ""
  echo "PROCESSING TIME $(echo "scale=1; $SECONDS / 60" | bc)  minutes" >> $logfile
  echo "PROCESSING TIME $(echo "scale=1; $SECONDS / 60" | bc)  minutes" 
  echo "TOTAL RUNGTIME $((($(date +%s)-$start)/3600)) hours" >> $logfile
  echo "TOTAL RUNGTIME $((($(date +%s)-$start)/3600)) hours" 

# SPECTRA are placed in the wrong place, so moving it here untill this is fixed in run_cascade
  if [ -d "$CASCADE_DATA_PATH/data/HST/WFC3/${planet}/SPECTRA" ]; then
    echo "Move extracted spectra"
    rm -r $CASCADE_DATA_PATH/data/HST/WFC3/${planet}/SPECTRA
    mv $CASCADE_DATA_PATH/data/WFC3/${planet}/SPECTRA $CASCADE_DATA_PATH/data/HST/WFC3/${planet}/
    rmdir $CASCADE_DATA_PATH/data/WFC3/${planet}
    rmdir $CASCADE_DATA_PATH/data/WFC3
  fi

  echo "logfile=$CASCADE_LOGS/run_cascade/transit_${planet}_$d.log"
  logfile=$CASCADE_LOGS/run_cascade/transit_${planet}_$d.log
  rm -f $logfile
  echo "$CASCADE_PATH/scripts/run_cascade.py $iof $itf -ip $ip &> $logfile"
  SECONDS=0
  $CASCADE_PATH/scripts/run_cascade.py $iof $itf -ip $ip &> $logfile
  ends=`date +%s`
  echo ""
  echo "PROCESSING TIME $(echo "scale=1; $SECONDS / 60" | bc)  minutes" >> $logfile
  echo "PROCESSING TIME $(echo "scale=1; $SECONDS / 60" | bc)  minutes" 

  echo "TOTAL RUNGTIME $((($(date +%s)-$start)/3600)) hours" >> $logfile
  echo "TOTAL RUNGTIME $((($(date +%s)-$start)/3600)) hours" 
  echo " "

done

