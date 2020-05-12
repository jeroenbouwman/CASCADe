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
function cleantmp {
    nrun=`ps -def | grep exo | grep run_cascade | grep python | wc -l`
    if [ $nrun -eq 0 ]; then
        rm -f /tmp/python*
        rm -rf /tmp/ray*
    fi
}


helpFunction()
{
   echo ""
   echo "Usage: $0 [-opts] [args]"
   echo ""
   echo "-a            Process all observations in (for now) \$CASCADE_DATA_PATH/data/HST/WFC3"
   echo "                  already processed observations are skipped unless -r is selected"
   echo "-c ncpu       Use this maximum number of cpu cores"
   echo "                  default set via ini file and code"
   echo "-f            Process failed observations"
   echo "                  observations with Error in the logfile"
   echo "-p names      Process listed planet observation(s)"
   echo "                  -p \"WASP-43b_iccz30 WASP-43b_iccz31\""
   echo "                  using -p implies -r"
   echo "-r            Rerun observations"
   echo "                  default only not yet reduced observations are processed"
   echo "-t            Transit calibration (explanet spectra extraction) only"
   echo ""
   echo "The data is processed by calling \$CASCADE_PATH/scripts/run_cascade.py"
   echo "for all requested planet observations."
   echo ""
   echo "Logfiles are saved in \$CASCADE_LOG_PATH/run_cascade/"
   echo "At each run the old logfiles are moved to \$CASCADE_LOG_PATH/run_cascade_old/"
   echo ""
   echo "A lockfile is put in \$CASCADE_DATA_PATH/locked/"
   echo "when present starting a run on the same object will skip the processing."
   echo ""
   exit 1 # Exit script after printing help
}
[ $# -eq 0 ] && helpFunction 
all=0
run_failed=0
rerun=0
order_by_size=0
transit_only=0
run_failed=0
ncpu=0
while getopts "afrstc:p:" opt
do
   case "$opt" in
      a ) all=1 ;;
      f ) run_failed=1 ;;
      r ) rerun=1 ;;
      s ) order_by_size=1 ;;
      t ) transit_only=1 ;;
      c ) ncpu="$OPTARG" ;;
      p ) planets="$OPTARG" ;;
      ? ) helpFunction ;; 
   esac
done

start=`date +%s`

[ -z ${CASCADE_HOME} ] && export CASCADE_HOME=${HOME}/CASCADe/
[ -z ${CASCADE_PATH} ] && export CASCADE_PATH=${CASCADE_HOME}/development/
[ -z ${CASCADE_LOG_PATH} ] && export CASCADE_LOG_PATH=${CASCADE_HOME}/logs/
[ -z ${CASCADE_DATA_PATH} ] && export CASCADE_DATA_PATH=${CASCADE_HOME}/observations/
[ -z ${CASCADE_SAVE_PATH} ] && export CASCADE_SAVE_PATH=${CASCADE_HOME}/observations/results/
[ -z ${CASCADE_INITIALIZATION_FILE_PATH} ] && export CASCADE_INITIALIZATION_FILE_PATH=${CASCADE_HOME}/init_files/

mkdir -p $CASCADE_LOG_PATH/run_cascade/
mkdir -p $CASCADE_LOG_PATH/run_cascade_old/

lockdir=$CASCADE_DATA_PATH/locked/
mkdir -p $lockdir
d=`date "+%Y-%m-%d"`

if [ $run_failed -eq 1 ]; then
    cwd=$PWD
    cd $CASCADE_LOG_PATH/run_cascade/
    if [ $transit_only -eq 1 ]; then
        planets=`grep Error transit_* | awk -F_ '{print $2"_"$3}' | uniq`
    else
        planets=`grep Error extract_* | awk -F_ '{print $2"_"$3}' | uniq`
    fi
    cd $cwd
    [ $transit_only -ne 1 ] && rerun=1
fi

# set order_size to run smallest observations first 
# nice to quicker get more results, good  when testing
if [ -n "$planets" ]; then
    [ $transit_only -ne 1 ] && rerun=1
else
    if [ $order_by_size -eq 1 ]; then
        pwd=$PWD
        cd $CASCADE_DATA_PATH/data/HST/WFC3/
        planets=`du -sm * | sort -n | awk '{print $2}'`
        cd $pwd
    else
        planets=`ls $CASCADE_DATA_PATH/data/HST/WFC3`
    fi
fi

echo ""
echo "Planet obervations being reduced"
echo $planets
echo ""

for planet in $planets; do
    echo $planet
    date
    lockfile=$lockdir/$planet
    if [ -f $lockfile -a $transit_only -eq 0 ]; then
        process=0
        echo "$planet already (being) processed"
        echo ""
    else 
        if [ -f $CASCADE_DATA_PATH/failed/$planet -a $rerun -eq 0 -a $transit_only -eq 0 ]; then
            echo "not processing failed observation $planet"
            echo ""
            process=0
        else
            process=1
            [ $transit_only -eq 0 ] && echo 1 > $lockfile
            if [ $rerun -eq 1 ]; then
                if [ $transit_only -eq 1 ]; then
                    echo "removing planet transit result for $planet"
                else
                    rm -rf $CASCADE_DATA_PATH/data/HST/WFC3/${planet}/SPECTRA
                    echo "removing SPECTRA and planet transit result for $planet"
                fi
                rm -rf $CASCADE_SAVE_PATH/$planet/
            fi
            iof=cascade_${planet}_object.ini
            ief=cascade_${planet}_extract_timeseries.ini
            itf=cascade_${planet}_calibrate_planet_spectrum.ini
            ip=$CASCADE_INITIALIZATION_FILE_PATH/HST/WFC3/${planet}
            iter=1
        fi
    fi
    processed=0

    while [ $process -eq 1 ]
    do
        process=0
        if [ -d "$CASCADE_DATA_PATH/data/HST/WFC3/${planet}/SPECTRA" ]; then
           n=`ls $CASCADE_DATA_PATH/data/HST/WFC3/${planet}/SPECTRA/*.fits | wc -l`
           if [ $n -eq 0 ]; then extract=1; else extract=0; fi
        else
           extract=1
        fi
        [ $transit_only -eq 1 ] && extract=0
        if [ $extract -eq 1 ]; then
            if [ $iter -eq 1 ]; then
                rm -rf $CASCADE_SAVE_PATH/$planet/
                mv $CASCADE_LOG_PATH/run_cascade/extract_${planet}_*.log $CASCADE_LOG_PATH/run_cascade_old/ 2>/dev/null
            fi
# set max cpu in init file
            if [ $ncpu -gt 0 ]; then
                cpu_str=`grep cascade_max_number_of_cpus $ip/$ief | awk -F: '{print $NF}'`
                sed s/"$cpu_str"/"cascade_max_number_of_cpus = $ncpu"/g $ip/$ief > tmp_$ief
                mv tmp_$ief $ip/$ief
            fi
            logfile=$CASCADE_LOG_PATH/run_cascade/extract_${planet}_$d.log
            echo $logfile
            date >> $logfile
            echo "$CASCADE_PATH/scripts/run_cascade.py $iof $ief -ip $ip >> $logfile"
            SECONDS=0
            cleantmp
            $CASCADE_PATH/scripts/run_cascade.py $iof $ief -ip $ip -nw &>> $logfile
            cleantmp
# reset max cpu in init file
            if [ $ncpu -gt 0 ]; then
                str=`grep cascade_max_number_of_cpus $ip/$ief | awk -F: '{print $NF}'`
                sed s/"$str"/"$cpu_str"/g $ip/$ief > tmp_$ief
                mv tmp_$ief $ip/$ief
            fi
            echo "" >> logfile
            echo "PROCESSING TIME $(echo "scale=1; $SECONDS / 60" | bc)  minutes" >> $logfile
            echo "TOTAL RUNGTIME $((($(date +%s)-$start)/3600)) hours" >> $logfile
            echo "" >> $logfile
            echo ""
            echo "PROCESSING TIME $(echo "scale=1; $SECONDS / 60" | bc)  minutes" 
            echo "TOTAL RUNGTIME $((($(date +%s)-$start)/3600)) hours" 
            echo ""
        fi
#
#= extracted SPECTRA present
#
        if [ -f "$CASCADE_SAVE_PATH/$planet/$planet_exoplanet_spectra.fits" ]; then
            error=0
        else
            logfile=$CASCADE_LOG_PATH/run_cascade/transit_${planet}_$d.log
            echo $logfile
            if [ $iter -eq 1 ]; then
                mv $CASCADE_LOG_PATH/run_cascade/transit_${planet}_*.log $CASCADE_LOG_PATH/run_cascade_old/ 2>/dev/null
                rm -f $CASCADE_SAVE_PATH/$planet/*exoplanet_spectra*
            else
                rm $logfile
            fi
            if [ -d "$CASCADE_DATA_PATH/data/HST/WFC3/${planet}/SPECTRA" ]; then
                date >> $logfile
                echo "$CASCADE_PATH/scripts/run_cascade.py $iof $itf -ip $ip -nw >> $logfile"
                SECONDS=0
# test a very small cpm_relative_sig_value_limit (was 0.05)
                #cpm_str=`grep cpm_relative_sig_value_limit $ip/$itf | awk -F: '{print $NF}'`
                #sed s/"$cpm_str"/"cpm_relative_sig_value_limit = 0.0005"/g $ip/$itf > tmp_$itf
                #mv tmp_$itf $ip/$itf
                cleantmp
                $CASCADE_PATH/scripts/run_cascade.py $iof $itf -ip $ip &>> $logfile
                cleantmp
#
# check if there is a broadcast errror, if so increase the rebin factor and repeat extraction and transit fit
#
                error=`grep Error $logfile | grep "operands could not be broadcast together with shapes" | wc -l`
                if [ $error -ne 0 ]; then
                    process=1 # rerun with modified rebin factor
                    errstr=`grep Error $logfile | grep "operands could not be broadcast together with shapes"`
                    echo errstr
                    echo errstr >> $logfile
                    rm -rf $CASCADE_DATA_PATH/data/HST/WFC3/${planet}/SPECTRA
                    rebin_str=`grep processing_rebin_factor_extract1d $ip/$ief`
                    rebin="$(cut -d'=' -f2 <<< $rebin_str)"
                    if [ $iter -eq 1 ]; then
                        err_str=`grep Error $logfile | grep "operands could not be broadcast together with shapes"`
                        s_one="$(cut -d'(' -f2 <<< $err_str)"
                        s_two="$(cut -d'(' -f3 <<< $err_str)"
                        s_one="$(cut -d',' -f1 <<< $s_one)"
                        s_two="$(cut -d',' -f1 <<< $s_two)"
                        rebin=`echo "$rebin*$s_one/$s_two" | bc -l`
                        iter=2
                    fi
                    rebin=`echo "1.05*$rebin" | bc -l`
                    echo "$planet rebin "$rebin >> $logfile
                    echo "$planet rebin "$rebin
                    sed s/"$rebin_str"/"processing_rebin_factor_extract1d = $rebin"/g $ip/$ief > tmp_$ief
                    mv tmp_$ief $ip/$ief
                    sed s/"processing_auto_adjust_rebin_factor_extract1d = True"/"processing_auto_adjust_rebin_factor_extract1d = False"/g $ip/$ief > tmp_$ief
                    mv tmp_$ief $ip/$ief
                else
                    process=0
                fi

                echo "" >> $logfile
                echo "PROCESSING TIME $(echo "scale=1; $SECONDS / 60" | bc)  minutes" >> $logfile
                echo "TOTAL RUNGTIME $((($(date +%s)-$start)/3600)) hours" >> $logfile
                echo "" >> $logfile
                echo " "
                echo "PROCESSING TIME $(echo "scale=1; $SECONDS / 60" | bc)  minutes" 
                echo "TOTAL RUNGTIME $((($(date +%s)-$start)/3600)) hours" 
                echo " "
            else
                echo "Error: No extracted spectra for $planet" > $logfile
            fi
        fi
        processed=1
    done
    [ $processed -eq 1 -a $transit_only -eq 0 ] && rm -f $lockfile
done

exit 0
