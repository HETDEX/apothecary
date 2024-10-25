#!/bin/bash

#. ../setup.sh

if [[ $# < 6 ]]; then
  echo "Too few parameters. Don't forget to add in the FP date and, optionally, GAIA vs SDSS"
  exit -1
fi

NIGHT=$1
SHOT=$2
RA=$3
DEC=$4
TRACK=$5
FP=$6

CFG=../${NIGHT}v${SHOT}.config
#LOG=${NIGHT}v${SHOT}.log
LOG=/dev/null

# Fall back to vdrp.config if no night/shot specific 
# configuration file exists.
#if [ ! -f $CFG ]; then
#    CFG=../vdrp.config
#fi
if [[ $# == 7 ]]; then
  if [[ $7 == 'GAIA' ]];then
    CFG=../vdrp.config.gaia
  elif [[ $7 == 'SDSS' ]]; then
    CFG=../vdrp.config.sdss
  else
    CFG=../vdrp.config.original
  fi
else
  if [ ! -f $CFG ]; then
    CFG=../vdrp.config.original
  fi
fi


echo Configuration file $CFG
vdrp_astrom --logfile $LOG -c $CFG $NIGHT $SHOT $RA $DEC $TRACK -t all --fplane_txt /scratch/03261/polonius/science_reductions/vdrp/fplane/fp$FP
rclean $NIGHT $SHOT
