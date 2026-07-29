#!/bin/bash

#copy a month of tars from corral, untar locally, and remove the local tar file
#usage copy_from_corral YYYYMM

CURRENT_USER=$(id -u -n)

if [ $CURRENT_USER == "hetdex" ]; then
  echo "Not okay. DO not run as hetdex"
  exit -1
else
  echo "Running as: $CURRENT_USER in 5 seconds ..."
  sleep 5
fi


#first copy the tars
for file in /corral/utexas/Hobby-Eberly-Telesco/het_raw/$1*.tar; do
  if [ -f "stop" ]; then
    echo "stop found. exiting ..."
    exit 1
  fi
  echo "copying $file"
  cp $file .
done

#next untar and remove the local copy 
for file in ./$1*.tar; do
  if [ -f "stop" ]; then
    echo "stop found. exiting ..."
    exit 1
  fi
  echo "untarring $file"
  tar -xvf $file
  mv $file old_tar
done

echo "Done"
