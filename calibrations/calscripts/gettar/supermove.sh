#!/bin/sh

TARLIST=$1
TYPE=$2

for i in `grep $TYPE $TARLIST | cut -b 1-27 | sort -g | uniq | awk '{print $1".tar"}'`
do
    NUMBER=`echo $i | cut -b 21-27`
    echo $NUMBER
    
    if [[ $NUMBER -gt 999 ]] && [[ $NUMBER -lt 7000 ]]; then
	echo $i
    fi
    #echo rm -f $i 
done


    
