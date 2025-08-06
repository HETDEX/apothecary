#runfn is the rcal_all_YYYYMM.run filename to use
runfn=$1
/work/03261/polonius/hetdex/science_reductions/myscripts/jobsplitter $runfn $(wc $runfn | awk '{print $1+1}') 1 "02:00:00" 1 64
