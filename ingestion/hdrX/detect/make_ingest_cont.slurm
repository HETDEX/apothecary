#runfn is the ingest_cont_YYYYMM.run filename to use
runfn=$1
/work/03261/polonius/hetdex/science_reductions/myscripts/jobsplitter $runfn $(wc $runfn | awk '{print $1+1}') 1 "00:15:00" 1 128
