#runfn is the ingest_line_YYYYMM.run filename to use, assumes 1 month of data
runfn=$1
/work/03261/polonius/hetdex/science_reductions/myscripts/jobsplitter $runfn $(wc $runfn | awk '{print $1+1}') 1 "00:45:00" 1 128
