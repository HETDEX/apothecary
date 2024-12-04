#spin up 10 (or change number) of calls to run the res, sres, rres calculations
#each one will see what is available to run (semaphore blocked), take the next and run
#
# to stop all of them, create a file called "stop"
#   at various points in their iterations, they check for the file and if it exits, they will exit
#
# files named *.open are available to be serviced
# files named *.run are currently being serviced by one of these instances OR were prematurely terminated (timeout or "stop" file)
# files named *.done have been completed
#
# any *.run files that need to be resumed, need to be renamed *.open (can use python reset_run.py)
# each instance checks for the next available *.open file and checks to see how much might already be complete (e.g. if it was a .run)
#   if some already done, it re-runs the most recent month with data (assumes it to be the one that was interrupted) and forward
for i in $(seq 1 10);
do
  python res_sres_rres.py &
done
