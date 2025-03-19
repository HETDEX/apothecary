Basic instructions:

>  rgetran <date> <shot>   or regtran_spr  makes the rmksim0c or rmksim0d call for skx and spr respectively
can just pump the output as:
> $(rgettran <data> <shot>)   s|t it makes the string and calls it as the command

IF the last line of rmksim0c/d is the rstep1 , then it will be immediately queued

Normally, just use launch_sims, which takes the queue and the number to launch as inputs:
e.g.: > launch_sims skx 40
will queue up 40 data+shots on the skx queue

on stampede3 , 40 on skx is the maximum ... anything extra goes to WaitQoS

-- currently spr is not working the way I want, so avoid it -- Note: spr is also 2x as expensive in SUs

To check the output:
make sure the HETDEX.oXXX file is free from errors and no timeout
make sure the ./output/*.find files are not all zeros
make sure the .sim file is populated

other??
