#!/usr/bin/env python

#commandline has one required parameter, the filename of the table to update

from astropy.table import Table
import os.path as op
import sys



cl_args = list(sys.argv)
#has to be there
fn_in = cl_args[1]
fn_out = cl_args[2]

T = Table.read(fn_in,format="ascii")

T.write(fn_out,format="fits",overwrite=False)



