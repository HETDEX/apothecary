#!/usr/bin/env python

#commandline has one required parameter, the filename of the table to update

from astropy.table import Table
from h5tools import amp_stats as AS
from hetdex_api.config import HDRconfig
import os.path as op
import sys



cl_args = list(sys.argv)
#has to be there
fn = cl_args[1]


AS.stats_make_simple_amp_flag_file(fn,outname="amp_flag",overwrite=False)



