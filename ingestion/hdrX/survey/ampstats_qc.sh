#!/usr/bin/env python

#commandline has one required parameter, the filename of the table to update

from astropy.table import Table
from h5tools import amp_stats as AS
from hetdex_api.config import HDRconfig
import os.path as op
import sys



cl_args = list(map(str.lower,sys.argv))
#has to be there
fn = cl_args[1]


T1 = Table.read(fn)
T1 = AS.stats_qc(T1,extend=True)

T1.write(fn,format="fits",overwrite=True)


