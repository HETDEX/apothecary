#!/usr/bin/env python

from h5tools import amp_stats as AS
#from hetdex_api.config import HDRconfig
import os.path as op
import numpy as np
from tqdm import tqdm

shotids = np.loadtxt("hdr4.shots",dtype=int)

sel = np.array(shotids > 20230101000)
sel = sel & np.array(shotids < 20230901000)

for shotid in tqdm(shotids[sel]):
	try:

		if op.exists(f"{shotid}_stats.pickle"):
			continue
		_ = AS.make_stats_for_shot(shotid=shotid)
	except Exception as e:
		print(f"{shotid}: exception {e}")
		

