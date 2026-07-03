import os
import glob
from tqdm import tqdm

BASEDIR = "." #do not end with /

#def rename(old_name,new_name):
#os.rename(old_name,new_name)

def get_old_names(pattern):
    """
    :return:
    """
    return sorted(glob.glob(f"{BASEDIR}/{pattern}"))

def generate_new_name(old_name):
    """
    this version takes res_XXX_YYY_ZZZ_AA.fits to resXXXAA.fits
    or rres_XXX_YYY_ZZZ_AA.fits to rresXXXAA.fits
    :param old_name:
    :return:
    """
    dirname = os.path.dirname(old_name)
    basename = os.path.basename(old_name)

    toks = basename.split("_")
    return os.path.join(dirname,toks[0]+toks[1]+toks[4])



old_names = get_old_names("res_*.fits")
for old_name in tqdm(old_names):
    new_name = generate_new_name(old_name)
    os.rename(old_name,new_name)

old_names = get_old_names("rres_*.fits")
for old_name in tqdm(old_names):
    new_name = generate_new_name(old_name)
    os.rename(old_name,new_name)