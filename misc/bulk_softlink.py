import glob
import os



linkpath = "../../original/"
dirs = sorted(glob.glob("sci*"))

for dir in dirs:
    fn = dir.replace("sci","")

    if not os.path.exists(f"{dir}/{fn}.h5"):
        cmd = f"ln -s {linkpath}{dir}/{fn}.h5 {dir}/{fn}.h5"
    os.system(f"{cmd}")

