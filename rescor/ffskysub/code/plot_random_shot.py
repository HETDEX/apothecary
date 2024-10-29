import os
import numpy as np

with open("all_shots.txt", "r") as ff:
    all_shots = ff.readlines()

all_shots = [x.strip() for x in all_shots]
all_shots = [int(x.split("/")[-1][:-3]) for x in all_shots]

this_shot = np.random.choice(all_shots)

exe_string = f"python plot_shot.py {this_shot}"
print(exe_string)
os.system(exe_string)