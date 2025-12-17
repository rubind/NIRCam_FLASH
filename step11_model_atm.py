import glob
from FileRead import readcol, file_to_fn, writecol
import Spectra
import numpy as np
import tqdm

[waves] = readcol("bosz2024_wave_r500.txt", 'f')
log_step = np.median(np.log(waves[1:]/waves[:-1]))

print("log_step", log_step)
assert np.all(np.abs(np.log(waves[1:]/waves[:-1]) - log_step) < 1e-3*log_step)

all_fls = np.sort(glob.glob("/home/drubin/koa_scratch/r500/m-0.50/bosz2024*txt")[::100])
all_Lsol_one_Rsol = []

for fl in tqdm.tqdm(all_fls):
    [H, continuum] = readcol(fl, 'ff')

    total_L_per_cm2 = sum(4 * np.pi * H * log_step*waves) * 1e-7 # ergs/s to Watts
    total_L = total_L_per_cm2 * 4*np.pi * (6.957e+10)**2.

    all_Lsol_one_Rsol.append(total_L/3.828e26)

writecol("model_atmosphere_grid.txt", [all_fls, all_Lsol_one_Rsol], headings = ["file", "all_Lsol_one_Rsol"])
