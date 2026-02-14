import glob
from FileRead import readcol, file_to_fn, writecol
import Spectra
import numpy as np
import tqdm
import matplotlib.pyplot as plt

def abs_AB_mag_one_Rsol(waves, dwaves, H, passband):
    one_Rsol_flux_10pc = 4 * np.pi * H
    one_Rsol_flux_10pc *= (6.957e+10/(10*3.086e+18))**2.

    return -2.5*np.log10(
        (passband * one_Rsol_flux_10pc * waves * dwaves).sum() /
        (passband * waves * dwaves * (0.10884806248/waves**2.)).sum()
        )

[waves] = readcol("bosz2024_wave_r500.txt", 'f')
log_step = np.median(np.log(waves[1:]/waves[:-1]))

print("log_step", log_step)
assert np.all(np.abs(np.log(waves[1:]/waves[:-1]) - log_step) < 1e-3*log_step)

dwaves = waves*log_step

all_fls = np.sort(glob.glob("/home/drubin/koa_scratch/r500/m-0.50/bosz2024*txt"))
all_Lsol_one_Rsol = []


filt_interp = {}

for filt in ["F090W_May2024_mean_system_throughput.txt", "F200W_May2024_mean_system_throughput.txt", "F335M_May2024_mean_system_throughput.txt", "F444W_May2024_mean_system_throughput.txt"]:
    filt_interp[filt.split("_")[0]] = file_to_fn("filters/" + filt, kind = 'linear', fill_value = 0., bounds_error = False)(waves/10000.)

filt_interp["GaiaG"] = file_to_fn("filters/passband.dat", kind = 'linear', fill_value = 0., bounds_error = False)(waves/10.)


for filt in glob.glob("filters/HST*dat"):
    filt_interp[filt.split(".")[-2]] = file_to_fn(filt, kind = 'linear', fill_value = 0., bounds_error = False)(waves)


abs_mags_one_Rsol = {}
for filt in filt_interp:
    abs_mags_one_Rsol[filt] = []

for filt in filt_interp:
    print(filt, np.sum(filt_interp[filt]*waves*dwaves)/np.sum(filt_interp[filt]*dwaves))

for fl in tqdm.tqdm(all_fls):
    [H, continuum] = readcol(fl, 'ff')

    total_L_per_cm2 = sum(4 * np.pi * H * dwaves) * 1e-7 # ergs/s to Watts
    total_L = total_L_per_cm2 * 4*np.pi * (6.957e+10)**2.

    all_Lsol_one_Rsol.append(total_L/3.828e26)

    for filt in filt_interp:
        abs_mags_one_Rsol[filt].append(abs_AB_mag_one_Rsol(waves = waves, dwaves = dwaves, H = H, passband = filt_interp[filt]))

    
to_write = [all_fls, all_Lsol_one_Rsol]
headings = ["file", "all_Lsol_one_Rsol"]

for filt in filt_interp:
    to_write.append(abs_mags_one_Rsol[filt])
    headings.append(filt + "_one_Rsol")
    
writecol("model_atmosphere_grid.txt", to_write, headings = headings)
