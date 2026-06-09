import pandas as pd
import numpy as np
from DavidsNM import miniNM_new
import matplotlib.pyplot as plt
from scipy.stats import scoreatpercentile
import subprocess
import tqdm
import sys

df = pd.read_csv("my_with_hst_fit.csv")

filt_name = sys.argv[1]
target = sys.argv[2]
blue_filt = "F150W"

print(df)
print(df.columns)


for key in ["chi2_SED_fit", "r_rsol", "mod_f475w", "F150W_count", "mod_f814w"]:
    df[key] = pd.to_numeric(df[key], errors="coerce")

inds = np.where((df["chi2_SED_fit"] < 200)*(df["F150W_count"] > 1)*(df["mod_f814w"] > 0))




log_R = np.log10(np.abs(np.array(df["r_rsol"])[inds]))
r_band = 0.5*(np.array(df["mod_f475w"])[inds] + np.array(df["mod_f814w"])[inds])
HSC_depth = 22.5

HSC_flux_ZP_equal_depth = 10.**(-0.4*(r_band - HSC_depth))

log_frac_unc = np.log10(0.2/HSC_flux_ZP_equal_depth)

bin_edges_log_R = np.linspace(np.nanmin(log_R), np.nanmax(log_R), int(10*(np.nanmax(log_R) - np.nanmin(log_R))) + 1)
bin_edges_log_frac_unc = np.linspace(np.nanmin(log_frac_unc), 0.0, int(10*(0.0 - np.nanmin(log_frac_unc))) + 1)

print("bin_edges_log_R", bin_edges_log_R)
print("bin_edges_log_frac_unc", bin_edges_log_frac_unc)


tot_star_hours = 0
plt_x = []
plt_y = []
plt_c = []

#log10_masses = np.linspace(-11, -9, 21)
log10_masses = np.arange(-11, -6 + 0.01, 0.1)

print("log10_masses", log10_masses)

subprocess.getoutput("rm -fr monte_carlo_results")
subprocess.getoutput("mkdir -p monte_carlo_results")

pwd = subprocess.getoutput("pwd")

jobs_by_filt = {filt_name: 0}

for i in tqdm.trange(len(bin_edges_log_frac_unc) - 1):
    for j in range(len(bin_edges_log_R) - 1):
        inds = np.where((log_frac_unc >= bin_edges_log_frac_unc[i])*(log_frac_unc < bin_edges_log_frac_unc[i+1])
                        *(log_R >= bin_edges_log_R[j])*(log_R < bin_edges_log_R[j+1])
                        )

        counts = sum(np.array(df[blue_filt + "_count"])[inds])


        if counts > 2:
            f = open("monte_carlo_results/tmp.sh", 'w')
            f.write("""#!/bin/bash
#SBATCH --job-name=mc
#SBATCH --partition=shared,kill-shared
#SBATCH --time=0-10:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G # Memory per node my job requires
#SBATCH --error=example-%A.err # %A - filled with jobid, where to write the stderr
#SBATCH --output=example-%A.out # %A - filled with jobid, wher to write the stdout
source ~/.bash_profile
""")

            median_log_R = np.median(log_R[inds])
            median_log_unc = np.median(log_frac_unc[inds])
            star_hours = counts*10.737*2./3600.

            print("filt_name", filt_name, "median_log_R", median_log_R, "median_log_unc", median_log_unc, "star_hours", star_hours)


            for log10_mass in log10_masses:
                f.write("cd " + pwd + "/monte_carlo_results/\n")
                f.write("echo 'median_log_R %f'\n" % median_log_R)
                f.write("echo 'star_hours %f'\n" % star_hours)
                f.write("echo 'filt_name %s'\n" % filt_name)
                f.write("echo 'log10_mass %f'\n" % log10_mass)
                f.write("python /home/drubin/NIRCam_ramp/step12_get_lens_count.py "  + str(10**median_log_R) + " " + str(star_hours) + " " + str(10**median_log_unc) + (" %.3g" % (10**log10_mass)) + " " + target + '\n')
                jobs_by_filt[filt_name] += 1

            f.write("echo 'done'\n")
            f.close()
            print(subprocess.getoutput("cd monte_carlo_results\n sbatch tmp.sh"))


f = open("jobs.txt", 'w')
f.write(str(jobs_by_filt))
f.close()

