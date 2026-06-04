import pandas as pd
import numpy as np
from DavidsNM import miniNM_new
import matplotlib.pyplot as plt
from scipy.stats import scoreatpercentile
import subprocess
import tqdm
import sys

df = pd.read_csv("my_with_hst_fit.csv")
inds = np.where((df["chi2_SED_fit"] < 200))

blue_filt = sys.argv[1]
red_filt = sys.argv[2]

print(df)

log_R = np.log10(np.abs(np.array(df["r_rsol"])[inds]))
log_frac_unc_blue = np.log10(np.array(df[blue_filt + "_unc"]*np.sqrt(df[blue_filt + "_count"])/df[blue_filt])[inds])
log_frac_unc_red = np.log10(np.array(df[red_filt + "_unc"]*np.sqrt(df[blue_filt + "_count"])/df[red_filt])[inds])

bin_edges_log_R = np.linspace(np.nanmin(log_R), np.nanmax(log_R), int(10*(np.nanmax(log_R) - np.nanmin(log_R))) + 1)
bin_edges_log_frac_unc_blue = np.linspace(np.nanmin(log_frac_unc_blue), 0.0, int(10*(0.0 - np.nanmin(log_frac_unc_blue))) + 1)
bin_edges_log_frac_unc_red = np.linspace(np.nanmin(log_frac_unc_red), 0.0, int(10*(0.0 - np.nanmin(log_frac_unc_red))) + 1)

print("bin_edges_log_R", bin_edges_log_R)
print("bin_edges_log_frac_unc_blue", bin_edges_log_frac_unc_blue)
print("bin_edges_log_frac_unc_red", bin_edges_log_frac_unc_red)

bin_edges_log_frac_unc = dict(blue_filt = bin_edges_log_frac_unc_blue, red_filt = bin_edges_log_frac_unc_red)
log_frac_unc = dict(blue_filt = log_frac_unc_blue, red_filt = log_frac_unc_red)


tot_star_hours = 0
plt_x = []
plt_y = []
plt_c = []

#log10_masses = np.linspace(-11, -9, 21)
log10_masses = np.arange(-11, -6 + 0.01, 0.1)

print("log10_masses", log10_masses)

subprocess.getoutput("rm -fr monte_carlo_results")

for log10_mass in log10_masses:
    subprocess.getoutput("mkdir -p monte_carlo_results/%.3f" % (-log10_mass))

pwd = subprocess.getoutput("pwd")

for filt_name in [blue_filt, red_filt]:
    for i in tqdm.trange(len(bin_edges_log_frac_unc[filt_name]) - 1):
        for j in range(len(bin_edges_log_R) - 1):
            inds = np.where((log_frac_unc[filt_name] >= bin_edges_log_frac_unc[filt_name][i])*(log_frac_unc[filt_name] < bin_edges_log_frac_unc[filt_name][i+1])
                            *(log_R >= bin_edges_log_R[j])*(log_R < bin_edges_log_R[j+1])
                            )

            counts = sum(np.array(df[filt_name + "_count"])[inds])

            
            if counts > 2:
                median_log_R = np.median(log_R[inds])
                median_log_unc = np.median(log_frac_unc[filt_name][inds])
                star_hours = counts*10.737*2./3600.
                
                print(filt_name, median_log_R, median_log_unc, star_hours)

                plt_x.append(median_log_R)
                plt_y.append(median_log_unc)
                plt_c.append(star_hours)

                for log10_mass in log10_masses:
                    f = open("monte_carlo_results/%.3f/tmp.sh" % (-log10_mass), 'w')
                    f.write("""#!/bin/bash
#SBATCH --job-name=mc
#SBATCH --partition=shared,kill-shared
#SBATCH --time=0-01:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G # Memory per node my job requires
#SBATCH --error=example-%A.err # %A - filled with jobid, where to write the stderr
#SBATCH --output=example-%A.out # %A - filled with jobid, wher to write the stdout
source ~/.bash_profile

cd """ + pwd + "/monte_carlo_results/%.3f\n" % (-log10_mass))
                    f.write("python ../../step12_get_lens_count.py "  + str(10**median_log_R) + " " + str(star_hours) + " " + str(10**median_log_unc) + " %.3g\n" % (10**log10_mass))
                
                    f.close()
                    print(subprocess.getoutput("cd monte_carlo_results/%.3f\n sbatch tmp.sh" % (-log10_mass)))
                
                if median_log_unc < -2.:
                    tot_star_hours += star_hours
print("tot_star_hours", tot_star_hours)

plt.scatter(plt_x, plt_y, c = plt_c)
plt.colorbar()
plt.savefig("binned_hours.pdf")
plt.close()
