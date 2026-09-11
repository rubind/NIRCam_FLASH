import pandas as pd
import numpy as np
from DavidsNM import miniNM_new
import matplotlib.pyplot as plt
from scipy.stats import scoreatpercentile
import subprocess
import tqdm
import sys

my_with_hst_fit_csv = sys.argv[1]
df = pd.read_csv(my_with_hst_fit_csv)

filt_name = sys.argv[2]
target = sys.argv[3]
HSC_depth = float(sys.argv[4])
hours = float(sys.argv[5])
cadence = float(sys.argv[6])


print(df)
print(df.columns)


for key in ["chi2_SED_fit", "r_rsol", "mod_f475w", "mod_f814w", "mod_f160w", "f160w_vega_hst"]:
    df[key] = pd.to_numeric(df[key], errors="coerce")

print("mod_f160w - f160w_vega_hst median", np.nanmedian(df["mod_f160w"] - df["f160w_vega_hst"]))

saturation_AB = 29.2843 - 2.5*np.log10(5e5/200)

print("saturation_AB", saturation_AB)

inds = np.where((df["chi2_SED_fit"] < 200)*(df["mod_f160w"] > saturation_AB))



log10_R = np.log10(np.abs(np.array(df["r_rsol"])[inds]))
r_band = np.array(df["mod_f160w"])[inds]

HSC_flux_ZP_equal_depth = 10.**(-0.4*(r_band - HSC_depth))

electrons_per_group = 10.**(-0.4*(r_band - 29.2843))*cadence

fractional_noise = np.sqrt((0.2/HSC_flux_ZP_equal_depth)**2. + 1./electrons_per_group)

log10_frac_unc = np.log10(fractional_noise)

plt.plot(r_band[::100], fractional_noise[::100], '.')
plt.yscale('log')
plt.savefig("fractional_noise_vs_mag.pdf", bbox_inches = 'tight')
plt.close()


bin_edges_log10_R = np.linspace(np.nanmin(log10_R), np.nanmax(log10_R), int(10*(np.nanmax(log10_R) - np.nanmin(log10_R))) + 1)
bin_edges_log10_frac_unc = np.linspace(np.nanmin(log10_frac_unc), 0.0, int(10*(0.0 - np.nanmin(log10_frac_unc))) + 1)

print("bin_edges_log10_R", bin_edges_log10_R)
print("bin_edges_log10_frac_unc", bin_edges_log10_frac_unc)


tot_star_hours = 0
plt_x = []
plt_y = []
plt_c = []

#log10_masses = np.linspace(-11, -9, 21)
log10_masses = np.arange(-12., -6 + 0.01, 0.1)

print("log10_masses", log10_masses)

subprocess.getoutput("rm -fr monte_carlo_results")
subprocess.getoutput("mkdir -p monte_carlo_results")

pwd = subprocess.getoutput("pwd")

jobs_by_filt = {filt_name: 0}

for i in tqdm.trange(len(bin_edges_log10_frac_unc) - 1):
    for j in range(len(bin_edges_log10_R) - 1):
        inds = np.where((log10_frac_unc >= bin_edges_log10_frac_unc[i])*(log10_frac_unc < bin_edges_log10_frac_unc[i+1])
                        *(log10_R >= bin_edges_log10_R[j])*(log10_R < bin_edges_log10_R[j+1])
                        )

        star_hours = len(log10_frac_unc[inds])*hours

        if star_hours > 0:
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

            median_log10_R = np.median(log10_R[inds])
            median_log10_unc = np.median(log10_frac_unc[inds])
            assert len(log10_frac_unc[inds]) == len(log10_frac_unc[inds])
        
            print("filt_name", filt_name, "median_log10_R", median_log10_R, "median_log10_unc", median_log10_unc, "star_hours", star_hours)
        
            
            for log10_mass in log10_masses:
                f.write("cd " + pwd + "/monte_carlo_results/\n")
                f.write("echo 'median_log_R %f'\n" % median_log10_R)
                f.write("echo 'star_hours %f'\n" % star_hours)
                f.write("echo 'filt_name %s'\n" % filt_name)
                f.write("echo 'log10_mass %f'\n" % log10_mass)
                f.write("python /home/drubin/NIRCam_ramp/step12_get_lens_count.py "  + str(10**median_log10_R) + " " + str(star_hours) + " " + str(10**median_log10_unc) + (" %.3g" % (10**log10_mass)) + " " + target + " " + str(cadence) + ' \n')
                jobs_by_filt[filt_name] += 1

            f.write("echo 'done'\n")
            f.close()
            print(subprocess.getoutput("cd monte_carlo_results\n sbatch tmp.sh"))

f = open("jobs.txt", 'w')
f.write(str(jobs_by_filt))
f.close()

