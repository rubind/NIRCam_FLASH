import pandas as pd
import numpy as np
from DavidsNM import miniNM_new
import matplotlib.pyplot as plt
from scipy.stats import scoreatpercentile
import subprocess

df = pd.read_csv("my_with_hst_fit.csv")
inds = np.where((df["chi2_SED_fit"] < 200))

print(df)

log_R = np.log10(np.abs(np.array(df["r_rsol"])[inds]))
log_frac_unc_090 = np.log10(np.array(df["F090W_unc"]/df["F090W"])[inds])
log_frac_unc_200 = np.log10(np.array(df["F200W_unc"]/df["F200W"])[inds])

bin_edges_log_R = np.linspace(np.nanmin(log_R), np.nanmax(log_R), int(10*(np.nanmax(log_R) - np.nanmin(log_R))) + 1)
bin_edges_log_frac_unc_090 = np.linspace(np.nanmin(log_frac_unc_090), 0.0, int(10*(0.0 - np.nanmin(log_frac_unc_090))) + 1)
bin_edges_log_frac_unc_200 = np.linspace(np.nanmin(log_frac_unc_200), 0.0, int(10*(0.0 - np.nanmin(log_frac_unc_200))) + 1)

print("bin_edges_log_R", bin_edges_log_R)
print("bin_edges_log_frac_unc_090", bin_edges_log_frac_unc_090)
print("bin_edges_log_frac_unc_200", bin_edges_log_frac_unc_200)

bin_edges_log_frac_unc = dict(F090W = bin_edges_log_frac_unc_090, F200W = bin_edges_log_frac_unc_200)
log_frac_unc = dict(F090W = log_frac_unc_090, F200W = log_frac_unc_200)


tot_star_hours = 0
plt_x = []
plt_y = []
plt_c = []

subprocess.getoutput("rm -fr monte_carlo_results")
subprocess.getoutput("mkdir monte_carlo_results")

for filt_name in ["F090W", "F200W"]:
    for i in range(len(bin_edges_log_frac_unc[filt_name]) - 1):
        for j in range(len(bin_edges_log_R) - 1):
            inds = np.where((log_frac_unc[filt_name] >= bin_edges_log_frac_unc[filt_name][i])*(log_frac_unc[filt_name] < bin_edges_log_frac_unc[filt_name][i+1])
                            *(log_R >= bin_edges_log_R[j])*(log_R < bin_edges_log_R[j+1])
                            )

            counts = sum(np.array(df[filt_name + "_count"])[inds])

            
            if counts > 0:
                median_log_R = np.median(log_R[inds])
                median_log_unc = np.median(log_frac_unc[filt_name][inds])
                star_hours = counts*10.737*2./3600.
                
                print(filt_name, median_log_R, median_log_unc, star_hours)

                plt_x.append(median_log_R)
                plt_y.append(median_log_unc)
                plt_c.append(star_hours)

                f = open("monte_carlo_results/tmp.sh", 'w')
                f.write("""
""")
                f.close()
                print(subprocess.getoutput("cd monte_carlo_results\n sbatch tmp.sh")
                
                if median_log_unc < -2.:
                    tot_star_hours += star_hours
print("tot_star_hours", tot_star_hours)

plt.scatter(plt_x, plt_y, c = plt_c)
plt.colorbar()
plt.savefig("binned_hours.pdf")
plt.close()
