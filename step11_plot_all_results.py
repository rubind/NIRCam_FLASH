import pandas as pd
import numpy as np
from DavidsNM import miniNM_new
import matplotlib.pyplot as plt



df = pd.read_csv("my_with_hst_fit.csv")

print(df)
print(df.columns)

plts = 15
sqrtn = int(np.ceil(np.sqrt(plts)))
ind = 1

plt.figure(figsize = (5*sqrtn, 4*sqrtn))

df["A_V+0.01"] = df["A_V"] + 0.01

for x_key, y_key in [("r_rsol", "chi2_SED_fit"),
                     ("r_rsol", "A_V+0.01"),
                     ("r_rsol", "Teff1000"),
                     ("r_rsol", "mod_f090w"),
                     ("chi2_SED_fit", "mod_f090w"),
                     ("Teff1000", "A_V+0.01")]:
    plt.subplot(sqrtn, sqrtn, ind)
    ind += 1
    
    inds = np.where((df["chi2_SED_fit"] < 200))
    plt.plot(  np.abs(np.array(df[x_key])[inds]), np.array(df[y_key])[inds], '.', color = 'b', alpha = 0.1)
    plt.title("chi2 < 200")
    plt.xlabel(x_key)
    plt.ylabel(y_key)

    if np.abs(np.array(df[x_key])[inds]).max() > np.abs(np.array(df[x_key])[inds]).min()*20:
        plt.xscale('log')

    if np.abs(np.array(df[y_key])[inds]).max() > np.abs(np.array(df[y_key])[inds]).min()*20:
        plt.yscale('log')
        


plt.subplot(sqrtn, sqrtn, ind)
ind += 1


inds = np.where((df["chi2_SED_fit"] < 200)*(np.abs(np.array(df["r_rsol"])) < 10))
plt.hist(np.abs(np.array(df["r_rsol"])[inds]), bins = 100, color = 'b')
plt.title("chi2 < 200, R < 10, count: %i" % len(inds[0]))
plt.xlabel("R/Rsol")


plt.subplot(sqrtn, sqrtn, ind)
ind += 1

inds = np.where((df["chi2_SED_fit"] < 200))    
plt.hist(np.abs(np.array(df["mod_f090w"])[inds]), bins = 100, color = 'b')
plt.title("chi2 < 200")
plt.xlabel("mod_f090w")

for filt in ["f090w", "f200w", "f335m", "f444w"]:
    plt.subplot(sqrtn, sqrtn, ind)
    ind += 1

    inds = np.where((df["chi2_SED_fit"] < 30))    

    ZP = np.nanmedian(df["mod_" + filt] + 2.5*np.log10(df[filt.upper()]))
    obs_mag = ZP - 2.5*np.log10(np.array(df[filt.upper()])[inds])
    
    plt.plot(np.array(df["mod_" + filt])[inds],
             obs_mag - np.array(df["mod_" + filt])[inds],
             '.', color = 'b')

    
    plt.title(filt + " " + str(ZP) + " 1s: " + str(ZP - 2.5*np.log10(10.737*2)))
    
plt.subplot(sqrtn, sqrtn, ind)
ind += 1
inds = np.where((df["chi2_SED_fit"] < 200)*(df["A_V"] < 10))    

plt.scatter(np.array(df["RA"])[inds],
            np.array(df["Dec"])[inds],
            c = np.array(df["A_V"])[inds], s=0.2)



plt.subplot(sqrtn, sqrtn, ind)
ind += 1
inds = np.where((df["chi2_SED_fit"] < 200))

plt.plot(np.array(df["r_rsol"])[inds], (2.5/np.log(10.))*np.array(df["F090W_unc"]/df["F090W"])[inds], '.', color = 'b')
plt.xlabel("r_rsol")
plt.ylabel("$F090W$ Magnitude Unc")
plt.xscale('log')
plt.yscale('log')

plt.tight_layout()
plt.savefig("SED_fit.pdf", bbox_inches = 'tight')
plt.close()

