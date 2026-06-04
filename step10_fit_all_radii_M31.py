import pandas as pd
import numpy as np
from DavidsNM import miniNM_new
from scipy.interpolate import LinearNDInterpolator
from FileRead import readcol
import matplotlib.pyplot as plt
import sys
import extinction
from tqdm.auto import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
from tqdm.auto import tqdm

def modelfn(P):
    #try:
    #    Teff, r_rsol, logg, meta = stellar_params_from_mass_age(P[0], P[1], feh=-0.5)
    #except:
    #    return np.array([100]*10)

    Teff1000 = P[0]
    r_rsol = P[1]
    logg = P[2]
    
    mags_1solR = np.array([
        ifns[filt](Teff1000, logg) for filt in HST_filt_list + JWST_filt_list])

    mags = mags_1solR - 2.5*np.log10(r_rsol**2.) + dist_mod

    
    mags += P[3]*A_lambda
    return mags


def chi2fn(P, passdata):
    obs_mags = passdata[0][0]
    obs_dmags = passdata[0][1]

    model = modelfn(P)
    pulls = (obs_mags - model)/obs_dmags
    if P[2] < 0:
        return 1e100

    if P[3] < 0:
        return 1e100

    
    chi2 = np.nansum(pulls**2.)
    #print(P, chi2)

    if np.isnan(chi2) or (chi2 == 0):
        return 1e100

    return chi2


def run_fit(miniscale, passdata, fit_extinction):
    bestF = 1e100


    bestP = np.sqrt([-1., -1., -1., -1.])
    
    for start_Teff1000 in [3., 5., 10., 20.]:
        P, F, Cmat = miniNM_new(ministart = [start_Teff1000, 1.0, 4.0, 0.3*fit_extinction],
                                miniscale = miniscale,
                                chi2fn = chi2fn,
                                passdata = passdata,
                                compute_Cmat = False)
        if F < bestF:
            bestP = P
            bestF = F
            
    return bestP, bestF


def fit_one_star(one_row):
    obs_mags_AB = []
    obs_dmags = []

    for filt, Vega_AB in zip(HST_filt_list, HST_Vega_10pc_AB):
        obs_mags_AB.append(one_row[filt.lower() + "_vega_hst"] + Vega_AB)
        obs_dmags.append(
            np.sqrt(one_row[filt.lower() + "_err_hst"]**2. + 0.05**2.)
        )

    for filt, AB_ZP in zip(JWST_filt_list, JWST_AB_ZPs):
        obs_mags_AB.append(AB_ZP - 2.5*np.log10(one_row[filt]/(10.73677*2.))
                           )
        obs_dmags.append(
            np.sqrt((1.0857*one_row[filt + "_unc"]/one_row[filt])**2. + 0.05**2)
            )

    
    passdata = (np.array(obs_mags_AB), np.array(obs_dmags))
    print("passdata", passdata)

    for fit_extinction in [1]:#[0, 1]:
        P, F = run_fit([1.0, 0.5, 1.0, 0.2*fit_extinction], passdata = passdata, fit_extinction = fit_extinction)

        mod = modelfn(P)

        return_dict = dict(A_V = P[3], logg = P[2], r_rsol = P[1], Teff1000 = P[0], chi2_SED_fit = F)

        for i, filt in enumerate(HST_filt_list + JWST_filt_list):
            return_dict["mod_" + filt.lower()] = mod[i]
            
        #plt.subplot(2,1,1+fit_extinction)
        #plt.errorbar(waves, passdata[0], yerr = np.clip(passdata[1], 0, 1), fmt = 'o')
        #plt.plot(waves, mod, '^')
        #plt.title(str(return_dict).replace(',', '\n'), size = 6)
    #plt.show()
    #plt.close()
    
    
    return return_dict

def load_model_atm():
    df_ma = pd.read_csv(model_atmosphere_grid, sep=r"\s+")
    print(df_ma)
    
    Teff1000 = [float(item.split("_t")[-1].split("_")[0])/1000. for item in df_ma["#file"]]
    logg = [float(item.split("_g")[-1].split("_")[0]) for item in df_ma["#file"]]

    print("Teff1000", Teff1000, len(Teff1000))
    print("logg", logg, len(logg))

    ifns = {}
    for filt in HST_filt_list + JWST_filt_list:
        ifns[filt] = LinearNDInterpolator(list(zip(Teff1000, logg)), df_ma[filt + "_one_Rsol"], fill_value=np.nan)

    return ifns

model_atmosphere_grid = sys.argv[1]

HST_filt_list =     ["F275W",                "F336W",            "F475W",              "F814W",                                    "F110W",            "F160W"]
# If Vega were at 10 parsecs, it would have these AB magnitudes. Absolute AB mags.
HST_Vega_10pc_AB = [1.5118909390959674, 1.1481593147969968,    -0.10883313322821817,  0.42976555918667475,               0.7767685103427203, 1.2742877498247838]

JWST_filt_list = ["F150W", "F277W"]
JWST_AB_ZPs = [27.8139, 27.8803]

waves = np.array([2750., 3360., 4750., 8140., 11000., 16000.] + [15000., 27500.])

A_lambda = extinction.fitzpatrick99(waves, a_v = 1.0, r_v = 3.1)

print("A_lambda", A_lambda)


dist_mod = 24.407 # 24.407 for M 31


ifns = load_model_atm()



job_index = int(sys.argv[2])
n_jobs = int(sys.argv[3])


df = pd.read_csv("my_with_hst.csv")
print(df)

# Keep only this job's rows: job_index, job_index+n_jobs, ...
df_sub = df.iloc[job_index::n_jobs].copy()

print(f"Running job {job_index} of {n_jobs}")
print(f"Fitting {len(df_sub)} rows out of {len(df)} total rows")

tqdm.pandas()

new_cols = df_sub.progress_apply(fit_one_star, axis=1, result_type="expand")
df_sub = pd.concat([df_sub, new_cols], axis=1)

outname = f"my_with_hst_fit_{job_index:02d}_of_{n_jobs:02d}.csv"
df_sub.to_csv(outname, sep=",", index=False)

print(f"Wrote {outname}")

#tqdm.pandas()

#new_cols = df.progress_apply(fit_one_star, axis=1, result_type="expand")
#df = pd.concat([df, new_cols], axis=1)

#df.to_csv("my_with_hst_fit.csv", sep=",", index=False)
