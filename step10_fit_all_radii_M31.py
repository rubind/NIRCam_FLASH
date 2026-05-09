import pandas as pd
import numpy as np
from DavidsNM import miniNM_new
from scipy.interpolate import LinearNDInterpolator
from FileRead import readcol
import matplotlib.pyplot as plt

def modelfn(P):
    #try:
    #    Teff, r_rsol, logg, meta = stellar_params_from_mass_age(P[0], P[1], feh=-0.5)
    #except:
    #    return np.array([100]*10)

    Teff1000 = P[0]
    r_rsol = P[1]
    logg = P[2]
    
    mags_1solR = np.array([
        ifns["F275W"](Teff1000, logg),
        ifns["F336W"](Teff1000, logg),
        ifns["F555W"](Teff1000, logg),
        ifns["F775W"](Teff1000, logg),
        ifns["F110W"](Teff1000, logg),
        ifns["F160W"](Teff1000, logg),
        
        ifns["F090W"](Teff1000, logg),
        ifns["F200W"](Teff1000, logg),
        ifns["F335M"](Teff1000, logg),
        ifns["F444W"](Teff1000, logg)])
    mags = mags_1solR - 2.5*np.log10(r_rsol**2.) + 18.477

    
    mags += P[3]*
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
        P, F, Cmat = miniNM_new(ministart = [start_Teff1000, 1.0, 4.0, 0.3*fit_extinction, 18.477],
                                miniscale = miniscale,
                                chi2fn = chi2fn,
                                passdata = passdata,
                                compute_Cmat = False)
        if F < bestF:
            bestP = P
            bestF = F
            
    return bestP, bestF


def fit_one_star(one_row):
    passdata = ([one_row["m_f275w_hst"] + 1.5118909390959674,
                 one_row["m_f336w_hst"] + 1.1481593147969968,
                 one_row["m_f555w_hst"] + -0.02779740708302566,
                 one_row["m_f775u_hst"] + 0.38789825812813267,
                 one_row["m_f110w_hst"] + 0.7767685103427203,
                 one_row["m_f160w_hst"] + 1.2742877498247838,
                 27.46090589142513 - 2.5*np.log10(one_row["F090W"]/(10.73677*2.)), #27.4519 - 2.5*np.log10(one_row["F090W"]/(10.73677*2.)), # These ZPs I worked out from integrating the passbands
                 28.09342044547302 - 2.5*np.log10(one_row["F200W"]/(10.73677*2.)), #27.9973 - 2.5*np.log10(one_row["F200W"]/(10.73677*2.)), # The uncommented ones are from the median of the residuals to the fits, i.e., they are from HST, should be more reliable because the star flats were arbitrarily normalized
                 27.376846702753912 - 2.5*np.log10(one_row["F335M"]/(10.73677*2.)), #27.2428 - 2.5*np.log10(one_row["F335M"]/(10.73677*2.)),
                 28.229016599400413 - 2.5*np.log10(one_row["F444W"]/(10.73677*2.))], #28.0491 - 2.5*np.log10(one_row["F444W"]/(10.73677*2.))],
                
                [np.sqrt(one_row["e_f275w_hst"]**2. + 0.05**2),
                 np.sqrt(one_row["e_f336w_hst"]**2. + 0.05**2),
                 np.sqrt(one_row["e_f555w_hst"]**2. + 0.05**2),
                 np.sqrt(one_row["e_f775u_hst"]**2. + 0.05**2),
                 np.sqrt(one_row["e_f110w_hst"]**2. + 0.05**2),
                 np.sqrt(one_row["e_f160w_hst"]**2. + 0.05**2),
                 np.sqrt((1.0857*one_row["F090W_unc"]/one_row["F090W"])**2. + 0.05**2),
                 np.sqrt((1.0857*one_row["F200W_unc"]/one_row["F200W"])**2. + 0.05**2),
                 np.sqrt((1.0857*one_row["F335M_unc"]/one_row["F335M"])**2. + 0.05**2),
                 np.sqrt((1.0857*one_row["F444W_unc"]/one_row["F444W"])**2. + 0.05**2)])

    print("passdata", passdata)

    for fit_extinction in [1]:#[0, 1]:
        P, F = run_fit([1.0, 0.5, 1.0, 0.2*fit_extinction], passdata = passdata, fit_extinction = fit_extinction)

        mod = modelfn(P)

        return_dict = dict(A_V = P[3], logg = P[2], r_rsol = P[1], Teff1000 = P[0], chi2_SED_fit = F,
                           mod_f275w = mod[0],
                           mod_f336w = mod[1],
                           mod_f555w = mod[2],
                           mod_f775w = mod[3],
                           mod_f110w = mod[4],
                           mod_f160w = mod[5],
                           mod_f090w = mod[6],
                           mod_f200w = mod[7],
                           mod_f335m = mod[8],
                           mod_f444w = mod[9])

        #plt.subplot(2,1,1+fit_extinction)
        #plt.errorbar([0.275, 0.336, 0.555, 0.775, 1.15, 1.55] + [0.9, 2.0, 3.35, 4.44], passdata[0], yerr = np.clip(passdata[1], 0, 1), fmt = 'o')
        #plt.plot([0.275, 0.336, 0.555, 0.775, 1.15, 1.55] + [0.9, 2.0, 3.35, 4.44], mod, '^')
        #plt.title(str(return_dict))
    #plt.show()
    #plt.close()
    
    
    return return_dict

def load_model_atm():
    [atm_fl, all_Lsol_one_Rsol, F090W_one_Rsol, F200W_one_Rsol, F335M_one_Rsol, F444W_one_Rsol, GaiaG_one_Rsol, F555W_one_Rsol, F336W_one_Rsol, F110W_one_Rsol, F775W_one_Rsol, F160W_one_Rsol, F275W_one_Rsol] = readcol("model_atmosphere_grid.txt", 'a,f,fffff,ffffff')

    atm_fl = np.array(atm_fl)

    good_mask = np.where(np.array([item.count("m-0.50_a+0.00_c+0.00") for item in atm_fl]))
    [atm_fl, all_Lsol_one_Rsol, F090W_one_Rsol, F200W_one_Rsol, F335M_one_Rsol, F444W_one_Rsol, GaiaG_one_Rsol, F555W_one_Rsol, F336W_one_Rsol, F110W_one_Rsol, F775W_one_Rsol, F160W_one_Rsol, F275W_one_Rsol] = [item[good_mask] for item in [atm_fl, all_Lsol_one_Rsol, F090W_one_Rsol, F200W_one_Rsol, F335M_one_Rsol, F444W_one_Rsol, GaiaG_one_Rsol, F555W_one_Rsol, F336W_one_Rsol, F110W_one_Rsol, F775W_one_Rsol, F160W_one_Rsol, F275W_one_Rsol]]


    
    Teff1000 = [float(item.split("_t")[-1].split("_")[0])/1000. for item in atm_fl]
    logg = [float(item.split("_g")[-1].split("_")[0]) for item in atm_fl]

    print("Teff1000", Teff1000, len(Teff1000))
    print("logg", logg, len(logg))

    ifns = {}
    for filt, idata in [
            ("F275W", F275W_one_Rsol),
            ("F336W", F336W_one_Rsol),
            ("F555W", F555W_one_Rsol),
            ("F775W", F775W_one_Rsol),
            ("F110W", F110W_one_Rsol),
            ("F160W", F160W_one_Rsol),
            
            ("F090W", F090W_one_Rsol),
            ("F200W", F200W_one_Rsol),
            ("F335M", F335M_one_Rsol),
            ("F444W", F444W_one_Rsol)]:
        ifns[filt] = LinearNDInterpolator(list(zip(Teff1000, logg)), idata, fill_value=np.nan)

    return ifns

ifns = load_model_atm()


df = pd.read_csv("my_with_hst.csv")

print(df)


from tqdm.auto import tqdm
tqdm.pandas()

new_cols = df.progress_apply(fit_one_star, axis=1, result_type="expand")
df = pd.concat([df, new_cols], axis=1)

df.to_csv("my_with_hst_fit.csv", sep=",", index=False)
