import sys
from DavidsNM import miniNM_new
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline
import subprocess
import tqdm

def NMAD(vals):
    return 1.4826*np.nanmedian(np.abs(vals - np.nanmedian(vals)))

def modelfn(P, ts):
    #P[0] is log10 radius ratio
    #P[1] is offset from center
    return np.array(
        [ifn(P[0], np.sqrt(
            ((t - P[2])*P[3])**2.
            + P[1]**2.)
             )[0,0] for t in ts])


def chi2fn(P, passdata):
    ts = passdata[0][0]
    ys_short = passdata[0][1]
    sigys_short = passdata[0][2]

    ys_long = passdata[0][3]
    sigys_long = passdata[0][4]

    mod = modelfn(P, ts)
    

    chi2 = np.nansum(((ys_short - mod)/sigys_short)**2.) + np.nansum(((ys_long - mod)/sigys_long)**2.)
    #print(P)
    
    return chi2


def load_LCs():
    f = fits.open("amplfication.fits")
    print(f.info())
    dat = f[0].data

    dim1 = f[2].data
    dim2 = f[1].data
    
    ifn = RectBivariateSpline(
        dim1,
        dim2,
        dat,
        kx = 2, ky = 2,
    )
    print("first dim", dim1)
    print("second dim", dim2)
    f.close()

    return ifn, dim1, dim2

def fit_data(ts, ys_short, sigys_short, ys_long, sigys_long, start_t, const_F, prior_results, chi2_threshold = 24.5021):
    assert np.all(np.array(ts[:-1]) <= np.array(ts[1:])), str(ts)
    
    P, F, NA = miniNM_new(ministart = [0.5, 0.0, start_t, 1.], miniscale = [0.1, 1.0, 1., 0.1],
                          chi2fn = chi2fn, passdata = [ts, ys_short, sigys_short, ys_long, sigys_long], compute_Cmat = False, verbose = False)
    best_mod = modelfn(P, ts)

    if F < const_F - chi2_threshold:
        peak_ampl = np.max(best_mod)
        first_ampl = best_mod[0]
        last_ampl = best_mod[-1]
        
        #two_sigma_increase_from_end = ((peak_ampl - first_ampl) > 2*sigys_short[0]) and ((peak_ampl - last_ampl) > 2*sigys_short[-1]):
        #double_amplification_from_end = ()

        #if two_sigma_increase_from_end or double_amplification_from_end

        weights_short = 1./sigys_short**2.
        weights_long = 1./sigys_long**2.
        
        ampl_short = np.nansum(weights_short*(ys_short - 1.)*(best_mod - 1.))/np.nansum(weights_short*(best_mod - 1.)**2.)
        ampl_long = np.nansum(weights_long*(ys_long - 1.)*(best_mod - 1.))/np.nansum(weights_long*(best_mod - 1.)**2.)

        unc_ampl_short = 1./np.sqrt(np.nansum(weights_short*(best_mod - 1.)**2.))
        unc_ampl_long = 1./np.sqrt(np.nansum(weights_long*(best_mod - 1.)**2.))

        ampl_string = "ampl_short %.4f +- %.4f ampl_long %.4f +- %.4f" % (ampl_short, unc_ampl_short, ampl_long, unc_ampl_long)
        print("ampl_string", ampl_string)
        
        if np.abs(ampl_short - ampl_long)/np.sqrt(unc_ampl_short**2. + unc_ampl_long**2.) < 2.5:
            # Good candidate!

            we_have_seen_before = 0
            for prior_result in prior_results:
                if (np.abs(prior_result[0] - F) < 0.2) and (np.abs(P[2] - prior_result[1]) < 0.5):
                    we_have_seen_before = 1


            inds = np.where((best_mod - 1) > (np.max(best_mod) - 1)*0.05)
            ts_in_range = ts[inds]
            
            fraction_good_short = sum(all_data["sigys_short"][inds] > 0)/float(len(all_data["sigys_short"][inds]))
            fraction_good_long = sum(all_data["sigys_long"][inds] > 0)/float(len(all_data["sigys_long"][inds]))

            print("fraction_good_short", fraction_good_short)
            print("fraction_good_long", fraction_good_long)
            
            if (we_have_seen_before == 0) and (fraction_good_short > 0.75) and (fraction_good_long > 0.75):
                plt.figure(figsize = (12, 24))
                for zoom in [0, 1, 2]:
                    plt.subplot(3,1,1+zoom)


                    
                
                    plt.plot(ts, best_mod)
                    xlim = (np.min(ts_in_range) - 5, np.max(ts_in_range) + 5)
                    
                    
                    assert len(all_data["fls"]) == len(ts)
                    for i in range(len(ts)):
                        if zoom == 0:
                            if i == 0:
                                plt.text(ts[i], 0.97, all_data["fls"][i], rotation = 90, fontsize = 4, color = 'k')
                            else:
                                if all_data["fls"][i] != all_data["fls"][i - 1]:
                                    plt.text(ts[i], 0.97, all_data["fls"][i], rotation = 90, fontsize = 4, color = 'k')

                        else:
                            if (ts[i] > xlim[0]) and (ts[i] < xlim[1]):
                                plt.text(ts[i], 0.97, all_data["fls"][i], rotation = 90, fontsize = 4, color = 'k')

                            
                    plt.title("chi2 " + str(F) + " decrease " + str(const_F - F) + " best fit " + str(P) + "\n" + ampl_string)
                    
                    plt.errorbar(ts, ys_short, yerr = sigys_short, fmt = '.', color = 'b')
                    plt.errorbar(ts + 0.1, ys_long, yerr = sigys_long, fmt = '.', color = 'r')

                    if zoom > 0:
                        plt.xlim(xlim)
                        
                        if zoom > 1:
                            y_range = np.max(best_mod) - np.min(best_mod)
                            plt.ylim(np.min(best_mod) - y_range, np.max(best_mod) + y_range)
                plt.title(sys.argv[1] + " " + sys.argv[2])
                pltname = "cand=%s_filt=%s_dchi2=%.2f_tmax=%.2f_vel=%.2g_ampl=%.3f_%s.pdf" % (cand_to_read, short_filt, const_F - F, P[2], P[3], peak_ampl, ampl_string.replace(" ", "_"))
                plt.savefig("candidate_plots/" + pltname, bbox_inches = 'tight')
                plt.close()
                
        return (F, P[2])
    return [-999999, -999999]


def read_data(cand_to_read, file_to_read):
    cmd = "grep ' %s ' %s | grep %s" % (cand_to_read, file_to_read, short_filt) # grep ' 24 ' ../photo_flattened_linear.txt

    print(cmd)
    lines_read = subprocess.getoutput(cmd)


    all_data = dict(ts = [],
                    fls = [],

                    ys_short = [],
                    sigys_short = [],
                    short_RMS = [],

                    ys_long = [],
                    sigys_long = [],
                    long_RMS = [],)

    # jw02729001001_02103_00001_nrca1_uncallin.fits:0 43545 84.65979430930318 -69.09633490080327 times 59732.93585596396 59732.93610450562 59732.93635304729 59732.93660158895 59732.93685013062 59732.937098672286 59732.937347213956 59732.93759575562 short_filt F090W 1236 185 0.22391294881990767 -0.7324882094491496 short_phot: 3128.933632232352 3348.080203180964 2957.0807376172215 3081.856176150162 2910.3612615652796 3239.5706550576815 3186.589395061583 3154.97468923033 short_RMS: 0.0251837691060483 0.037125630734881156 0.041270264365757967 0.02718021102519018 0.03360368767288839 0.03388760941899862 0.03129939907389596 0.019034987414355553 short_uncs: 106.02008455176548 108.42765363164573 104.84326415285726 105.54090804966056 104.13676597522833 107.29114078839399 106.37848824343389 106.4754073008059 F335M 581 70 0.041244334067704336 -0.02928854937332942 long_phot: 1265.8623862128609 820.5394548580801 1033.8384973268792 983.8265661074148 878.2334024673046 899.9517964159764 1260.1505683980033 937.1954965324775 long_RMS: 0.04440373921966271 0.04701686220764472 0.048333798646865905 0.06334216090299061 0.060390343082170936 0.04878571089172337 0.04619503347952911 0.04406195187427845 long_uncs: 110.91810299181091 107.93771235420843 109.73492725101804 110.19203070251358 108.6635545698453 108.20401374181277 113.60316814860768 108.82106617545661


    short_filts = []
    
    for line in lines_read.split('\n'):
        parsed = line.split(None)


        if parsed[1] == cand_to_read:            
            fl = parsed[0]
            times_start = parsed.index("times")
            times_end = parsed.index("short_filt")
            short_filts.append(parsed[times_end + 1])

            
            short_phot_start = parsed.index("short_phot:")
            short_RMS_start = parsed.index("short_RMS:")
            short_uncs_start = parsed.index("short_uncs:")
            
            long_phot_start = parsed.index("long_phot:")
            long_RMS_start = parsed.index("long_RMS:")
            long_uncs_start = parsed.index("long_uncs:")
            
            n_obs = times_end - times_start - 1

            for i in range(n_obs):
                all_data["fls"].append(fl)
                all_data["ts"].append(float(parsed[times_start + i + 1]))
                
                all_data["short_RMS"].append(   float(parsed[short_RMS_start + i + 1])   )
                all_data["ys_short"].append(   float(parsed[short_phot_start + i + 1])   )
                all_data["sigys_short"].append(   float(parsed[short_uncs_start + i + 1])   )
                
                all_data["long_RMS"].append(   float(parsed[long_RMS_start + i + 1])   )
                all_data["ys_long"].append(   float(parsed[long_phot_start + i + 1])   )
                all_data["sigys_long"].append(   float(parsed[long_uncs_start + i + 1])   )

    assert len(set(short_filts)) == 1, str(short_filts)
                
    for key in all_data:
        all_data[key] = np.array(all_data[key])

    for short_long in ["short", "long"]:
        median_RMS = np.nanmedian(all_data[short_long + "_RMS"])
        NMAD_RMS = NMAD(all_data[short_long + "_RMS"])

        RMS_threshold = median_RMS + 5*NMAD_RMS
        assert np.isnan(RMS_threshold) == 0
        RMS_threshold = np.clip(RMS_threshold, 0.01, 0.2)

        print("RMS_threshold", RMS_threshold)

        inds = np.where((all_data[short_long + "_RMS"] > RMS_threshold) + np.isnan(all_data[short_long + "_RMS"]))
        all_data["sigys_" + short_long][inds] = np.nan

        med_flux = np.nanmedian(all_data["ys_" + short_long])
        all_data["ys_" + short_long] /= med_flux
        all_data["sigys_" + short_long] /= med_flux
        
    all_data["ts"] = (all_data["ts"] - np.min(all_data["ts"]))*86400./(10.737*2.)
    return all_data

    
ifn, dim1, dim2 = load_LCs()

cand_to_read = sys.argv[1]
file_to_read = sys.argv[2]
short_filt = sys.argv[3]


subprocess.getoutput("mkdir candidate_plots")

all_data = read_data(cand_to_read = cand_to_read, file_to_read = file_to_read)

print(all_data)


if 0:
    ts = np.arange(20, dtype=np.float64)
    ys = np.exp(-1*(ts - 10)**2.)*0.05 + 1.
    sigys = np.ones(20)*0.01
    
    sigys /= np.median(ys)
    ys /= np.median(ys)
    
    const_F = 2*sum(   ((ys - 1.)/sigys)**2.   )
    
    fit_data(ts = ts, start_t = 8, const_F = const_F,
             ys_short = ys + np.random.normal(size = len(ys))*sigys, sigys_short = sigys,
             ys_long = ys + np.random.normal(size = len(ys))*sigys*2, sigys_long = sigys*2)
if 1:
    
    const_F = np.nansum(   ((all_data["ys_short"] - 1.)/all_data["sigys_short"])**2.   ) + np.nansum(   ((all_data["ys_long"] - 1.)/all_data["sigys_long"])**2.   )

    print("const_F", const_F)
    prior_results = [] # F, P[2] which is t_max

    

    for this_start_t in all_data["ts"]:
        fit_results = fit_data(ts = all_data["ts"], start_t = this_start_t, const_F = const_F,
                               ys_short = all_data["ys_short"], sigys_short = all_data["sigys_short"],
                               ys_long = all_data["ys_long"], sigys_long = all_data["sigys_long"], prior_results = prior_results)

        if fit_results[0] > 0:
            prior_results.append(fit_results)

