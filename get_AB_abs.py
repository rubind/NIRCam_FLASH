from FileRead import readcol, file_to_fn
from scipy.interpolate import interp1d
import numpy as np
import sys
import glob

def get_AB_mag(sed, filt, lambs):
    return -2.5*np.log10(
        (lambs*filt*sed).sum() /
        (0.10884806248*filt/lambs).sum()
        )

sed_fl = sys.argv[1]
distance_pc = float(sys.argv[2])

lambs = np.arange(2000., 50000.)

sed = file_to_fn(sed_fl, kind = 'linear')(lambs)

for fl in np.sort(glob.glob("nircam_throughputs/mean_throughputs/F*W_May*txt") + ["nircam_throughputs/mean_throughputs/F335M_May2024_mean_system_throughput.txt"] + glob.glob("nircam_throughputs/HST*dat")):
    [x, y] = readcol(fl, 'ff')

    if fl.count("HST") == 0:
        x *= 10000.
    else:
        print("HST found, not scaling")
    
    filt = interp1d(x, y, kind = 'linear', fill_value = 0., bounds_error = False)(lambs)
    
    AB_mag = get_AB_mag(sed, filt, lambs)
    print("AB_mag absolute", fl.split("/")[-1].split("_")[0], sum(filt*lambs)/sum(filt), AB_mag - 5*np.log10(distance_pc/10.))
    


"""
python get_AB_abs.py ~/Dropbox/SCP_Stuff/calspec/alpha_lyr_mod_004.ascii 10 | grep AB
AB_mag absolute HST 11623.77313876173 0.7767685103427203
AB_mag absolute HST 15392.322533082977 1.2742877498247838
AB_mag absolute HST 2708.0293217479802 1.5118909390959674
AB_mag absolute HST 3358.349944847046 1.1481593147969968
AB_mag absolute HST 5334.314169724697 -0.02779740708302566
AB_mag absolute HST 7659.937010704895 0.38789825812813267

Add to Vega, get AB
"""
