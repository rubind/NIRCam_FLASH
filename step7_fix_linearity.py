import numpy as np
import matplotlib.pyplot as plt
import tqdm
import sys


# jw02729001001_02103_00001_nrca1_uncallin.fits:0  42540  84.6691509440011  -69.09478505652189  times  59732.93585596396  59732.93610450562  59732.93635304729  59732.93660158895  59732.93685013062  59732.937098672286  59732.937347213956  59732.93759575562  short_filt  F090W  1652  109  0.21169471059977285  -0.5980001439610282  short_phot:  89755.0517853141  89708.1518904868  89145.60936671353  88886.33740864879  88659.26734903519  98734.58346348564  103986.62487594178  102015.1799725641  short_RMS:  0.012119567413669685  0.010217574570189418  0.012325207514663445  0.012664128409625193  0.015881244321152968  1.0  1.0  1.0  short_uncs:  409.78569585433246  408.9274638458859  407.21456345317034  406.4690248120894  405.51550404520543  -1.0  -1.0  -1.0  F335M  786  32  0.5275937151912111  0.24807903101063355  long_phot:  27871.11665091134  27712.748437159324  27306.745534768204  27435.222882108443  27562.409399594435  27496.601986955175  27079.246635789805  26902.169818424587  long_RMS:  0.023133097136320575  0.02010421516108168  0.021063772914666445  0.016929242087888964  0.024067144981814507  0.024955015338508014  0.019824746861178712  0.016987721497903008  long_uncs:  239.19720551362118  238.49701767534935  237.0403893882927  237.16065375574422  237.87221568006885  237.91525194452507  236.08709400591965  236.1535784985878



filt = sys.argv[1]
short_or_long = sys.argv[2]
frame = int(sys.argv[3])

all_data_by_star = {}

f = open("photo_unflat.txt", 'r')
lines = f.read().split('\n')
f.close()

unique_filts = []

xvals = []
yvals = []
detectors = []

subpix_x = []
subpix_y = []

for line in tqdm.tqdm(lines[::]):
    if line.count(".fits") and line.count(filt):
        parsed = line.split(None)

        

        short_phot_start = parsed.index("short_phot:")
        short_phot_end = parsed.index("short_RMS:")
        short_RMS_end = parsed.index("short_uncs:")

        n_points = short_phot_end - short_phot_start - 1
        
        short_phot_vals = []
        short_unc_vals = []
        
        for i, j, k in zip(range(short_phot_start + 1, short_phot_end), range(short_phot_end + 1, short_RMS_end), range(short_RMS_end + 1, short_RMS_end + 1 + n_points)):
            if float(parsed[j]) < 0.2:
                short_phot_vals.append(float(parsed[i]))
                short_unc_vals.append(float(parsed[k]))
                
        short_filt = parsed[parsed.index("short_filt") + 1]
        
        long_phot_start = parsed.index("long_phot:")
        long_phot_end = parsed.index("long_RMS:")
        long_RMS_end = parsed.index("long_uncs:")
        long_phot_vals = []
        long_unc_vals = []

        
        for i, j, k in zip(range(long_phot_start + 1, long_phot_end), range(long_phot_end + 1, long_RMS_end), range(long_RMS_end + 1, long_RMS_end + 1 + n_points)):
            if float(parsed[j]) < 0.2:
                long_phot_vals.append(float(parsed[i]))
                long_unc_vals.append(float(parsed[k]))
                
        long_filt = parsed[long_phot_start - 3]

        if short_or_long == "short":
            phot_to_use = np.array(short_phot_vals)
        else:
            phot_to_use = np.array(long_phot_vals)
            
        if len(phot_to_use) == n_points:
            #plt.plot(short_phot_vals[0], (short_phot_vals[-1] - short_phot_vals[0])/short_phot_vals[0], '.', color = 'b')
            detectors.append(parsed[0].split("_uncallin")[0].split("_")[-1])
            assert detectors[-1].count("nrc") == 1
            xvals.append(phot_to_use[frame])#phot_to_use[0])
            yvals.append(np.log(phot_to_use[frame]/np.nanmedian(phot_to_use)))#(phot_to_use[-1] - np.median(phot_to_use))/np.median(phot_to_use))

            subpix_ind = parsed.index(short_or_long + "_phot:") - 2
            subpix_x.append(float(parsed[subpix_ind]) % 1.)
            subpix_y.append(float(parsed[subpix_ind+1]) % 1.)
            

bin_edges = 10**np.linspace(3, np.log10(np.max(xvals)), 25)
print("bin_edges", bin_edges)

xvals = np.array(xvals)
yvals = np.array(yvals)
subpix_x = np.array(subpix_x)
subpix_y = np.array(subpix_y)

print("subpix_x", subpix_x)
print("subpix_y", subpix_y)

detectors = np.array(detectors)

plt.subplot(1,2,1)
plt.hist(subpix_x, bins = 100) #hexbin(subpix_x, subpix_y, gridsize = 100)
plt.subplot(1,2,2)
plt.hist(subpix_y, bins = 100) #hexbin(subpix_x, subpix_y, gridsize = 100)
plt.savefig("subpix_xy_" + filt + ".pdf", bbox_inches = 'tight')
plt.close()


for detector in np.unique(detectors):
    print(detector)
    
    binx = []
    biny = []
    for i in range(len(bin_edges) - 1):
        inds = np.where((xvals > bin_edges[i])*(xvals <= bin_edges[i+1])*(detectors == detector))
        binx.append(np.nanmedian(xvals[inds]))
        biny.append(np.nanmedian(yvals[inds]))

    plt.plot(binx, biny, label = detector)
plt.axhline(0)
plt.title(filt)
plt.legend(loc = 'best')

plt.xscale('log')
plt.savefig("bin_by_detector_" + filt + "_nbins=%02i_frame=%02i.pdf" % (len(bin_edges) - 1, frame))
plt.close()

subpix_bins = np.linspace(0, 1, 4)
for j in range(len(subpix_bins) - 1):
    for k in range(len(subpix_bins) - 1):

    
        binx = []
        biny = []
        for i in range(len(bin_edges) - 1):
            inds = np.where((xvals > bin_edges[i])*(xvals <= bin_edges[i+1])
                            *(subpix_x > subpix_bins[j])*(subpix_x <= subpix_bins[j+1])
                            *(subpix_y > subpix_bins[k])*(subpix_y <= subpix_bins[k+1]))
            binx.append(np.nanmedian(xvals[inds]))
            biny.append(np.nanmedian(yvals[inds]))
            
        just_subpix_inds = np.where((subpix_x > subpix_bins[j])*(subpix_x <= subpix_bins[j+1])
                                    *(subpix_y > subpix_bins[k])*(subpix_y <= subpix_bins[k+1]))
            
        subpix_pos_x = np.nanmedian(subpix_x[just_subpix_inds])
        subpix_pos_y = np.nanmedian(subpix_y[just_subpix_inds])
        
        plt.plot(binx, biny, label = "%.2f %.2f" % (subpix_pos_x, subpix_pos_y))
plt.axhline(0)
plt.title(filt)
plt.legend(loc = 'best')
plt.xscale('log')
plt.savefig("bin_by_subpix_" + filt + "_nbins=%02i_frame=%02i.pdf" % (len(bin_edges) - 1, frame))
plt.close()
