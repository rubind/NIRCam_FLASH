import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import tqdm
import sys
from scipy.interpolate import interp1d
from DavidsNM import miniLM_new



# jw02729001001_02103_00001_nrca1_uncallin.fits:0  42540  84.6691509440011  -69.09478505652189  times  59732.93585596396  59732.93610450562  59732.93635304729  59732.93660158895  59732.93685013062  59732.937098672286  59732.937347213956  59732.93759575562  short_filt  F090W  1652  109  0.21169471059977285  -0.5980001439610282  short_phot:  89755.0517853141  89708.1518904868  89145.60936671353  88886.33740864879  88659.26734903519  98734.58346348564  103986.62487594178  102015.1799725641  short_RMS:  0.012119567413669685  0.010217574570189418  0.012325207514663445  0.012664128409625193  0.015881244321152968  1.0  1.0  1.0  short_uncs:  409.78569585433246  408.9274638458859  407.21456345317034  406.4690248120894  405.51550404520543  -1.0  -1.0  -1.0  F335M  786  32  0.5275937151912111  0.24807903101063355  long_phot:  27871.11665091134  27712.748437159324  27306.745534768204  27435.222882108443  27562.409399594435  27496.601986955175  27079.246635789805  26902.169818424587  long_RMS:  0.023133097136320575  0.02010421516108168  0.021063772914666445  0.016929242087888964  0.024067144981814507  0.024955015338508014  0.019824746861178712  0.016987721497903008  long_uncs:  239.19720551362118  238.49701767534935  237.0403893882927  237.16065375574422  237.87221568006885  237.91525194452507  236.08709400591965  236.1535784985878



def read_data():
    f = open("photo_unflat.txt", 'r')
    lines = f.read().split('\n')
    f.close()

    data_to_fit = {}
    
    for line in tqdm.tqdm(lines[::]):
        if line.count(".fits"):
            parsed = line.split(None)

            short_detector = parsed[0].split("_")[3]
            long_detector = short_detector[:-1] + "long"

            short_filt = parsed[parsed.index("short_filt") + 1]


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


            long_phot_start = parsed.index("long_phot:")
            long_phot_end = parsed.index("long_RMS:")
            long_RMS_end = parsed.index("long_uncs:")
            long_phot_vals = []
            long_unc_vals = []
            long_filt = parsed[long_phot_start - 5]

            
            #print("short_detector", short_filt, short_detector)
            #print("long_detector", long_filt, long_detector)


            for i, j, k in zip(range(long_phot_start + 1, long_phot_end), range(long_phot_end + 1, long_RMS_end), range(long_RMS_end + 1, long_RMS_end + 1 + n_points)):
                if float(parsed[j]) < 0.2:
                    long_phot_vals.append(float(parsed[i]))
                    long_unc_vals.append(float(parsed[k]))


            if (len(short_phot_vals) == n_points) and (np.min(short_phot_vals) > 1000):
                for frame in range(n_points):
                    short_data_key = (short_detector, short_filt, frame)
                    if short_data_key not in data_to_fit:
                        data_to_fit[short_data_key] = dict(frame_phot = [], log_frame_resid = [], unc_log_frame_resid = [])
                        
                    data_to_fit[short_data_key]["frame_phot"].append(short_phot_vals[frame])
                    data_to_fit[short_data_key]["log_frame_resid"].append(np.log(  short_phot_vals[frame]/np.nanmedian(short_phot_vals)  ))
                    data_to_fit[short_data_key]["unc_log_frame_resid"].append(np.abs(  short_unc_vals[frame]/np.nanmedian(short_phot_vals)  ))
                    
            if (len(long_phot_vals) == n_points) and (np.min(long_phot_vals) > 1000):
                for frame in range(n_points):
                    long_data_key = (long_detector, long_filt, frame)
                    if long_data_key not in data_to_fit:
                        data_to_fit[long_data_key] = dict(frame_phot = [], log_frame_resid = [], unc_log_frame_resid = [])
                        
                    data_to_fit[long_data_key]["frame_phot"].append(long_phot_vals[frame])
                    data_to_fit[long_data_key]["log_frame_resid"].append(np.log(  long_phot_vals[frame]/np.nanmedian(long_phot_vals)  ))
                    data_to_fit[long_data_key]["unc_log_frame_resid"].append(np.abs(  long_unc_vals[frame]/np.nanmedian(long_phot_vals)  ))
    return data_to_fit

def get_ifn(nodes, P):
    ifn = interp1d(nodes, P, kind = 'cubic', bounds_error = False, fill_value = P[0])
    return ifn

def residfn(P, passdata):
    nodes, these_data_to_fit = passdata[0]

    ifn = get_ifn(nodes, P)
    model = ifn(these_data_to_fit["frame_phot"])
    pulls = (these_data_to_fit["log_frame_resid"] - model)/these_data_to_fit["unc_log_frame_resid"]

    return np.sign(pulls)*np.sqrt(np.abs(pulls))


def fit_data(data_to_fit):
    ifns = {}
    pdf = PdfPages("all_splines.pdf")

    
    for key in tqdm.tqdm(data_to_fit):
        nodes = np.exp(np.linspace(np.log(1000.) - 0.0001,
                                   np.log(np.max(data_to_fit[key]["frame_phot"])) + 0.0001,
                                   11))
        print(key, nodes)
        P, F, Cmat = miniLM_new(ministart = np.zeros(len(nodes), dtype=np.float64),
                                miniscale = np.ones(len(nodes), dtype=np.float64),
                                residfn = residfn,
                                passdata = (nodes, data_to_fit[key]), verbose=True)
        ifn = get_ifn(nodes, P)
        ifns[key] = ifn
        
        fig = plt.figure()
        plt_x = np.exp(np.linspace(np.log(100.) - 0.0001,
                                   np.log(np.max(data_to_fit[key]["frame_phot"])) + 0.0001,
                                   101))
        plt_y = ifn(plt_x)
        plt.plot(plt_x, plt_y)
        plt.errorbar(nodes, P, yerr= np.sqrt(np.diag(Cmat)), fmt = 'o')
        plt.title(str(key) + '\n' + str(F) + " npts: %i" % len(data_to_fit[key]["frame_phot"]))
        plt.xscale('log')
        
        pdf.savefig(fig)
        plt.close(fig)
    pdf.close()
    return ifns


def fix_linearity(ifns):
    f = open("photo_unflat.txt", 'r')
    lines = f.read().split('\n')
    f.close()

    new_lines = []
    
    # jw02729001001_02103_00001_nrca1_uncallin.fits:0  42540  84.6691509440011  -69.09478505652189  times  59732.93585596396  59732.93610450562  59732.93635304729  59732.93660158895  59732.93685013062  59732.937098672286  59732.937347213956  59732.93759575562  short_filt  F090W  1652  109  0.21169471059977285  -0.5980001439610282  short_phot:  89755.0517853141  89708.1518904868  89145.60936671353  88886.33740864879  88659.26734903519  98734.58346348564  103986.62487594178  102015.1799725641  short_RMS:  0.012119567413669685  0.010217574570189418  0.012325207514663445  0.012664128409625193  0.015881244321152968  1.0  1.0  1.0  short_uncs:  409.78569585433246  408.9274638458859  407.21456345317034  406.4690248120894  405.51550404520543  -1.0  -1.0  -1.0  F335M  786  32  0.5275937151912111  0.24807903101063355  long_phot:  27871.11665091134  27712.748437159324  27306.745534768204  27435.222882108443  27562.409399594435  27496.601986955175  27079.246635789805  26902.169818424587  long_RMS:  0.023133097136320575  0.02010421516108168  0.021063772914666445  0.016929242087888964  0.024067144981814507  0.024955015338508014  0.019824746861178712  0.016987721497903008  long_uncs:  239.19720551362118  238.49701767534935  237.0403893882927  237.16065375574422  237.87221568006885  237.91525194452507  236.08709400591965  236.1535784985878



    for line in tqdm.tqdm(lines):
        if line.count("short_phot:") == 0:
            new_lines.append(line)
        else:
            parsed = line.split(None)

            short_detector = parsed[0].split("_")[3]
            long_detector = short_detector[:-1] + "long"

            
            short_start = parsed.index("short_phot:")
            short_end = parsed.index("short_RMS:")

            n_frames = short_end - short_start - 1

            short_uncs_start = parsed.index("short_uncs:")


            long_start = parsed.index("long_phot:")
            long_end = parsed.index("long_RMS:")
            assert long_end - long_start - 1 == n_frames

            long_uncs_start = parsed.index("long_uncs:")


            short_filt = parsed[short_start - 5]
            long_filt = parsed[long_start - 5]

            for i in range(n_frames):
                key = (short_detector, short_filt, i)
                parsed[short_start + 1 + i] = str(float(parsed[short_start + 1 + i])*np.exp(-   ifns[key](float(parsed[short_start + 1 + i]))   ))
                parsed[short_uncs_start + 1 + i] = str(float(parsed[short_uncs_start + 1 + i])*np.exp(-   ifns[key](float(parsed[short_uncs_start + 1 + i]))   ))

            for i in range(n_frames):
                key = (long_detector, long_filt, i)
                parsed[long_start + 1 + i] = str(float(parsed[long_start + 1 + i])*np.exp(-   ifns[key](float(parsed[long_start + 1 + i]))   ))
                parsed[long_uncs_start + 1 + i] = str(float(parsed[long_uncs_start + 1 + i])*np.exp(-   ifns[key](float(parsed[long_uncs_start + 1 + i]))   ))

            new_lines.append(" ".join(parsed))    

    f = open("photo_unflat_linear.txt", 'w')
    f.write('\n'.join(new_lines))
    f.close()

        
data_to_fit = read_data()
print(data_to_fit.keys())

ifns = fit_data(data_to_fit)

fix_linearity(ifns)
