from FileRead import readcol
import numpy as np
import sys
from DavidsNM import miniLM_new
import tqdm
import json
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.interpolate import RectBivariateSpline


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)

    
def spatial_median(x, y, c, R = 50):
    coords = np.column_stack([x, y])
    tree = cKDTree(coords)


    # For each point, find neighbors within R
    idxs = tree.query_ball_point(coords, r=R)

    c_smooth = np.empty_like(c)
    for i, neighbors in enumerate(idxs):
        if len(neighbors) > 0:
            c_smooth[i] = np.median(c[neighbors])
        else:
            c_smooth[i] = c[i]  # fallback if isolated point

    plt.scatter(x, y, c=c_smooth, s=scatter_size, cmap='viridis')
    plt.colorbar()



stan_code = """
data {
    int n_star;
    int n_detector;
    int n_obs;

    vector [n_obs] phots;
    int <lower = 0, upper = n_detector - 1> detector_inds[n_obs];
    int <lower = 0, upper = n_detector - 1> star_inds[n_obs];
}

parameters {
    vector [n_detector - 1] rel_sens;
    vector [n_star] true_fluxes;
}

model {
    for (i in 1:n_obs) {
        phots[i] ~ laplace();
    }
}
"""

def parseP(P):
    parsed = {}
    ind = 0

    parsed["rel_sens"] = P[ind:ind+all_data["n_detectors"]-1]
    parsed["rel_sens"] = np.append(parsed["rel_sens"], 1.)
    ind += all_data["n_detectors"]-1

    for term in terms:
        parsed["rel_sens_" + term] = P[ind:ind+all_data["n_detectors"]]
        ind += all_data["n_detectors"]

    parsed["true_fluxes"] = P[ind:ind+all_data["n_stars"]]
    ind += all_data["n_stars"]
    return parsed
    

def modelfn(P):
    parsed = parseP(P)

    #the_mod = np.zeros(len(all_data["phots"]), dtype=np.float64)

    
    rel_sens_term = parsed["rel_sens"][all_data["detector_inds"]]

    for term in terms:
        rel_sens_term += (parsed["rel_sens_" + term][all_data["detector_inds"]])*all_data["term_" + term]

    the_mod = parsed["true_fluxes"][all_data["star_inds"]] * rel_sens_term
    
    return the_mod

def residfn(P, NA, save_plot = ""):
    parsed = parseP(P)
    the_mod = modelfn(P)

    the_resid = all_data["phots"] - the_mod
    if save_plot != "":
        print("median abs(resid)/model ", np.median(np.abs(the_resid)/the_mod))
        print(parsed)

        plt.figure(figsize = (12, 12))
        for i in range(all_data["n_detectors"]):
            plt.subplot(3,3,i+1)
            inds = np.where(all_data["detector_inds"] == i)
            plt.scatter(all_data["xs"][inds], all_data["ys"][inds], c = the_resid[inds]/the_mod[inds], vmin = -0.05, vmax = 0.05, s = scatter_size)
            plt.title(all_data["unique_detectors"][i])
        plt.tight_layout()
        plt.savefig(save_plot, bbox_inches = 'tight')
        plt.close()

        plt.figure(figsize = (12, 12))
        for i in range(all_data["n_detectors"]):
            plt.subplot(3,3,i+1)
            inds = np.where(all_data["detector_inds"] == i)
            spatial_median(x = all_data["xs"][inds], y = all_data["ys"][inds], c = the_resid[inds]/the_mod[inds], R = 100.)
            plt.title(all_data["unique_detectors"][i])
        plt.tight_layout()
        plt.savefig(save_plot.replace(".pdf", "_median.pdf"), bbox_inches = 'tight')
        plt.close()

        plt.figure(figsize = (12, 12))
        for i in range(all_data["n_detectors"]):
            plt.subplot(3,3,i+1)
            inds = np.where(all_data["detector_inds"] == i)
            plt.scatter(all_data["xs"][inds], the_mod[inds], c = the_resid[inds]/the_mod[inds], vmin = -0.05, vmax = 0.05, s = scatter_size)
            plt.title(all_data["unique_detectors"][i])
            plt.yscale('log')
        plt.tight_layout()
        plt.savefig(save_plot.replace(".pdf", "_vs_x.pdf"), bbox_inches = 'tight')
        plt.close()

        
    the_resid = np.sign(the_resid)*np.sqrt(np.abs(the_resid)/all_data["phots"])

    return the_resid


def read_data(phot_file, the_filt, min_flux):
    if int(the_filt.split("F")[-1].split("W")[0].split("M")[0]) < 210:
        print("short filter!")
        short_or_long = "short"
    else:
        print("long filter!")
        short_or_long = "long"
    
    
    f = open(phot_file, 'r')
    lines = f.read().split('\n')
    f.close()
    
    lines = [item for item in lines if item.count(the_filt) == 1]
    
    print("len(lines)", len(lines))

    all_data = dict(file_names = [], star_ids = [], xs = [], ys = [], phots = [])
    # jw02729001001_02103_00001_nrca1_uncallin.fits:0 42675 84.67559442487747 -69.09434245895841 times 59732.93585596396 59732.93610450562 59732.93635304729 59732.93660158895 59732.93685013062 59732.937098672286 59732.937347213956 59732.93759575562 short_filt F090W 1920 125 0.12782874275499426 -0.3970284832298483 short_phot: 152306.0199830377 150466.43500187932 542606.2706632036 236334.27375341379 212080.35905619944 205343.62184735062 202917.46598281458 195538.86125442907 short_RMS: 0.08163929671212547 0.08345615947833676 1.0 1.0 1.0 1.0 1.0 1.0 short_uncs: 567.8273376732021 564.3262569780543 -1.031831911070209 -1.0286161952102475 -1.032009875792325 -1.0343184100768454 -1.028933625336535 -1.030080128043523 F335M 919 40 0.6054811398401565 -0.07892624559492395 long_phot: 40059.70805165472 41233.85534541359 40290.384811619675 41128.7155941244 41222.279680248655 41262.91974549036 40979.92289193649 41228.45538038781 long_RMS: 0.06746241813984886 0.0682877260981806 0.07104524531340808 0.07606528086040933 0.07366265018681913 0.07671533610407383 0.0816865383586041 0.08708138531034623 long_uncs: 284.1567623876938 286.5305291607755 284.0160399513766 286.2041381815019 284.17797734085156 284.69832614835826 284.6893634189089 281.988601121711

    
    for line in tqdm.tqdm(lines):
        parsed = line.split(None)
        phot = parsed[parsed.index(short_or_long + "_phot:") + 1: parsed.index(short_or_long + "_RMS:")]
        RMS = parsed[parsed.index(short_or_long + "_RMS:") + 1: parsed.index(short_or_long + "_RMS:") + 1 + len(phot)]
        
        assert parsed[parsed.index(short_or_long + "_RMS:") + 1 + len(phot)] == short_or_long + "_uncs:"
        
        phot = [float(item) for item in phot]
        RMS = [float(item) for item in RMS]

        assert len(RMS) == len(phot)
        ind = parsed.index(the_filt)

        for i in range(len(phot))[::-1]:
            if RMS[i] > 0.2:
                del phot[i]
                del RMS[i]
                

        if (len(phot) > 3) and (min_flux > 0): # Just read in bright stars for fitting the model
            med_phot = np.median(phot)
            if med_phot > min_flux:
                
                all_data["file_names"].append(parsed[0].split(":")[0])
                all_data["star_ids"].append(parsed[1])
                
                all_data["xs"].append(float(parsed[ind+1]))
                all_data["ys"].append(float(parsed[ind+2]))
                
                all_data["phots"].append(med_phot)
        if min_flux <= 0: # If we need all stars, which we do when constructing the model!
            all_data["file_names"].append(parsed[0].split(":")[0])
            all_data["star_ids"].append(parsed[1])
            all_data["xs"].append(float(parsed[ind+1]))
            all_data["ys"].append(float(parsed[ind+2]))


    if short_or_long == "short":
        all_data["detectors"] = [item.split("_")[3] for item in all_data["file_names"]]
    else:
        all_data["detectors"] = [item.split("_")[3][:4] + "long" for item in all_data["file_names"]]

    all_data["unique_detectors"] = list(np.sort(np.unique(all_data["detectors"])))

    all_data["xs"] = np.array(all_data["xs"])
    all_data["ys"] = np.array(all_data["ys"])
    
    
    
    all_data["unique_stars"] = list(set(all_data["star_ids"]))
    print("len(unique_stars)", len(all_data["unique_stars"]))
    if min_flux > 0:
        all_data["star_inds"] = np.array([all_data["unique_stars"].index(item) for item in all_data["star_ids"]])
    else:
        all_data["star_inds"] = []
        
    all_data["n_stars"] = len(all_data["unique_stars"])
    all_data["n_detectors"] = len(all_data["unique_detectors"])
    all_data["detector_inds"] = np.array([all_data["unique_detectors"].index(item) for item in all_data["detectors"]])

    print("len(unique_detectors)", len(all_data["unique_detectors"]))

    all_data["phots"] = np.array(all_data["phots"])
    all_data["star_inds"] = np.array(all_data["star_inds"])





    try:
        int(sys.argv[4])
        poly_terms = 1
    except:
        poly_terms = 0

    if poly_terms:
        terms = ["x", "y", "xx", "yy", "xy"][:int(sys.argv[4])]

        for term in terms:
            if len(term) == 1:
                all_data["term_" + term] = (all_data[term + "s"] - 1024.)/1024.
            elif len(term) == 2:
                all_data["term_" + term] = (all_data[term[0] + "s"] - 1024.)*(all_data[term[1] + "s"] - 1024.)/ (1024**2.)
    else:
        assert sys.argv[4][0].upper() == "S"
        grid_size = int(sys.argv[4][1:])
        mid_point = int(np.floor(grid_size/2.))
        print("grid_size", grid_size, "mid_point", mid_point)

        x01 = all_data["xs"]/2048. # Scaled between 0 and 1
        y01 = all_data["ys"]/2048. # Scaled between 0 and 1

        x_nodes = np.linspace(0, 1., grid_size)
        terms = []

        plt.figure(figsize = (3*grid_size, 3*grid_size))
        for i in tqdm.trange(grid_size):
            for j in range(grid_size):
                if (i == mid_point) and (j == mid_point):
                    pass
                else:
                    coeffs = np.zeros([grid_size]*2, dtype=np.float64)
                    coeffs[i, j] = 1.
                    coeffs -= np.mean(coeffs)

                    ifn = RectBivariateSpline(x_nodes, x_nodes, coeffs, kx = 2, ky = 2)
                    term_name = "S" + str(i) + str(j)

                    all_data["term_" + term_name] = ifn(x01, y01, grid=False)
                    terms.append(term_name)
                    if min_flux > 0:
                        plt.subplot(grid_size, grid_size, i*grid_size + j + 1)
                        plt.scatter(x01, y01, c = all_data["term_" + term_name])
        plt.savefig("spline_coeffs.pdf", bbox_inches = 'tight')
        plt.close()
    
    return all_data, terms

scatter_size = 3
all_data, terms = read_data(phot_file = sys.argv[1],
                            the_filt = sys.argv[2], # E.g., F090W
                            min_flux = float(sys.argv[3]) # E.g., 10000
                            )


prefix = sys.argv[1].split(".")[0] + "_" + sys.argv[2] + "_minflux=" + sys.argv[3] + "_terms=" + "+".join(terms)

ministart = np.concatenate((
    np.ones(all_data["n_detectors"]-1, dtype=np.float64),
    np.zeros(all_data["n_detectors"]*len(terms), dtype=np.float64),
    [np.median(all_data["phots"][np.where(all_data["star_inds"] == i)]) for i in range(all_data["n_stars"])]
    ))

miniscale = np.ones(len(ministart), dtype=np.float64)

residfn(ministart, None, save_plot = prefix + "_resid_initial.pdf")

P, NA, NA = miniLM_new(ministart = ministart,
                       miniscale = miniscale,
                       passdata = None,
                       residfn = residfn, verbose = True, return_Cmat = False)#, maxiter = 20)


parsed = parseP(P)


residfn(P, None, save_plot = prefix + "_resid_final.pdf")


all_data, NA = read_data(phot_file = sys.argv[1],
                         the_filt = sys.argv[2], # E.g., F090W
                         min_flux = -1 # E.g., 10000
                         )

assert NA == terms

rel_sens_term = parsed["rel_sens"][all_data["detector_inds"]]

for term in terms:
    rel_sens_term += (parsed["rel_sens_" + term][all_data["detector_inds"]])*all_data["term_" + term]


with open(prefix + "_fit_params.json", "w") as f:
    json.dump([parsed, rel_sens_term], f, indent=2, cls=NumpyEncoder)
