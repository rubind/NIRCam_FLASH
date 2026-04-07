import sys
import glob
import numpy as np
from astropy.io import fits
from DavidsNM import miniNM_new, miniLM_new, save_img, save_patches
from FileRead import readcol
from scipy.interpolate import RectBivariateSpline
import astropy.wcs as wcs
import tqdm
import subprocess
from astropy.time import Time




def weighted_rms(x, w, subtract_mean=True):
    """
    Compute the weighted root-mean-square of x.

    Parameters
    ----------
    x : array_like
        Data values.
    w : array_like, optional
        Weights. If None, equal weights are assumed.
    subtract_mean : bool, default=False
        If True, compute RMS relative to the weighted mean (like stdev).

    Returns
    -------
    wrms : float
        Weighted RMS value.
    """
    x = np.asarray(x, dtype=float)
    if w is None:
        w = np.ones_like(x)
    else:
        w = np.asarray(w, dtype=float)

    if subtract_mean:
        mu = np.nansum(x*w)/np.nansum(w)
        x = x - mu

    return np.sqrt(np.nansum(w*x**2)/np.nansum(w))



def modelfn(patch, n_samps, P, psf_FN, data_cube):
    """Parameter vector P is x0, y0, ampl0, ampl1, ..."""
    xs_1d = np.arange(patch, dtype=np.float64)
    xs_1d -= np.mean(xs_1d)

    model = []
    for i in range(n_samps):
        model.append(psf_FN(xs_1d - P[0], xs_1d - P[1])*P[2+i])
        model[-1] += np.nanmedian(data_cube[i] - model[-1])
        
    return np.array(model)

def residfn(P, passdata):
    data_cube, psf_FN = passdata[0]

    model = modelfn(patch = half_patch*2 + 1,
                    n_samps = len(data_cube),
                    P = P,
                    psf_FN = psf_FN, data_cube = data_cube)

    resid = data_cube - model
    resid = resid.flatten()
    inds = np.where(1 - np.isnan(resid))
    return resid[inds]


def chi2fn(P, passdata):
    if np.abs(P[0]) > 1.5:
        return 1e100
    if np.abs(P[1]) > 1.5:
        return 1e100

    
    resid = residfn(P, passdata)
    return np.dot(resid, resid)

def do_phot(data_cube, psf_FN, read_noise, save_result = ""):
    median_data_cube = np.array([np.nanmedian(data_cube, axis = 0)])
    assert median_data_cube.shape == (1, half_patch*2 + 1, half_patch*2 + 1), str(median_data_cube.shape)

    median_P = [0., 0., 500.]
    for i in range(2):
        median_P, NA, NA = miniLM_new(ministart = median_P, miniscale = [0., 0., median_P[2]/10.],
                                      residfn = residfn, passdata = [median_data_cube, psf_FN])

        median_P, NA, NA = miniNM_new(ministart = median_P, miniscale = [0.5, 0.5, 0.0],
                                      chi2fn = chi2fn, passdata = [median_data_cube, psf_FN], compute_Cmat = False)

    print("median_P", median_P)

    P = np.concatenate((median_P[:2], [median_P[2]]*len(data_cube)))
    
    for i in range(2):
        P, NA, NA = miniLM_new(ministart = P,
                               miniscale = [0., 0.] + [median_P[2]/5.]*len(data_cube),
                               residfn = residfn, passdata = [data_cube, psf_FN])

        scale_in_flux = median_P[2]/5.
        scale_in_flux = max(scale_in_flux, -2*min(P[2:]))
        
        P, NA, NA = miniNM_new(ministart = P, miniscale = [0., 0.] + [scale_in_flux]*len(data_cube),
                               chi2fn = chi2fn, passdata = [data_cube, psf_FN], compute_Cmat = False)
        
        print("iteration", i, P)

    model = modelfn(patch = half_patch*2 + 1,
                    n_samps = len(data_cube),
                    P = P,
                    psf_FN = psf_FN, data_cube = data_cube)

    RMSs = []
    uncs = []
    xs_1d = np.arange(half_patch*2 + 1, dtype=np.float64)
    xs_1d -= np.mean(xs_1d)

    for i in range(len(data_cube)):
        tmp_PSF = psf_FN(xs_1d - P[0], xs_1d - P[1])
        tmp_PSF = tmp_PSF[np.where(1 - np.isnan(data_cube[i]))]
        tmp_PSF_sum = np.sum(tmp_PSF)
        
        if tmp_PSF_sum > 0.5:
            RMSs.append(weighted_rms(data_cube[i] - model[i], w = psf_FN(xs_1d - P[0], xs_1d - P[1]))
                        /(model[i].max())
                        )

            data_cube_variance = np.abs(data_cube[i]) + read_noise**2.

            tmp_PSF = psf_FN(xs_1d - P[0], xs_1d - P[1])
            
            
            uncs.append(
                np.sqrt(
                np.nansum(   tmp_PSF**2. * data_cube_variance   ) /
                (np.nansum(   tmp_PSF**2.  )**2.)
                )
                )
            
        else:
            RMSs.append(1.)
            uncs.append(-1.)
            
        
    if len(save_result) > 0:
        model = modelfn(patch = half_patch*2 + 1,
                        n_samps = len(data_cube),
                        P = P,
                        psf_FN = psf_FN, data_cube = data_cube)
        save_patches(np.concatenate((   data_cube, model, data_cube - model, psf_FN(xs_1d - P[0], xs_1d - P[1])   )), save_result)
        
    return P[2:], RMSs, model, uncs, P[0], P[1]


def read_PSF(fl):
    f = fits.open(fl)
    dat = f[0].data
    f.close()

    assert fl.count("10x_") == 1

    xs1d = np.arange(len(dat), dtype=np.float64)*0.1
    xs1d -= np.mean(xs1d)

    return RectBivariateSpline(xs1d, xs1d, dat, kx=2, ky=2)

    

    
    

half_patch = 5

print("python step6_do_phot.py WD_jw02559001001_02101_nrca2.txt 0 0 1 2 3 ... the first 0 is for writing out a fits file")

if __name__ == "__main__":
    input_fl = sys.argv[1]
    short_wave_fl = sys.argv[2]
    #[source_ids, x_sw, y_sw, NA, NA, ras, decs] = readcol(input_fl, 'f,ffff,ff')
    [x_pix, y_pix, ras, decs] = readcol(input_fl, 'ff,ff') # E.g., WD_jw02729001003_02105_nrca3.txt

    write_fits = int(sys.argv[3])
    use_model_PSF = int(sys.argv[4])

    prefix = "photo_subset_" + sys.argv[1].split(".")[0] + "_" + short_wave_fl.split(".")[0] + "_" + use_model_PSF*"modelPSF" + (1 - use_model_PSF)*("empiricalPSF")

    if len(sys.argv) > 5:
        f_phot = open(prefix + "_" + sys.argv[5] + "--" + sys.argv[-1] + ".txt", 'w')
        i_range = sys.argv[5:]
    else:
        f_phot = open(prefix + ".txt", 'w')
        i_range = range(len(ras))





    whoami = subprocess.getoutput("whoami")


    psf_FNs = {}
    if use_model_PSF:
        for fl in glob.glob("/home/" + whoami + "/NIRCam_ramp/F*10x_PSF.fits"):
            filt_name = fl.split("/")[-1].split("_")[0]
            print("Reading ", filt_name, fl)
            psf_FNs[filt_name] = read_PSF(fl)
    else:
        # WD_jw02729001001_02103_nrca1.txt
        # PSF_10x_jw02729001004_02105_nrca3_cal.fits

        for fl in glob.glob("PSF_10x_" + short_wave_fl.split("_")[3] + "_*fits"):
            print("fl", fl)
            psf_FNs[fl.split("_")[-1].split(".fits")[0]] = read_PSF(fl)

        for fl in glob.glob("PSF_10x_" + short_wave_fl.split("_")[3][:-1] + "long_*fits"):
            print("fl", fl)
            psf_FNs[fl.split("_")[-1].split(".fits")[0]] = read_PSF(fl)

        for key in psf_FNs:
            print(key)
            print(psf_FNs[key](0., 0.))


    short_data_cubes = []
    long_data_cubes = []

    short_read_noise = []
    long_read_noise = []

    short_xys = []
    long_xys = []
    fls = []
    short_filts = []
    long_filts = []
    mjds = []


    tmp_fls = [short_wave_fl]

    for fl in tqdm.tqdm(tmp_fls):
        print(fl)

        f = fits.open(fl)
        print(f.info())
        n_integrations = len(f["SCI"].data)
        print("n_integrations", n_integrations)

        f.close()




        for int_ind in range(n_integrations):
            fls.append(fl + ":" + str(int_ind))

            f = fits.open(fl)
            print(f.info())
            mjds.append(f["INT_TIMES"].data["int_start_MJD_UTC"][int_ind] + np.arange(f[0].header["NGROUPS"])*f[0].header["TGROUP"]/86400.)
            dat = f["SCI"].data[int_ind]*1.
            sat_mask = (f["GROUPDQ"].data[int_ind] & 2) != 0    # shape: (nint, ngroup, ny, nx)
            dat[sat_mask] = np.nan

            short_data_cubes.append(dat)
            short_filts.append(f[0].header["FILTER"])
            short_read_noise.append(f["RN"].data*1.)
            f.close()

            f = fits.open(fl.replace("_uncallin", "_tweakreg"))
            print(f.info())
            w = wcs.WCS(f["SCI"].header, f).celestial
            short_xys.append(
                np.array(np.around(w.all_world2pix(np.array([ras, decs]).T, 1, quiet = True)), dtype=np.int32)
            )
            f.close()

            print(short_xys[-1].shape)

            #f = open("ds9.reg", 'w')
            #f.write("""# Region file format: DS9 version 4.1
            #global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
            #image
            #""")
            #for i in range(len(short_xys[-1])):
            #    f.write("circle(%f,%f,5)\n" % (short_xys[-1][i][0], short_xys[-1][i][1]))
            #f.close()

            f = fits.open(fl.split("_uncallin")[0][:-1] + "long_uncallin.fits")
            print(f.info())
            dat = f["SCI"].data[int_ind]*1.
            sat_mask = (f["GROUPDQ"].data[int_ind] & 2) != 0    # shape: (nint, ngroup, ny, nx)
            dat[sat_mask] = np.nan
            long_data_cubes.append(dat)
            long_filts.append(f[0].header["FILTER"])
            long_read_noise.append(f["RN"].data*1.)

            f.close()


            f = fits.open(fl.split("_uncallin")[0][:-1] + "long_tweakreg.fits")
            w = wcs.WCS(f["SCI"].header, f).celestial
            long_xys.append(
                np.array(np.around(w.all_world2pix(np.array([ras, decs]).T, 1, quiet = True)), dtype=np.int32)
            )
            f.close()




    f_phot.write("#WD_fl  star_ind shortx shorty short_phot short_RMS longx longy long_phot long_RMS\n")

    all_ims = []

    for istr in tqdm.tqdm(i_range):
        i = int(istr)


        for j in range(len(short_xys)): # Each image
            short_xy = short_xys[j][i]
            long_xy = long_xys[j][i]

            #print("short_xy", short_xy)
            #print("long_xy", long_xy)


            if short_xy[0] > half_patch + 1 and short_xy[0] < 2048 - half_patch and short_xy[1] > half_patch + 1 and short_xy[1] < 2048 - half_patch:
                if long_xy[0] > half_patch + 1 and long_xy[0] < 2048 - half_patch and long_xy[1] > half_patch + 1 and long_xy[1] < 2048 - half_patch:

                    print(short_xy, long_xy)
                    print("short_data_cubes[j]", short_data_cubes[j].shape)

                    short_cutout = short_data_cubes[j][:, short_xy[1] - half_patch - 1: short_xy[1] + half_patch,
                                                       short_xy[0] - half_patch - 1: short_xy[0] + half_patch]

                    long_cutout = long_data_cubes[j][:, long_xy[1] - half_patch - 1: long_xy[1] + half_patch,
                                                       long_xy[0] - half_patch - 1: long_xy[0] + half_patch]

                    short_cutout = short_cutout[1:] - short_cutout[:-1]
                    long_cutout = long_cutout[1:] - long_cutout[:-1]

                    short_read_noise_cutout = short_read_noise[j][short_xy[1] - half_patch - 1: short_xy[1] + half_patch,
                                                                  short_xy[0] - half_patch - 1: short_xy[0] + half_patch]

                    long_read_noise_cutout = long_read_noise[j][long_xy[1] - half_patch - 1: long_xy[1] + half_patch,
                                                                long_xy[0] - half_patch - 1: long_xy[0] + half_patch]


                    short_PSF_key = short_filts[j]
                    long_PSF_key = long_filts[j]


                    short_phot, short_RMSs, short_model, short_uncs, cent_x_short, cent_y_short = do_phot(short_cutout, psf_FNs[short_PSF_key], read_noise = short_read_noise_cutout)
                    #, save_result = "short_" + fls[j].split(".")[0] + ".fits")
                    long_phot, long_RMSs, long_model, long_uncs, cent_x_long, cent_y_long = do_phot(long_cutout, psf_FNs[long_PSF_key], read_noise = long_read_noise_cutout)
                    #, save_result = "long_" + fls[j].split(".")[0] + ".fits")

                    print("short_cutout", short_cutout.shape)
                    print("long_cutout", long_cutout.shape)
                    print("short_model", short_model.shape)

                    print("short_RMSs", short_RMSs, "long_RMSs", long_RMSs)

                    this_im = []
                    for t in range(len(short_cutout)):
                        this_im.append(np.concatenate((short_cutout[t], short_model[t], short_cutout[t] - short_model[t],
                                                       long_cutout[t], long_model[t], long_cutout[t] - long_model[t])))
                    this_im = np.concatenate(tuple([item.T for item in this_im]))

                    to_write = [fls[j], i, ras[i], decs[i], "times"] + list(mjds[j][1:]) + [
                        "short_filt", short_filts[j], short_xy[0], short_xy[1], cent_x_short, cent_y_short, "short_phot:"] + list(short_phot) + [
                            "short_RMS:"] + list(short_RMSs) + [
                                "short_uncs:"] + list(short_uncs) + [long_filts[j], long_xy[0], long_xy[1], cent_x_long, cent_y_long, "long_phot:"] + list(long_phot) + [
                                    "long_RMS:"] + list(long_RMSs) + [
                                        "long_uncs:"] + list(long_uncs)

                    to_write = [str(item) for item in to_write]

                    f_phot.write("  ".join(to_write) + '\n')
                    all_ims.append(this_im)
                else:
                    print("Out of range!")
            else:
                print("Out of range!")

    all_ims = np.concatenate((tuple([item for item in all_ims])))
    if write_fits:
        save_img(all_ims, "cutout_" + prefix + "_" +  istr + ".fits")

    f_phot.close()

