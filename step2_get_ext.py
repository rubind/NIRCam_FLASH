import glob
from astropy.io import fits
import tqdm

for fl in tqdm.tqdm(glob.glob("F*resample*fits")):
    # jw02559001003_03101_nrcblong
    output = "_".join(fl.split("_")[:3])

    for ext in ["SCI", "ERR"]:
        f = fits.open(fl)
        print(f.info())
        for i in range(1, len(f))[::-1]:
            if f[i].name != ext:
                del f[i]

        f.writeto(output + "_" + ext + ".fits")
        f.close()


