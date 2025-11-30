import glob
import subprocess
from astropy.io import fits
import tqdm


fls = glob.glob("j*_tweakreg.fits")

#fls = [item for item in fls if item.count("long") == 0]

filts = []
for fl in tqdm.tqdm(fls):
    f = fits.open(fl)
    filts.append(f[0].header["FILTER"])
    f.close()



for filt in set(filts):
    outputfl = filt + "_stacked.fits"

    these_fls = []
    for i in range(len(fls)):
        if filts[i] == filt:
            these_fls.append(fls[i])


    tmpfl = "tmp_" + filt + ".sh"
    f = open(tmpfl, 'w')
    f.write("""#!/bin/bash
#SBATCH --job-name=resamp
#SBATCH --partition=shared,kill-shared
#SBATCH --time=0-02:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G # Memory per node my job requires
#SBATCH --error=example-%A.err # %A - filled with jobid, where to write the stderr
#SBATCH --output=example-%A.out # %A - filled with jobid, wher to write the stdout
source ~/.bash_profile

source activate jwst
pip install jwst

python ~/NIRCam_ramp/resamp.py """ + outputfl + " " + " ".join(these_fls))
    f.close()

    print(outputfl)

    if len(glob.glob(outputfl[:-5] + "*")) == 1:
        print("Already did ", outputfl)
    else:
        print(subprocess.getoutput("sbatch " + tmpfl))
