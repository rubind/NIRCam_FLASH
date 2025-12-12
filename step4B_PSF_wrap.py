import glob
import subprocess
from astropy.io import fits
import sys
import tqdm



fls = glob.glob("j*_tweakreg.fits")

print(fls)


filts = []

for fl in tqdm.tqdm(fls):
    f = fits.open(fl)
    filts.append(f[0].header["FILTER"])
    f.close()
    

chip_filters = []

for i in tqdm.trange(len(fls)):
    chip_filters.append(fls[i].split("_")[3] + "_" + filts[i])


assert len(chip_filters) == len(fls)

for chip_filter in set(chip_filters):
    print("chip_filter", chip_filter)

    these_ims = []
    for i in range(len(fls)):
        if chip_filters[i] == chip_filter:
            these_ims.append(fls[i])
            
    
    f = open("tmp.sh", 'w')
    f.write("""#!/bin/bash
#SBATCH --job-name=phot
#SBATCH --partition=shared
#SBATCH --time=1-12:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G # Memory per node my job requires
#SBATCH --error=example-%A.err # %A - filled with jobid, where to write the stderr
#SBATCH --output=example-%A.out # %A - filled with jobid, wher to write the stdout
source ~/.bash_profile

python ~/NIRCam_ramp/PSF_builder_mine.py """ + sys.argv[1] + " 0 " + chip_filter + " " + " ".join(these_ims))
    f.close()

    print(subprocess.getoutput("sbatch tmp.sh"))
