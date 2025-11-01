import glob
import subprocess


fls = glob.glob("j*_cal.fits")

pointings = [(item.split("_")[0], item.split("_")[1], item.split("_")[3]) for item in fls]

print(pointings)

for pointing in set(pointings):
    these_fls = glob.glob(pointing[0] + "_" + pointing[1] + "*" + pointing[2] + "*_cal.fits")
    assert len(these_fls) < 17

    outputfl = "_".join(pointing) + "_stacked.fits"

    
    f = open("tmp.sh", 'w')
    f.write("""#!/bin/bash
#SBATCH --job-name=resamp
#SBATCH --partition=shared,kill-shared
#SBATCH --time=0-01:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G # Memory per node my job requires
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
        print(subprocess.getoutput("sbatch tmp.sh"))
