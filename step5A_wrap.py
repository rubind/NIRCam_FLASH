import glob
import subprocess


fls = glob.glob("*_uncal.fits")

print(fls)



f = open("tmp.sh", 'w')
f.write("""#!/bin/bash
#SBATCH --job-name=phot
#SBATCH --partition=shared,kill-shared
#SBATCH --time=0-01:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G # Memory per node my job requires
#SBATCH --error=example-%A.err # %A - filled with jobid, where to write the stderr
#SBATCH --output=example-%A.out # %A - filled with jobid, wher to write the stdout
source ~/.bash_profile

source activate jwst
pip install jwst
""")

chips_done = []
for fl in fls:
    chip = fl.split("_")[3]
    assert fl.split("_")[4] == "uncal.fits"

    if chips_done.count(chip) == 0:
        chips_done.append(chip)
        
        f.write("python step5_nonlin.py " + fl + '\n')
f.close()
assert len(chips_done) == 10

print(subprocess.getoutput("sbatch tmp.sh"))
