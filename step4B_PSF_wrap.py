import glob
import subprocess


fls = glob.glob("WD*txt")

print(fls)



for fl in fls:

    
    f = open("tmp.sh", 'w')
    f.write("""#!/bin/bash
#SBATCH --job-name=phot
#SBATCH --partition=shared
#SBATCH --time=0-12:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G # Memory per node my job requires
#SBATCH --error=example-%A.err # %A - filled with jobid, where to write the stderr
#SBATCH --output=example-%A.out # %A - filled with jobid, wher to write the stdout
source ~/.bash_profile

python ~/NIRCam_ramp/PSF_builder_mine.py """ + fl + " 0 ")
    f.close()

    print(subprocess.getoutput("sbatch tmp.sh"))
