import glob
import subprocess
import sys

fls_per_job = int(sys.argv[1])


fls = glob.glob("*_uncal.fits")
done_fls = glob.glob("*_uncallin.fits")

print(fls)



while len(fls) > 0:
    #F150W2_SCI = sys.argv[1]
    #F150W2_ERR = sys.argv[2]
    #F322W2_SCI = sys.argv[3]
    #F322W2_ERR = sys.argv[4]
    #OUT_ECSV   = sys.argv[5]

    
    f = open("tmp.sh", 'w')
    f.write("""#!/bin/bash
#SBATCH --job-name=phot
#SBATCH --partition=shared,kill-shared
#SBATCH --time=0-02:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G # Memory per node my job requires
#SBATCH --error=example-%A.err # %A - filled with jobid, where to write the stderr
#SBATCH --output=example-%A.out # %A - filled with jobid, wher to write the stdout
source ~/.bash_profile

source activate jwst
pip install jwst

""")

    for i in range(fls_per_job):
        if len(fls) > 0:
            if done_fls.count(fls[-1].replace("_uncal.fits", "_uncallin.fits")) == 0:
                f.write("python ~/NIRCam_ramp/step5_nonlin.py " + fls[-1] + '\n')
            else:
                print("Already done")
            del fls[-1]
    f.close()

    print(subprocess.getoutput("sbatch tmp.sh"))
