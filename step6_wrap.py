import glob
import subprocess
import sys
import numpy as np
import time

def qcount():
    the_count = int(subprocess.getoutput("squeue | grep drubin | wc").split(None)[0])
    print("count", the_count)
    return the_count


use_model_PSF = int(sys.argv[1])
n_jobs_per_fl = int(sys.argv[2])
n_jobs_max = int(sys.argv[3])

redo_old_jobs = int(sys.argv[4])

if redo_old_jobs:
    print("Ready to remove old photometry")
    input()
    print(subprocess.getoutput("rm -fv photo_subset_*"))


wd_fl = glob.glob("WD*csv")[0]

f = open(wd_fl, 'r')
lines = f.read().split('\n')
f.close()

good_stars = 0

for line in lines[1:]:
    parsed = line.split(None)

    
    try:
        float(parsed[1])
        good_line = 1
    except:
        good_line = 0
        print("skipping", parsed)
        

    if good_line:
        good_stars += 1




fls = glob.glob("j*nrc??_uncallin.fits")

print(fls, len(fls))
print("star_inds", good_stars)

assert good_stars > 100

pwd = subprocess.getoutput("pwd")

star_inds = np.arange(good_stars)
star_inds = [str(item) for item in star_inds]


existing_phot_files = glob.glob("photo_subset*txt")

for fl in fls:
    for i in range(n_jobs_per_fl):
        f = open("tmp.sh", 'w')
        f.write("""#!/bin/bash
#SBATCH --job-name=phot
#SBATCH --partition=shared,kill-shared
#SBATCH --time=0-20:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G # Memory per node my job requires
#SBATCH --error=example-%A.err # %A - filled with jobid, where to write the stderr
#SBATCH --output=example-%A.out # %A - filled with jobid, wher to write the stdout
source ~/.bash_profile

cd """ + pwd + '\n')


        #output_name = "photo_subset_" + sys.argv[1].split(".")[0] + "_" + short_wave_fl.split(".")[0] + "_" + use_model_PSF*"modelPSF" + (1 - use_model_PSF)*("empiricalPSF") + "_" + sys.argv[5] + "--" + sys.argv[-1] + ".txt"

        sys_argv = [None, wd_fl, fl, 0, use_model_PSF, " ".join(star_inds[i::n_jobs_per_fl])]
        
        output_name = "photo_subset_" + sys_argv[1].split(".")[0] + "_" + sys_argv[2].split(".")[0] + "_" + sys_argv[4]*"modelPSF" + (1 - sys_argv[4])*("empiricalPSF") + "_" + sys_argv[5] + "--" + sys_argv[-1] + ".txt"

        print("output_name", output_name)
        
        if existing_phot_files.count(output_name) == 0:
            f.write("python ~/NIRCam_ramp/step6_do_phot.py " + " ".join([str(item) for item in sys_argv[1:]]) + '\n')

        f.close()

        while qcount() > n_jobs_max:
            print("Waiting...")
            time.sleep(60)
            
        print(subprocess.getoutput("sbatch tmp.sh"))

