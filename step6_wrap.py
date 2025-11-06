import glob
import subprocess
import sys

n_per_job = 50

print("Ready to remove old photometry")
input()

print(subprocess.getoutput("rm -fv photo_subset_WD_*"))


fls = glob.glob("WD*txt")

print(fls)

pwd = subprocess.getoutput("pwd")


for fl in fls:
    f = open(fl, 'r')
    lines = f.read().split('\n')
    f.close()

    line_count = 0
    for line in lines:
        parsed = line.split(None)
        try:
            float(parsed[1])
            line_count += 1
        except:
            pass
    n_jobs = int(line_count/float(n_per_job))

    print("n_jobs", n_jobs)

    star_inds = [str(ind) for ind in range(line_count)]
    
    for i in range(n_jobs):
        f = open("tmp.sh", 'w')
        f.write("""#!/bin/bash
#SBATCH --job-name=phot
#SBATCH --partition=shared
#SBATCH --time=0-20:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G # Memory per node my job requires
#SBATCH --error=example-%A.err # %A - filled with jobid, where to write the stderr
#SBATCH --output=example-%A.out # %A - filled with jobid, wher to write the stdout
source ~/.bash_profile



cd """ + pwd + '\n')


        f.write("python ~/NIRCam_ramp/step6_do_phot.py " + fl + " 0 " + " ".join(star_inds[i::n_jobs]) + '\n')

        f.close()
        
        print(subprocess.getoutput("sbatch tmp.sh"))

