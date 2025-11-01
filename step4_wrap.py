import glob
import subprocess


fls = glob.glob("j*csv")

print(fls)



for fl in fls:
    #F150W2_SCI = sys.argv[1]
    #F150W2_ERR = sys.argv[2]
    #F322W2_SCI = sys.argv[3]
    #F322W2_ERR = sys.argv[4]
    #OUT_ECSV   = sys.argv[5]

    
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

python ~/NIRCam_ramp/step4_find_WD.py """ + fl)
    f.close()

    print(subprocess.getoutput("sbatch tmp.sh"))
