import glob
import subprocess


fls = glob.glob("*_nrc??_SCI.fits")

print(fls)



#for fl in fls:
    #F150W2_SCI = sys.argv[1]
    #F150W2_ERR = sys.argv[2]
    #F322W2_SCI = sys.argv[3]
    #F322W2_ERR = sys.argv[4]
    #OUT_ECSV   = sys.argv[5]

    #long_sci = fl[:len("jw02559001002_02101_nrca")] + "long_SCI.fits"

    #print(fl, long_sci)


fl = "F090W_stacked_resample.fits_SCI.fits"
long_sci = "F200W_stacked_resample.fits_SCI.fits"

#fl = "F200W_stacked_resample.fits_SCI.fits"
#long_sci = "F444W_stacked_resample.fits_SCI.fits"

    
f = open("tmp.sh", 'w')
f.write("""#!/bin/bash
#SBATCH --job-name=phot
#SBATCH --partition=shared,kill-shared
#SBATCH --time=0-02:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G # Memory per node my job requires
#SBATCH --error=example-%A.err # %A - filled with jobid, where to write the stderr
#SBATCH --output=example-%A.out # %A - filled with jobid, wher to write the stdout
source ~/.bash_profile

source activate jwst
pip install jwst

python ~/NIRCam_ramp/step3_photometry.py """ + fl + " " + fl.replace("SCI", "ERR") + " " + long_sci + " " + long_sci.replace("SCI", "ERR") + " " + fl.split("_SCI")[0] + ".csv")
f.close()

print(subprocess.getoutput("sbatch tmp.sh"))
