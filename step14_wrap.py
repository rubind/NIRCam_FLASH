import subprocess
import sys
import time

n_jobs_max = int(sys.argv[1])
fl_to_read = sys.argv[2]
n_jobs = 2000

assert fl_to_read.count("photo_flattened_linear.txt") == 1

subprocess.getoutput("rm -f candidate_plots/*")


f = open(fl_to_read, 'r')
lines = f.read().split('\n')
f.close()

cand_filt = []

pwd = subprocess.getoutput("pwd")

# jw02729001001_02103_00001_nrca1_uncallin.fits:0 42765 84.67433880042816 -69.09451901518868 times 59732.93585596396 59732.93610450562 59732.93635304729 59732.93660158895 59732.93685013062 59732.937098672286 59732.937347213956 59732.93759575562 short_filt F090W 1865 132 -0.005733687781277149 0.2988610672299714 short_phot: 749.7492432194513 912.9939744813416 681.9511397456146 689.5893147069232 666.8361807539247 874.548396550577 644.2872762646092 816.0915522607518 short_RMS: 0.06280033091727942 0.06211021119469391 0.05700267527881943 0.06429440017079321 0.05922172344897474 0.06697926647705094 0.095320370807359 0.0722433165175593 short_uncs: 77.58967946246801 78.44982395939222 76.47519112005148 76.37263119780641 75.89305199957796 78.69846949137442 75.75512706003889 78.00942422084056 F335M 892 44 -0.22829237574048822 -0.7936367540164886 long_phot: 665.4048843228056 1151.9336344521175 926.1793600613605 871.3833163387774 864.9636226656806 655.0733978949836 1023.0692452535848 885.3975311723932 long_RMS: 0.06890916639979111 0.07188587719375107 0.06908721017662826 0.05774687492442679 0.07040118566675428 0.050663296744133134 0.05751647425000646 0.06560656739665044 long_uncs: 119.12184173444639 123.12117978845818 120.97324039775218 122.51270653510066 120.34274404886268 120.50478655524762 123.15475505756675 122.67703698117494


for line in lines:
    parsed = line.split(None)
    if parsed.count("short_filt") == 1:
        ind = parsed.index("short_filt")
        cand = parsed[1]
        short_filt = parsed[ind + 1]

        cand_filt.append((cand, short_filt))

print(len(cand_filt), "read")
cand_filt = list(set(cand_filt))

print(len(cand_filt), "unique")


for i in range(n_jobs):
    print("i", i, n_jobs)

    while int(subprocess.getoutput("squeue | grep drubin | wc").split(None)[0]) > n_jobs_max:
        time.sleep(200)
        
    f = open("tmp.sh", 'w')
    f.write("""#!/bin/bash
#SBATCH --job-name=phot
#SBATCH --partition=shared,kill-shared
#SBATCH --time=1-12:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G # Memory per node my job requires
#SBATCH --error=example-%A.err # %A - filled with jobid, where to write the stderr
#SBATCH --output=example-%A.out # %A - filled with jobid, wher to write the stdout
source ~/.bash_profile

cd """ + pwd + '\n')

    for this_cand_filt in cand_filt[i::n_jobs]:
        f.write("python /home/drubin/NIRCam_ramp/step14_find_best_candidates.py %s %s %s\n" % (this_cand_filt[0], sys.argv[2], this_cand_filt[1]))
    f.close()

    print(subprocess.getoutput("sbatch tmp.sh"))
