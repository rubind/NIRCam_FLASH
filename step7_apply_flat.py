import glob
import numpy as np
import json
import sys
import tqdm


fit_results = {}
all_rel_sens_term = {}

for fl in glob.glob("photo_unflat_*terms=S00+S01+S02+S03+S10*json"):
    with open(fl) as f:
        parsed, rel_sens_term = json.load(f)
    filt = fl.split("_")[2]
    
    fit_results[filt] = parsed
    all_rel_sens_term[filt] = rel_sens_term

f = open(sys.argv[1], 'r')
lines = f.read().split('\n')
f.close()

new_lines = []

# jw02729001001_02103_00001_nrca1_uncallin.fits:0  42060  84.67204340672123  -69.09417751846111  times  59732.93585596396  59732.93610450562  59732.93635304729  59732.93660158895  59732.93685013062  59732.937098672286  59732.937347213956  59732.93759575562  short_filt  F090W  1784  71  short_phot:  99038.25159382843  98475.24968160014  98357.0322967695  99019.88003624743  116098.24143723339  101435.43218628655  102204.82377358645  107005.21612752581  short_RMS:  0.012010876400377507  0.012324904433704367  0.01647370475676274  0.011793001168715388  1.0  1.0  1.0  1.0  short_uncs:  436.4979471730728  435.67193584858677  434.65219278655957  436.74201257497333  -1.0  -1.0  -1.0  -1.0  F335M  852  14  long_phot:  21085.859856858573  20887.207250437947  21414.73431012046  21200.739580474255  20943.169486532046  20992.79387136198  20708.128127690543  21047.12900772877  long_RMS:  0.025800818270554284  0.012528602957625243  0.025509751343246127  0.014756844018666037  0.017642565100633976  0.014920407171486834  0.014773035924090362  0.020059470543729627  long_uncs:  208.87083743312962  207.1498396007937  209.84605016659322  209.0289244866579  207.28554633951163  207.61539309549167  206.21911971920196  206.74326388838986


line_count = dict()

for line in tqdm.tqdm(lines):
    if line.count("short_phot:") == 0:
        new_lines.append(line)
    else:
        parsed = line.split(None)
        short_start = parsed.index("short_phot:")
        short_end = parsed.index("short_RMS:")

        short_uncs_start = parsed.index("short_uncs:")
        short_uncs_end = parsed.index("short_uncs:")
        

        long_start = parsed.index("long_phot:")
        long_end = parsed.index("long_RMS:")

        long_uncs_start = parsed.index("long_uncs:")

        
        short_filt = parsed[short_start - 3]
        long_filt = parsed[long_start - 3]

        if short_filt not in line_count:
            line_count[short_filt] = 0

        if long_filt not in line_count:
            line_count[long_filt] = 0
        
        for i in range(short_start+1, short_end):
            parsed[i] = str(float(parsed[i])/all_rel_sens_term[short_filt][line_count[short_filt]])

        for i in range(short_uncs_start+1, short_uncs_end):
            parsed[i] = str(float(parsed[i])/all_rel_sens_term[short_filt][line_count[short_filt]])

            
        for i in range(long_start+1, long_end):
            parsed[i] = str(float(parsed[i])/all_rel_sens_term[long_filt][line_count[long_filt]])

        for i in range(long_uncs_start+1, len(parsed)):
            parsed[i] = str(float(parsed[i])/all_rel_sens_term[long_filt][line_count[long_filt]])

        new_lines.append(" ".join(parsed))

        line_count[short_filt] += 1
        line_count[long_filt] += 1
        
print(line_count)


for key in line_count:
    assert line_count[key] == len(all_rel_sens_term[key])
    



f = open("photo_flattened.txt", 'w')
f.write('\n'.join(new_lines))
f.close()
