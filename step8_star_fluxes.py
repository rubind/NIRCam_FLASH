import numpy as np
import matplotlib.pyplot as plt
import tqdm
import sys

# jw02729001001_02103_00001_nrca1_uncallin.fits:0  42060  84.67204340672123  -69.09417751846111  times  59732.93585596396  59732.93610450562  59732.93635304729  59732.93660158895  59732.93685013062  59732.937098672286  59732.937347213956  59732.93759575562  short_filt  F090W  1784  71  short_phot:  99038.25159382843  98475.24968160014  98357.0322967695  99019.88003624743  116098.24143723339  101435.43218628655  102204.82377358645  107005.21612752581  short_RMS:  0.012010876400377507  0.012324904433704367  0.01647370475676274  0.011793001168715388  1.0  1.0  1.0  1.0  short_uncs:  436.4979471730728  435.67193584858677  434.65219278655957  436.74201257497333  -1.0  -1.0  -1.0  -1.0  F335M  852  14  long_phot:  21085.859856858573  20887.207250437947  21414.73431012046  21200.739580474255  20943.169486532046  20992.79387136198  20708.128127690543  21047.12900772877  long_RMS:  0.025800818270554284  0.012528602957625243  0.025509751343246127  0.014756844018666037  0.017642565100633976  0.014920407171486834  0.014773035924090362  0.020059470543729627  long_uncs:  208.87083743312962  207.1498396007937  209.84605016659322  209.0289244866579  207.28554633951163  207.61539309549167  206.21911971920196  206.74326388838986

all_data_by_star = {}

f = open("photo_flattened.txt", 'r')
lines = f.read().split('\n')
f.close()

unique_filts = []


for line in tqdm.tqdm(lines):
    if line.count(".fits"):
        parsed = line.split(None)

        short_phot_start = parsed.index("short_phot:")
        short_phot_end = parsed.index("short_RMS:")
        short_RMS_end = parsed.index("short_uncs:")

        n_points = short_phot_end - short_phot_start - 1
        
        short_phot_vals = []
        short_unc_vals = []
        
        for i, j, k in zip(range(short_phot_start + 1, short_phot_end), range(short_phot_end + 1, short_RMS_end), range(short_RMS_end + 1, short_RMS_end + 1 + n_points)):
            if float(parsed[j]) < 0.2:
                short_phot_vals.append(float(parsed[i]))
                short_unc_vals.append(float(parsed[k]))
                
        short_filt = parsed[parsed.index("short_filt") + 1]
        
        long_phot_start = parsed.index("long_phot:")
        long_phot_end = parsed.index("long_RMS:")
        long_RMS_end = parsed.index("long_uncs:")
        long_phot_vals = []
        long_unc_vals = []
        
        for i, j, k in zip(range(long_phot_start + 1, long_phot_end), range(long_phot_end + 1, long_RMS_end), range(long_RMS_end + 1, long_RMS_end + 1 + n_points)):
            if float(parsed[j]) < 0.2:
                long_phot_vals.append(float(parsed[i]))
                long_unc_vals.append(float(parsed[k]))
                
        long_filt = parsed[long_phot_start - 3]

        star_ind = int(parsed[1])

        unique_filts.append(short_filt)
        unique_filts.append(long_filt)

        if (star_ind in all_data_by_star) and (short_filt in all_data_by_star[star_ind]):
            all_data_by_star[star_ind][short_filt].extend(short_phot_vals)
            all_data_by_star[star_ind][short_filt + "_unc"].extend(short_unc_vals)
        else:
            if star_ind not in all_data_by_star:
                all_data_by_star[star_ind] = {}
            all_data_by_star[star_ind][short_filt] = short_phot_vals
            all_data_by_star[star_ind][short_filt + "_unc"] = short_unc_vals

        if (star_ind in all_data_by_star) and (long_filt in all_data_by_star[star_ind]):
            all_data_by_star[star_ind][long_filt].extend(long_phot_vals)
            all_data_by_star[star_ind][long_filt + "_unc"].extend(long_unc_vals)
        else:
            if star_ind not in all_data_by_star:
                all_data_by_star[star_ind] = {}
            all_data_by_star[star_ind][long_filt] = long_phot_vals
            all_data_by_star[star_ind][long_filt + "_unc"] = long_unc_vals

        all_data_by_star[star_ind]["RA"] = float(parsed[2])
        all_data_by_star[star_ind]["Dec"] = float(parsed[3])


        
plt.figure(figsize = (16, 12))
unique_filts = list(set(unique_filts))
unique_filts.sort()

for star_ind in tqdm.tqdm(all_data_by_star):
    for filt in unique_filts:
        if filt in all_data_by_star[star_ind]:
            assert len(all_data_by_star[star_ind][filt]) == len(all_data_by_star[star_ind][filt + "_unc"])

            all_data_by_star[star_ind][filt + "_count"] = len(all_data_by_star[star_ind][filt])
            all_data_by_star[star_ind][filt] = np.median(all_data_by_star[star_ind][filt])
            all_data_by_star[star_ind][filt + "_unc"] = np.median(all_data_by_star[star_ind][filt + "_unc"])


            
    #if "F090W" in all_data_by_star[star_ind] and "F200W" in all_data_by_star[star_ind]:
    #    plt.plot(-2.5*np.log10(all_data_by_star[star_ind]["F090W"]) - -2.5*np.log10(all_data_by_star[star_ind]["F200W"]),
    #             -2.5*np.log10(all_data_by_star[star_ind]["F200W"]), '.', color = 'b', alpha = 0.05)
#plt.savefig("color-mag.png", bbox_inches = 'tight')
#plt.close()

f = open("star_fluxes.txt", 'w')
f.write("#ID RA Dec ")
for unique_filt in unique_filts:
    f.write(unique_filt + " " + unique_filt + "_count " + unique_filt + "_unc ")
f.write('\n')

for star_ind in all_data_by_star:
    to_write = [star_ind, all_data_by_star[star_ind]["RA"], all_data_by_star[star_ind]["Dec"]]

    for filt in unique_filts:
        if filt in all_data_by_star[star_ind]:
            to_write.append(all_data_by_star[star_ind][filt])
            to_write.append(all_data_by_star[star_ind][filt + "_count"])
            to_write.append(all_data_by_star[star_ind][filt + "_unc"])
        else:
            to_write.append(-1)
            to_write.append(-1)
            to_write.append(-1)
    to_write = [str(item) for item in to_write]
    f.write(" ".join(to_write) + '\n')
f.close()
    
