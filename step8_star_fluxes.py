import numpy as np
import matplotlib.pyplot as plt
import tqdm
import sys

# jw02609001001_02101_00001_nrca1_uncallin.fits:0 112725 11.516891413359724 42.1028667173646 times 59810.65731655358 59810.65756509524 59810.65781363691 59810.658062178576 short_filt F150W 1692 84 -0.2545770336301516 0.3069138938737631 short_phot: 384.7877272883814 254.77103534903657 301.19988419486225 326.7866035828571 short_RMS: 0.20643282206613464 0.3521488953995088 0.24753677149282702 0.21870731338693064 short_uncs: 71.99318204740538 70.68577165384343 69.80280239704348 72.1128288893627 F277W 806 20 1.4808303990829403 1.332920846475241 long_phot: -19.469660006994427 -56.0165567174156 27.148463982767908 8.792403239409223 long_RMS: 2.058382546919319 1.0809718482030204 9.692690537157398 5.207217927878371 long_uncs: 54.025734576355006 57.28714790484584 56.50117230452 54.53941578494179


all_data_by_star = {}

f = open("photo_flattened_linear.txt", 'r')
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
                
        long_filt = parsed[long_phot_start - 5]

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
            all_data_by_star[star_ind][filt + "_unc"] = np.median(all_data_by_star[star_ind][filt + "_unc"])/np.sqrt(len(all_data_by_star[star_ind][filt + "_unc"]))


            
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
