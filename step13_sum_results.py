import matplotlib.pyplot as plt
import numpy as np
from FileRead import readcol, writecol
import sys

scale_by = float(sys.argv[1]) # E.g., 6000 for 30 hours, not 1 hour, all PHAT, not 1%, all M31, not half

for fl in sys.argv[2:]:
    f = open(fl, 'r')
    lines = f.read().split('\n')
    f.close()


    all_counts = {}

    for line in lines:

        if line.count("filt_name") == 1:
            current_filt = line.split(None)[1]

        if line.count("log10_mass") == 1:
            current_log10_mass = line.split(None)[-1]

            cur_key = (current_log10_mass, current_filt)
            if cur_key not in all_counts:
                all_counts[cur_key] = 0.000001
                print("setting", cur_key, "to zero")



        if line.count("N_exp") == 1:
            if line.count("Milky Way") == 1:
                all_counts[cur_key] += float(line.split("N_exp=")[1])*scale_by



    save_x = []
    save_y = []
    
    for key in all_counts.keys():
        print(key, all_counts[key])

        plt.subplot(2,1,1)
        plt.plot(10**(float(key[0])), all_counts[key], '.', color = 'b'*key.count("F150W") + 'r'*key.count("F277W") + 'k'*key.count("r") + 'b'*key.count("H"))
        plt.xscale('log')

        plt.subplot(2,1,2)
        plt.plot(10**(float(key[0])), 3./all_counts[key], '.', color = 'b'*key.count("F150W") + 'r'*key.count("F277W") + 'k'*key.count("r") + 'b'*key.count("H"))
        plt.xscale('log')
        plt.yscale('log')

        if key.count("F150W"):
            save_x.append(10**(float(key[0])))
            save_y.append(3./all_counts[key])

    writecol("constraint_sum_F150W_" + str(scale_by) + "x.txt", [save_x, save_y])
        
    plt.subplot(2,1,1)
    plt.ylim(0, plt.ylim()[1])
plt.show()
    
