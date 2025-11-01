import glob
import subprocess

print(subprocess.getoutput("rm -fv photo*merge*"))

fls =  glob.glob("photo_subset_WD_j*--*.txt")

prefixes = ["_".join(item.split("_")[:6]) for item in fls]

unique_prefixes = list(set(prefixes))

print("unique_prefixes", unique_prefixes)

for unique_prefix in unique_prefixes:
    these_fls = []
    for i in range(len(fls)):
        if prefixes[i] == unique_prefix:
            these_fls.append(fls[i])

    print("unique_prefix", unique_prefix)
    assert len(these_fls) > 0
    print(len(these_fls))
    
    print(subprocess.getoutput("cat " + " ".join(these_fls) + " > " + unique_prefix + "_merge.txt"))
    
print(subprocess.getoutput("tar -cvzf photo.tar.gz photo*merge*"))
