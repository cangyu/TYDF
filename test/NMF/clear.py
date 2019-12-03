import os

target = ["report.txt", "map_blessed.nmf"]

for f in os.listdir('.'):
    if os.path.isdir(f):
        print("Cleaning \"" + f + "\" ...")
        os.chdir(f)
        for ff in os.listdir('.'):
            if ff in target:
                os.remove(ff)
        os.chdir("..")
print("Finished.")
