import os

target = ["report.txt", "xyz_blessed.fmt"]

for f in os.listdir('.'):
    if os.path.isdir(f):
        print("Cleaning \"" + f + "\" ...")
        os.chdir(f)
        for ff in os.listdir('.'):
            if ff in target:
                os.remove(ff)
        os.chdir("..")
print("Finished.")
