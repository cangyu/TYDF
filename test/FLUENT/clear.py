import os

target = ["report.txt", "blessed.msh"]

for f in os.listdir('.'):
    if os.path.isdir(f):
        print("Cleaning \"" + f + "\" ...")
        os.chdir(f)
        for ff in os.listdir('.'):
            if ff in target:
                os.remove(ff)
        os.chdir("..")
print("Finished.")
