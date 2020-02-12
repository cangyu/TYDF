import os

fluent_target = ["report.txt", "blessed.msh"]
nmf_target = ["report.txt", "map_blessed.nmf"]
p3d_target = ["report.txt", "xyz_blessed.fmt"]

for f in os.listdir('.'):

    # Loop over all cases
    if os.path.isdir(f):
        print("Cleaning case \"" + f + "\" ...")
        os.chdir(f)
        
        # FLUENT MESH FILE
        FLUENT_DIR = os.path.join('.', "FLUENT")
        if os.path.isdir(FLUENT_DIR):
            print("  FLUENT MESH ...")
            os.chdir(FLUENT_DIR)
            for ff in os.listdir('.'):
                if ff in fluent_target:
                    os.remove(ff)
            os.chdir("..")
        
        # NEUTRAL MAP FILE     
        NMF_DIR = os.path.join('.', "NMF")
        if os.path.isdir(NMF_DIR):
            print("  NEUTRAL MAP FILE ...")
            os.chdir(NMF_DIR)
            for ff in os.listdir('.'):
                if ff in nmf_target:
                    os.remove(ff)
            os.chdir("..")

        # PLOT3D GRID FILE
        P3D_DIR = os.path.join('.', "PLOT3D")
        if os.path.isdir(P3D_DIR):
            print("  PLOT3D GRID ...")
            os.chdir(P3D_DIR)
            for ff in os.listdir('.'):
                if ff in p3d_target:
                    os.remove(ff)
            os.chdir("..")
        
        print("  DONE!")
        os.chdir("..")

print("Finished!")
