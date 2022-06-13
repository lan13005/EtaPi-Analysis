#!/usr/bin/python3


import os
import multiprocessing
import shutil
import subprocess

workingDir=os.getcwd()

#ts=["010016", "016021", "021027", "027034", "034042", "042051", "051061", "061072", "072085", "085100"]
#ts=["010020", "0200325", "0325050", "050075", "075100"]
ts=["kmatrix_fit_results"]


#Ls="S_D_pD"
Ls="S_D"

def draw(t):
    if shutil.which("etapi_plotter") is None:
        print("HEY DUMMY! etapi_plotter program does not exist. You probably forgot to source your environment...")
        exit()

    os.chdir(t)
    waves=[
#        ############### WITH TMD D-PRIME
#        "S0+-_S0++_D1--_D0+-_D1+-_D0++_D1++_D2++_pD1--_pD0+-_pD1+-_pD0++_pD1++_pD2++",
#        "S0+-;S0++", # individual S waves
#        "D1--;D0+-;D1+-;D0++;D1++;D2++", # individual D waves  
#        "D1--_pD1--;D0+-_pD0+-;D1+-_pD1+-;D0++_pD0++;D1++_pD1++;D2++_pD2++", # Merge a2/a2prime with same M and reflectivity
#        "pD1--;pD0+-;pD1+-;pD0++;pD1++;pD2++", # individual D prime
#        "S0+-_S0++", #sum S waves
#        "D1--_D0+-_D1+-_D0++_D1++_D2++", #sum D waves
#        "pD1--_pD0+-_pD1+-_pD0++_pD1++_pD2++", #sum D prime waves
#        "D1--_D0+-_D1+-_D0++_D1++_D2++_pD1--_pD0+-_pD1+-_pD0++_pD1++_pD2++", # sum ALL D waves
#        "S0++_D0++_D1++_D2++_pD0++_pD1++_pD2++", # all + ref
#        "S0+-_D1--_D0+-_D1+-_pD1--_pD0+-_pD1+-" # all - ref
#        ############### NO TMD D-PRIME
#        "S0+-_S0++_D1--_D0+-_D1+-_D0++_D1++_D2++",
#        "S0+-;S0++", # individual S waves
#        "D1--;D0+-;D1+-;D0++;D1++;D2++", # individual D waves  
#        "S0+-_S0++", #sum S waves
#        "D1--_D0+-_D1+-_D0++_D1++_D2++", #sum D waves
#        "S0++_D0++_D1++_D2++", # all + ref
#        "S0+-_D1--_D0+-_D1+-" # all - ref
#        ############## KMATRIX 
        "S0+-_S0++_D0+-_D0++_D2+-_D2++",
        "S0+-;S0++", # individual S waves
        "D0+-;D0++;D2+-;D2++", # individual D waves  
        "S0+-_S0++", #sum S waves
        "D0+-_D0++_D2+-_D2++", #sum D waves
        "S0++_D0++_D2++", # all + ref
        "S0+-_D0+-_D2+-" # all - ref
    ]
    waves=";".join(waves)

    files=os.listdir()
    fitFiles=[afile for afile in files if "etapi0_SD_TMD_piecewise_update" in afile]
    if "etapi0_SD_TMD_piecewise_update.fit" not in fitFiles:
        fitFile=[int(afile.split(".")[0].split("etapi0_SD_TMD_piecewise_update")[1]) for afile in fitFiles]
        fitFile="etapi0_SD_TMD_piecewise_update"+str(max(fitFile))+".fit" 
    else:
        fitFile="etapi0_SD_TMD_piecewise_update.fit"

    cmd="sed -i 's/folder=\".*\"/folder=\"%s\"/' %s" % (workingDir+"/"+t,workingDir+"/overlayBins.C")
    subprocess.Popen(cmd.split(" "), stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()

    cmd="python3 ../overlayBins.py 2 '"+waves+"' '"+fitFile+"' '..' "+Ls
    print(cmd)
    os.system(cmd)
    os.chdir("..")

    return 0

with multiprocessing.Pool(len(ts)) as p:
    p.map(draw,ts)
