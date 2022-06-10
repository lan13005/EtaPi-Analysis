#!/usr/bin/python


import os
import multiprocessing
import shutil

#ts=["010016", "016021", "021027", "027034", "034042", "042051", "051061", "061072", "072085", "085100"]
ts=["010020", "0200325", "0325050", "050075", "075100"]
#ts=["075100"]


Ls="S_D_pD"
#Ls="S_D"

def draw(t):
    if shutil.which("etapi_plotter") is None:
        print("HEY DUMMY! etapi_plotter program does not exist. You probably forgot to source your environment...")
        exit()

    os.chdir(t)
    waves=[
        "S0+-_S0++_D1--_D0+-_D1+-_D0++_D1++_D2++_pD1--_pD0+-_pD1+-_pD0++_pD1++_pD2++",
        "S0+-;S0++", # individual S waves
        "D1--;D0+-;D1+-;D0++;D1++;D2++", # individual D waves  
        "D1--_pD1--;D0+-_pD0+-;D1+-_pD1+-;D0++_pD0++;D1++_pD1++;D2++_pD2++", # Merge a2/a2prime with same M and reflectivity
        "pD1--;pD0+-;pD1+-;pD0++;pD1++;pD2++", # individual D prime
        "S0+-_S0++", #sum S waves
        "D1--_D0+-_D1+-_D0++_D1++_D2++", #sum D waves
        "pD1--_pD0+-_pD1+-_pD0++_pD1++_pD2++", #sum D prime waves
        "D1--_D0+-_D1+-_D0++_D1++_D2++_pD1--_pD0+-_pD1+-_pD0++_pD1++_pD2++", # sum ALL D waves
        "S0++_D0++_D1++_D2++_pD0++_pD1++_pD2++", # all + ref
        "S0+-_D1--_D0+-_D1+-_pD1--_pD0+-_pD1+-" # all - ref
#        ############### NO D-PRIME
#        "S0+-_S0++_D1--_D0+-_D1+-_D0++_D1++_D2++",
#        "S0+-;S0++", # individual S waves
#        "D1--;D0+-;D1+-;D0++;D1++;D2++", # individual D waves  
#        "S0+-_S0++", #sum S waves
#        "D1--_D0+-_D1+-_D0++_D1++_D2++", #sum D waves
#        "S0++_D0++_D1++_D2++", # all + ref
#        "S0+-_D1--_D0+-_D1+-" # all - ref
#        ############## REDUCED WAVESET
#        "S0+-_S0++_D2++",
#        "S0+-;S0++;D2++",
    ]
    waves=";".join(waves)

    files=os.listdir()
    fitFile=[int(afile.split(".")[0].split("etapi0_SD_TMD_piecewise_update")[1]) for afile in files if "etapi0_SD_TMD_piecewise_update" in afile]
    fitFile="etapi0_SD_TMD_piecewise_update"+str(max(fitFile))+".fit" 
    #fitFile="etapi0_SD_TMD_piecewise_update.fit" 

    cmd="python3 ../overlayBins.py 2 '"+waves+"' '"+fitFile+"' '..' "+Ls
    print(cmd)
    os.system(cmd)
    os.chdir("..")

    return 0

with multiprocessing.Pool(len(ts)) as p:
    p.map(draw,ts)
