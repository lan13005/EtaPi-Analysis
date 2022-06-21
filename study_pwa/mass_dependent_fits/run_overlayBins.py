#!/usr/bin/python3


import os
import multiprocessing
import shutil
import subprocess

workingDir=os.getcwd()

#ts=["010016", "016021", "021027", "027034", "034042", "042051", "051061", "061072", "072085", "085100"]
ts=["010020", "0200325", "0325050", "050075", "075100"]


doAccCorr="false"
Ls="S_D_pD"
#Ls="S_D"

waves=[
    ############### TMD WITH D-PRIME
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
    ############### TMD NO D-PRIME
#    "S0+-_S0++_D1--_D0+-_D1+-_D0++_D1++_D2++",
#    "S0+-;S0++", # individual S waves
#    "D1--;D0+-;D1+-;D0++;D1++;D2++", # individual D waves  
#    "S0+-_S0++", #sum S waves
#    "D1--_D0+-_D1+-_D0++_D1++_D2++", #sum D waves
#    "S0++_D0++_D1++_D2++", # all + ref
#    "S0+-_D1--_D0+-_D1+-" # all - ref
    ############## KMATRIX 
    #"S0+-_S0++_D0+-_D0++_D2+-_D2++",
    #"S0+-;S0++", # individual S waves
    #"D0+-;D0++;D2+-;D2++", # individual D waves  
    #"S0+-_S0++", #sum S waves
    #"D0+-_D0++_D2+-_D2++", #sum D waves
    #"S0++_D0++_D2++", # all + ref
    #"S0+-_D0+-_D2+-" # all - ref
]
waves=";".join(waves)

def draw(t):
    if shutil.which("etapi_plotter") is None:
        print("HEY DUMMY! etapi_plotter program does not exist. You probably forgot to source your environment...")
        exit()

    os.chdir(t)

    files=os.listdir()
    fitFiles=[afile for afile in files if "etapi_result" in afile]
    if "etapi_result.fit" not in fitFiles:
        fitFile=[int(afile.split(".")[0].split("etapi_result")[1]) for afile in fitFiles]
        fitFile="etapi_result"+str(max(fitFile))+".fit" 
    else:
        fitFile="etapi_result.fit"

    folder=workingDir+"/"+t
    cmd="python3 ../overlayBins.py 1 '"+waves+"' '"+fitFile+"' '"+workingDir+"' '"+Ls+"' '"+doAccCorr+"' '"+folder+"'"
    print(cmd)
    os.system(cmd)
    os.chdir("..")

    return 0

with multiprocessing.Pool(len(ts)) as p:
    ### CHANGING THE WAVESET TO LOOP OVER
    groupVec=waves.split(";")
    groupsWithQuotes=['"_'+tmp+'"' if tmp!="" else '""' for tmp in groupVec]
    repStr="vector<string> groups="
    vecStr="{"
    vecStr+=",".join(groupsWithQuotes)
    vecStr+="}"
    sedArgs=["sed","-i",'s@'+repStr+'.*;@'+repStr+vecStr+';@g',workingDir+"/overlayBins.C"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()

    p.map(draw,ts)
