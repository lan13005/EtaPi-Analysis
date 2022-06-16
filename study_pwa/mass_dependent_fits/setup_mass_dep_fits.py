#!/usr/bin/python

import subprocess
import os
import random


baseDir="/d/grid17/ln16/dselector_v3/"

def replaceStr(search,replace,fileName):
    print("replace str: "+replace)
    sedArgs=["sed","-i",'s@'+search+'@'+replace+'@g',fileName]
    #print("replacing: "+search+" with "+replace)
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()

fileName=baseDir+"study_pwa/mass_dependent_fits/config_files/etapi_hybrid.cfg"

filePrefix=fileName.split(".")[0]
fileAffix=fileName.split(".")[1]
newFileName=filePrefix+"-copy."+fileAffix
print("copying "+fileName+" to "+newFileName)
os.system("cp "+fileName+" "+newFileName)

for pol in ["000","045","090","135"]:
    baseLoc=baseDir+"phase1_selected/t010020_m104180/"
    if not os.path.exists(baseLoc):
        raise ValueError("YOU ARE REQUESTING FOR A FOLDER THAT DOES NOT EXIST! FIX IN SETUP FIT SCRIPT")

    search="DATAFILE_"+pol
    fileLoc="pol"+pol+"_t010020_m104180_DTOT_selected_data_flat.root"
    replace=baseLoc+fileLoc
    replaceStr(search,replace,newFileName)

    search="BKGNDFILE_"+pol
    fileLoc="pol"+pol+"_t010020_m104180_DTOT_selected_bkgnd_flat.root"
    replace=baseLoc+fileLoc
    replaceStr(search,replace,newFileName)

    search="ACCMCFILE_"+pol
    fileLoc="polALL_t010020_m104180_FTOT_selected_acc_flat.root"
    replace=baseLoc+fileLoc
    replaceStr(search,replace,newFileName)

    search="GENMCFILE_"+pol
    fileLoc="polALL_t010020_m104180_FTOT_gen_data_flat.root"
    replace=baseLoc+fileLoc
    replaceStr(search,replace,newFileName)
    

waves=[
        "D0+-", "D0++", "D1+-", "D1++", "D2++", "D1--",
        "pD0+-", "pD0++", "pD1+-", "pD1++", "pD2++", "pD1--",
        ] # TMD waveset
#waves=["D2++","pD2++"]
refs=["Negative","Positive"]
parts=["Re","Im"]

def reinitWave(wave,anchor):
    for j,ref in enumerate(refs): 
        for i,part in enumerate(parts):
            refpart=ref+part
            if i==0:
                rsample=random.uniform(-500,500)
                isample=random.uniform(-500,500)
            search="initialize LOOPREAC::"+refpart+"::"+wave
            if anchor:
                replace=search+" cartesian "+str(rsample)+" 0 real"
            else:
                replace=search+" cartesian "+str(rsample)+" "+str(isample)
            replaceStr(search+".*",replace,newFileName)

print("\n------------------------------------------------\n")
print("reintializing production amplitudes")
print("------------------------------------------------\n")
#reinitWave("S0++",True)
#reinitWave("S0+-",True)
for wave in waves:
    reinitWave(wave,False)

searchStr="parameter pcwsBin"
constrainedParMap={}
print("\n------------------------------------------------\n")
print("reintializing piecewise production parameters between [-100,100]")
print("------------------------------------------------\n")
with open(newFileName) as newFile:
    for line in newFile:
        if line.startswith(searchStr):
            line=line.rstrip()
            parType , parName, parVal = line.split(" ")
            #if "SKIP" in parName:
            if "Bin_4Im" in parName:
                sample=0.0
                replaceStr(line,parType+" "+parName+" "+str(sample)+" fixed",newFileName)
            else:
                sample=random.uniform(0,1000)
                replaceStr(line,parType+" "+parName+" "+str(sample),newFileName)
             

commentWavesNotInThisList=waves+["S0++","S0+-"]
listAllWaves=[]
excludeLines=[]
with open(newFileName) as newFile:
    for line in newFile:
        if line.startswith("initialize"):
            listAllWaves.append(line.split("::")[-1].split(" ")[0])
    excludeList=list(set(listAllWaves)-set(commentWavesNotInThisList))
    newFile.seek(0)
    for linenum,line in enumerate(newFile):
        if any(exclude in line for exclude in excludeList):
            excludeLines.append(linenum)
            os.system("sed -i '"+str(linenum+1)+"s/^/#/' "+newFileName) 









            
