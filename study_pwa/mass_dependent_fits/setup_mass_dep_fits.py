#!/usr/bin/python3

import subprocess
import os
import random
import sys
import argparse

parser = argparse.ArgumentParser(description='Make a config file',formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('fname', type=str, help='config file location')
parser.add_argument('seed', type=int, help='seed to randomize the fits with')
parser.add_argument('pcwsBins', type=int, help='number of piecewise mass bins')
parser.add_argument('pcwsMassMin', type=float, help='piecewise mass minimum')
parser.add_argument('pcwsMassMax', type=float, help='piecewise mass maximum')
parser.add_argument('tmin', type=float, help='t miniimum')
parser.add_argument('tmax', type=float, help='t maximum')

args = parser.parse_args()

fileName=args.fname
seed=args.seed
pcwsBins=args.pcwsBins
pcwsMassMin=args.pcwsMassMin
pcwsMassMax=args.pcwsMassMax
tmin=args.tmin
tmax=args.tmax

if seed!=-1:
    random.seed(seed)

#baseDir="/d/grid17/ln16/dselector_v3/phase1_selected/"
#baseDir="/scratch-fits/"
baseDir="rootFiles/"


def replaceStr(search,replace,fileName):
    print("replace str: "+replace)
    sedArgs=["sed","-i",'s@'+search+'@'+replace+'@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()


filePrefix=fileName.split(".")[0].split("/")[-1]
fileAffix=fileName.split(".")[1]
newFileName=filePrefix+"-copy."+fileAffix
print("copying "+fileName+" to "+newFileName)
os.system("cp "+fileName+" "+newFileName)


t="010020"
m="104180" #"104180"
extraTag="_vh" 
for pol in ["000","045","090","135"]:
    baseLoc=baseDir+"t"+t+"_m"+m+extraTag+"/"
    if not os.path.exists(baseLoc):
        raise ValueError("YOU ARE REQUESTING FOR A FOLDER THAT DOES NOT EXIST! FIX IN SETUP FIT SCRIPT")

    search="DATAFILE_"+pol
    fileLoc="pol"+pol+"_t"+t+"_m"+m+extraTag+"_DTOT_selected_data_flat.root"
    replace=baseLoc+fileLoc
    replaceStr(search,replace,newFileName)

    search="BKGNDFILE_"+pol
    fileLoc="pol"+pol+"_t"+t+"_m"+m+extraTag+"_DTOT_selected_bkgnd_flat.root"
    replace=baseLoc+fileLoc
    replaceStr(search,replace,newFileName)

    search="ACCMCFILE_"+pol
    fileLoc="polALL_t"+t+"_m"+m+extraTag+"_FTOT_selected_acc_flat.root"
    #fileLoc="pol"+pol+"_t"+t+"_m"+m+"_FTOT_selected_acc_flat.root"
    replace=baseLoc+fileLoc
    replaceStr(search,replace,newFileName)

    search="GENMCFILE_"+pol
    fileLoc="polALL_t"+t+"_m"+m+extraTag+"_FTOT_gen_data_flat.root"
    #fileLoc="pol"+pol+"_t"+t+"_m"+m+"_FTOT_gen_data_flat.root"
    replace=baseLoc+fileLoc
    replaceStr(search,replace,newFileName)
    

waves=[
        "D0+-", "D0++", "D1+-", "D1++", "D2++", "D1--",
        "pD0+-", "pD0++", "pD1+-", "pD1++", "pD2++", "pD1--",
        ] # TMD waveset
#waves=["D2++","D2+-","D0++","D0+-"] # KMATRIX waveset
refs=["Negative","Positive"]
parts=["Re","Im"]

reaction="LOOPREAC"
#reaction="EtaPi0_000"

def reinitWave(wave,anchor):
    for j,ref in enumerate(refs): 
        for i,part in enumerate(parts):
            refpart=ref+part
            #scale=0 if wave[0]=="p" else 100 # start with null hypothesis where a2p does not exist
            scale=100
            if i==0:
                rsample=random.uniform(-1*scale,scale)
                isample=random.uniform(-1*scale,scale)
            search="initialize "+reaction+"::"+refpart+"::"+wave
            if anchor:
                replace=search+f" cartesian {rsample:0.5f} 0 real"
            else:
                replace=search+f" cartesian {rsample:0.5f} {isample:0.5f}"
            replaceStr(search+".*",replace,newFileName)

print("\n------------------------------------------------\n")
print("reintializing production amplitudes")
print("------------------------------------------------\n")
for wave in waves:
    reinitWave(wave,False)


print("\n------------------------------------------------\n")
print("intializing piecewise production parameters")
print("------------------------------------------------\n")
with open(newFileName) as cfg:
    lines=[line.rstrip() for line in cfg.readlines() if "ROOTDataReaderFilter" in line]
    for line in lines:
        accReplace=f" Mpi0eta {pcwsMassMin} {pcwsMassMax}" 
        genReplace=f" Mpi0eta_thrown {pcwsMassMin} {pcwsMassMax} mandelstam_t_thrown {tmin} {tmax}"
        #accReplace=accReplace+genReplace if "accmc" in line else accReplace
        replace=line+accReplace if any([ftype in line for ftype in ["accmc","data","bkgnd"]]) else line+genReplace
        replaceStr(line,replace,newFileName)


print("\n------------------------------------------------\n")
print("intializing piecewise production parameters")
print("------------------------------------------------\n")
realBin=4
pcwsMin=-200
pcwsMax=200
searchStrs=['PIECEWISE PARAMETER DEFINITIONS', 'PIECEWISE PARSCAN DEFINITIONS', 'PIECEWISE AMPLITUDE DEFINITIONS']
with open(newFileName) as cfg:
    lines=cfg.readlines()
    lineNums=[i+1 for i,line in enumerate(lines) if any([searchStr in line for searchStr in searchStrs])] 

parLines=''
for ref in ["Neg","Pos"]:
    for binNum in range(pcwsBins):
        for part in ["Re","Im"]:
            sample=random.uniform(pcwsMin,pcwsMax)
            if binNum==realBin and part=="Im":
                sample=0
            parLines+=f'parameter pcwsBin_{binNum}{part}{ref} {sample:0.5f}\n'
    if ref=="Neg":
        parLines+="\n"

parScanLines=""
parScanLines+=f'define uplimit {pcwsMax}\n'
parScanLines+=f'define lowlimit {pcwsMin}\n'
for ref in ["Neg","Pos"]:
    for binNum in range(pcwsBins):
        for part in ["Re","Im"]:
            if binNum==realBin and part=="Im":
                parScanLines+="#"
            parScanLines+=f'parRange pcwsBin_{binNum}{part}{ref} lowlimit uplimit\n'
    if ref=="Neg":
        parScanLines+="\n"

ampLines=""
for ref,ref2,ref3 in zip(["Neg","Pos"],["Negative","Positive"],["-","+"]):
    for ampPart in ["Re","Im"]:
        ampLines+=f'amplitude {reaction}::{ref2}{ampPart}::S0+{ref3} Piecewise {pcwsMassMin} {pcwsMassMax} {pcwsBins} 23 {ref} ReIm '
        for binNum in range(pcwsBins):
            for part in ["Re","Im"]:
                ampLines+=f'[pcwsBin_{binNum}{part}{ref}] '
        ampLines=ampLines[:-1]+"\n" # remove the extra space in the end and include a newline

nNewLines=0
for lineNum, newLines in zip(lineNums, [parLines, parScanLines, ampLines]):
    lines.insert(lineNum+nNewLines, newLines)
    nNewLines+=1 # the number of new lines added is not dependent on \n characters. instead the entire string counts as 1 line

with open(newFileName, "w") as f:
    lines = "".join(lines)
    f.write(lines)



print("\n------------------------------------------------\n")
print("COMMENT OUT WAVES WE ARE NOT USING")
print("------------------------------------------------\n")
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









            
