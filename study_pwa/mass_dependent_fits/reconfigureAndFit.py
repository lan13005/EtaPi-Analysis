#!/usr/bin/python3

import os
import re
import sys
import subprocess
import argparse
import numpy as np
import time
import random
from rotateProductionCoeffs import rotateProductionCoeffs

'''
The purpose of this code is to reconfigure a .fit results file and seed the 
  functionally equivalent cfg file at the bottom of the fit file with the best fit results 
  From here we can answer questions like how the modification of the data, i.e. by including a Delta cut, influences the amplitudes
'''

def checkFits(floc):
    searchMinimumStr="bestMinimum"
    searchMinuitStatusStr="lastMinuitCommandStatus"
    searchErrMatrixStatus="eMatrixStatus"
    with open(floc) as infile:
        lines=infile.readlines()
        for line in lines:
            if searchMinimumStr in line:
                NLL=float(line.split(" ")[-1].split("\t")[1].rstrip().lstrip())
            if searchMinuitStatusStr in line:
                minuitStatus=float(line.split(" ")[-1].split("\t")[1].rstrip().lstrip())
            if searchErrMatrixStatus in line:
                eMatrixStatus=float(line.split(" ")[-1].split("\t")[1].rstrip().lstrip())
        convergenceStatus=NLL<0 and np.isfinite(NLL) and minuitStatus==0 and eMatrixStatus==3
    return convergenceStatus

def replace(cfg,old,new):
    '''
    Replace the value of varName to varValue in the file called fileName. Depending on the value type we can include quotes or not
    '''
    sedArgs=["sed","-i","s@"+old+"@"+new+"@g",cfg]
    print(" ".join(sedArgs))
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()

def printMap(map1,header):
    print(header)
    for wave, v in map1.items():
        print(f'{wave}: {v}')
        
def zeroComplex(x):
    approxZero=1e-13
    real = 0.0 if abs(x.real)<approxZero else x.real
    imag = 0.0 if abs(x.imag)<approxZero else x.imag
    return real+1j*imag

def updateProductionAmpsAndRotate(parameterMap,posAngle,negAngle):
    '''
    This function will rotate all production coefficients to anchor another piecewise bin
    '''
    ###########
    ## Get the complex production amplitudes
    ###########
    waves_used={}
    for k,v in parameterMap.items():
        if k.startswith("EtaPi0"):
            wave, part=k.split("::")[-1].split("_")
            value = float(v)
            if wave not in waves_used:
                waves_used[wave]=complex()
            if part=="re":
                waves_used[wave]=complex(value,waves_used[wave].imag)
            else:
                waves_used[wave]=complex(waves_used[wave].real,value)

    ###########
    ## Rotate
    ###########
    #printMap(waves_used,"OLD")
    for wave,value in waves_used.items():
        angle=posAngle if wave.endswith('+') else negAngle
        #print(f'{wave}: {np.angle(waves_used[wave])}')
        # anchoring a wave to be real = rotating to nearest real line which could be be 0 or pi 
        # np.angle outputs on [-pi,pi]
        # the actual calculation is kind of messy but it turns out that if we wish to rotate by the shortest amount to the real line
        #     then we can multiply by exp(-ia) and by an additional -1 if greater than absolutle value is greater than pi/2
        waves_used[wave]*=np.cos(-angle)+1j*np.sin(-angle) # cosine is even so we could technically use angle instead of -angle
        if abs(angle)>np.pi/2:
            waves_used[wave]*=-1
        waves_used[wave]=zeroComplex(waves_used[wave])
        #print(f' {wave}: {np.angle(waves_used[wave])}')
    #printMap(waves_used,"NEW")

    ###########
    ## Update
    ###########
    #printMap(parameterMap,"OLD")
    for k,v in parameterMap.items():
        if k.startswith("EtaPi0"):
            wave, part=k.split("::")[-1].split("_")
            if part=="re":
                parameterMap[k]=waves_used[wave].real
            else:
                parameterMap[k]=waves_used[wave].imag
    #printMap(parameterMap,"NEW")

    return parameterMap

def reconfigure(fname,otag,anchor=-1,perturbScale=0,usePolar=True):
    '''
    This script creates a new configuration file based off the converged fit values in a .fit file
    The intention is to not change the actual values that are fitted for but you can modify the input data
        so that we can see the systemtatic changes things 
        1. Gauss constraint mean/mins
        2. Datasets
        3. The RootDataReaderFilter scheme

    anchor=-1 means we do not modify the anchor bin
    perturbScale determines how much we should perturb the production amplitudes in the perturbMe list
        perturbScale=X will multiply the real and imaginary parts by (1+X) 
        By default the most intense production parameter will be perturbed by perturbScale
        "fixed" production parameters will not be perturbed
    usePolar: bool to decide whether we hould convert production parameters to polar coordinates
    '''

    perturbMe=[] # list of amplitudes to perturb by the perturbScale variable

    with open(fname) as f:
        lines=f.readlines()
    
    cfgFirstLine=[i for i,v in enumerate(lines) if "## FIT CONFIGURATION ###" in v][0]
    fitResultLine=[i for i,v in enumerate(lines) if "+++ Parameter Values and Errors +++" in v][0]
    cfg=lines[cfgFirstLine+1:]
    fitresults=lines[fitResultLine:cfgFirstLine+1]
    fitresults=[re.sub(r" +","",line).rstrip().lstrip() for line in fitresults]
    # make sure the line is a parameterName and value that is separated by a tab
    fitresults=[line.split("\t") for line in fitresults if (len(line.split("\t"))==2)&(~line[0].isdigit())] 
    fitMap={k:float(v) for k,v in fitresults}
    #print("FITMAP BELOW")
    #printMap(fitMap,'')
    #print("FITMAP ABOVE")


    ###### EXTRACT PARTIAL WAVE AMPLITUDES, INTENSITIES #######
    ###### Originally developed so we can extract the most intense D-wave, recall piecewise terms contains the scale of the S-wave production amp #####
    ###### We can randomly perturb the most intense D-wave ####
    productionAmplitudes={}
    # You can also choose the nth most intense wave also. Dont have to settle for the most intense
    #nth_most_intense=1 
    nth_most_intense=random.randint(1,3) # for the waveset syst scan you might actually be removing the nth most intense wave. We can randomize for flexibility
    for k,v in fitMap.items():
        if "::" in k:
            amp, part = k.split("::")[-1].split("_")
            if amp not in productionAmplitudes.keys():
                productionAmplitudes[amp]=complex()
            if part=="re":
                productionAmplitudes[amp]=complex(v,productionAmplitudes[amp].imag)
            else:
                productionAmplitudes[amp]=complex(productionAmplitudes[amp].real,v)
    for k,v in productionAmplitudes.items():
        productionAmplitudes[k]=abs(v)
    for ith_most_intense in range(nth_most_intense):
        mostIntenseWave=max(productionAmplitudes, key=productionAmplitudes.get)
        productionAmplitudes.pop(mostIntenseWave) # pop the most intense, can cause problems if we wish to use productionAmplitudes
    perturbMe.append(mostIntenseWave)
    print(f"Perturbing {nth_most_intense}th most intense prodution amplitude: {mostIntenseWave} in {fname}")
    if perturbScale!=0:
        print(f"Perturbing waves: {perturbMe} with scale: {perturbScale}")

    old_anchor={}
    if anchor!=-1:
        ######################
        ## 1. find old anchor
        ######################
        old_anchor={i:l for i,l in enumerate(cfg) if l.startswith('parameter pcwsBin') and 'fixed' in l}
        ######################
        ## 2. update map with rotated production amplitudes 
        ######################
        ## np.angle returns a value on [-pi,pi]
        pos_anchor_angle=np.angle(fitMap[f'pcwsBin_{anchor}RePos']+1j*fitMap[f'pcwsBin_{anchor}ImPos'])
        neg_anchor_angle=np.angle(fitMap[f'pcwsBin_{anchor}ReNeg']+1j*fitMap[f'pcwsBin_{anchor}ImNeg'])
        fitMap=updateProductionAmpsAndRotate(fitMap,pos_anchor_angle,neg_anchor_angle)
        #print("".join(cfg))
    
    ######################
    ## reinitialize cfg file with best fit parameters
    ######################
    for i, line in enumerate(cfg):
        line=line.rstrip().lstrip()
        fields=line.split(" ")
        if fields[0]=="parameter": 
            # parameter lines have >=3 fields respectively
            if len(fields)==3:
                keyword, name, value = line.split() 
                cfg[i]=f"{keyword} {name} {fitMap[name]}\n"
            elif len(fields)==4:
                # if the line has 4 fields then the 4th field is likely to be "fixed" so we dont need to change anything anyways
                pass
            elif len(fields)==6:
                keyword, name, value, constraint, constraintVar1, constraintVar2 = line.split()
                cfg[i]=f"{keyword} {name} {fitMap[name]} {constraint} {constraintVar1} {constraintVar2}\n"
        if fields[0]=="initialize": ## These are the production parameters
            if len(fields)==5 or (len(fields)==6 and fields[-1]=='real'): 
                keyword, amp, coordinate, val1, val2 = line.split() 
                perturbation=perturbScale if amp.split("::")[-1].split("_")[0] in perturbMe else 0
                ## Cant move rePart and imPart out! Will result in bug!
                rePart=fitMap[amp+"_re"] if amp+"_re" in fitMap.keys() else 0
                imPart=fitMap[amp+"_im"] if amp+"_im" in fitMap.keys() else 0
                cfg[i]=f'{keyword} {amp} {coordinate}  {rePart*(1+perturbation)} {imPart*(1+perturbation)}\n'

    if anchor!=1:
        ######################
        ## 3. swap anchors - fixed keyword
        ######################
        for i,l in old_anchor.items():
            cfg[i]=" ".join(l.split(" ")[:-1])+"\n" # unanchor old bin
        cfg=[l.rstrip()+" fixed\n" if l.startswith(f'parameter pcwsBin_{anchor}Im') else l for l in cfg]
    
    fout=fname.split(".")[0].split("/")[1]
    fout=fout+"_"+otag+".cfg"

    with open(fout,"w") as newcfg:
        newcfg.write("".join(cfg))

    return fout

#############################
## Previous implementation ##
#############################
#if __name__ == "__main__":
#    parser = argparse.ArgumentParser()
#    parser.add_argument('otag', type=str, nargs="?", default='vary',  help='output tag to append to file names')
#    parser.add_argument('anchor', type=int, nargs="?", default=-1, help='piecewise anchor bin')
#    parser.add_argument('--method', type=int, nargs="?", default=0, help='method for rotateProductionCoeffs function')
#    parser.add_argument('--ignore_waves', type=int, nargs="?", default=[], help='ignore_waves for rotateProductionCoeffs function')
#    parser.add_argument('--tbins', type=int, nargs="+")
#    args = parser.parse_args()

def reconfigureAndFit(otag='vary',anchor=-1,method=0,ignore_waves=["NOTHING_TO_MATCH"],waveMatch='',tbinsChosen=[]):
    ts=["010020","0200325","0325050","050075","075100"]
    startPerturbScales=[0.0, 0.0, 0.0, 0.0, 0.0]

    cfgs=[
            # choose which config files you want to run over
            f"{t}/etapi_result_src.fit" 
            #f"{t}/etapi_result.fit" 
             for t in ts
            ]
    nprocesses=9
    niters=1

    if tbinsChosen==[]: # empty list will actually be a full list. Why would we want to run over no t-bins?
        tbinsChosen=[0,1,2,3,4]
    ts=list(np.array(ts)[tbinsChosen])
    startPerturbScales=list(np.array(startPerturbScales)[tbinsChosen])
    cfgs=list(np.array(cfgs)[tbinsChosen])
    for t,cfgFile,startPerturbScale in zip(ts,cfgs,startPerturbScales):
        ofolder=cfgFile.split("/")
        ofolder=ofolder[len(ofolder)-2:len(ofolder)-1][0]
        ofolder+="_"+otag

        print("Starting fits")
        perturbScales=[startPerturbScale]
        for i in range(15):
            perturbScales.append(perturbScales[0]+(i+1)*0.02)
            perturbScales.append(perturbScales[0]-(i+1)*0.02)
        for perturbScale in perturbScales:
            newCfgLoc=reconfigure(cfgFile,otag,anchor=anchor,perturbScale=perturbScale)
            rotateProductionCoeffs(newCfgLoc,ignore_waves,waveMatch,method) # modifies the file to incorporate incorporate a phase parameter to Zlm amplitudes
            #cmd="mpirun -np "+str(nprocesses)+" fitMPI -c "+newCfgLoc+" -r "+str(niters)+" -m 1000000 -t 1.0 -x 1 -f 0.15 -H 1" 
            cmd="mpirun -np "+str(nprocesses)+" fitMPI -c "+newCfgLoc+" -m 1000000 -t 0.1 -x 1" 
            pipeCmd=' > fit.log'
            print(cmd+pipeCmd)
            if True: # can quickly turn off if we wish to not run the fit and just generate a file
                os.system(cmd+pipeCmd)
                fitConverged=checkFits("etapi_result.fit")
#                fitConverged=True
                subdir='' if fitConverged else f'failed_perturbScale_{perturbScale:0.2f}'

                print(f"moving output to {ofolder}/{subdir}")
                os.system(f"mkdir -p {ofolder}/{subdir}")
                os.system(f"mv -f etapi_result.fit {ofolder}/{subdir}/etapi_result_0.fit")# 2>/dev/null")
                os.system(f"mv -f {newCfgLoc} {ofolder}/{subdir}")# 2>/dev/null")
                os.system(f"mv -f *.log {ofolder}/{subdir}")# 2>/dev/null")
                os.system(f"mv -f *.ni {ofolder}/{subdir}")# 2>/dev/null")
                
                # For w.e. reason running the same config file twice can lead one to reach MACHINE ACCURACY and one to converge nicely...
                ### IDK what is going on under the hood but perhaps letting things reset is a good idea
                print("Sleeping for 10 seconds to hopefully let things reset completely...")
                time.sleep(10)

                if fitConverged:
                    break

    return 0



