#!/usr/bin/python3

import os
import re
import sys
import subprocess
import argparse
import numpy as np

'''
The purpose of this code is to reconfigure a .fit results file and seed the 
  functionall equivalent cfg file at the bottom of the fit file with the best fit results 
  From here we can answer questions like how the modification of the data, i.e. by including a Delta cut, influences the amplitudes
'''

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
    ###########
    ## Get the complex production amplitudes
    ###########
    waves_used={}
    for k,v in parameterMap.items():
        if k.startswith("EtaPi0"):
            wave, part=k.split("::")[-1].split("_")
            value = float(v)
            if wave in waves_used: 
                if part=="re":
                    waves_used[wave]=value+1j*waves_used[wave].imag
                else:
                    waves_used[wave]=waves_used[wave].real+1j*value
            else:
                waves_used[wave]=1+1j*1 # arbitrarily set to 1+1j
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
                parameterMap[k]=f'{waves_used[wave].real}'
            else:
                parameterMap[k]=f'{waves_used[wave].imag}'
    #printMap(parameterMap,"NEW")

    return parameterMap

def reconfigure(fname,otag,anchor=-1):
    '''
    This script creates a new configuration file based off the converged fit values in a .fit file
    The intention is to not change the actual values that are fitted for but you can modify the input data
        so that we can see the systemtatic changes things 
        1. Gauss constraint mean/mins
        2. Datasets
        3. The RootDataReaderFilter scheme

    anchor=-1 means we do no modify the anchor bin
    '''

    with open(fname) as f:
        lines=f.readlines()
    
    cfgFirstLine=[i for i,v in enumerate(lines) if "## FIT CONFIGURATION ###" in v][0]
    fitResultLine=[i for i,v in enumerate(lines) if "+++ Parameter Values and Errors +++" in v][0]
    cfg=lines[cfgFirstLine+1:]
    fitresults=lines[fitResultLine:cfgFirstLine+1]
    fitresults=[re.sub(r" +","",line).rstrip().lstrip() for line in fitresults]
    # make sure the line is a parameterName and value that is separated by a tab
    fitresults=[line.split("\t") for line in fitresults if (len(line.split("\t"))==2)&(~line[0].isdigit())] 
    fitMap={k:v for k,v in fitresults}

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
        pos_anchor_angle=np.angle(float(fitMap[f'pcwsBin_{anchor}RePos'])+1j*float(fitMap[f'pcwsBin_{anchor}ImPos']))
        neg_anchor_angle=np.angle(float(fitMap[f'pcwsBin_{anchor}ReNeg'])+1j*float(fitMap[f'pcwsBin_{anchor}ImNeg']))
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
            # if the line has 4 fields then the 4th field is likely to be "fixed" so we dont need to change anything anyways
            if len(fields)==3:
                keyword, name, value = line.split() 
                cfg[i]=f"{keyword} {name} {fitMap[name]}\n"
            if len(fields)==6:
                keyword, name, value, constraint, constraintVar1, constraintVar2 = line.split()
                cfg[i]=f"{keyword} {name} {fitMap[name]} {constraint} {constraintVar1} {constraintVar2}\n"
        # Sometimes it will have 6 where the 6th field is "fixed", so we wont have to change it anyways
        if fields[0]=="initialize":
            if len(fields)==5: 
                keyword, amp, coordinate, val1, val2 = line.split() 
                cfg[i]=f'{keyword} {amp} {coordinate}  {fitMap[amp+"_re"]} {fitMap[amp+"_im"]}\n'
            if len(fields)==6: 
                keyword, amp, coordinate, val1, val2, fixed = line.split() 
                cfg[i]=f'{keyword} {amp} {coordinate}  {fitMap[amp+"_re"]} {fitMap[amp+"_im"]} fixed\n'

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('otag', type=str, nargs="?", default='vary',  help='output tag to append to file names')
    parser.add_argument('anchor', type=int, nargs="?", default=-1, help='piecewise anchor bin')
    args = parser.parse_args()
    otag=args.otag
    anchor=args.anchor

    ts=["010020","0200325","0325050","050075","075100"]
    #ts=["0200325","0325050","050075","075100"]
    cfgs=[
            # choose which config files you want to run over
            f"{t}/etapi_result_src.fit" 
            #f"{t}/etapi_result.fit" 
             for t in ts
            ]
    
    nprocesses=9
    niters=1
    for t,cfgFile in zip(ts,cfgs):
        newCfgLoc=reconfigure(cfgFile,otag,anchor)

        ofolder=cfgFile.split("/")
        ofolder=ofolder[len(ofolder)-2:len(ofolder)-1][0]
        ofolder+="_"+otag
    
        print("Starting fits")
        #cmd="mpirun -np "+str(nprocesses)+" fitMPI -c "+newCfgLoc+" -r "+str(niters)+" -m 1000000 -t 1.0 -x 1 -f 0.15" 
        cmd="mpirun -np "+str(nprocesses)+" fitMPI -c "+newCfgLoc+" -m 1000000 -t 1.0 -x 0" 
        pipeCmd=' > fit.log'
        print(cmd+pipeCmd)
        os.system(cmd+pipeCmd)

        print(f"moving output to {ofolder}")
        os.system(f"mkdir -p {ofolder}")
        os.system("mv -f etapi_result.fit "+ofolder+"/etapi_result_0.fit 2>/dev/null")
        os.system("mv -f "+newCfgLoc+" "+ofolder+" 2>/dev/null")
        os.system("mv -f *.log "+ofolder+" 2>/dev/null")
        os.system("mv -f *.ni "+ofolder+" 2>/dev/null")



