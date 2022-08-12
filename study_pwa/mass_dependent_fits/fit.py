#!/usr/bin/python3

import os
import sys
import numpy as np
import subprocess
import uproot3 as uproot
import pandas as pd
import random
import time

def checkParLimits(fitFile):
    '''
    Check to see if the fit parameters are at the limits. If so, the fit did not 
    converge properly
    '''
    boundedPars={}
    with open(fitFile) as infile:
        for line in infile:
            if "parameter" in line and "bounded" in line:
                _, parName, parVal, _, parMin, parMax = line.split(" ")
                boundedPars[parName]=[float(parMin),float(parMax.rstrip())]
    with open(fitFile) as infile:
        parNotAtLimit=True
        infile.seek(0)
        for line in infile:
            for key in boundedPars.keys():
                if line.lstrip().startswith(key+"\t"):
                    parName, parVal=[ele.rstrip().lstrip() for ele in line.split("\t")]
                    parVal = float(parVal)
                    minVal, maxVal = boundedPars[parName]
                    shiftedRatio=(parVal-minVal)/(maxVal-minVal)
                    if shiftedRatio<percent/100 or shiftedRatio>(1-percent/100):
                        parNotAtLimit *= False
                        print(parName+" is within "+str(percent)+"% of limits! "+str(parVal)+" bounded on ["+str(minVal)+","+str(maxVal)+"]")
    return parNotAtLimit

def has_numbers(inputString):
    return any(char.isdigit() for char in inputString)

def checkFits(folder,returnGoodFilesInstead=False):
    '''
    Do not care about fit status, only care that the likelihood is reasonable (negative and finite)
    0 Not calculated at all
    ## Error matrix status?
    1 Diagonal approximation only, not accurate
    2 Full matrix, but forced positive-definite
    3 Full accurate covariance matrix (After MIGRAD, this is the indication of normal convergence.)

    ## Minuit status?
    *       -1 = undefined status
    *        0 = normal
    *        1 = blank command
    *        2 = unreadable command
    *        3 = unkown command
    *        4 = abnormal termination (e.g., MIGRAD not converged)
    '''
    searchMinimumStr="bestMinimum"
    searchMinuitStatusStr="lastMinuitCommandStatus"
    
    files=[folder+"/"+f for f in os.listdir(folder) if ".fit" in f and has_numbers(f)]
    convergenceStatuses=[]
    for fitFile in files:
        with open(fitFile) as infile:
            NLL=0
            for line in infile:
                if searchMinimumStr in line:
                    NLL=float(line.split(" ")[-1].split("\t")[1].rstrip().lstrip())
                if searchMinuitStatusStr in line:
                    minuitStatus=float(line.split(" ")[-1].split("\t")[1].rstrip().lstrip())
            if NLL==0:
                raise ValueError("NLL not found in "+fitFile+"! Terminating")
            #parNotAtLimit=checkParLimits(fitFile)
            parNotAtLimit=True # ignore if we wish
            convergenceStatus=NLL<0 and np.isfinite(NLL) and minuitStatus==0 and parNotAtLimit # Require a negative log likelihood and minuit status = 0
            convergenceStatuses.append(convergenceStatus)
    if not returnGoodFilesInstead:
        return sum(convergenceStatuses) # how much results files converged properly
    else:
        return np.array(files)[convergenceStatuses]


def spawnProcessChangeSetting(old,new):
    '''
    Replace the value of varName to varValue in the file called fileName. Depending on the value type we can include quotes or not
    '''
    sedArgs=["sed","-i","s@"+old+"@"+new+"@g","setup_mass_dep_fits.py"]
    print(" ".join(sedArgs))
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()

def checkCfgFileProperMassLimits(cfgLoc):
    searchKeyWord="Mpi0eta"
    with open(cfgLoc) as cfg:
        lines=cfg.readlines()
        lines=[line for line in lines if not line.startswith("#")]
        isLooping=len([line for line in lines if "loop" in line])>0
        if isLooping:
            rootFileLines=[line for line in lines if ".root" in line]
        piecewise=[line for line in lines if "Piecewise" in line]
    
    oneSetOfFiles=rootFileLines[0].rstrip().lstrip().split(" ")[2:]
    oneFile=oneSetOfFiles[0]
    df=uproot.open(oneFile)["kin"].arrays(["Mpi0eta"],outputtype=pd.DataFrame)
    minVal=df.Mpi0eta.min()
    maxVal=df.Mpi0eta.max()
    
    piecewise=piecewise[0].split(" ")
    minPW=float(piecewise[3])
    maxPW=float(piecewise[4])
    
    if abs(minVal-minPW)>0.1 or abs(maxVal-maxPW)>0.1:
        print("\n****************************")
        print("THERE IS A MISMATCH BETWEEN YOUR MASS RANGE IN YOUR ROOT FILE AND THE RANGE SPECIFIED IN YOUR PIECEWISE DEFINITION!")
        print("Min/Max value in root file: {}/{}".format(minVal,maxVal))
        print("Min/Max value in piecewise: {}/{}".format(minPW,maxPW))
        print("TERMINATING THE PROGRAM!")
        print("****************************\n")
        exit()

def replaceStr(search,replace,fileName):
    print("replace str: "+replace)
    sedArgs=["sed","-i",'s@'+search+'@'+replace+'@g',fileName]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()

def main(arg):
    ##########################################
    # Basic setup
    ##########################################
    nprocesses=9
    fitFileName="etapi_result.fit"
    percent=3.0 # parameters must not be within percent of the defined parameter limits
    niters=1
    workingDir=os.getcwd()

    # create a seed to sample another seed value that is input to setup_mass_dep_fits for reproducible series of fits
    # Set seed to -1 to not use a seed
    seed=1992 
    if seed!=-1:
        random.seed(seed)

    ts=["010020","0200325","0325050","050075","075100"]
    tmins=[0.1,0.2,0.325,0.5,0.75]
    tmaxs=[0.2,0.325,0.5,0.75,1.0]
    ms=["104156"]
    mmins=[1.04]
    mmaxs=[1.56]

    baseCfgFile="config_files/etapi_hybrid.cfg"

    runFits=False
    getSummary=False
    if arg==0:
        runFits=True
    elif arg==1:
        getSummary=True
    elif arg==2:
        runFits=True
        getSummary=True
    else:
        raise ValueError("Argument [0/1/2] to runFit, getSummary or both")
    
    if runFits:
        for j, tmin, tmax in zip(range(len(ts)),tmins,tmaxs):
            t=ts[j]
        
            os.system("rm -f "+fitFileName)
            os.system("rm -f fit*log")

            for ic,m,mmin,mmax in zip(range(len(ms)),ms,mmins,mmaxs):
                cfgFile=baseCfgFile.split("/")[-1].split(".")[0]+"-copy"
    
                os.system("mkdir -p overlayPlots")
        
                setup_seed=random.randint(0,99999999) if seed!=-1 else -1
                cmd="python setup_mass_dep_fits.py "+baseCfgFile+" "+str(setup_seed)
                print(cmd)
                os.system(cmd) # reinitialize
                # checkCfgFileProperMassLimits(cfgFile+".cfg") ## DO NOT CHECK ANYMORE SINCE WE USE FILTERING DATA READERS

                with open(cfgFile+".cfg") as cfg:
                    lines=[line.rstrip() for line in cfg.readlines() if "ROOTDataReaderFilter" in line]
                    lines=[line for line in lines if 'mc' in line]
                    for line in lines:
                        extraReplace=" Mpi0eta {:0.3f} {:0.3f}".format(mmin,mmax) if "accmc" in line else ""
                        replace=line+extraReplace+" Mpi0eta_thrown {:0.3f} {:0.3f} mandelstam_t_thrown {:0.3f} {:0.3f}".format(
                                mmin-0.1,mmax+0.1,tmin-0.1,tmax+0.1) 
                        replaceStr(line,replace,cfgFile+".cfg")

                print("Starting fits")
                cmd="mpirun -np "+str(nprocesses)+" fitMPI -c "+cfgFile+".cfg -r "+str(niters)+" -m 1000000 -t 1.0 -x 0 -f 0.15" 
                pipeCmd=' > fit.log'
                os.system(cmd+pipeCmd)
                
                # Move results to the desired folder
                ofolder=t+"_cfg"+str(ic)
                os.system("mkdir -p "+ofolder)
                os.system("mv -f etapi_result*.fit "+ofolder+" 2>/dev/null")
                os.system("mv -f "+cfgFile+"*.cfg "+ofolder+" 2>/dev/null")
                os.system("mv -f *.log "+ofolder+" 2>/dev/null")
                os.system("mv -f *.ni "+ofolder+" 2>/dev/null")
                os.system("mv -f overlayPlots "+ofolder)
                #os.system("ln -snfr rootFiles "+ofolder+"/rootFiles")

            _j=j+1 if j!=len(ts)-1 else 0
            spawnProcessChangeSetting(ts[j],ts[_j]) # prepare for new t bin
    
    ################################################
    ### Write out a summary of the fits
    ################################################
    if getSummary:
        wall_times_allt=[]
        for t in ts:
            for ic,mmin,mmax in zip(range(len(ms)),mmins,mmaxs):
                print("===================================================")
                print(" -----------   "+t+" * mass ["+str(mmin)+","+str(mmax)+"]   ------------")
                print("===================================================")
        
                interestingLines=["bestMinimum","lastMinuitCommandStatus","eMatrixStatus"]
                fitFiles=[t+"_cfg"+str(ic)+"/"+f for f in os.listdir(t+"_cfg"+str(ic)) if ".fit" in f and has_numbers(f)]
                fitFiles.sort(key=lambda x: int(x.split("_")[-1].split(".")[0]))
                fits=[]
                firstMin=0
                percConverged=0
                percApproxConverged=0
                for j,fitFile in enumerate(fitFiles):
                    with open(fitFile) as fit:
                        lines = [ line.rstrip() for line in fit if any(ele for ele in interestingLines if ele in line) ]
                        bestMinimum = [float(line.split("\t")[-1]) for line in lines if "bestMinimum" in line][0]
                        if j==0:
                            firstMin=bestMinimum
                        status = [line.split("\t")[-1] for line in lines if "lastMinuitCommandStatus" in line][0]
                        estatus = [line.split("\t")[-1] for line in lines if "eMatrixStatus" in line][0]
                        if status=="0" and estatus=="3":
                            print(f"{fitFile} shifted minimum: {bestMinimum-firstMin} status: {status} estatus: {estatus} - converged!")
                            percConverged+=1
                            percApproxConverged+=1
                        elif status=="0" and estatus=="1":
                            print(f"{fitFile} shifted minimum: {bestMinimum-firstMin} status: {status} estatus: {estatus} - approx. converged!")
                            percApproxConverged+=1
                        else:
                            print(f"{fitFile} shifted minimum: {bestMinimum-firstMin} status: {status} estatus: {estatus} - failed!")
                percConverged/=len(fitFiles)
                percApproxConverged/=len(fitFiles)
                print(f"{percConverged*100:0.2f}({percApproxConverged*100:0.2f}) percent (>approximately) converged...")


if __name__ == "__main__":
    start_time=time.time()
    argc=len(sys.argv)
    if argc!=2:
        raise ValueError("Argument [0/1/2] to runFit, getSummary or both")
    arg=int(sys.argv[1])
    main(arg)
    stop_time=time.time()
    print(f"Total elapsed time: {stop_time-start_time}s or {(stop_time-start_time)/60}m")
    
