#!/usr/bin/python

import os
import sys
import numpy as np
import subprocess

nprocesses=9
cfgFile="etapi_hybrid-copy"
fitFileName="etapi_result.fit"
percent=3.0 # parameters must not be within percent of the defined parameter limits
workingDir=os.getcwd()

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


def checkFit(fitFile):
    '''
    Do not care about fit status, only care that the likelihood is reasonable (negative and finite)
    0 Not calculated at all
    ## Error matrix status?
    1 Diagonal approximation only, not accurate
    2 Full matrix, but forced positive-definite
    3 Full accurate covariance matrix (After MIGRAD, this is the indication of normal convergence.)

    ## Minuit status?
    status = 1    : Covariance was made pos defined
    status = 2    : Hesse is invalid
    status = 3    : Edm is above max 
    status = 4    : Reached call limit
    status = 5    : Any other failure 
    '''
    searchMinimumStr="bestMinimum"
    searchMinuitStatusStr="lastMinuitCommandStatus"
    if not os.path.exists(fitFile):
        return False
    else:
        with open(fitFile) as infile:
            for line in infile:
                if searchMinimumStr in line:
                    NLL=float(line.split(" ")[-1].split("\t")[1].rstrip().lstrip())
                if searchMinuitStatusStr in line:
                    minuitStatus=float(line.split(" ")[-1].split("\t")[1].rstrip().lstrip())
        #parNotAtLimit=checkParLimits(fitFile)
        parNotAtLimit=True # ignore if we wish
    return NLL<0 and np.isfinite(NLL) and minuitStatus==0 and parNotAtLimit # Require a negative log likelihood and minuit status = 0


def spawnProcessChangeSetting(old,new):
    '''
    Replace the value of varName to varValue in the file called fileName. Depending on the value type we can include quotes or not
    '''
    sedArgs=["sed","-i","s@"+old+"@"+new+"@g","setup_mass_dep_fits.py"]
    print(" ".join(sedArgs))
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()



#ts=["010016", "016021", "021027", "027034", "034042", "042051", "051061", "061072", "072085", "085100","010016"]
ts=["010020","0200325","0325050","050075","075100","010020"]
#ts=["010020","010020"]

argc=len(sys.argv)
if argc!=2:
    raise ValueError("Argument [0/1/2] to runFit, getSummary or both")
arg=int(sys.argv[1])

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
    for j in range(len(ts)-1):
        t=ts[j]
    
        os.system("mkdir -p overlayPlots")
    
        if os.path.exists(fitFileName):
            print("Fit file already exists... Deleting it to reset fitting")
            os.system("rm "+fitFileName)
            os.system("rm fitAttempt*log")
        
        i=0
        newFitFile="" # initialize fitFile
        nPassedCheck=0
        while not checkFit(newFitFile):
#        while nPassedCheck<1:
#            if checkFit(newFitFile):
#                nPassedCheck+=1
            os.system("python setup_mass_dep_fits.py") # reinitialize
            print("Starting a new fit attempt...")
            cmd="mpirun -np "+str(nprocesses)+" fitMPI -c "+cfgFile+".cfg -m 80000 -t 0.1"
            pipeCmd=' > fitAttempt'+str(i)+'.log'
            os.system(cmd+pipeCmd)
            newFitFile="etapi_result"+str(i)+".fit"
            os.system("mv etapi_result.fit "+newFitFile)
            os.system("mv "+cfgFile+".cfg "+cfgFile+str(i)+".cfg")
            i+=1
    
        # Move results to the desired folder
        os.system("mkdir -p "+t)
        #os.system("mv -f *.root "+t)
        os.system("mv -f etapi_result*.fit "+t)
        os.system("mv -f "+cfgFile+"*.cfg "+t)
        os.system("mv -f *.log "+t)
        os.system("mv -f *.ni "+t)
        os.system("mv -f overlayPlots "+t)
        spawnProcessChangeSetting(ts[j],ts[j+1]) # prepare for new t bin


################################################
### Write out a summary of the fits
################################################
if getSummary:
    wall_times_allt=[]
    for t in ts[:-1]:
        print("===================================================")
        print(" ----------------    "+t+"    -----------------")
        print("===================================================")
    
        files=os.listdir(t)
        iterations=[int(afile.split("fitAttempt")[1].split(".log")[0]) for afile in files if "fitAttempt" in afile]
        iterations.sort()
        
        interestingLines=["bestMinimum","lastMinuitCommandStatus","eMatrixStatus"]
        wall_times=[]
        for iteration in iterations:
            print("******** ITERATION "+str(iteration)+" *********")
            with open("./"+t+"/fitAttempt"+str(iteration)+".log") as log:
                for line in log:
                    if "time" in line:
                        print(line.rstrip().lstrip())
                        if "wall" in line:
                            wall_time=float(line.split(" ")[-2])
                            wall_times.append(wall_time)
            fitFile="./"+t+"/etapi_result"+str(iteration)+".fit"
            with open(fitFile) as fit:
                for line in fit:
                    if any(ele for ele in interestingLines if ele in line):
                        print(line.rstrip().lstrip())
            checkParLimits(fitFile)
        wall_times_allt.append(wall_times)
    
    print("\n\n-----------------------------------------------")
    print("SUMMARY")
    print("----------------------------------------------")
    for it,t in enumerate(ts[:-1]):
        print("tbin {} required {} fits to find one converging".format(t,len(wall_times_allt[it])))
        print("  Total wall time in t bin: {}s".format(sum(wall_times_allt[it])))
    total_wall_time=sum([sum(wall_times) for wall_times in wall_times_allt[:-1]])
    print("Total Wall Time: {}s or {}m or {}h".format(total_wall_time,total_wall_time/60,total_wall_time/3600))


