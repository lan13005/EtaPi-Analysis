#!/usr/bin/python

import os
import numpy as np
import pandas as pd

wall_times=[]
total_iterations=[]
eMatrixCodes=[]
minuitCodes=[]
likelihoods=[]

eMatrixCodeMap={
        0:"ematrix not calculated",
        1:"approximation only, not accurate",
        2:"full matrix, forced pos-def",
        3:"full accurate covariance matrix"
        }
minuitCodeMap={
        0:"converged",
        1:"IDK",
        2:"IDK",
        3:"IDK",
        4:"Lawrence's Call limit code",
        5:"General Minuit failure" # Needed to manually split 4/5 since they were giving the same code
        }

files=os.listdir(".")
iterations=[int(afile.split("fitAttempt")[1].split(".log")[0]) for afile in files if "fitAttempt" in afile]
iterations.sort()

def checkFileForString(infile, searchStr):
    with open(infile) as f:
        for line in f:
            if searchStr in line:
                return True
        return False


interestingLines=["bestMinimum","lastMinuitCommandStatus","eMatrixStatus","a2mass\t","a2width\t"]
for iteration in iterations:
    print("******** ITERATION "+str(iteration)+" *********")
    with open("fitAttempt"+str(iteration)+".log") as log:
        for line in log:
            if "time" in line:
                print(line.rstrip().lstrip())
                if "MIGRAD" in line:
                    wall_times.append(float(line.split(" ")[7].rstrip().lstrip()))
    with open("etapi0_SD_TMD_piecewise_update"+str(iteration)+".fit") as fit:
        for line in fit:
            if any(ele for ele in interestingLines if ele in line):
                newline=line.rstrip().lstrip()
                print(newline)
                if "lastMinuitCommandStatus" in newline:
                    minuitCode=int(newline.split("\t")[-1])
                    if minuitCode==4:
                        minuitCode = 4 if checkFileForString("fitAttempt"+str(iteration)+".log","CALL LIMIT EXCEEDED IN MIGRAD") else 5
                    minuitCodes.append(minuitCode)
                if "eMatrixStatus" in newline:
                    eMatrixCodes.append(int(newline.split("\t")[-1]))
                if "bestMinimum" in newline:
                    likelihoods.append(float(newline.split("\t")[-1]))

total_iterations.append(len(iterations))

print("\n---------------------------------------")
print("total iterations: {}".format(sum(total_iterations)))
print("TOTAL WALL TIME: {:0.0f}s  min={:0.0f}s, max={:0.0f}s, std={:0.0f}s".format(1.0*sum(wall_times),min(wall_times),max(wall_times),np.std(wall_times)))
print("TOTAL WALL TIME: {:0.0f}m".format(1.0*sum(wall_times)/60))
print("TOTAL WALL TIME: {:0.2f}h".format(1.0*sum(wall_times)/60/60))
print("---------------------------------------")


def prettyPrintSeries(series,tag):
    for k,v in zip(series.index,series.values):
        print("Number of {} with '{}' = {}".format(tag,k,v))

m=minuitCodes
e=eMatrixCodes
minuitCodes=pd.Series(minuitCodes).value_counts()
eMatrixCodes=pd.Series(eMatrixCodes).value_counts()
minuitCodes.index=minuitCodes.index.map(minuitCodeMap)
eMatrixCodes.index=eMatrixCodes.index.map(eMatrixCodeMap)
prettyPrintSeries(minuitCodes,"minuitCode")
prettyPrintSeries(eMatrixCodes,"eMatrixCode")
#print("minuit codes for each fit {}".format(m))
#print("eMatrixCodes for each fit {}".format(e))
print("")


goodFitBool = [True if m[i]==0 and e[i]==3 else False for i in range(len(m))]
goodFitBool=np.array(goodFitBool)

likelihoods=np.array(likelihoods).astype(int)
likelihoods=likelihoods[goodFitBool]
print("LIKELIHOODS FOR CONVERGED FITS")
print("likelihoods for each fit {}".format(likelihoods))
print("minimum shifted likelihoods {}".format(likelihoods-likelihoods.min()))
wall_times=np.array(wall_times)
print("wall times: {}".format(wall_times[goodFitBool]))
print("Average time per converged fit: {}s".format(sum(wall_times[goodFitBool])/sum(goodFitBool)))
print("")








