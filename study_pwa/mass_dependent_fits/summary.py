import os
import numpy as np
import pandas as pd

wall_times=[]
total_iterations=[]
eMatrixCodes=[]
minuitCodes=[]

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
        4:"Could be call limit or general failure"
        }


ts=["010020", "0200325", "0325050", "050075", "075100"]
for t in ts:
    print("===================================================")
    print(" ----------------    "+t+"    -----------------")
    print("===================================================")

    files=os.listdir(t)
    iterations=[int(afile.split("fitAttempt")[1].split(".log")[0]) for afile in files if "fitAttempt" in afile]
    iterations.sort()
    
    interestingLines=["bestMinimum","lastMinuitCommandStatus","eMatrixStatus","a2mass\t","a2width\t"]
    for iteration in iterations:
        print("******** ITERATION "+str(iteration)+" *********")
        with open("./"+t+"/fitAttempt"+str(iteration)+".log") as log:
            for line in log:
                if "time" in line:
                    print(line.rstrip().lstrip())
                    if "MIGRAD" in line:
                        wall_times.append(float(line.split(" ")[7].rstrip().lstrip()))
        with open("./"+t+"/etapi0_SD_TMD_piecewise_update"+str(iteration)+".fit") as fit:
            for line in fit:
                if any(ele for ele in interestingLines if ele in line):
                    newline=line.rstrip().lstrip()
                    print(newline)
                    if "lastMinuitCommandStatus" in newline:
                        minuitCodes.append(int(newline.split("\t")[-1]))
                    if "eMatrixStatus" in newline:
                        eMatrixCodes.append(int(newline.split("\t")[-1]))

    total_iterations.append(len(iterations))

print("\n---------------------------------------")
print("t bins: {}".format(ts))
print("total iterations: {}".format(sum(total_iterations)))
print("TOTAL WALL TIME: {:0.0f}s  min={:0.0f}s, max={:0.0f}s, std={:0.0f}s".format(1.0*sum(wall_times),min(wall_times),max(wall_times),np.std(wall_times)))
print("TOTAL WALL TIME: {:0.0f}m".format(1.0*sum(wall_times)/60))
print("TOTAL WALL TIME: {:0.2f}h".format(1.0*sum(wall_times)/60/60))
print("---------------------------------------")

print(minuitCodes)
print(eMatrixCodes)
minuitCodes=pd.Series(minuitCodes).value_counts()
eMatrixCodes=pd.Series(eMatrixCodes).value_counts()
minuitCodes.index=minuitCodes.index.map(minuitCodeMap)
eMatrixCodes.index=eMatrixCodes.index.map(eMatrixCodeMap)
print("")
print(minuitCodes)
print("")
print(eMatrixCodes)








