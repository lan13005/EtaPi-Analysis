import os
import numpy as np


def has_numbers(inputString):
    return any(char.isdigit() for char in inputString)

def checkFits(folder):
    '''
    ## Do not care about fit status, only care that the likelihood is reasonable (negative and finite)

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
    searchErrMatrixStatus="eMatrixStatus"

    files=[folder+"/"+f for f in os.listdir(folder) if ".fit" in f and has_numbers(f)]
    
    convergenceStatuses=[]
    minuitStatuses=[]
    eMatrixStatuses=[]
    for fitFile in files:
        with open(fitFile) as infile:
            NLL=0
            for line in infile:
                if searchMinimumStr in line:
                    NLL=float(line.split(" ")[-1].split("\t")[1].rstrip().lstrip())
                if searchMinuitStatusStr in line:
                    minuitStatus=float(line.split(" ")[-1].split("\t")[1].rstrip().lstrip())
                if searchErrMatrixStatus in line:
                    eMatrixStatus=float(line.split(" ")[-1].split("\t")[1].rstrip().lstrip())
            if NLL==0:
                raise ValueError("NLL not found in "+fitFile+"! Terminating")
            parNotAtLimit=True # ignore if we wish
            convergenceStatus=NLL<0 and np.isfinite(NLL)# and minuitStatus==0 and eMatrixStatus==3 and parNotAtLimit # Require a negative log likelihood and minuit status = 0
            convergenceStatuses.append(convergenceStatus)
            minuitStatuses.append(minuitStatus)
            eMatrixStatuses.append(eMatrixStatus)

    #print("file\tminuitStatus\teMatrixStatus\tstatus")
    mapStatus={}
    for f,m,e,s in zip(files,minuitStatuses,eMatrixStatuses,convergenceStatuses):
        print(f"{f}\t{m:0.0f}\t{e:0.0f}\t{s}")
    return np.array(files)[convergenceStatuses]




