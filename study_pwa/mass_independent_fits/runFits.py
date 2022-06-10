#!/usr/bin/python

import os
import numpy as np
import shutil
import math
import glob
import random
import sys
import subprocess
import time
import fileinput
from multiprocessing import Pool
from generate_cfg import writeCfg, constructOutputFileName
from determineAmbiguities import executeFinders
import pandas as pd
import itertools
import operator as op
from functools import reduce

fitName="EtaPi_fit" # location of the folder containing the inputs that were created from divideData.pl
seedFileTag="param_init" # seedFile name. Should also match the variable from divideData.pl 
seedAmpInit=9183 # choose a seed to randomly start to sample from
verbose=False
doAccCorr="true" # has to be a string input
factorSample=[False,5] # (first arg)should we bootstrap the acc mc with a sampling factor(second argument)
keep_logs_fit_seed=False

start = time.time()
workingDir = os.getcwd()
print("\n")
print("current working directory: %s" % (workingDir))
fitDir = workingDir+"/"+fitName
print("fit directory: %s" % (fitDir))

def reorderWaveset(waveset):
    tmp=waveset.split("_")
    tmp.sort()
    return "_".join(tmp)

def getAmplitudesInBin(params):
    binNum,lmes,j=params

    entriesDF=pd.read_csv(fitDir+"/yields.txt",delimiter=" ")
    mapBinToEntries={k:v for k,v in zip(entriesDF["bin"],entriesDF["entries"])}
    entries=mapBinToEntries[binNum]

    waveset=constructOutputFileName(lmes)
    waveset=reorderWaveset(waveset)
    if os.path.exists(workingDir+"/finalAmps"):
        resultsFolder=os.listdir(workingDir+"/finalAmps")
        if any([os.path.exists(workingDir+"/finalAmps/"+f) for f in resultsFolder if set(f.split("_"))==set(waveset.split("_"))]):
            print("skipping fit since the waveset {} already exists...".format(waveset))
            return 1
    
    seedFile=seedFileTag+"_"+str(j)+".cfg"
    os.chdir("bin_"+str(binNum))

    binCfgSrc = "bin_"+str(binNum)+"-full.cfg"
    binCfgDest, pols=writeCfg(lmes,binCfgSrc,seedAmpInit,j)
    replaceStr="fit {}({})".format(waveset,j)
    searchStr="^fit .*";
    sedArgs=["sed","-i",'s@'+searchStr+'@'+replaceStr+'@g',binCfgDest]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
        
    if keep_logs_fit_seed:
        os.system("touch "+seedFile)

    if factorSample[0]:
        with open(binCfgDest,'r') as cfgFile:
            for line in cfgFile.readlines():
                if line.startswith("accmc"):
                    searchStr=line.rstrip().lstrip()
        prefix=" ".join(searchStr.split(" ")[:2])
        files=searchStr.split(" ")[3]
        args=" -s "+str(random.randrange(9999999))+" -n "+str(factorSample[1])
        replaceStr=prefix+" ROOTDataReaderBootstrap "+files+args
        sedArgs=["sed","-i",'s@'+searchStr+'@'+replaceStr+'@g',binCfgDest]
        subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    
    haveHeader=False
    if keep_logs_fit_seed:
        logFile=open(logDir+"_"+waveset+"/amplitudeFits"+str(j)+".log","w+")
    with open(logDir+"_"+waveset+"/amplitude"+str(j)+".txt","w") as outFile:
        if keep_logs_fit_seed:
            callFit = "fit -c "+binCfgDest+" -s "+seedFile
        else:
            callFit = "fit -c "+binCfgDest
        if verbose:
            print(("({0:.1f}s)(Bin {1})(Iteration {2}): "+callFit).format(time.time()-start,binNum,j))
    
        # We needed to create another fit results file so we can run things in parallel
        resultsFile="{}({}).fit".format(waveset,j)
        if os.path.exists(resultsFile):
            os.remove(resultsFile)

        try:
            output=subprocess.check_output(callFit.split(" "))#, stdout=logFile, stderr=logFile)
        except subprocess.CalledProcessError as err:
            print("*** ABOVE CALL FAILED TO COMPLETE - TRY TO RUN THE FIT MANUALLY IN THE CORRESPONDING BIN FOLDER AS FOLLOWS ***\n*** cd {0}".format(fitDir+"/bin_"+str(binNum))+" ***\n*** "+callFit+" ***\n*** IF RUNFITS.PY DOES NOT EXIT BY ITSELF YOU SHOULD KILL IT MANUALLY NOW ***")
            exit(0)

        os.rename(binCfgDest,logDir+"_"+waveset+"/"+binCfgDest)

        if keep_logs_fit_seed:
            logFile.write("ITERATION: "+str(j)+"\n")
            logFile.write(output)
        if "STATUS=CONVERGED" in output:
            status="C" # (C)onverged
        elif "STATUS=FAILED" in output:
            if "ERR MATRIX NOT POS-DEF" in [line.lstrip().rstrip() for line in output.split("\n")]:
                status="H" # (H)essian failed
            else:
                status="F" # (F)ailed 
        elif "STATUS=CALL LIMIT" in output:
            status="L" # (L)imit call 
        else:
            status="U" # (U)ncertain / unsure / unqualified
        if verbose:
            print("Status: "+status)

        resultsFilePath=logDir+"_"+waveset+"/"+resultsFile
        if verbose:
            print("Moving fit results to: "+os.getcwd()+"/"+resultsFilePath)
        if os.path.exists(resultsFile) and os.stat(resultsFile).st_size!=0: # if the fit files is not empty then we will try and use it
            shutil.move(resultsFile,resultsFilePath)
            if os.path.exists(seedFile) and os.stat(seedFile).st_size!=0: # param_init.cfg only exists if the fit converged
                shutil.move(seedFile,os.getcwd()+"/"+logDir+"_"+waveset+"/param_init_"+str(j)+".cfg")
            getAmplitudeCmd='getAmpsInBin "'+binCfgDest+'" "'+resultsFilePath+'" "'+pols+'" "'+str(j)+'" "'+doAccCorr+'" "'+str(entries)
            if verbose:
                print(getAmplitudeCmd)
            getAmplitudeCmd=getAmplitudeCmd.split(" ")
            getAmplitudeCmd=[cmd.replace('"','') for cmd in getAmplitudeCmd]
            output=subprocess.check_output(getAmplitudeCmd)
            output=output.split("\n")
            for out in output:
                if len(out.split("\t"))>1:
                    if haveHeader==False and out.split("\t")[-1]=="iteration": #not out[0].isdigit():
                        outFile.write("status\t")
                        outFile.write("solution\t") # In here solution is always O. When determining ambiguities solution would be like A1, A2 ... 
                        outFile.write(out+"\n")
                        haveHeader=True
                    if out[0].isdigit():
                        outFile.write(status+"\t")
                        outFile.write("O\t") 
                        outFile.write(out+"\n")
        else:
            print("fit file does not exist or is empty! The fit program failed to complete on bin {}".format(binNum))
            print("  fit file exists? {}".format(os.path.exists(resultsFile)))
            print("  fit file has non-zero size? {}".format(os.stat(resultsFile).st_size!=0))

    os.chdir("..")
    if keep_logs_fit_seed:
        logFile.close()

    return 0

def makeLogFolders(binNum,lmess):
    '''
    Remove all the log/ directories in each bin folder
    '''
    print("Cleaning relevant log folder of bin"+str(binNum))
    for lmes in lmess: 
        waveset=constructOutputFileName(lmes)
        waveset=reorderWaveset(waveset)
        if not os.path.exists("bin_"+str(binNum)+"/"+logDir+"_"+waveset):
            print("making directory bin_"+str(binNum)+"/"+logDir+"_"+waveset)
            os.mkdir("bin_"+str(binNum)+"/"+logDir+"_"+waveset)

def gatherResults(binNum,lmess):
    '''
    Gather all the fit results into one file
    '''
    if verbose:
        print("Grabbing all the results")
    for lmes in lmess:
        waveset=constructOutputFileName(lmes)
        waveset=reorderWaveset(waveset)
        if not os.path.exists(workingDir+"/"+"finalAmps/"+waveset):
            print("making directory "+workingDir+"/"+"finalAmps/"+waveset)
            os.system("mkdir -p "+workingDir+"/"+"finalAmps/"+waveset)
        binName="bin_"+str(binNum)
        os.chdir(binName)
        files=[]
        listFiles=os.listdir(logDir+"_"+waveset)
        # sort the directories by length so names with "ambig" tags are last. They do not have a header for amplitude.txt
        listFiles=sorted(listFiles, key=lambda x: (len(x), x)) 
        for afile in listFiles:
            fileTag=afile.split(".")[0]
            if "amplitude" in fileTag:
                if "Fit" not in fileTag:
                    files.append(logDir+"_"+waveset+"/"+afile) 
                    if verbose:
                        print(binName+"/"+logDir+"_"+waveset+"/"+afile)
            if len(files)>0:
                os.system("cat "+files[0]+" > amplitudes.txt")
            if len(files)>1:
                os.system("tail -q -n 1 "+" ".join(files[1:])+" >> amplitudes.txt")
        os.system("cp amplitudes.txt "+workingDir+"/finalAmps/"+waveset+"/amplitudes-binNum"+str(binNum)+".txt")
        os.chdir("..")

def gatherMomentResults(lmes,verbose):
    '''
    Gather all the fit results into one file
    '''
    # Need to grab the mass binning to input to project_moments_polarized
    waveset=constructOutputFileName(lmes)
    waveset=reorderWaveset(waveset)
    with open(workingDir+"/divideData.pl") as f:
        for line in f:
            if line.startswith("$lowMass"):
                lowMass=line.split("=")[-1].split(";")[0].rstrip().lstrip()
            if line.startswith("$highMass"):
                highMass=line.split("=")[-1].split(";")[0].rstrip().lstrip()
            if line.startswith("$nBins"):
                nBins=line.split("=")[-1].split(";")[0].rstrip().lstrip()
            if line.startswith("$fitName"):
                fitName=line.split("=")[-1].split(";")[0].rstrip().lstrip()
    print("Grabbing all the moments results")
    os.chdir(workingDir)
    for binNum in range(startBin,endBin):
        outfile="moments-binNum"+str(binNum)+".txt"
        cmd="project_moments_polarized -o "+workingDir+"/finalAmps/"+waveset+"/"+outfile+" -w "+waveset+" -imax "+str(numIters)+" -b "+str(binNum)
        cmd+=" -mmin "+lowMass+" -mmax "+highMass+" -mbins "+nBins+" -fitdir "+fitName+"/bin_"+str(binNum)+"/logs_"+waveset+" -v "+str(verbose)
        print("running: "+cmd)
        os.system(cmd)
        os.chdir("..")


def getVectorOfPotentialLMEs(fitDir,startBin,endBin,numIters):
    '''
    THIS FUNCTION WILL BE USED TO CREATE A LIST OF ALL POTENTIAL WAVESETS OF ALL POSSIBLE SIZES
        GIVEN SOME BASE WAVESET AND THE SET OF POSSIBLE WAVES
    LOOPING THROUGH ALL POTENTIAL WAVE SETS WILL ALLOW USE TO DETERMINE SYSTEMATICS RELATED
        TO WAVESET AND DETERMINE HOW LEAKAGE OCCURS. THIS IS AN EXHAUSTIVE SEARCH
    '''
    def comb(n, r):
        ''' Formula for n choose r '''
        r = min(r, n-r)
        numer = reduce(op.mul, range(n, n-r, -1), 1)
        denom = reduce(op.mul, range(1, r+1), 1)
        return numer / denom  # or / in Python 2
    def constructSeed(seed):
        '''
        Is a nested fucntion, will only be used as nested
        Function to be used for growing waveset
        seed wavesets will also be our anchors. In most cases it should be the two S-waves, +/- reflectivity
        '''
        waveset=[]
        for wave in seed:
            waveset.append(list(mapSpectToVectNotation[wave])+[True])
        return waveset

    def determineBestBranch(potential_spects_perBin,seed_waveset,bins):
        '''
        Is a nested fucntion, will only be used as nested
        Function to be used for growing waveset
        Determine best branch to grow the waveset in the direction of 
        '''
        resultsFolder=os.listdir(workingDir+"/finalAmps")
        aics=[]
        startBin, endBin = bins
        best_lmess_perBin=[]
        best_spect_perBin=[]
        best_aics=[]
        for ibin in range(startBin, endBin):
            potential_spects=potential_spects_perBin[ibin-startBin]
            aics_across_wavesets=[]
            print(ibin,potential_spects)
            for potential_spect in potential_spects:
                srcFolder=["finalAmps/"+f for f in resultsFolder if set(f.split("_"))==set(potential_spect)]
                print(potential_spect)
                print(srcFolder)
                if len(srcFolder)!=1:
                    raise ValueError("Results folder does not exist for {}! Fix me...".format("_".join(potential_spect)))
                amplitudeDF=pd.read_csv(srcFolder[0]+"/amplitudes-binNum"+str(ibin)+".txt",delimiter='\t')
                min_aic=amplitudeDF[amplitudeDF.status=="C"]['aic'].min()
                aics_across_wavesets.append(min_aic)
            nans=np.isnan(aics_across_wavesets)
            if nans.any():
                print("AIC array has a nan value! Will zero it aics={}".format(aics_across_wavesets))
            aics_across_wavesets=np.nan_to_num(aics_across_wavesets)
            print("AICS FOR CURRENT WAVESET SIZE={}".format(len(potential_spects[0])))
            print("(bin{}: {}".format(ibin,aics_across_wavesets))
            imin=np.argmin(aics_across_wavesets)
            best_lmess=[]
            for spect in potential_spects[imin]:
                if spect in seed_waveset:
                    best_lmess.append(list(mapSpectToVectNotation[spect])+[True])
                else:
                    best_lmess.append(list(mapSpectToVectNotation[spect])+[False])
            best_lmess_perBin.append(best_lmess)
            best_spect_perBin.append(potential_spects[imin])
            best_aics.append(aics_across_wavesets[imin])
        return best_lmess_perBin, best_spect_perBin, best_aics


    ##################################
    # What are the potential waves? 
    #   Determine all the potential waves given the Ls and es. 
    #   Basically filling in all the m-projections
    ##################################
    potential_Ls=["S","D"]
    es=["+","-"]
    seed_waveset = ["S0+","S0-"]

    all_potential_spect=[]
    all_potential_vect=[]
    spectroscopicMapping={"S":0, "P":1, "D":2}
    for L in potential_Ls:
        Lmax=spectroscopicMapping[L]
        for l in range(-1*Lmax,Lmax+1): # inclusive of Lmax
            for e in es:
                all_potential_spect.append(L+str(l)+e)
                all_potential_vect.append((spectroscopicMapping[L],l,e))
    mapSpectToVectNotation={spect:vect for spect,vect in zip(all_potential_spect,all_potential_vect)}
    mapVectToSpectNotation={vect:spect for vect,spect in zip(all_potential_vect,all_potential_spect)}
    all_potential_spect=all_potential_spect#[:4]

    if processes > (endBin-startBin)*numIters*(len(all_potential_vect)-len(seed_waveset)):
        print("You are trying to spawn more processes than I expect you can use currently")
        print(" choose better")
        exit();
    print("Potential waves: {}".format(all_potential_spect))
    print("----------------------------------")

    ###################################
    #### ALL POTENTIAL WAVESETS ####
    ###################################
    unused_waves=list(set(all_potential_spect)-set(seed_waveset))
    growIters=len(unused_waves)
    print("Total number jobs to complete: {}".format((endBin-startBin)*numIters*sum([comb(growIters,j+1) for j in range(growIters)])))
    print("** Beginning in 2 seconds **")
    print("----------------------------------")
    time.sleep(2)
    potential_vects=[]
    potential_spects=[]
    for growIter in range(growIters):
        for combination in list(itertools.combinations(unused_waves, growIter+1)):
            current_spect=seed_waveset+list(combination)
            current_vect=constructSeed(seed_waveset)+[list(mapSpectToVectNotation[wave])+[False] for wave in combination]
            potential_vects.append(current_vect)
            potential_spects.append(current_spect)

    return potential_vects

def mapYields(startBin,endBin,fitDir):
    '''
    FUNCTION TO LOAD THE WEIGHTED YIELDS IN ALL MASS BINS FOR A GIVEN FIT DIRECTORY
    THE RESULT WILL BE USED IN THE CALCULATION OF THE BAYESIAN INFORMATION CRITERIA
    '''
    print("\nGrabbing yields in each mass bin to be used for BIC calculation...")
    print("    this surprisingly takes a long time\n")
    datareader="ROOTDataReader"
    if os.path.exists(fitDir+"/yields.txt"):
        return 1 
    with open(fitDir+"/yields.txt","w") as outfile:
        outfile.write("bin entries\n")
        for i in range(startBin,endBin):
            dir1=fitDir+"/bin_"+str(i)+"/"
            mapToFile={}
            ## Load the cfg file to grab the {data, bkgnd, accmc, genmc} if applicaple
            with open(dir1+"bin_"+str(i)+"-full.cfg") as f:
                lines=f.readlines()
                for line in lines:
                    if datareader in line:
                        l=line.rstrip().lstrip().split(" ")
                        ## If there is a single root file for the ROOTDataReader we can keep its location
                        ##    if there were multiple then looping must have happened and we can rescan
                        ##    the config file for the correct lines
                        if ".root" in l[-1]:
                            mapToFile[l[0]]=[l[-1].split("/")[-1]]
                        else:
                            for line2 in lines:
                                if l[-1] in line2 and datareader not in line2:
                                    if l[0][0]!="#":
                                        mapToFile[l[0]]=[ele.split("/")[-1] for ele in line2.rstrip().lstrip().split(" ")[2:]]
            
            bkgndfiles=mapToFile["bkgnd"] if "bkgnd" in mapToFile else []
            datafiles=mapToFile["data"]
            print("bkgnd files: {}".format(bkgndfiles))
            print("data files: {}".format(datafiles))
            print("")
            
            yields=[]
            for files in [datafiles, bkgndfiles]:
                total=0
                for f in files:
                    cmd=["getIntegral.py", dir1+f, "0"]
                    cmd=" ".join(cmd)
                    total+=float(os.popen(cmd).read())
                yields.append(total)
            # the entries that will be used for the likelihood fit is the signal yield
            #     we do not include the yields of the accmc - malte
            entries=yields[0]-yields[1]
            outfile.write(str(i)+" "+str(entries)+"\n")
    return 0 # return output is never used


### CHOOSE BIN NUMBER
if __name__ == '__main__':
    os.chdir(fitDir)
    startBin=0
    endBin=25
    numIters=50 # number of iterations to randomly sample and try to fit. No guarantees any of them will converge
    # EACH BIN SHARES THE SAME SEED FOR A GIVEN ITERATION
    seeds=[random.randint(1,100000) for _ in range(numIters)]
    processes=50 # number of process to spawn to do the fits
    logDir="logs"

    for i in range(startBin,endBin):
        os.system("rm -rf bin_{}/logs_*".format(i))
    os.system("rm -rf "+workingDir+"/finalAmps")

    ##############################################################
    # DETERMINE ALL THE WAVESETS WE WANT TO RUN OVER
    #      + MAKE SOME FOLDERS
    ##############################################################
    #potential_vects=getVectorOfPotentialLMEs(fitDir,startBin,endBin,numIters)
    potential_vects=[
            #### S + TMD
            #[
            #[0,0,"+",True],
            #[0,0,"-",True],
            #[2,-1,"-",False],
            #[2,0,"+",False],
            #[2,0,"-",False],
            #[2,1,"+",False],
            #[2,1,"-",False],
            #[2,2,"+",False]
            #],
            ### K-MATRIX
            [
            [0,0,"+",True],
            [0,0,"-",True],
            [2,0,"+",False],
            [2,0,"-",False],
            [2,2,"-",False],
            [2,2,"+",False]
            ],
    ]
    os.chdir(fitDir)
    for ibin in range(startBin,endBin):
        makeLogFolders(ibin,potential_vects)
    
    ##############################################################
    ## NEED TO EXTRACT THE WEIGHTED YIELDS (DATA-BKGND) SO WE CAN 
    ##      DETERMINE BAYESIAN INFORMATION CRITERA
    ##############################################################
    os.chdir(workingDir)
    mapYields(startBin,endBin,fitDir)

    ##############################################################
    ## BEGIN FITTING AND GETTING AMPLITUDE RESULTS
    ##############################################################
    print("begin fitting...")
    os.chdir(fitDir)
    params=[(i,potential_vect,j) for i in range(startBin,endBin) for potential_vect in potential_vects for j in range(numIters)]
    p=Pool(processes)
    p.map(getAmplitudesInBin, params)
    p.terminate()
    
    ##############################################################
    ## GATHER ALL THE RESULTS INTO A SINGLE LOCATION
    ##############################################################
    for ibin in range(startBin,endBin):
        os.chdir(fitDir)
        gatherResults(ibin,potential_vects)
        

stop = time.time()
print("\nExecution time in seconds: %s" % (stop-start))











