#!/usr/bin/python

import os
import re #regex
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
from generate_cfg import writeCfg, constructOutputFileName, getPreamble
from determineAmbiguities import executeFinders
import pandas as pd
import itertools
import operator as op
from functools import reduce

logDirBaseName="logs" # creates a folder in each bin_x folder with the name logDirBaseName+"_"+waveset
fitName="EtaPi_fit" # location of the folder containing the inputs that were created from divideData.pl
seedFileTag="param_init" # seedFile name. Should also match the variable from divideData.pl 
writeCfgSeed=-1 # choose a seed to randomly start to sample from
verbose=False
doAccCorr="false" # has to be a string input

############################
# [True, ["data", "bkgnd"], [1, 1]] to bootstrap data 
# [True, ["acc"], [1]] to bootstrap acceptance mc 
# [True, ["data", "bkgnd", "accmc"], [1,1,1]] to bootstrap data and acceptance mc together (I think this is the proper final uncertainty) 
# [True, ["acc", "genmc"], [N,N]] to under or oversample the MC that will be used for acceptance correction. "genmc" needs to be here also or else
#      the acceptance changes
bootstrapSettings=[False,["data","bkgnd"],[1,1]] 
bsFolderTag="_bs_"+"_".join([sample+str(factor)+"x" for sample,factor in zip(bootstrapSettings[1],bootstrapSettings[2])]) # appends a tag to the logs folder in each bin subfolder
bsFolderTag+="" # starts with underscore: include another tag for more folder separation
forceDataBkngdSameSeed=True # data and bkgnd trees can be read with the same seed. If tree same size then it would grab same set of indicies
############################

start = time.time()
workingDir = os.getcwd()
print("\n")
print("current working directory: %s" % (workingDir))
fitDir = workingDir+"/"+fitName
print("fit directory: %s" % (fitDir))

def reorderWaveset(waveset):
    '''
    We want a invariant way to denote a string of underscore separated waves
      i.e. S0++_S0+- = S0+-_S0++
    '''
    tmp=waveset.split("_")
    tmp.sort() # We can apply a sort to order things
    return "_".join(tmp)

def determineBestFit(srcLogDir,binNum):
    '''
    When running fits we can save a file whose name contains the best fit iteration
       This file name can then be grabbed by the bootstrap function to construct
       a new config file that bootstraps around this best fit
    '''
    bestIteration=9999
    bestMinimum=0
    os.system("rm -f FIT_ITER_USED*")
    for i in range(numIters):
        fname=srcLogDir+"/bin_"+str(binNum)+"-"+str(i)+".fit"
        if not os.path.exists(fname):
            raise ValueError("There should be a fit file at this location but there isnt! {}".format(fname))
        with open(fname) as f:
            for line in f:
                ## While reading each file from top to bottom we will first encounter lastMinuitCommandStatus, eMatrixStatus, and then bestMinimum
                if "lastMinuitCommandStatus" in line:
                    #print(line.rstrip())
                    #print(int(line.split("\t")[1]))
                    if int(line.split("\t")[1])!=0:
                        break
                if "eMatrixStatus" in line:
                    #print(line.rstrip())
                    #print(int(line.split("\t")[1]))
                    if int(line.split("\t")[1])!=3:
                        break
                if "bestMinimum" in line:
                    NLL=float(line.rstrip().lstrip().split("\t")[1])
                    if NLL<bestMinimum:
                        bestIteration=i
                        bestMinimum=NLL
                    #print(bestIteration,bestMinimum)
    if bestIteration==9999:
        raise ValueError("(Bin "+str(binNum)+")Doesn't seem like any of the nominal fits "\
                +"reached a good convergence = [lastMinuitCommandStatus, eMatrixStatus] = [0,3]")
    os.system("touch FIT_ITER_USED_FOR_BOOTSTRAP-"+str(bestIteration))

    #### CAN MANUALLY OVERWRITE OR SET bestIterations HERE


def bootstrapCfgGenerator(srcLogDir,binNum,waveset,j):
    '''
    Generates a new config file by loading the best fit parameters from a 
        set of nominal fits (without bootstrapping)
    '''
    bestIteration=9999
    binDir=fitDir+"/bin_"+str(binNum)
    fs=[f for f in os.listdir(binDir) if f.startswith("FIT_ITER")]
    if len(fs)!=1:
        raise ValueError("bootstrapCfgGenerator error: not able to locate file that holds the best fit iteration. Did you run a nominal pre-fit?")
    bestIteration=int(fs[0].split("-")[1])

    ##########################
    # Load the best solution from seedFileTag
    ##########################
    mapInitializeLines={}
    mapParameters={}
    f=srcLogDir+"/"+seedFileTag+"_"+str(bestIteration)+".cfg"
    if not os.path.exists(f):
        raise ValueError("bootstrapCfgGenerator could not find the seed file: {}".format(f))
    with open(srcLogDir+"/"+seedFileTag+"_"+str(bestIteration)+".cfg") as param:
        lines=param.readlines()
        lines=[line.rstrip().lstrip() for line in lines]
        for line in lines:
            if line.startswith("initialize"):
                mapInitializeLines[line.split("::")[2].split(" ")[0]] = " ".join(line.split("::")[2].split(" ")[2:])
            if line.startswith("parameter"):
                mapParameters[line.split(" ")[1]] = line.split(" ")[2]
    #    print(mapInitializeLines)
    #    print(mapParameters)
    
    
    ##########################
    # Load the best solutions configuration file and manually update the cfg file 
    # with the seed file with the best results
    ##########################
    f=srcLogDir+"/"+waveset+"("+str(bestIteration)+").cfg"
    if not os.path.exists(f):
        raise ValueError("bootstrapCfgGenerator could not find the config file: {}".format(f))
    with open(srcLogDir+"/"+waveset+"("+str(bestIteration)+").cfg","r") as cfg:
        lines=cfg.readlines()
        for i,line in enumerate(lines):
            if line.startswith("initialize") and "fixed" not in line:
                wave=line.split("::")[2].split(" ")[0]
                replacement=" ".join(line.split(" ")[:3])+" "+mapInitializeLines[wave]+"\n"
                lines[i]=replacement
            if line.startswith("parameter") and "fixed" not in line:
                parName=line.split(" ")[1]
                replacement=" ".join(line.split(" ")[:2])+" "+mapParameters[parName]+"\n"
                lines[i]=replacement
            if line.startswith("fit"):
                lines[i]="fit bin_"+str(binNum)+"-"+str(j) 
    
    outputCfg=waveset+"("+str(j)+").cfg"
    with open(outputCfg,"w") as cfg:
        cfg.write("".join(lines))
    return outputCfg


def modifyCfgForBootstrapping(input_cfg,bootstrapSettings,binNum,j):
    '''
    Modifes input_cfg such that dataset(i.e. accmc/data/genmc) will use ROOTDataReaderBootstrap with sampling factor given by bootstrapSettings[1]
    '''
    #print("modifying cfg for bootstrapping")
    if bootstrapSettings[0]:
        with open(input_cfg,'r') as cfgFile:
            for line in cfgFile.readlines():
                if any(line.startswith(sample) for sample in bootstrapSettings[1]):
                    line=re.sub('\s+',' ',line) # standardize input by converting multispaces into single spaces
                    if forceDataBkngdSameSeed:
                        random.seed(seeds[binNum][j%numIters]) # rolling seeds. With mod(20): when j=0,20,40 the seed(and thus cfg) will be the same
                    searchStr=line.rstrip().lstrip()
                    prefix=" ".join(searchStr.split(" ")[:2])
                    files=searchStr.split(" ")[3]
                    index=bootstrapSettings[1].index(searchStr.split(" ")[0])
                    args=" -s "+str(random.randrange(9999999))+" -n "+str(bootstrapSettings[2][index])
                    replaceStr=prefix+" ROOTDataReaderBootstrap "+files+args
                    sedArgs=["sed","-i",'s@'+searchStr+'@'+replaceStr+'@g',input_cfg]
                    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    #print("  modified!")


def getAmplitudesInBin(params):
    binNum,lmes,j,logDir=params
    waveset=constructOutputFileName(lmes)
    waveset=reorderWaveset(waveset)
    os.chdir("bin_"+str(binNum))
    seedFile=seedFileTag+"_"+str(j)+".cfg"
    os.system("touch "+seedFile)

    binCfgSrc = "bin_"+str(binNum)+"-full.cfg"

    if os.path.exists(logDir+"_"+waveset):
        fs=os.listdir(logDir+"_"+waveset)
        if len([f for f in fs if f.endswith(".fit")])==numIters:
            print("skipping fits for waveset {} in bin {} since we have all the expected number of fit files already...".format(waveset,binNum))
            os.chdir("..")
            return 0

    ## Load the total {data-bkngd} yield. Used in BIC calculation
    entriesDF=pd.read_csv(fitDir+"/yields.txt",delimiter=" ")
    mapBinToEntries={k:v for k,v in zip(entriesDF["bin"],entriesDF["entries"])}
    entries=mapBinToEntries[binNum]

    if bootstrapSettings[0]:
        if verbose:
            print("Making configuration file for bootstrapping")
        _, pols=getPreamble(binCfgSrc,j) # use it just to get the polarizations
        binCfgDest = bootstrapCfgGenerator(logDirBaseName+"_"+waveset, binNum, waveset, j) # should point to the log base directory where the nominal fits should exist
        modifyCfgForBootstrapping(binCfgDest,bootstrapSettings,binNum,j)
    else:
        # Config generator with random initialization
        if verbose:
            print("Making configuration file for nominal fits")
        binCfgDest, pols=writeCfg(lmes,binCfgSrc,writeCfgSeed,j)
        # The following is left over code from a merge - no use I think
        # replaceStr="fit {}({})".format(waveset,j)
        # searchStr="^fit .*";
        # sedArgs=["sed","-i",'s@'+searchStr+'@'+replaceStr+'@g',binCfgDest]
        # subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
        
    haveHeader=False
    with open(logDir+"_"+waveset+"/amplitude"+str(j)+".txt","w") as outFile, open(logDir+"_"+waveset+"/amplitudeFits"+str(j)+".log","w+") as logFile:
        callFit = "fit -c "+binCfgDest+" -s "+seedFile
        if verbose:
            print(("({0:.1f}s)(Bin {1})(Iteration {2}): "+callFit).format(time.time()-start,binNum,j))
    
        # We needed to create another fit results file so we can run things in parallel
        #resultsFile="{}({}).fit".format(waveset,j)
        resultsFile="bin_"+str(binNum)+"-"+str(j)+".fit"
        if os.path.exists(resultsFile):
            os.remove(resultsFile)

        try:
            output=subprocess.check_output(callFit.split(" "))#, stdout=logFile, stderr=logFile)
        except subprocess.CalledProcessError as err:
            print("*** ABOVE CALL FAILED TO COMPLETE - TRY TO RUN THE FIT MANUALLY IN THE CORRESPONDING BIN FOLDER AS FOLLOWS ***\n*** cd {0}".format(fitDir+"/bin_"+str(binNum))+" ***\n*** "+callFit+" ***\n*** IF RUNFITS.PY DOES NOT EXIT BY ITSELF YOU SHOULD KILL IT MANUALLY NOW ***")
            exit(0)

        os.rename(binCfgDest,logDir+"_"+waveset+"/"+binCfgDest)

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
            if os.path.exists(seedFile) and os.stat(seedFile).st_size!=0: # seedFile only exists if the fit converged
                shutil.move(seedFile,os.getcwd()+"/"+logDir+"_"+waveset+"/"+seedFileTag+"_"+str(j)+".cfg")
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

    return 0

def makeLogFolders(logDir,binNum,lmess):
    '''
    Remove all the log/ directories in each bin folder
    '''
    for lmes in lmess: 
        waveset=constructOutputFileName(lmes)
        waveset=reorderWaveset(waveset)
        fname="bin_"+str(binNum)+"/"+logDir+"_"+waveset
        if not os.path.exists(fname):
            print("making directory "+fname)
            os.mkdir(fname)

def gatherResults(params):
    '''
    Gather all the fit results into one file
    '''
    logDir,fitDir,finalAmpsFolder,binNum,lmes=params
    os.chdir(fitDir)
    if verbose:
        print("---------------------------")
        print("Spawning a process to grab all the amplitudes")
    waveset=constructOutputFileName(lmes)
    waveset=reorderWaveset(waveset)
    if not os.path.exists(workingDir+"/"+finalAmpsFolder+"/"+waveset):
        print("making directory "+workingDir+"/"+finalAmpsFolder+"/"+waveset)
        os.system("mkdir -p "+workingDir+"/"+finalAmpsFolder+"/"+waveset)
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
                    print("READING {}".format(binName+"/"+logDir+"_"+waveset+"/"+afile))
        if len(files)>0:
            os.system("cat "+files[0]+" > amplitudes.txt")
        if len(files)>1:
            os.system("tail -q -n 1 "+" ".join(files[1:])+" >> amplitudes.txt")
    os.system("cp amplitudes.txt "+workingDir+"/"+finalAmpsFolder+"/"+waveset+"/amplitudes-binNum"+str(binNum)+".txt")
    
    if not bootstrapSettings[0]: # If we are not bootstrapping then we will save the best fit
        # Include an extra step here to grab the best fit result
        determineBestFit(logDir+"_"+waveset, binNum) # using logDir as the source directory to load desired cfg
    os.chdir("..")

def gatherMomentResults(params):
    '''
    Gather all the fit results into one file
    '''
    if verbose:
        print("---------------------------")
        print("Spawning a process to grab all the moments")
    logDir,fitDir,finalAmpsFolder,lmes=params
    os.chdir(fitDir)
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
    outfile="moments-binNum"+str(binNum)+".txt"
    cmd="project_moments_polarized -o "+workingDir+"/"+finalAmpsFolder+"/"+waveset+"/"+outfile+" -w "+waveset+" -imax "+str(numIters)+" -b "+str(binNum)
    cmd+=" -mmin "+lowMass+" -mmax "+highMass+" -mbins "+nBins+" -fitdir "+fitName+"/bin_"+str(binNum)+"/"+logDir+"_"+waveset+" -v "+str(verbose)
    print("running: "+cmd)
    os.system(cmd)
    os.chdir("..")


def getVectorOfPotentialLMEs(fitDir,startBin,endBin,numIters):
    '''
    THIS FUNCTION WILL BE USED TO CREATE A LIST OF ALL POTENTIAL WAVESETS OF ALL POSSIBLE SIZES
        GIVEN SOME BASE WAVESET AND THE SET OF POSSIBLE WAVES
    LOOPING THROUGH ALL POTENTIAL WAVE SETS WILL ALLOW US TO DETERMINE SYSTEMATICS RELATED
        TO WAVESET AND DETERMINE HOW LEAKAGE OCCURS. THIS IS AN EXHAUSTIVE SEARCH - EXPENSIVE
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
        seed wavesets are our anchors. In most cases it should be the two S-waves, +/- reflectivity
        '''
        waveset=[]
        for wave in seed:
            waveset.append(list(mapSpectToVectNotation[wave])+[True])
        return waveset

    def determineBestBranch(potential_spects_perBin,seed_waveset,bins):
        '''
        Is a nested fucntion, will only be used as nested
        Function to be used for growing a waveset
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
    if verbose:
        print("---------------------------")
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
            if verbose:
                print("bkgnd files: {}".format(bkgndfiles))
                print("data files: {}\n".format(datafiles))
            
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
    endBin=13
    numIters=50 # number of iterations to randomly sample and try to fit. No guarantees any of them will converge 
    # (int) number of seeds such that the seed used for iteration j is j%nRollingSeeds. 
    #     This gives us a way to reinitialize the parameters in a fit but keep the same bootstrap seed
    nRollingSeeds=numIters 
    # EACH BIN SHARES THE SAME SEED FOR A GIVEN ITERATION
    seeds=[[random.randint(1,999999) for _ in range(nRollingSeeds)] for _ in range(endBin)]
    processes=50 # number of process to spawn to do the fits
    
    finalAmpsFolder="finalAmps"
    logDir=logDirBaseName
    if bootstrapSettings[0]:
        ## Modify logs + finalAmps folders for bootstrapping
        finalAmpsFolder+=bsFolderTag
        logDir+=bsFolderTag 

    ## We have to be careful what we want to delete. Clearly deleting
    ##   all logs_* folders might be unwanted...
    #for i in range(startBin,endBin):
    #    os.system("rm -rf bin_{}/logs_*".format(i))
    #os.system("rm -rf "+workingDir+"/finalAmps")

    ##############################################################
    # DETERMINE ALL THE WAVESETS WE WANT TO RUN OVER
    #      + MAKE SOME FOLDERS
    # MIGHT BE CONFUSING NOTATION I USE {SPECT, VECT}
    #   SPECT IS THE WAVE REPRESENTATION AS A STRING, I.E. S0++
    #   VECT IS A CORRESPONDING VECTOR NOTATION, I.E. [0,0,"+"] = S0+ in [L,M,e] format
    ##############################################################
    #potential_vects=getVectorOfPotentialLMEs(fitDir,startBin,endBin,numIters)
    potential_vects=[
            #### S + TMD
            [
            [0,0,"+",True],
            [0,0,"-",True],
            [2,-1,"-",False],
            [2,0,"+",False],
            [2,0,"-",False],
            [2,1,"+",False],
            [2,1,"-",False],
            [2,2,"+",False]
            ],
            ### K-MATRIX
#            [
#            [0,0,"+",True],
#            [0,0,"-",True],
#            [2,0,"+",False],
#            [2,0,"-",False],
#            [2,2,"-",False],
#            [2,2,"+",False]
#            ],
    ]
    os.chdir(fitDir)
    for ibin in range(startBin,endBin):
        makeLogFolders(logDir,ibin,potential_vects)
    
    ##############################################################
    ## NEED TO EXTRACT THE WEIGHTED YIELDS (DATA-BKGND) SO WE CAN 
    ##      DETERMINE BAYESIAN INFORMATION CRITERA
    ##############################################################
    os.chdir(workingDir)
    mapYields(startBin,endBin,fitDir)

    ##############################################################
    ## BEGIN FITTING AND GETTING AMPLITUDE RESULTS
    ##############################################################
    print("Total number jobs to complete: {}".format((endBin-startBin)*numIters*len(potential_vects)))
    print("** Beginning in 2 seconds **")
    print("----------------------------------\n\n")
    os.chdir(fitDir)
    params=[(i,potential_vect,j,logDir) for i in range(startBin,endBin) for potential_vect in potential_vects for j in range(numIters)]
    p=Pool(processes)
    p.map(getAmplitudesInBin, params)
    p.terminate()
    
    ##############################################################
    ## GATHER ALL THE RESULTS INTO A SINGLE LOCATION
    ##############################################################
    p=Pool(len(potential_vects)*(endBin-startBin))
    params=[(logDir,fitDir,finalAmpsFolder,i,potential_vect) for i in range(startBin,endBin) for potential_vect in potential_vects]
    p.map(gatherResults, params)
    #p.map(gatherMomentResults, params)
    #gatherResults(logDir,fitDir,finalAmpsFolder,ibin,potential_vect)
    #gatherMomentResults(logDir,fitDir,finalAmpsFolder,potential_vect)
    p.terminate()
        

stop = time.time()
print("\nExecution time in seconds: %s" % (stop-start))











