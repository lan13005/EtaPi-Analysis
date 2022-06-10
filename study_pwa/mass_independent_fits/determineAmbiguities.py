#!/usr/bin/python3
import glob
import os
import numpy as np
import subprocess
from multiprocessing import Pool
import shutil

def runcommand(cmd):
    print("running: "+cmd)
    os.system(cmd)

def constructMapping(paramFile):
    '''
    Creates a mapping of values: i.e.
    {S0+_re: 1, S0+_im: 0, ... }
    '''
    with open(paramFile) as fitFile:
        lines=fitFile.readlines()
        lines=[line for line in lines if line.startswith("initialize")]
        lines=[line for line in lines if searchStrForPols in line]
        lines=[line for line in lines if sum([1 for searchStr in searchStrForAmps if searchStr in line])>0]
        amps=[line.split(" ")[1].split("::")[-1] for line in lines] 
        reals=np.array([float(line.split(" ")[3].rstrip().lstrip()) for line in lines])
        imags=np.array([float(line.split(" ")[4].rstrip().lstrip()) for line in lines]) 
        complexs=reals+1j*imags
        mapping=dict(zip(amps,complexs))
    return mapping

def modifyOriginalLines(ampAmbig,paramFile):
    '''
    Modify original values with the ambiguous amplitude solutions
    Include a fake delimiter (TMP) so we can split and rejoin lines
    '''
    with open(paramFile) as fitFile:
        modifiedLines=fitFile.readlines()
    initLines=[line for line in modifiedLines if line.startswith("initialize")]
    modifiedInitLines=[]
    for key,val in ampAmbig.items():
        lineSelect=[line for line in initLines if key in line][0].rstrip().split(" ")
        lineSelect[3]=str(val.real)
        if str(val.imag)=="-0.0":
            lineSelect[4]=str(-1*val.imag)
        else:
            lineSelect[4]=str(val.imag)
        lineSelect=" ".join(lineSelect)
        modifiedInitLines.append(lineSelect+"\n")
    modifiedLines=[line for line in modifiedLines if not line.startswith("initialize")]
    # Also include a modification to fit output name
    fitSearchLine=[line for line in modifiedLines if line.startswith("fit")][0].rstrip()
    fitReprLine=fitSearchLine+'-ambig'+str(ambigIdx)
    modifiedLines=modifiedLines+["\n\n"]+modifiedInitLines
    modifiedLines="".join(modifiedLines)
    modifiedLines=modifiedLines.replace(fitSearchLine,fitReprLine)
    newFitName=fitReprLine.split(" ")[1]
    return modifiedLines, newFitName

def dumpToFile(cfgFile, modifiedLines):
    '''
    Dump the modified lines to a new file
    '''
    filePrefix=cfgFile.split(".")[0]
    fileType=cfgFile.split(".")[1]
    newCfgName=filePrefix+"-ambig"+str(ambigIdx)+"."+fileType
    if verbose:
        print("\tDumping new cfg file to: "+newCfgName)
    with open(newCfgName,"w") as newCfg:
        newCfg.write(modifiedLines)
    return newCfgName


def a2(amps):
    S=amps["S0+"]
    D0=amps["D0+"]
    D1=amps["D1+"]
    return S+np.sqrt(5)*D0

def delta(amps):
    S=amps["S0+"]
    D0=amps["D0+"]
    D1=amps["D1+"]
#     return 30.*D1**2-8*(S**2-5*D0**2) # CASE 1, assuming Vincent's first iteration eqn is correct
    return 30.*D1**2-4*(S+np.sqrt(5)*D0)*(4*S-2*np.sqrt(5)*D0) # CASE 2, Me and Gabriels Solution

def v(root,amps):
    S=amps["S0+"]
    D0=amps["D0+"]
    D1=amps["D1+"]
    if root==1:
        return 1./(2*a2(amps))*(np.sqrt(30)*D1+np.sqrt(delta(amps)))
    if root==2:
        return 1./(2*a2(amps))*(np.sqrt(30)*D1-np.sqrt(delta(amps)))

def S(xv1,xv2,xa2):
#     return xa2/6*(2+xv1*xv2)
#     return xa2/4*(2+xv1*xv2) # Assuming Case 1
    return xa2/3*(xv1*xv2/2+1) # Assuming Case 2

def D0(xv1,xv2,xa2):
#     return xa2/(6*np.sqrt(5))*(4-xv1*xv2)
#     return xa2/(4*np.sqrt(5))*(2-xv1*xv2) # Assuming Case 1
    return xa2/(3*np.sqrt(5))*(2.-xv1*xv2/2) # Assuming Case 2

def D1(xv1,xv2,xa2):
    return xa2/np.sqrt(30)*(xv1+xv2)

def alignSwave(amp):
    '''
    We typically align the S-wave to the real axis, acting as an anchor
    '''
    amp_copy=amp.copy()
    sangle=np.angle(amp["S0+"])
    srotation=np.exp(-1j*sangle)
    for key in amp_copy.keys():
        amp_copy[key] *= srotation
        amp_copy[key] = np.round(amp_copy[key],5)
    return amp_copy

def determineAmbiguity(ampData):
    '''
    Combines all the steps which takes as input a dictionary of complex amplitudes
      and returns another dictionary of complex amplitudes which are the calculated ambiguities
    '''
    # Determine a2, v1, v2, Delta
    xa2=a2(ampData)
    xv1=v(1,ampData)
    xv2=v(2,ampData)
    # Conjugate v2 and calculate new S, D0, D1
    xv2star=np.conjugate(xv2)
    Salt=S(xv1,xv2star,xa2)
    D0alt=D0(xv1,xv2star,xa2)
    D1alt=D1(xv1,xv2star,xa2)
    ampAmbig={
            "S0+":Salt,
            "D0+":D0alt,
            "D1+":D1alt
    }
    # Align the S-wave again
    ampAmbig=alignSwave(ampAmbig)
    return ampAmbig


def executeFinder(params):
    binLocation, iteration = params

    workingDir=binLocation+'/'+logFolder+'/'
    os.chdir(workingDir)
    filesInBin=glob.glob('*')
    if verbose:
        print("cwd: {}".format(workingDir))
        print("running over: {}".format(workingDir,params))
        print("LOOKING IN: "+workingDir+'*')
        #print(" - Found: {}".format(filesInBin))

    anyCfgFile=[afile.split('/')[-1] for afile in filesInBin if afile.endswith(".cfg")] # this could give param_init.cfg or the amptools input cfg
    aSampleCfgFile=[afile for afile in anyCfgFile if not afile.startswith("param_init") and not "ambig" in afile][0]
    cfgFilePrefix = aSampleCfgFile.split("(")[0]
    cfgFileSuffix = aSampleCfgFile.split(")")[1]
    cfgFile = cfgFilePrefix+"("+str(iteration)+")"+cfgFileSuffix

    paramFile = "param_init_"+str(iteration)+".cfg"
    ampData=constructMapping(paramFile)
    ampAmbig=determineAmbiguity(ampData)
    if verbose:
        print("\tAMPLITUDES IN FIT FILE: {}".format(ampData))
        print("\tCALCULATED AMBIGUITIY: {}".format(ampAmbig))

    modifiedLines, newFitName=modifyOriginalLines(ampAmbig,cfgFile)
    newCfgName=dumpToFile(cfgFile,modifiedLines)

    refittedSeed="param_init_"+str(iteration)+"-ambig"+str(ambigIdx)+".cfg"
    os.system("touch "+refittedSeed)
    
    #######################
    # Run the fit command and check for status
    #######################
    #fitCmd='fit -c "'+newCfgName+'" -s "'+workingDir+refittedSeed+'"' # > /dev/null 2>&1' # we arent saving `fit` program output, send to dev null
    fitCmd='fit -c '+newCfgName+' -s '+workingDir+refittedSeed # > /dev/null 2>&1' # we arent saving `fit` program output, send to dev null
    if verbose:
        print('\tFITTING AMBIGUITY CMD: {}'.format(fitCmd))
    output=subprocess.check_output(fitCmd.split(" "))#, stdout=logFile, stderr=logFile)
    try:
        output=subprocess.check_output(fitCmd.split(" "))#, stdout=logFile, stderr=logFile)
        #os.system(fitCmd)
    except subprocess.CalledProcessError as err:
        print("error in completing fit function in ambiguity finder")
        exit()
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
    print("Status: "+status)
    
    #######################
    # Run getAmpsInBin
    #######################
    with open("amplitude"+str(iteration)+"-ambig"+str(ambigIdx)+".txt","w") as outFile:
        ambigFitFile=newFitName+'.fit'
        getAmplitudeCmd='getAmpsInBin "'+newCfgName+'" "'+ambigFitFile+'" "'+'_'.join(polarizations)+'" "'+str(iteration)+'"'
        if verbose:
            print("\tRUNNING GETAMPINBIN: "+getAmplitudeCmd)
        getAmplitudeCmd=getAmplitudeCmd.split(" ")
        getAmplitudeCmd=[cmd.replace('"','') for cmd in getAmplitudeCmd]
        output=subprocess.check_output(getAmplitudeCmd).decode("utf-8") 
        output=output.split("\n")
        for out in output:
            if len(out.split("\t"))>1:
                if out[0].isdigit():
                    outFile.write(status+"\t")
                    outFile.write("A"+str(ambigIdx)+"\t")
                    outFile.write(out+"\n")

def executeFinders(binLocations,numIters,processes,waveset,pols,xsearchStrForAmps,xsearchStrForPols,xverbose):
    global logFolder
    global verbose
    global searchStrForAmps
    global searchStrForPols
    global polarizations
    global ambigIdx
    logFolder="logs_"+waveset
    verbose=True#xverbose
    searchStrForAmps=xsearchStrForAmps
    searchStrForPols=xsearchStrForPols
    polarizations=pols
    ambigIdx=1

    params=[(i,j) for i in binLocations for j in range(numIters)]
    p=Pool(processes)
    p.map(executeFinder,params)
    p.terminate()

##################################################
## REQUIRING YOUR INPUT
##################################################
#baseDir="/d/grid17/ln16/myDSelector/amptools/"
#fitName="EtaPi_fit"
#newSource=baseDir+fitName+"_logs"
#startBin=0
#endBin=45
#numIters=30
#processes=10
#
## What is your waveset. lmes define your logFolder
#lmes=["S0+","D0+","D1+"]
#logFolder="logs_"+"_".join(lmes)
#
## What polarization string to give to getAmpsInBin
#polarizations=["000","045","090","135","AMO"]
#
## Make a copy of the fit folder to keep things clean
#def remakeFolder():
#    runcommand("rm -rf "+newSource)
#    runcommand("rsync -av --exclude '*root' /dev/shm/"+fitName+" "+newSource)
#    runcommand("mv "+newSource+"/"+fitName+"/* "+newSource)
#    runcommand("rm -r "+newSource+"/"+fitName)
#
##binLocations=glob.glob(newSource+'/*')
#binLocations=[newSource+"/bin_"+str(binNum) for binNum in range(startBin,endBin)]
#print(binLocations)
#print("\n")
#
##remakeFolder()
#params=[(i,j) for i in binLocations for j in range(numIters)]
##with Pool(processes) as p:    
##    p.map(executeFinder,params)
#gatherResults(cwd,newSource,"_".join(lmes))

















