import os
import sys
import subprocess
import time
import shutil

def printHelp():
    print("\nUsage: \n-----------------")
    print("Running etapi_plotter and plotting their results into a giant pdf")
    print("python overlayBins.py v1 v2 v3 v4 v5 v6 v7")
    print("where v1=0 -- for all bins run etapi_plotter")
    print("      v1=1 -- gather all results from etapi_plotter into overlayPlots folder")
    print("      v1=2 -- do both")
    print('where v2 is a string -- Set to "" to plot everything' )
    print("      v2 is a string containing _ separated amplitudes to plot that are again ; separated to group another set")
    print('      i.e "S0+_D0+;S0+_D0+_P1+" will return plots of the S0+ and D0+ contributions only.')
    print("      The second will include P1+ contribution also")
    print("where v3 is the .fit file")
    print("where v4 is the folder where overlayBins.C is in")
    print("where v5 is a semicolon separated list to gather fit fractions for (used in the etapi_plotter program)")
    print('where v6 is either "T" or "F" to do or not apply acceptance correction')
    print('where v7 is the folder to look for the fit file') 
    
args=sys.argv
nargs=len(args)
if nargs!=8:
    printHelp()
    exit()
else:
    option=int(args[1])
    ampString=str(args[2])
    fitFile=str(args[3])
    cFolder=str(args[4])
    Ls=str(args[5])
    doAccCorr=str(args[6])
    folder=str(args[7])
    groupVec=ampString.split(";")
    groups=[tmp if tmp!="" else tmp for tmp in groupVec]
    groupTags=["_"+tmp if tmp!="" else tmp for tmp in groupVec]
    groupsWithQuotes=['"_'+tmp+'"' if tmp!="" else '""' for tmp in groupVec]
    print("GROUPS: {}".format(groups))

baseDir=os.getcwd()

verbose=True

def runEtaPiPlotterForAllBins():
    print("   running etapi_plotter...")
    #for igroup in range(len(groups)):
    #    cmd=["etapi_plotter",fitFile,"-o","etapi_plot"+groupTags[igroup]+".root","-s",groups[igroup]]
    cmd=["etapi_plotter",fitFile,"-s",ampString,"-a",doAccCorr,"-F",Ls]
    print("calling: "+" ".join(cmd))
    if verbose:
        try:
            subprocess.check_call(cmd)
        except:
            print("  ERROR OCCURED IN ABOVE CALL. Not really sure what the root cause is. Will continue anyways")
            pass
    else:
        with open(os.devnull, 'wb') as devnull:
            try:
                subprocess.check_call(cmd, stdout=devnull, stderr=subprocess.STDOUT)
            except:
                print("  ERROR OCCURED IN ABOVE CALL. Not really sure what the root cause is. Will continue anyways")
                pass

def gatherPlotResultsIntoPDFs():
    cmd=["root","-l","-b","-q","\'"+cFolder+'/overlayBins.C("'+folder+'"'+")\'"]
    print(cmd)
    print(" ".join(cmd))
    os.system(" ".join(cmd)+" 1> /dev/null 2> /dev/null")
    #if verbose:
    #    subprocess.check_call(cmd)
    #else:
    #    with open(os.devnull, 'wb') as devnull:
    #        subprocess.check_call(cmd, stdout=devnull, stderr=subprocess.STDOUT)

if option==0 or option==2:
    print("RUNNING etapi_plotter FOR ALL BINS...")
    runEtaPiPlotterForAllBins()
if option==1 or option==2:
    print("Recreating overlayPlots folder...")
    os.chdir(baseDir)
    if os.path.exists("overlayPlots") and os.path.isdir("overlayPlots"):
        shutil.rmtree("overlayPlots")
    os.mkdir("overlayPlots")
    print("GATHERING ALL THE RESULTS INTO A PDF...")
    gatherPlotResultsIntoPDFs()

print("DONE! Dumped in overlayPlots folder...")
