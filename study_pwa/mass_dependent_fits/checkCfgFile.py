import collections
import os

def checkCfgFile(cfgLoc):
    #######################################################################################################
    ######### CHECKS IF THE MASS LIMITS OF THE DATAREADERFILTER AND PIECEWISE LIMITS MATCH ################
    #######################################################################################################
    searchKeyWord="Mpi0eta"
    print(cfgLoc)
    print(os.getcwd())
    with open(cfgLoc) as cfg:
        lines=cfg.readlines()
        lines=[line for line in lines if not line.startswith("#")]
        piecewise=[line for line in lines if "Piecewise" in line][0]
        massFilter=[line for line in lines if "ROOTDataReaderFilter" in line][0]
    
    piecewise=piecewise.split(" ")
    minPW=float(piecewise[3])
    maxPW=float(piecewise[4])

    bMassNotMatching=False
    bFilterVarsNotDefined=False
    if searchKeyWord not in set(massFilter.split(" ")):
        bMassNotMatching=True
        bFilterVarsNotDefined=True
    else: # if search word exists
        if len(massFilter.split(" "))>=7:
            kwIdx = massFilter.split(" ").index(searchKeyWord)
            minFilter, maxFilter=massFilter.split(" ")[kwIdx+1:kwIdx+3]
            minFilter, maxFilter=float(minFilter), float(maxFilter)
            if abs(minFilter-minPW)>0.0001 or abs(maxFilter-maxPW)>0.0001:
                bMassNotMatching=True

    ###############################################################################
    ################# CHECK IF THERE ARE ANY DUPLICATED PARAMETERS ################
    ###############################################################################
    parLines=[line.split(" ")[1] for line in lines if line.startswith("parameter")]
    parScanLines=[line.split(" ")[1] for line in lines if line.startswith("parScan")]
    
    bParDuplicated=False
    bParScanDuplicated=False
    duplicatedPars=[]
    if len(parLines)!=len(set(parLines)):
        duplicatedPars=([item for item, count in collections.Counter(parLines).items() if count > 1])
        bParDuplicated=True
    if len(parScanLines)!=len(set(parScanLines)):
        duplicatedParScans=([item for item, count in collections.Counter(parScanLines).items() if count > 1])
        bParScanDuplicated=True

    ###############################################################################
    ############### CHECK IF THERE ARE ENOUGH PIECEWISE PARAMETERS ################
    ###############################################################################
    piecewiseLines=[line for line in lines if line.startswith("amplitude")]
    piecewiseLines=[line.rstrip().lstrip() for line in lines if "Piecewise" in line]
    nbins=piecewiseLines[0].split(" ")[5]
    ampLines=[line.split(" ")[9:] for line in piecewiseLines]
    correctForm=[all([amp.startswith("[")&amp.endswith("]") for amp in amps]) for amps in ampLines]
    bPiecewiseBracketMissing=False
    if not all(correctForm):
        bPiecewiseBracketMissing=True
        incorrectForms=[v for (v,b) in zip(piecewiseLines,correctForm) if not b]
    bPiecewiseParNotDefined=False

    bParsInitAndUsedNotMatching=False
    ampsUsed=set([x[1:-1] for xs in ampLines for x in xs])
    ampsInit=set(parLines)
    if len(ampsUsed-ampsInit)!=0:
        bParsInitAndUsedNotMatching=True

    bNumParsUsedNotMatching=False
    correctNumPars=[len(set(x))==nbins*2 for x in ampLines] 
    if not all(correctNumPars):
        bNumParsUsedNotMatching=True
    incorrectNumPars=[v for (v,b) in zip(piecewiseLines,correctNumPars) if not b]

    print("\n===========================================")
    if bMassNotMatching:
        print("\n****************************")
        print("THERE IS A MISMATCH BETWEEN YOUR MASS RANGE IN YOUR FILTER DATAREADER AND THE RANGE SPECIFIED IN YOUR PIECEWISE DEFINITION!")
        if bFilterVarsNotDefined:
            print("No filter variables defined in ROOTDataReaderFilter!")
        else:
            print("Min/Max value in root file: {}/{}".format(minFilter,maxFilter))
            print("Min/Max value in piecewise: {}/{}".format(minPW,maxPW))
        print("TERMINATING THE PROGRAM!")
        print("****************************")
    if bParDuplicated:
        print("\n****************************")
        print("There are duplicated 'parameter' lines which might lead to errors!")
        print(duplicatedPars)
        print("****************************")
    if bParScanDuplicated:
        print("\n****************************")
        print("There are duplicated 'parScan' lines which might lead to errors!")
        print(duplicatedParScans)
        print("****************************")
    if bPiecewiseBracketMissing:
        print("\n****************************")
        print("piecewise amplitude parts missing brackets detected:")
        for incorrectForm in incorrectForms:
            print(incorrectForm)
        print("****************************")
    if bParsInitAndUsedNotMatching:
        print("\n****************************")
        print("Difference is piecewise bins initialized vs used! Used parameters must be surrounded by brackets so it might be flagged above also")
        print("amps used - initialized: {}".format(ampsUsed-ampsInit))
        print("amps initialized - used: {}".format(ampsInit-ampsUsed))
        print("****************************")
    if bParsInitAndUsedNotMatching:
        print("\n****************************")
        print("The following Piecewise lines has an unexpected number of arguments that does not match with the bins requested")
        print(incorrectNumPars)
        print("****************************")
    print("===========================================\n")
    if any([bMassNotMatching, bParDuplicated, bParScanDuplicated]):
        return False
    else:
        return True

