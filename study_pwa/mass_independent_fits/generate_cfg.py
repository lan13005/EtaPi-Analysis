import random
import sys
import re
import os

'''
The goal of this code is to generate config files with any type of waveset
Part of the reference config file will be taken to form the preamble of the output
config file. This preamble will contain the dataset names, fit name, scales, etc

For each wave we want to include we need to do the following
    1. Initialize the wave
    2. Constrain the wave
    3. Scale the waves
'''

verbose=False
delimiterForOutFileName="_" # The output config file would have a name like 000__045__090

cfgFileLoops=True # DO NOT MODIFY THIS ON YOUR OWN
reactionName="LOOPREAC" # DO NOT MODIFY THIS ON YOUR OWN
reactionAngle="LOOPPOLANG" # DO NOT MODIFY THIS ON YOUR OWN
reactionPol="LOOPPOLVAL" # DO NOT MODIFY THIS ON YOUR OWN

def constructWave(l,m,e):
    '''
    constructs the wave in the typical notation
    '''
    mapSpinToSpectroscopic={0:"S",1:"P",2:"D"}
    waveString=mapSpinToSpectroscopic[l]+str(m)
    if e=="+":
        waveString+="+"
    else:
        waveString+="-"
    return waveString

def constructWaveString(l,m,e,c,prefix):
    '''
    l = Spin {0, 1, 2 ...}
    m = Spin projection {..., -2, -1, 0, 1+, 2+, ... }
    e = reflectivity {+/-}
    c = component (Re)al or (Im)aginary part
    prefix = string containing the dataset name
    '''
    assert l <= 2 and l >= 0  # Just do up to D waves for now
    assert abs(m) <= l
    assert e in ["+","-"]

    waveString=prefix+"::"
    if e=="+":
        waveString+="Positive"
    else:
        waveString+="Negative"
    waveString+=c+"::"
    waveString+=constructWave(l,m,e)
    return waveString
    
def defineWave(l,m,e):
    '''
    l = Spin {0, 1, 2 ...}
    m = Spin projection {..., -2, -1, 0, 1+, 2+, ... }
    e = reflectivity {+/-}
    '''
    assert l <= 2 and l >= 0  # Just do up to D waves for now
    assert abs(m) <= l
    assert e in ["+","-"]

    prefix="amplitude "
    outputStrs=[]
    for c in ["Re","Im"]:
        outputStr = prefix
        outputStr += constructWaveString(l,m,e,c,reactionName)
        outputStr+=" Zlm "
        waveValues=[str(l), str(m)]
        if e=="+" and c=="Re":
            waveValues+=["+1","+1"]
        if e=="+" and c=="Im":
            waveValues+=["-1","-1"]
        if e=="-" and c=="Re":
            waveValues+=["+1","-1"]
        if e=="-" and c=="Im":
            waveValues+=["-1","+1"]
        outputStr += " ".join(waveValues)
        outputStr += " "+reactionAngle+" "+reactionPol
        outputStrs.append(outputStr)
    return "\n".join(outputStrs)

def initializeWave(l,m,e,anchor):
    '''
    l = Spin {0, 1, 2 ...}
    m = Spin projection {..., -2, -1, 0, 1+, 2+, ... }
    e = reflectivity {+/-}
    anchor = boolean to set this wave as the anchor, anchor wave requires wave to be positive
    '''
    if verbose:
        print("initializing wave lme={0}{1}{2}".format(l,m,e))
    prefix="initialize "
    c="Re" # since we constrain the (Re)al and (Im)aginary parts of the waves to be the same we only need to initialize one part
    outputStr = prefix
    outputStr += constructWaveString(l,m,e,c,reactionName)
    outputStr += " cartesian "
    if anchor:
        outputStr += str(random.uniform(-100,100)) + " 0.0 real"
    else:
        outputStr += str(random.uniform(-100,100)) + " " + str(random.uniform(-100,100)); 
    return outputStr

def constrainWave(l,m,e,preamble):
    '''
    l = Spin {0, 1, 2 ...}
    m = Spin projection {..., -2, -1, 0, 1+, 2+, ... }
    e = reflectivity {+/-}
    preamble = data copied over from the reference config, will be used to the fit name
    '''
    if verbose:
        print("constraining wave lme={0}{1}{2}".format(l,m,e))
    prefix="constrain "
    outputStr = prefix
    # First we constrain the Re and Im parts of a given wave
    if cfgFileLoops:
        outputStr += constructWaveString(l,m,e,"Re","LOOPREAC") + " " + constructWaveString(l,m,e,"Im","LOOPREAC")
        outputStr += "\n"
        dataset_name = [line for line in preamble.split("\n") if line.startswith("loop LOOPREAC")]
        dataset_name = dataset_name[0].split(" ")[2]
        # Second we constrain the Re/Im parts of a given wave for a given dataset to the Re/Im parts of the same wave for a different waveset
        outputStr += "constrain " + constructWaveString(l,m,e,"Re",dataset_name) + " " + constructWaveString(l,m,e,"Re","LOOPREAC")
        outputStr += "\n"
        outputStr += "constrain " + constructWaveString(l,m,e,"Im",dataset_name) + " " + constructWaveString(l,m,e,"Im","LOOPREAC")
    else:
        outputStr += constructWaveString(l,m,e,"Re",reactionName) + " " + constructWaveString(l,m,e,"Im",reactionName)
    return outputStr

def scaleWave(l,m,e):
    '''
    l = Spin {0, 1, 2 ...}
    m = Spin projection {..., -2, -1, 0, 1+, 2+, ... }
    e = reflectivity {+/-}
    '''
    if verbose:
        print("scaling wave lme={0}{1}{2}".format(l,m,e))
    if cfgFileLoops:
        prefix="scale "
        outputStrs=[]
        for c in ["Re","Im"]:
            outputStr = prefix
            outputStr += constructWaveString(l,m,e,c,"LOOPREAC") + " LOOPSCALE"
            outputStrs.append(outputStr)
        return "\n".join(outputStrs)
    else:
        return ""
    
def writeWave(l,m,e,anchor,preamble):
    '''
    l = Spin {0, 1, 2 ...}
    m = Spin projection {..., -2, -1, 0, 1+, 2+, ... }
    e = reflectivity {+/-}
    anchor = boolean to set this wave as the anchor, anchor wave requires wave to be positive
    '''
    if verbose:
        print("writing wave lme={0}{1}{2}".format(l,m,e))
    outputList=[
            defineWave(l,m,e),
            initializeWave(l,m,e,anchor),
            constrainWave(l,m,e,preamble),
            scaleWave(l,m,e)
            ]
    return "\n".join(outputList)

def constructOutputFileName(lmes,i=-1):
    mapLtoSpect={0:"S",1:"P",2:"D"};
    names=[mapLtoSpect[lme[0]]+str(lme[1])+lme[2] for lme in lmes]
    cfgFileName=delimiterForOutFileName.join(names)
    cfgFileName=cfgFileName.split("_")
    cfgFileName.sort() # We can apply a sort to order things
    cfgFileName="_".join(cfgFileName)
    if i==-1:
        return cfgFileName # would basically be the waveset
    else:
        return cfgFileName+"("+str(i)+").cfg"

def getPreamble(reference_file,i):
    '''
    reference_file: we will get our preamble from here
    i: iteration number, we should set this when doing multiple fits with random initializations
    '''
    # First need to update the global variables otherwise we are creating local versions
    global cfgFileLoops
    global reactionAngle
    global reactionPol
    global reactionName
    with open(reference_file,"r") as ref:
        '''
        We will be very specific on what we write to the new config file. We will
        remove all commented lines and remove all lines that deal with setting up
        the partial waves. We will also remove repeated new lines
        '''
        preamble=ref.readlines()
        preamble=[line for line in preamble if not line.startswith("#")]
        preamble=[line for line in preamble if not line.startswith("amplitude")]
        preamble=[line for line in preamble if not line.startswith("initialize")]
        preamble=[line for line in preamble if not line.startswith("constrain")]
        preamble=[line for line in preamble if not line.startswith("scale")]

        preamble=[line.rstrip().lstrip()+"-"+str(i)+"\n" if line.startswith("fit") else line for line in preamble]

        reactionLine=[line for line in preamble if line.startswith("reaction")]

        assert len(reactionLine)==1
        fitName=reactionLine[0].split(" ")[1].rstrip().lstrip()
        cfgFileLoops = any([line.startswith("loop") for line in preamble])
        if cfgFileLoops:
            loopLines=[line for line in preamble if line.startswith("loop")]
            loopLineKeys=[line.split(" ")[1] for line in loopLines]
            if "LOOPREAC" not in loopLineKeys or "LOOPPOLANG" not in loopLineKeys or "LOOPPOLVAL" not in loopLineKeys:
                raise ValueError("LOOPREAC or LOOPPOLANG or LOOPPOLVAL not defined/found. Make sure you use something like: loop LOOPREAC reac1 reac2 ...\n   in your reference config file or else the program becomes unsure what the name of the parameter should be")
            reactionName="LOOPREAC"
            reactionAngle="LOOPPOLANG"
            reactionPol="LOOPPOLVAL"
            pols = [line for line in preamble if line.startswith("loop LOOPREAC")]
            pols = pols[0].rstrip().lstrip().split(" ")[2:]
            pols = [pol.split("_")[1] for pol in pols]
            pols = "_".join(pols)
        else:
            defineLines=[line for line in preamble if line.startswith("define")]
            defineLineKeys=[line.split(" ")[1] for line in defineLines]
            if "polAngle" not in defineLineKeys or "polVal" not in defineLineKeys:
                raise ValueError("polAngle or polVal not defined/found. Make sure you use something like: define polAngle 0\n   in your reference config file or else the program becomes unsure what the name of the parameter should be")
            reactionName=fitName
            reactionAngle="polAngle"
            reactionPol="polVal"
            pols = fitName.split("_")[1].rstrip().lstrip()
            if pols not in ["000","045","090","135","AMO"]:
                raise ValueError("fitName not expected. Should be in the format myFitName_pol where myFitName is your choice and pol = {000,045,090,135,AMO}")
            pols = "_".join([pols])

        if verbose:
            print("FITNAME: "+fitName)
            print("POLARIZATIONS USED: ")
            print(pols)

        preamble="".join(preamble)
        preamble=re.sub(r'\n+','\n',preamble).strip()

    return preamble,pols

    
def writeCfg(lmes,reference_file,seed,i):
    '''
    lmes: List of lists. Each sublist is in the [l,m,e,anchor] format
    reference_file: we will get our preamble from here
    seed: set the random seed we will sample from to initialize our waveset
    i: iteration number, we should set this when doing multiple fits with random initializations
    '''
    if seed!=-1:
        random.seed(seed)
        
    preamble, pols = getPreamble(reference_file,i)

    waveStrings=[preamble]
    for lme in lmes:
        waveStrings.append(writeWave(*lme,preamble=preamble))
    
    # output config file name?
    cfgFileName=constructOutputFileName(lmes,i)
    #mapLtoSpect={0:"S",1:"P",2:"D"};
    #names=[mapLtoSpect[lme[0]]+str(lme[1])+lme[2] for lme in lmes]
    #cfgFileName="_".join(names)+"("+str(i)+").cfg"
    
    with open(cfgFileName,"w") as cfgFile:
        outputString="\n\n".join(waveStrings)
        outputString=outputString.replace("\r\n", os.linesep)
        cfgFile.writelines(outputString)
    return cfgFileName, pols
