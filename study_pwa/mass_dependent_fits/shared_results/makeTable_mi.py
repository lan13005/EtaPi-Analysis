#!/usr/bin/python3

import os
import pandas as pd
import numpy as np
import sys
import argparse

'''
The goal of this program is to extract a table of central values and associated errors from a set of mass independent fits
    if bootstrap exists it will also calculate it
'''


def loadMassIndependentData(baseFolder,waveset,masses,factor=1.0,bsTag="",keep=[True,True,True]):
    '''
    [baseFolder]: location to look for the finalAmps+[bsTag] which should contain a folder [waveset]
    [factor]: should we keep likelihoods within this factor of the best?
    [keep]: what data should we keep? only [converged, best likelihood, smallest err]
    '''
    keepConvergedOnly, keepBestLikelihoodOnly, keepSmallestErrorOnly = keep
    
    width=masses[1]-masses[0]
    masses=(masses)[:-1]+width/2
    binToMass={binNum:mass for binNum,mass in enumerate(masses)}

    datasets=[]
    totalFits=0
    filteredFits=0
    if not os.path.exists(baseFolder):
        print("folder does not exist! skipping: "+baseFolder)
        return None
    base=baseFolder+"finalAmps"+bsTag+"/"+waveset
    print("loading results in: "+base)
    for binNum in range(0,len(masses)):
        fullData=pd.read_csv(base+"/amplitudes-binNum"+str(binNum)+".txt",delimiter="\t")#,nrows=50)
        fullData["binNum"]=binNum
        totalFits+=len(fullData)

        # Convert datatypes
        ampCols=[tmp for tmp in fullData.columns if tmp not in ["status","likelihood","binNum","iteration","solution"]]
        fullData[ampCols]=fullData[ampCols].astype("float")

        # Filter results
        if keepConvergedOnly:
            fullData=fullData[fullData.status=="C"]
        if keepBestLikelihoodOnly:
            bestLikelihoods=fullData.groupby("binNum").likelihood.min()
            fullData=fullData[fullData.likelihood<=bestLikelihoods.values[0]*factor] # likelihoods within factor of the best is kept
        fullData.columns=[rearrangeWaveNotation[col] if col in rearrangeWaveNotation.keys() else col for col in fullData.columns]
        if keepSmallestErrorOnly:
            fullData=fullData[fullData["S0++_err"]==fullData["S0++_err"].min()]

        filteredFits+=len(fullData)
        datasets.append(fullData)

    print("\tt - total entries pre filtering {0}".format(totalFits))
    print("\tt - total entries post filtering {0}".format(filteredFits))
    datasets=pd.concat(datasets)
    datasets["mass"]=datasets.binNum.map(binToMass)
    datasets=datasets.reset_index(drop=True)

    return datasets


#######################################
###### Make things look prettier #######
#######################################

prettyWave={
    "S0++":"S_{0}^{+}",
    "D2-+":"D_{-2}^{+}",
    "D1-+":"D_{-1}^{+}",
    "D0++":"D_{0}^{+}",
    "D1++":"D_{1}^{+}",
    "D2++":"D_{2}^{+}",
    "S0+-":"S_{0}^{-}",
    "D2--":"D_{-2}^{-}",
    "D1--":"D_{-1}^{-}",
    "D0+-":"D_{0}^{-}",
    "D1+-":"D_{1}^{-}",
    "D2+-":"D_{2}^{-}",
    "D2-+_pD2-+":"D_{-2}^{+}",
    "D1-+_pD1-+":"D_{-1}^{+}",
    "D0++_pD0++":"D_{0}^{+}",
    "D1++_pD1++":"D_{1}^{+}",
    "D2++_pD2++":"D_{2}^{+}",
    "S0+-_pS0+-":"S_{0}^{-}",
    "D2--_pD2--":"D_{-2}^{-}",
    "D1--_pD1--":"D_{-1}^{-}",
    "D0+-_pD0+-":"D_{0}^{-}",
    "D1+-_pD1+-":"D_{1}^{-}",
    "D2+-_pD2+-":"D_{2}^{-}",
}

rearrangeWaveNotation={
    "S0+":"S0++",
    "S0-":"S0+-",
    "D-2+":"D2-+",
    "D-2-":"D2--",
    "D-1+":"D1-+",
    "D-1-":"D1--",
    "D0+":"D0++",
    "D0-":"D0+-",
    "D1+":"D1++",
    "D1-":"D1+-",
    "D2+":"D2++",
    "D2-":"D2+-",
}
rearrangeWaveNotationErr={k+"_err":v+"_err" for k,v in rearrangeWaveNotation.items()}
rearrangeWaveNotation={**rearrangeWaveNotation,**rearrangeWaveNotationErr}

#################
# Basic Setup
#################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Construct a csv output of the mass independent fit results')
    parser.add_argument('baseFolder', type=str, help='base folder to scour')
    parser.add_argument('mi_lowMass', type=float, help='lower mass')
    parser.add_argument('mi_upMass', type=float, help='upper mass')
    parser.add_argument('mi_nbins', type=int, help='number of bins')
    parser.add_argument('useBSerr_MI', type=bool, help='does bootstrap folder exist? if so we can integrate the percentiles into the output')
    args = parser.parse_args()
    mi_floc=args.baseFolder
    mi_lowMass=args.mi_lowMass
    mi_upMass=args.mi_upMass
    mi_nbins=args.mi_nbins
    useBSerr_MI=args.useBSerr_MI
    
    mi_waveset="D-1-_D0+_D0-_D1+_D1-_D2+_S0+_S0-"
    mi_waveset_vec=mi_waveset.split("_")
    mi_waveset_vec=[rearrangeWaveNotation[ele] for ele in mi_waveset_vec]
    
    ts=["010020","0200325","0325050","050075","075100"]
    midts=[0.15,0.2625,0.4125,0.625,0.875]
    tLabels=[f"0.1<t<0.2 $GeV^2$", f"0.2<t<0.325 $GeV^2$", f"0.325<t<0.5 $GeV^2$", f"0.5<t<0.75 $GeV^2$", f"0.75<t<1.0 $GeV^2$"]
    bsTag="_bs_data1x_bkgnd1x"
    
    mis=[]
    for t in ts:
        #################
        # Load Nominal results
        #################
        mi=loadMassIndependentData(mi_floc+"/"+t+"/",mi_waveset,
                                   np.linspace(mi_lowMass,mi_upMass,mi_nbins+1)
                                  ) # by default only the first best NLL with smallest errors and converged fits are saved
        mi['t']=t
    
        #################
        # Load Bootstrap results
        #################
        if useBSerr_MI:
            mi_bs=loadMassIndependentData(mi_floc+"/"+t+"/",mi_waveset,
                                          np.linspace(mi_lowMass,mi_upMass,mi_nbins+1),bsTag=bsTag,keep=[True,False,False])
            # Extract 16th, 50th, 84th percentiles and standard deviation
            lstds=mi_bs.groupby("mass")[mi_waveset_vec].quantile(0.16).reset_index()
            mstds=mi_bs.groupby("mass")[mi_waveset_vec].quantile(0.50).reset_index()
            ustds=mi_bs.groupby("mass")[mi_waveset_vec].quantile(0.84).reset_index()
            stds=mi_bs.groupby("mass")[mi_waveset_vec].std().reset_index()
            ustds[mi_waveset_vec] = ustds[mi_waveset_vec]-mstds[mi_waveset_vec]
            lstds[mi_waveset_vec] = mstds[mi_waveset_vec]-lstds[mi_waveset_vec]
            lstds=lstds.rename(columns={ele:ele+"_err"+"_bsl" for ele in mi_waveset_vec}).drop("mass",axis=1)
            ustds=ustds.rename(columns={ele:ele+"_err"+"_bsu" for ele in mi_waveset_vec}).drop("mass",axis=1)
            stds=stds.rename(columns={ele:ele+"_err"+"_bss" for ele in mi_waveset_vec}).drop("mass",axis=1)
            mi=pd.concat([mi,lstds,ustds,stds],axis=1)
    
        mis.append(mi)
    mis=pd.concat(mis)
    
    mis.to_csv("mi_results.csv",sep=" ",float_format="%.4f",index=False)


















