#!/usr/bin/python3

import os
import glob
import sys
import argparse
import pandas as pd
import uproot3 as uproot
import numpy as np

def loadValue(floc,search="bestMinimum"):
    with open(floc) as f:
        lines=[line.rstrip().lstrip() for line in f.readlines()]
        line=[float(line.split("\t")[1]) for line in lines if line.split("\t")[0]==search][0]
    return line

def loadMergedPols(fileName,histBaseName,pols):
    '''
    We use amptools' plotter program to make histograms that were output into separate root files
    Use this function to load the root file and grab the histogram you want ~ integrated over some set of polarizations
    '''
    totalValues=0
    for pol in pols:
        hist=uproot.open(fileName)["EtaPi0_"+pol+"_"+histBaseName]
        edges=hist.edges
        width=edges[1]-edges[0]
        value=hist.values
        totalValues+=value
    return totalValues, edges, width

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Scour folder for PlotGen output root files')
    parser.add_argument('baseFolders', nargs="+", type=str, help='base folder to scour')
    args = parser.parse_args()
    baseFolders=args.baseFolders

    deepFolders=[]
    ## Get a list of subfolders that contain fit files
    for baseFolder in baseFolders:
        deepFolders+=glob.glob(baseFolder+"/**/*.fit",recursive=True)
    deepFolders=sorted(list(set(["/".join(deepFolder.split("/")[:-1]) for deepFolder in deepFolders])))
    deepFolders=[f for f in deepFolders if len(glob.glob(f+"/normInt*"))!=0] # the deepest folder should have normInt files in them
    print(f"Listing results for: {deepFolders}")
    
    mapMItoMDcols={
            "D1--":"D1--_pD1--",
            "D0+-":"D0+-_pD0+-",
            "D1+-":"D1+-_pD1+-",
            "D0++":"D0++_pD0++",
            "D1++":"D1++_pD1++",
            "D2++":"D2++_pD2++",
            "S0++":"S0++",
            "S0+-":"S0+-"
                    }

    ts=['010020']#,'0200325','0325050','050075','075100']
    folderHasSpecificForm=set([f.split("/")[1] for f in deepFolders])==set(ts)
    if folderHasSpecificForm:
        flocs=glob.glob(baseFolder+"/**/etapi_plot_*.root",recursive=True)
        flocs=list(set(["/".join(f.split("/")[:-1]) for f in flocs]))
        print(flocs)

        for ofile, tag in zip(["md_results.csv","md_results_finerBins.csv"],["_40MeVBin",""]):
            md_df={k:[] for k in mapMItoMDcols.keys()}
            md_df["t"]=[]
            md_df["mass"]=[]
            md_df["data"]=[]
            md_df["iteration"]=[]
            md_df["nll"]=[]
            md_df["minuitstatus"]=[]
            md_df["ematrixstatus"]=[]
            for floc in flocs:
                t, iteration=floc.split("/")[-1].split("_")
                ###################
                ## Testing purposes
                #if t!="010020" or iteration!="12": 
                #    continue
                ###################
                dataExists=False
                for k,v in mapMItoMDcols.items():
                    histacc, edgesacc, widthacc = loadMergedPols(floc+"/etapi_plot_"+v+".root","Metapi"+tag+"acc",["000","045","090","135"])
                    centers = list(np.around(edgesacc[:-1]+(edgesacc[1]-edgesacc[0])/2,5))
                    md_df[k]+=list(histacc)
                    if not dataExists:
                        md_df['mass']+=centers
                        histdat, edges, width = loadMergedPols(floc+"/etapi_plot_"+v+".root","Metapi"+tag+"dat",["000","045","090","135"])
                        histdat = list(histdat)
                        halfwidth = (edges[1]-edges[0])/2
                        dataExists=True
                md_df["data"]+=histdat
                md_df["nll"]+=[loadValue(floc+"/etapi_result.fit",search="bestMinimum")]*len(histdat)
                md_df["minuitstatus"]+=[int(loadValue(floc+"/etapi_result.fit",search="lastMinuitCommandStatus"))]*len(histdat)
                md_df["ematrixstatus"]+=[int(loadValue(floc+"/etapi_result.fit",search="eMatrixStatus"))]*len(histdat)
                md_df["t"]+=[t]*len(histdat)
                md_df["iteration"]+=[int(iteration)]*len(histdat)
            md_df=pd.DataFrame(md_df)
            md_df.to_csv(ofile,index=False,sep=" ",float_format="%.4f")



