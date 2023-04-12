#!/usr/bin/python3

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import uproot3 as uproot
import argparse
from multiprocessing import Pool
import itertools
import time

import matplotlib as mpl
import mplhep as hep
mpl.rcParams.update(mpl.rcParamsDefault)
plt.style.use([hep.styles.ATLAS])

SMALL_SIZE = 16
MEDIUM_SIZE = 20
BIGGER_SIZE = 24

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


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


def plotWaves(floc,variable,selectRef,resonances,histtype,ftype,variationsMap):
    print(f" ****** DRAWING PLOTS IN {floc} / {variable} ******* ")
    selectRef=selectRef if selectRef else [0,1]
    resonances=resonances if resonances else ["","p"]
    histtype=histtype if histtype else "fill"
    variable=variable if variable else "mass"

    assert( variable in ["mass","cos"] )
    assert( histtype in ["fill","step"] )
    assert( resonances in [["","p"],[""],["p"],["both"]] )
    assert( selectRef in [[0],[1],[0,1]] )
    assert( ftype in ["png","pdf"] )

    if variable=="mass":
        xvar="Metapi" 
    else:
        xvar="cosTheta"

    floc=floc+"/" if floc[-1]!="/" else floc
    variation=variationsMap[int(floc.split('_')[-1].rstrip('/'))]

    wavesets=np.array([["S0+-","D1--","D0+-","D1+-"],["S0++","D0++","D1++","D2++"]])
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

    fileName="etapi_plot_"
    wavesets=wavesets[selectRef]
    refs=np.array(["negativeRef","positiveRef"])[selectRef]
    naturalities=np.array(["Unnatural\nProduction","Natural\nProduction"])[selectRef]
    if variable=="mass":
        xlabel=r"$M(\eta\pi)$ $GeV^2$"
    elif variable=="cos":
        xlabel=r"$cos\theta$"
    ylabel=r"Entries / 40 MeV"
            
#    ########################################################
#    ########################################################
    #completeWaveset="S0+-_S0++_D1--_D0+-_D1+-_D0++_D1++_D2++_pD1--_pD0+-_pD1+-_pD0++_pD1++_pD2++"
    completeWaveset="all"
    fig,axes=plt.subplots(3,3,figsize=(16,12))
    axes=axes.flatten()
    xvars=["Metapi_40MeVBin","cosTheta","phi","Phi","psi","t","Mpi0proton","Metaproton"]
    xlabels=[r"$M(\eta\pi) GeV$",r"$cos\theta$",r"$\phi$ rad",r"$\Phi$ rad",r"$\psi$ rad",r"$-t$ $GeV^2$",r"$M(\pi p) GeV$",r"$M(\eta p) GeV$"]
    for ax,_xvar,_xlabel in zip(axes,xvars,xlabels):
        histdat, edges, width = loadMergedPols(floc+"/"+fileName+completeWaveset+".root",_xvar+"dat",["000","045","090","135"])
        histbkg, _, _ = loadMergedPols(floc+"/"+fileName+completeWaveset+".root",_xvar+"bkg",["000","045","090","135"])
        histacc, edgesacc, widthacc = loadMergedPols(floc+"/"+fileName+completeWaveset+".root",_xvar+"acc",["000","045","090","135"])
        values=histdat-histbkg
        binScaleFactor=width/widthacc
        valuesacc=histacc*binScaleFactor

        xmin=edges[:-1][histdat>0][0] # edges has 1 more dimension than the counts, need to account for that
        xmax=edges[1:][histdat>0][-1]
        hep.histplot(values,edges,c='black',ax=ax)
        hep.histplot(valuesacc,edgesacc,color="forestgreen",ax=ax,alpha=0.7,linewidth=2,histtype="fill")

        ax.set_xlim(xmin,xmax)
        ax.set_ylim(0)
        ax.set_xlabel(_xlabel,size=24)
        ax.set_ylabel(f"Entries / {width:0.3f}",size=24)

    plt.tight_layout()
    fig.subplots_adjust(top=0.90)
    st = fig.suptitle(f'Variation: {variation}', fontsize=30)
    st.set_y(0.95)
    plt.savefig(floc+"/summed_waves."+ftype)

    ########################################################
    # Begin drawing plots we have selected
    ########################################################
    md_curves=[]
    for iws,waveset,ref,naturality in zip(range(len(wavesets)),wavesets,refs,naturalities):
        fig,axes=plt.subplots(2,2,figsize=(14,10),sharex=True,sharey=True)
        axes=axes.flatten()
        for prime in resonances:
            for iw,wave in enumerate(waveset):
                ############### Load the data to be plotted - ok to overwrite histdat and histacc now
                if wave[0]=="D":
                    selectWave=prime+wave if prime in ["","p"] else wave+"_p"+wave
                    if prime=="":
                        histLabel=f"$a_2(1320)$ ${prettyWave[wave]}$"
                    elif prime=="p":
                        histLabel=f"$a_2(1700)$ ${prettyWave[wave]}$"
                    else:
                        histLabel=r"$a_2/a_2(1700)$"+f" ${prettyWave[wave]}$"
                else: #All non-D waves
                    selectWave=prime+wave if prime in [""] else wave
                    histLabel=f"${prettyWave[wave]}$"
                    
                rootFileName=fileName+selectWave+".root"

                if not os.path.exists(floc+rootFileName):
                    print(f"file not found {floc+rootFileName}")
                    continue
                _xvar = xvar+"_40MeVBin" if variable=="mass" else xvar
                histdat,edges,width=loadMergedPols(floc+rootFileName,_xvar+"dat",["000","045","090","135"])
                histbkg,_,_=loadMergedPols(floc+rootFileName,_xvar+"bkg",["000","045","090","135"])
                values=histdat-histbkg
                histacc,edgesacc,widthacc=loadMergedPols(floc+rootFileName,xvar+"acc",["000","045","090","135"])
                binScaleFactor=width/widthacc
                valuesacc=histacc*binScaleFactor

                xmin=edges[:-1][histdat>0][0] # edges has 1 more dimension than the counts, need to account for that
                xmax=edges[1:][histdat>0][-1]

                hep.histplot(values,edges,c='black',ax=axes[iw])
                    
                if (prime=="" or prime=="both")*(wave[0]=="D"):
                    color="orange"
                elif (prime=="p")*(wave[0]=="D"):
                    color="royalblue"
                else:
                    color='darkgray'    
                
                _label = histLabel
                _alpha=0.8
                _linewidth=2
                hep.histplot(valuesacc,edgesacc,color=color,ax=axes[iw],alpha=_alpha,label=_label,linewidth=_linewidth,histtype=histtype)

                axes[iw].legend(prop={"size":18},loc=1)#,bbox_to_anchor=(1,0.95))
                axes[iw].set_xlim(xmin,xmax)
                axes[iw].set_ylim(0,axes[iw].get_ylim()[1]*1.02)

        axes[2].set_xlabel(xlabel,size=30)
        axes[0].set_ylabel(ylabel,size=30)
        axes[0].yaxis.set_label_coords(-0.15, 0.1)
        axes[2].xaxis.set_label_coords(1.05, -0.1)
        
        plt.text(0.75, 0.62, naturality, fontsize=18, transform=axes[0].transAxes, color="red", weight="bold")
        plt.tight_layout()
        fig.subplots_adjust(top=0.85,wspace=0, hspace=0)
        st = fig.suptitle(f'Variation: {variation}', fontsize=30)
        st.set_y(0.95)
        plt.savefig(floc+"/"+variable+"_"+ref+"."+ftype)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Draw diagnostic plots')
    parser.add_argument('floc', help='folder location containing the etapi_plotter root files')

    refs=[0,1] # 'where 0,1 represents -,+ reflectivies to plot. i.e. [0,1] will plot both [1] will just plot + ref
    resonances=["","p"] # where ""=a2 or "p"=a2prime or set resonances=["both"] to merge (for mass-indep comparison)
    histtype='fill' # fill or step determines how to draw histograms. Step is a normal histogram with no fill
    recurse=True # recursively descend through folder and draw all plots. Searches folder for any folder that contains atleast 1 root file with right format
    ftype="png" # ["png","pdf"] to save the final images as. pdf gives an unimportant warning.
    maxthreads=50

    args = parser.parse_args()
    floc=args.floc

    time.sleep(2)

    variationsList=[]
    with open(f'{floc}/variationsList.log','r') as variationsFile:
        for line in variationsFile.readlines():
            if line[0]!='#':
                variationsList.append(line) 
    nvariations=len(variationsList)
    print(f"Variations in variationsList.log:\n{variationsList}")
    variationsMap={i:variation for i,variation in enumerate(variationsList)}

    if recurse:
        def star_plotWaves(args):
            return plotWaves(*args) # gotta expand the arguments if using map

        flocs=glob.glob(floc+"/**/etapi_plot_*.root",recursive=True)
        flocs=list(set(["/".join(f.split("/")[:-1]) for f in flocs]))
        print(len(flocs))
        print(nvariations)
        args=list(itertools.product(flocs,["mass","cos"])) 
        args=[list(arg)+[refs,resonances,histtype,ftype,variationsMap] for arg in args]
        nthreads=len(args) if len(args)<maxthreads else maxthreads
        print(flocs)
        print(f"\nSpawning {nthreads} threads to draw all {len(flocs)} the plots faster...\n")
        with Pool(nthreads) as p:
            p.map(star_plotWaves,args)
        os.system(f"mkdir -p {floc}/gathered_results")
        flocs=sorted(flocs)
        filessToCopy=[glob.glob(f+"/*."+ftype) for f in flocs]

        # will try and minimize the suffixes needed but if there are any duplicates we will take a longer version which should have no duplicates
        #suffixes=["_".join(f.split("/")[1:-1]) for f in flocs] # LONG NAMING SUFFIX 
        suffixes=["_".join(f.split("/")[-1:]) for f in flocs] # SHORT NAMING SUFFIX. The t-bin is already included in the name. Dont include in the suffix also
        if len(suffixes)!=len(set(suffixes)): # CHECK IF WE ARE LOSING ANY FILES DUE TO NAMING COLLISIONS
            print("LOSING FILES DUE TO NAMING COLLISIONS! RESORTING TO USING ~FULL PATH")
            suffixes=["_".join(f.split("/")[1:]) for f in flocs] # We can include basically the full path to the file
        assert(len(suffixes)==len(set(suffixes)))

        for suffix, filesToCopy in zip(suffixes, filessToCopy):
            for fileToCopy in filesToCopy: 
                in_f=fileToCopy
                out_f=floc+"/gathered_results/"+fileToCopy.split("/")[-1].split(".")[0]+"_"+suffix+"."+ftype
                #print(f"cp {in_f} {out_f}")
                os.system(f"cp -f {in_f} {out_f}")

    else:
        plotWaves(floc,"mass",refs,resonances,histtype,ftype)
        plotWaves(floc,"cos",refs,resonances,histtype,ftype)

    print(" ---------------------------------------------------------------")
    print('NOTE! Only folders with etapi_plotter root files will have been found!')







