#!/usr/bin/python3

import os
import glob
import sys
import argparse
import pandas as pd
import uproot3 as uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 200
import mplhep as hep

mpl.rcParams.update(mpl.rcParamsDefault)
plt.style.use([hep.styles.ATLAS])

SMALL_SIZE = 20
MEDIUM_SIZE = 22
BIGGER_SIZE = 24

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=17)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare MassIndep and specific iterations of MassDep fits, draw with finer binning if possible')
    parser.add_argument('mi_csv', type=str, help='csv file containing the mass independent fit results')
    parser.add_argument('md_csv', type=str, help='csv file containing the mass dependent fit results')
    parser.add_argument('md_finer_csv', type=str, help='csv file containing the mass dependent fit results with a finer mass binning')
    parser.add_argument('md_iterations', type=str, help='underscore separated iterations of md fits that we want to draw, semicolon separated t-bins')
    args = parser.parse_args()
    mi_csv=args.mi_csv
    md_csv=args.md_csv
    md_finer_csv=args.md_finer_csv
    iterationss=args.md_iterations
    iterationss=iterationss.split(";")
    iterationss=[iterations.split("_") if "_" in iterations else [iterations] for iterations in iterationss]
    print(iterationss)

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
    }


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
    
    ts=['010020','0200325','0325050','050075','075100']
    mi_exists=mi_csv!="none"
    md_exists=md_csv!="none"

    if mi_exists and md_exists:
        # Load Mass indep fits
        mi_df_full=pd.read_csv(mi_csv,delimiter=' ',dtype={"t":str})
        nBins=mi_df_full.binNum.max()+1
        masses=mi_df_full.mass.unique()
        halfwidth=(masses[1]-masses[0])/2

        # Load Mass dep fits
        md_df_full=pd.read_csv(md_csv,delimiter=' ',dtype={"t":str})
        md_finer_df_full=pd.read_csv(md_finer_csv,delimiter=' ',dtype={"t":str})

        fig,axes=plt.subplots(5,8,figsize=(14,8),sharex=True)
        for it, t in enumerate(ts):
            drawnData={k:False for k in mapMItoMDcols.keys()}
            for iteration in iterationss[it]:
                iteration=int(iteration)
                md_df=md_df_full.loc[(md_df_full.t==t)&(md_df_full.iteration==iteration)].reset_index(drop=True)
                md_finer_df=md_finer_df_full.loc[(md_finer_df_full.t==t)&(md_finer_df_full.iteration==iteration)].reset_index(drop=True)
                binRatio=(md_df.mass[1]-md_df.mass[0])/(md_finer_df.mass[1]-md_finer_df.mass[0])
                mi_df=mi_df_full.loc[(mi_df_full.t==t)].reset_index(drop=True)

                for ik,k in enumerate(mapMItoMDcols.keys()):
                    axes[it,ik].plot(md_finer_df.mass,md_finer_df[k]*binRatio)#,label='MD')
                    axes[it,ik].axhline(0,c='black',linewidth=1)
                    
                    if not drawnData[k]:
                        #axes[it,ik].errorbar(mi_df.mass,mi_df["D2++"],yerr=(mi_df["D2++_err_bsl"],mi_df["D2++_err_bsu"]),c='black',label='MI',fmt='o')
                        axes[it,ik].errorbar(mi_df.mass,mi_df[k],yerr=mi_df[k+"_err_bss"],c='black',label='MI',fmt='o',markersize=2)
                        histdat=md_df["data"].reset_index(drop=True)
                        nonzero_idxs=np.where(histdat!=0)[0]
                        data_edges=np.array(list((md_df.mass-halfwidth).values)+[md_df.mass.max()+halfwidth])
                        hep.histplot((histdat,data_edges), ax=axes[it,ik], c='black')
                        histdat=histdat[nonzero_idxs]
                        data_edges=data_edges[nonzero_idxs+1]
                        drawnData[k]=True
                    axes[it,ik].set_xticklabels([])
                    axes[it,ik].set_yticklabels([])
                    axes[it,ik].set_ylim(0,max(histdat)*1.1)

        for ik,k in enumerate(mapMItoMDcols.keys()):
            axes[0,ik].set_title(f"${prettyWave[k]}$",size=14)

        axes[0,0].set_xlim(data_edges.min(),data_edges.max())
        #fig.text(0.5, 0.04, r"$M(\eta\pi)$ GeV", ha='center')
        #fig.text(0.04, 0.5, f"Events / {2*halfwidth:0.3f} GeV", va='center', rotation='vertical')
        axes[it,ik].legend(loc=1,prop={"size":14}) # should be the last plot
        plt.tight_layout()
        plt.subplots_adjust(wspace=0, hspace=0)
        plt.savefig("compareMIMD.pdf")
        #################
    
    
