
#!/usr/bin/python3

import glob
import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['figure.dpi'] = 200
from scipy.optimize import curve_fit
import scipy.stats as stats

import hiplot as hip
import mplhep as hep
mpl.rcParams.update(mpl.rcParamsDefault)
plt.style.use([hep.styles.ATLAS])

SMALL_SIZE = 20
MEDIUM_SIZE = 24
BIGGER_SIZE = 28

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

waves=['D+','pD+','D-','pD-','S+','S-']
labels=[r'$D^+_{a_2(1320)}$',r'$D^+_{a_2(1700)}$',r'$D^-_{a_2(1320)}$',r'$D^-_{a_2(1700)}$',r'$S^{+}$',r'$S^{-}$']
ts=['010020','0200325','0325050','050075','075100']
tlabels=['0.1<t<0.2 GeV','0.2<t<0.325 GeV','0.325<t<0.5 GeV','0.5<t<0.75 GeV','0.75<t<1.0 GeV']


for ofolder, fracOrYield in zip(['bootstrap_results_fitFrac', 'bootstrap_results'], ['frac','yield']):
    figAx_stdVsSample=list(plt.subplots(5,1,figsize=(18,24)))
    figAx_stdVsSample[1]=figAx_stdVsSample[1].flatten()
    figAx_yield=list(plt.subplots(5,1,figsize=(18,24)))
    figAx_yield[1]=figAx_yield[1].flatten()
    figAx_scaleFactor=list(plt.subplots(2,3,figsize=(18,12)))
    figAx_scaleFactor[1]=figAx_scaleFactor[1].flatten()

    os.system(f'mkdir -p {ofolder}')
    
    maxX=-1
    for it,t,tlabel in zip(range(len(ts)),ts,tlabels):
        ####################################
        #### LOAD BOOTSTRAP FIT RESULTS ####
        ####################################
        folders=glob.glob(f'results_t{t}/results_t{t}*')
        prefix=''
        value=''
        prefixes=[]
        values=[]
        for i,folder in enumerate(folders):
            fname=f'{folder}/etapi_plotter_output.log'
            print(f'loading: {fname}')
            with open(fname) as f:
                lines=f.readlines()
                assert(len(lines)==138) # this is very specific to my current setup!
                for line in lines:
                    line=line.split(' ')
                    prefix, value = ' '.join(line[:-3]), float(line[-3])
                    if i==0:
                        prefixes.append(prefix)
                    if prefix.startswith("TOTAL EVENTS"):
                        total=float(line[3])
                        values.append(value)
                    else:
                        if fracOrYield=='yield':
                            values.append(value*total)
                        else:
                            values.append(value)
        values=np.array(values)
        if len(values)>len(prefixes): 
            values=values.reshape(-1,len(prefixes))
        means=values.mean(axis=0)
        uncertainties=values.std(axis=0)
        
        lines=[] 
        bootstrap_uncerts={}
        for prefix, mean, uncertainty in zip(prefixes,means,uncertainties):
            lines.append(f'{prefix} {mean:0.8f} +- {uncertainty:0.8f}')
            for wave in waves:
                if f' {wave} ' in prefix:
                    bootstrap_uncerts[wave]=uncertainty
        lines='\n'.join(lines)
        
        with open(f'{ofolder}/bootstrap_result_mean_std_t{t}.log','w') as f:
            f.write(lines)
            f.close()
    
        ####################################
        #### LOAD NOMINAL FIT RESULTS (TO COMPARE MINUIT/BOOTSTRAP)
        ####################################
        nominal_floc=f'/d/grid17/ln16/dselector_v3/study_pwa/mass_dependent_fits/shared_results/systematic_nominal_v9/nominal/{t}/{t}_0/etapi_plotter_output.log'
        nominal_uncerts={wave:-1 for wave in waves}
        with open(nominal_floc) as f:
            lines=f.readlines()
            assert(len(lines)==138)
            for line in lines:
                if "TOTAL EVENTS" in line:
                    total=float(line.split(' ')[3])
                for wave in waves:
                    if f' {wave} ' in line:
                        if fracOrYield=='yield':
                            uncert=float(line.split(' ')[-1])*total
                        else:
                            uncert=float(line.split(' ')[-1])
                        nominal_uncerts[wave]=uncert
        print('\n***************************')
        print(f"nominal fit: {nominal_floc}")
        for k,v in nominal_uncerts.items():
            print(f't{t}: {k} has bootstrapped/nominal uncertainity {bootstrap_uncerts[k]:0.8f}/{v:0.8f} = {bootstrap_uncerts[k]/v:0.8f}')
        print('***************************\n')
    
        ####################################
        #### Diagnostic plots as a function of Number of samples
        ####################################
        ax=figAx_stdVsSample[1][it]
        samples=range(2,len(values))
        for wave,label in zip(waves,labels):
            wave=f' {wave} ' # need spaces to distinguish it from other prefixes
            uncert_vs_sample=[]
            j=[i for i,v in enumerate(prefixes) if wave in v][0]
            for i in samples:
                uncert_vs_sample.append(values[:i,j].std(axis=0))
            ax.plot(samples,uncert_vs_sample,label=label)
            if maxX<len(values):
                maxX=len(values)
        ax.set_ylabel(f"Std() of {fracOrYield}")
        ax.set_xlabel("Number of Samples")
        if it==0:
            ax.legend()
        if it==len(ts)-1:
            for ax in figAx_stdVsSample[1]:
                ax.set_xlim(0,maxX*1.15)
            figAx_stdVsSample[0].tight_layout()
            figAx_stdVsSample[0].savefig(f'{ofolder}/nsamples_diagnostic.png')
    
        ####################################
        #### Histogram of all bootstrapped samples
        ####################################
        #fig,axes=plt.subplots(2,3,figsize=(18,12))
        #axes=axes.flatten()
        #for ax,wave in zip(axes,waves):
        #    wave=f' {wave} ' # need spaces to distinguish it from other prefixes
        #    j=[i for i,v in enumerate(prefixes) if wave in v][0]
        #    hep.histplot(np.histogram(values[:,j],bins=30),c='black',density=True,ax=ax)
        #    ax.set_title(f'{wave} std={values[:,j].std():0.5f}')
        #[axes[i].set_ylabel("Density") for i in [0,3]]
        #[axes[i].set_xlabel("Yield") for i in [3,4,5]]
    

        ax=figAx_yield[1][it]
        for wave,label in zip(waves,labels):
            wave=f' {wave} ' # need spaces to distinguish it from other prefixes
            j=[i for i,v in enumerate(prefixes) if wave in v][0]
            hep.histplot(np.histogram(values[:,j],bins=30),ax=ax, label=label, histtype='fill', alpha=0.6)
            hep.histplot(np.histogram(values[:,j],bins=30),ax=ax)
            if fracOrYield=='yield':
                ax.set_xlabel('Acceptance Corrected Yields')
            else:
                ax.set_xlabel('Fit Fraction')
            ax.set_ylabel('Density')
        ax.set_ylim(0)
        ax.set_title(tlabel)
        if it==0:
            ax.legend()
        if it==len(ts)-1:
            figAx_yield[0].tight_layout()
            figAx_yield[0].savefig(f'{ofolder}/bs_distribution_diagnostic.png')
    
        ax=figAx_scaleFactor[1][it]
        ####################################
        #### Histogram of all bootstrapped samples
        ####################################
        scaleFactors=[]
        for wave in waves:
            scaleFactors.append(bootstrap_uncerts[wave]/nominal_uncerts[wave])
        ax.scatter(range(len(waves)),scaleFactors,c='black',marker='o')
        ax.set_title(tlabel)
        ax.set_xticks(range(len(waves)))
        ax.set_xticklabels(labels, rotation=45) 
        ax.set_xlim(-0.5,len(waves)-0.5)
        ax.set_ylabel("Scale Factor Bootstrap/Minuit",size=16)
        ax.set_xlabel("Coherent Sum")
        if it==len(ts)-1:
            figAx_scaleFactor[0].tight_layout()
            figAx_scaleFactor[0].savefig(f'{ofolder}/bs_minuit_scaleFactor.png')

    








