
#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib as mpl
import uproot3 as uproot
import os
import glob
import seaborn as sns

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


def getValues(floc):
    with open(floc) as f:
        lines=f.readlines()
        bestMin=-1
        minuitStatus=-1
        ematrixStatus=-1
        a2mass=-1
        a2width=-1
        potentialSeeds=[]
        for line in lines:
            line=line.lstrip()
            if line.startswith("bestMinimum\t"):
                bestMin=float(line.split('\t')[-1].rstrip())
            if line.startswith("eMatrixStatus\t"):
                ematrixStatus=int(line.split('\t')[-1].rstrip())
            if line.startswith("lastMinuitCommandStatus\t"):
                minuitStatus=int(line.split('\t')[-1].rstrip())
            if line.startswith("a2mass\t"):
                a2mass=float(line.split('\t')[-1].rstrip())
            if line.startswith("a2width\t"):
                a2width=float(line.split('\t')[-1].rstrip())
            for src in ['data','bkgnd','accmc']:
                if line.startswith(f"{src} ") and "ROOTDataReaderBootstrap" in line:
                    potentialSeeds.append(int(line.split(" ")[4].rstrip()))
        potentialSeeds=list(set(potentialSeeds))
        assert(len(potentialSeeds)==1)
    return bestMin, ematrixStatus, minuitStatus, a2mass, a2width, potentialSeeds[0]


def loadFractions(flocs,wave,baseDir,tag=''):
    values=[]
    errors=[]
    for floc in flocs:
        iteration = floc.split("_")[-1].split('.')[0]
        floc=f'{baseDir}/results_{tag}{iteration}/etapi_plotter_output.log'
        if not os.path.exists(floc):
            continue
        with open(floc) as f:
            lines=f.readlines()
            lines=[line for line in lines if line.startswith(f"FIT FRACTION {wave}")]
            values.append(float(lines[0].split(" ")[4].rstrip()))
            errors.append(float(lines[0].split(" ")[6].rstrip()))
    return np.array(values), np.array(errors)

def loadNominal(floc):
    # i = t-bin
    waves=[
        'S0++','S0+-',
        'D0++','D1++','D2++','D1--','D0+-','D1+-',
        'pD0++','pD1++','pD2++','pD1--','pD0+-','pD1+-',
        'D+','pD+','D-','pD-',
        'S+','S-'
    ]
    with open(floc) as f:
        lines=f.readlines()
    df={}
    for wave in waves:
        line=[line for line in lines if line.startswith(f"FIT FRACTION {wave}")][0]
        value=float(line.split(" ")[4].rstrip())
        err=float(line.split(" ")[6].rstrip())
        df[wave]=[value]
        df[wave+'_err']=[err]
    yields=[line for line in lines if line.startswith("TOTAL EVENTS")][0].split(" ")
    totalYield=float(yields[3])
    totalYieldErr=float(yields[5])
    df['totalYield']=totalYield
    df['totalYieldErr']=totalYieldErr
    return pd.DataFrame(df)


def loadCorrectedYields(flocs,baseDir,tag):
    values=[]
    errors=[]
    for floc in flocs:
        iteration = floc.split("_")[-1].split('.')[0]
        floc=f'{baseDir}/results_{tag}{iteration}/etapi_plotter_output.log'
        if not os.path.exists(floc):
            continue
        with open(floc) as f:
            lines=f.readlines()
            lines=[line for line in lines if line.startswith(f"TOTAL EVENT")]
            values.append(float(lines[0].split(" ")[3].rstrip()))
            errors.append(float(lines[0].split(" ")[5].rstrip()))
    return np.array(values), np.array(errors)


def loadInfo(baseDir,tag):
    # i = t-bin
    waves=[
        'S0++','S0+-',
        'D0++','D1++','D2++','D1--','D0+-','D1+-',
        'pD0++','pD1++','pD2++','pD1--','pD0+-','pD1+-',
        'D+','pD+','D-','pD-',
        'S+','S-'
    ]

    minimas=[]
    ematrixStatuses=[]
    minuitStatuses=[]
    convergences=[]
    seeds=[]
    a2masses=[]
    a2widths=[]
    files=[f'{baseDir}/etapi_result_{i}.fit' for i in range(100)]
    for i,f in enumerate(files):
        if not os.path.exists(f):
            continue
        values=getValues(f)
        minimas.append(values[0])
        ematrixStatuses.append(values[1])
        minuitStatuses.append(values[2])
        convergences.append(values[1]==3 and values[2]==0)
        a2masses.append(values[3])
        a2widths.append(values[4])
        seeds.append(values[5])

    minimas=np.array(minimas)
    ematrixStatuses=np.array(ematrixStatuses)
    minuitStatuses=np.array(minuitStatuses)
    convergences=np.array(convergences)
    a2masses=np.array(a2masses)
    a2widths=np.array(a2widths)
    files=np.array(files)

    df={}
    df['minuit']=minuitStatuses
    df['ematrix']=ematrixStatuses
    df['converged']=convergences
    df['NLL']=minimas#[convergences]
    df['a2mass']=a2masses
    df['a2width']=a2widths
    df['seed']=seeds

    values,errors=loadCorrectedYields(files,baseDir,tag)
    df['totalYield']=values
    df['totalYieldErr']=errors

    for wave in waves:
        value,err=loadFractions(files,wave,baseDir,tag)#[convergences],wave)
        df[wave]=value
        df[wave+'_err']=err

    return pd.DataFrame(df)

def drawDiagnostics(df1,df2,nominalDelta,additionalTag=''):
    '''
    df1 should be alternative worse fit
    df2 should be the better nominal fit
    '''
    x=df1
    y=df2
    print(f'There are {len(set(x.seed.values))} unique seeds')
    assert(np.prod(x.seed.values==y.seed.values))
    x=x['NLL'].values
    y=y['NLL'].values
    
    #### DIFFERENCE PLOT
    plt.figure(figsize=(10,8))
    fig,ax=plt.subplots(1,1,figsize=(12,6))
    hep.histplot(np.histogram(x-y,bins=50),c='black')
    plt.xlabel(r"$\Delta(NLL)$ paired on seed")
    plt.ylabel("Counts")
    plt.axvline((x-y).mean(),c='red',linestyle='--',label='mean')
    fracLT0=100.0*sum((x-y)<0)/len(x)
    #plt.title(f'mean(std)[%<0]:{(x-y).mean():0.2f}({(x-y).std():0.2f})[{fracLT0:0.0f}%]',size=30)
    plt.title(f'mean(std):{(x-y).mean():0.2f}({(x-y).std():0.2f})',size=30)
    plt.legend(prop={'size':30})
    plt.savefig(f'diffPlot_delta{nominalDelta}{additionalTag}.pdf')
    
    #############
    
    #### CUMULATIVE DISTRIBUTION SHIFTED
    fractions=[]
    thresholds=np.linspace(0,100,200)
    for threshold in thresholds:
        fractions.append(sum(x-y-(x-y).mean()+threshold<0)/len(x))
    fractions=np.array(fractions)
    
    fig,ax=plt.subplots(1,1,figsize=(12,6))
    plt.plot(thresholds,fractions,c='black')
    plt.ylim(0)
    plt.axhline(0.05,c='red',linestyle='--',label='Fraction=0.05')
    plt.axvline(nominalDelta,c='blue',linestyle='--',label=r'Nominal $\Delta$NLL')
    plt.legend(prop={'size':30})
    plt.ylabel("Fraction Less Than Zero")
    plt.xlabel(r"$\Delta(NLL)$ For Hypothetical Nominal Fits")
    plt.savefig(f'cumulative_delta{nominalDelta}.pdf')

ifolder='/d/grid17/ln16/dselector_v3/study_pwa/mass_dependent_fits/shared_results/systematic_nominal_v9/precision_of_NLL/'

#####################
#### PILOT STUDY ####
#####################

## Technically we only bootstrap {data,bkgnd} here
#baseDir=f'{ifolder}/bootstrap_test/results_data_bkgnd/'
#df_v0_data=loadInfo(baseDir,'')
#
## Technically we only bootstrap {data,bkgnd} here
#baseDir=f'{ifolder}/bootstrap_test_altMinimum/results_data_bkgnd/'
#df_v0_altMin_data=loadInfo(baseDir,'')
#
## Technically we only bootstrap {data,bkgnd} AND {accmc} here
#baseDir=f'{ifolder}/bootstrap_t010020_altMin/results_t010020/'
#df_v1_altMin=loadInfo(baseDir,'t010020_')
#
## Technically we only bootstrap {data,bkgnd} AND {accmc} here
#baseDir=f'{ifolder}/results_t010020/'
#df_v1=loadInfo(baseDir,'t010020_')

#drawDiagnostics(df_v0_data,df_v0_altMin_data,100)
#drawDiagnostics(df_v1,df_v1_altMin,5)

#####################
#### V1 STUDY ####
#####################
baseDir=f'{ifolder}/v1/010020_delta5NLL/bootstrap_altMin/results/'
df_010020_altMin=loadInfo(baseDir,'')
baseDir=f'{ifolder}v1/010020_delta5NLL/bootstrap_best/results/'
df_010020_best=loadInfo(baseDir,'')

baseDir=f'{ifolder}/v1/075100_delta9NLL/bootstrap_altMin/results/'
df_075100_altMin=loadInfo(baseDir,'')
baseDir=f'{ifolder}v1/075100_delta9NLL/bootstrap_best/results/'
df_075100_best=loadInfo(baseDir,'')

drawDiagnostics(df_010020_altMin,df_010020_best,5,'_t010020')
drawDiagnostics(df_075100_altMin,df_075100_best,9,'_t075100')














