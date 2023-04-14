#!/usr/bin/python3

from helper import *
import re
import math
import numpy as np
import time


'''
DID NOT COMPLETE THIS PROGRAM! USE THE JUPYTER NOTEBOOK
'''


#################
# SETUP - LOAD DATA
#################
#parser = argparse.ArgumentParser(description='Draw diagnostic plots')
#parser.add_argument('folder', help='Location of fit results file to grab production coefficients from')
#args = parser.parse_args()
#floc=args.floc

folder='systematic_nominal_v9/precision_of_NLL/v1/nominal/'
waves=["D2-", "D1-", "D0+", "D1+", "D2+","pD2-", "pD1-", "pD0+", "pD1+", "pD2+"] # will tack on reflectivity later


prettyMap={
        "D2--":r'$D_{-2}$',
        "D1--":r'$D_{-1}$',
        "D0+-":r'$D_{0}$',
        "D1+-":r'$D_{1}$',
        "D2+-":r'$D_{2}$',
        "pD2--":r'$D_{-2}$',
        "pD1--":r'$D_{-1}$',
        "pD0+-":r'$D_{0}$',
        "pD1+-":r'$D_{1}$',
        "pD2+-":r'$D_{2}$',
        "D2-+":r'$D_{-2}$',
        "D1-+":r'$D_{-1}$',
        "D0++":r'$D_{0}$',
        "D1++":r'$D_{1}$',
        "D2++":r'$D_{2}$',
        "pD2-+":r'$D_{-2}$',
        "pD1-+":r'$D_{-1}$',
        "pD0++":r'$D_{0}$',
        "pD1++":r'$D_{1}$',
        "pD2++":r'$D_{2}$',
        }

input_ts=[]
input_fits=[]
for it,t in enumerate(ts):
    for fit in glob.glob(f'{folder}/{t}/etapi_result_*.fit'):
        input_fits.append(fit)
        input_ts.append(midts[it])

print(" **** SOURCE FILES ****")
for fit,t in zip(input_fits,input_ts):
    print(f't={t}: {fit}')
print()

def loadResults(fitFile,dump=False):
    with open(fitFile) as f:
        lines=f.readlines()
        cfgFirstLine=[i for i,v in enumerate(lines) if "## FIT CONFIGURATION ###" in v][0]
        fitResultLine=[i for i,v in enumerate(lines) if "+++ Parameter Values and Errors +++" in v][0]
        cfg=lines[cfgFirstLine+1:]
        fitresults=lines[fitResultLine:cfgFirstLine+1]
        fitresults=[re.sub(r" +","",line).rstrip().lstrip() for line in fitresults]
        # make sure the line is a parameterName and value that is separated by a tab
        fitresults=[line.split("\t") for line in fitresults if (len(line.split("\t"))==2)&(~line[0].isdigit())] 
        fitMap={k:float(v) for k,v in fitresults}
        if dump:
            for k,v in fitMap.items():
                print(f'{k} {v}')
        return fitMap

def extractRTheta(fitMap,dump=False):
    # prototype to search for: EtaPi0_000::NegativeRe::D0+-_re
    reaction="EtaPi0_000"
    part="Re"
    refls=["Positive","Negative"]
    phaseMap={}
    magMap={}
    for refl in refls:
        e='+' if refl=='Positive' else '-'
        for wave in waves:
            prefix=f'{reaction}::{refl}{part}::{wave}{e}_'
            if prefix+'re' not in fitMap.keys():
                magMap[wave+e]=0
                phaseMap[wave+e]=0
            else:
                cartesian=fitMap[prefix+'re']+1j*fitMap[prefix+'im']
                r, theta = abs(cartesian), np.angle(cartesian)
                magMap[wave+e]=r
                phaseMap[wave+e]=theta
        if dump:
            for map1, label in zip([magMap,phaseMap],['Magnitudes','Phases']):
                print(f' ***** {label} *****')
                for k,v in map1.items():
                    print(f'{k}: {v}')
                print()
    return magMap, phaseMap


def mergeMaps(listOfMaps):
    for i in range(1,len(listOfMaps)):
        assert( listOfMaps[0].keys()==listOfMaps[i].keys() )
    merged={}
    for map1 in listOfMaps:
        for k,v in map1.items():
            if k not in merged.keys():
                merged[k]=[v]
            else:
                merged[k].append(v)
    return merged

def plotData(x,y,ax,c):
    x=np.array(x)
    y=np.array(y)

    # i refers to the starting/reference point
    for i in range(math.ceil(len(x)/2)):
        # j is the candidate nearest point(forward/backward)
        ## LOOP 1: Goes in the forward direction
        for j in range(i,math.ceil(len(x)/2)):
            if x[j]>x[i]:
                ax.plot([x[i],x[j]],[y[i],y[j]], linestyle='-', c=c)
                break
        ## LOOP 2: Goes in the backward direction.
        for j in range(i,math.ceil(len(x)/2)):
            ii=len(x)-i-1
            jj=len(x)-j-1
            if x[ii]>x[jj]:
                ax.plot([x[ii],x[jj]],[y[ii],y[jj]], linestyle='-', c=c)
                break

magMaps=[]
phaseMaps=[]
for i in range(len(input_fits)):
    fitMap=loadResults(input_fits[i])
    magMap, phaseMap = extractRTheta(fitMap)
    magMaps.append(magMap)
    phaseMaps.append(phaseMap)

merged=mergeMaps(phaseMaps)

fig,ax=plt.subplots(2,2,figsize=(20,15))

cs=plt.rcParams['axes.prop_cycle'].by_key()['color']
cmap={}
ic=0
fitIterations=[int(fit.split('.')[0].split('_')[-1]) for fit in input_fits]
for k,v in merged.items():
    reducedAmp=k.lstrip('p')[:3]
    if reducedAmp not in cmap.keys():
        cmap[reducedAmp]=cs[ic]
    icol = 1 if k.startswith('p') else 0
    irow = 1 if k.endswith('-') else 0
    plotData(input_ts,v,ax[irow,icol],c=cmap[reducedAmp])

    kwarg_map={0:{'marker':'o','facecolor':cmap[reducedAmp]}, 1:{'marker':'o','facecolor':'none'}}
    for x,y,i,p in zip(input_ts,v,fitIterations,range(len(input_ts))):
        label=prettyMap[k] if p==0 else ''
        ax[irow,icol].scatter(x,y,edgecolor=cmap[reducedAmp],linewidths=3,label=label,s=80, **kwarg_map[i])
    ic+=1

for irow, rlabel in enumerate(['Pos. Ref.','Neg. Ref.']):
    for icol, clabel in enumerate([r'$a_2(1320)$',r'$a_2(1700)$']):
        ax[irow,icol].set_xlabel(r'$-t$ $GeV^2$',size=24)
        ax[irow,icol].set_ylabel('Production Coef. Phase',size=24)
        ax[irow,icol].set_title(f'{clabel} {rlabel}',size=30)
ax[0,1].legend(prop={'size':36},bbox_to_anchor=(1.05,1))

plt.tight_layout()
plt.savefig('drawProductionCoeffs.png')




















