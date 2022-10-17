#!/usr/bin/python3

import os
import re
import numpy as np
from getYieldsFromConfig import getYield

veryNegNum=-999
veryPosNum=999

class getNominal:
    def __init__(self,isgen):
        self.isgen=isgen
        self.nominal={
                'Mpi0eta':[1.04,1.72],
                'pVH':[0.5,veryPosNum],
                'unusedEnergy':[veryNegNum,0.01],
                'chiSq':[veryNegNum,13.277],
                '!photonTheta1':[veryNegNum,2.5,10.3,11.9],
                '!photonTheta2':[veryNegNum,2.5,10.3,11.9],
                '!photonTheta3':[veryNegNum,2.5,10.3,11.9],
                '!photonTheta4':[veryNegNum,2.5,10.3,11.9],
                'photonE1':[0.1,veryPosNum],
                'photonE2':[0.1,veryPosNum],
                'photonE3':[0.1,veryPosNum],
                'photonE4':[0.1,veryPosNum],
                'proton_momentum':[0.3,veryPosNum],
                'proton_z':[52,78],
                'mmsq':[-0.05,0.05],
                }
        self.nominal_gen={
                'Mpi0eta_thrown':[1.04,1.72],
                }

        self.resetToNominal()

    def resetToNominal(self):
        if self.isgen:
            self.selections=self.nominal_gen
        else:
            self.selections=self.nominal

    def vary(self,vs):
        for v in vs:
            self.selections[v[0]]=v[1:]
            # if the selection is chosen to take the maximum range defined by the very large numbers
            #    that implies we do not wish to make a selection... Delete the key
            if v[1]==veryNegNum and v[2]==veryPosNum and v[0] in self.selections:
                del self.selections[v[0]]

    def getCondition(self):
        condition=''
        for k,v in self.selections.items():
            nrepeats=int(len(v)/2)
            assert len(v)%2==0
            for i in range(nrepeats):
                condition+=f" {k} {v[2*i]} {v[2*i+1]}"
        return condition


#############################################################
#### THIS PROGRAM ONLY CHANGES THE ROOTDATREADERFILTER LINE 
#############################################################
def setVariation(vs):
    ts=["010020"]#,"0200325","0325050","050075","075100"]
    
    ofiles=[]
    for t in ts:
        fname=f"{t}/etapi_result.fit"
        with open(fname,'r') as f:
            srclines=f.readlines()
            filterIdxs={l:i for i,l in enumerate(srclines) if "ROOTDataReaderFilter" in l}
        
        ################################################## 
        ### set the variation, i.e. nominal or some other systematic
        ################################################## 
        lines=srclines
        for k,v in filterIdxs.items():
            if any([data in k for data in ['data','bkgnd','accmc']]):
                manager=getNominal(False)
                manager.vary(vs)
                condition=manager.getCondition()
                lines[v] = ' '.join(k.split(' ')[:4])+condition+'\n'
            if 'gen' in k:
                manager=getNominal(True)
                #manager.vary(vs) # we generally will not vary genmc selections/cuts
                condition=manager.getCondition()
                lines[v] = ' '.join(k.split(' ')[:4])+condition+'\n'

        lines="".join(lines)

        ofile=f"{t}/etapi_result_src.fit"
        with open(ofile,"w") as out:
            #print(f"writing to {ofile}")
            out.write(lines)
        ofiles.append(ofile)
    return ofiles

###########################
#### SET OUTPUT FOLDER ####
ofolder="shared_results/anchorNone"
if os.path.exists(ofolder):
    goAhead=input(f'{ofolder} exists. Do you want to overwrite? (y/n)')
    if goAhead=='y':
        os.system(f'rm -r {ofolder}')
    else:
        print('exiting...')
        exit()
print()
os.system(f'mkdir -p {ofolder}')
###########################


###############################################################
#####  SYSTEMATICALLY VARY FORMAT 1: THE EVENT SELECTIONS IN THE ROOTDataReaderFilter
###############################################################
variations={}
variations['nominal']=[]
#variations['ueL1']=[['unusedEnergy',veryNegNum,0.25]]
#variations['ueL2']=[['unusedEnergy',veryNegNum,0.40]]
#variations['chiT']=[['chiSq',veryNegNum,10]]
#variations['chiL']=[['chiSq',veryNegNum,18]]
#variations['thetaBeamL']=[
#            ['!photonTheta1',veryNegNum,2.0,10.3,11.9],
#            ['!photonTheta2',veryNegNum,2.0,10.3,11.9],
#            ['!photonTheta3',veryNegNum,2.0,10.3,11.9],
#            ['!photonTheta4',veryNegNum,2.0,10.3,11.9]]
#variations['thetaBeamT']=[
#            ['!photonTheta1',veryNegNum,3.0,10.3,11.9],
#            ['!photonTheta2',veryNegNum,3.0,10.3,11.9],
#            ['!photonTheta3',veryNegNum,3.0,10.3,11.9],
#            ['!photonTheta4',veryNegNum,3.0,10.3,11.9]]
#variations['thetaTransL']=[ # this looser selection basically removes the cut on the transition region
#            ['!photonTheta1',veryNegNum,2.5],
#            ['!photonTheta2',veryNegNum,2.5],
#            ['!photonTheta3',veryNegNum,2.5],
#            ['!photonTheta4',veryNegNum,2.5]]
#variations['thetaTransT']=[
#            ['!photonTheta1',veryNegNum,2.5,10.0,12.5],
#            ['!photonTheta2',veryNegNum,2.5,10.0,12.5],
#            ['!photonTheta3',veryNegNum,2.5,10.0,12.5],
#            ['!photonTheta4',veryNegNum,2.5,10.0,12.5]]
#variations['eT1']=[
#            ['photonE1',0.12,veryPosNum],
#            ['photonE2',0.12,veryPosNum],
#            ['photonE3',0.12,veryPosNum],
#            ['photonE4',0.12,veryPosNum]]
#variations['eT2']=[
#            ['photonE1',0.13,veryPosNum],
#            ['photonE2',0.13,veryPosNum],
#            ['photonE3',0.13,veryPosNum],
#            ['photonE4',0.13,veryPosNum]]
#variations['pMomT1']=[['proton_momentum',0.325,veryPosNum]]
#variations['pMomT2']=[['proton_momentum',0.35,veryPosNum]] # DR study is 0.35 and 0.40. This is directly related to the t binning. So we have to be careful here
#variations['pZL']=[['proton_z',50,80]]
#variations['pZT']=[['proton_z',54,76]]
#variations['mmsqT1']=[['mmsq',-0.025,0.025]]
#variations['mmsqT2']=[['mmsq',-0.020,0.020]]

### SCANS OF OTHER EVENT SELECTIONS ###
#for s,v in zip(['15','16','17','18','19','20'],[1.5,1.6,1.7,1.8,1.9,2.0]):
#    variations[f'Mpi0pGT{s}']=[['Mpi0p',v,veryPosNum],['pVH',veryNegNum,veryPosNum]]
#for s,v in zip(['00','02','04','06','08'],[0.0,0.2,0.4,0.6,0.8]):
#    variations[f'cosThetaLower{s}']=[['cosTheta_eta_hel',veryNegNum,v]]#,['pVH',veryNegNum,veryPosNum]]
#for s,v in zip(['00','02','04','06','08'],[0.0,-0.2,-0.4,-0.6,-0.8]):
#    variations[f'cosThetaUpper{s}']=[['cosTheta_eta_hel',v,veryPosNum]]#,['pVH',veryNegNum,veryPosNum]]


print('variation: t-bin yields')
#i=0
#for k,vs in variations.items():
#    ofiles=setVariation(vs)
#    os.system(f"./reconfigureAndFit.py {k}")
#    os.system(f'mv 0*_{k} {ofolder}') 
#    yields=np.array([getYield(ofile) for ofile in ofiles])
#    if i==0:
#        initial_yields=yields.copy()
#    yields/=initial_yields
#    yields=[f'{x:0.2f}' for x in yields]
#    yields=' '.join(yields)
#    print(f"{k:<12}: {yields}")
#    i+=1

    

################################################################
######  SYSTEMATICALLY VARY FORMAT 2: SEARCH AND REPLACE VALUES
################################################################
## [polMag, uncertainty] of the phase1 dataset as calculated in 
## /d/grid17/ln16/dselector_v3/getPolMags/makePolValsV9/computeEtaPiPolMags.py 
#variations={}
#polmag000=[0.35062,0.00397] 
#polmag045=[0.34230,0.00412]
#polmag090=[0.34460,0.00404]
#polmag135=[0.35582,0.00420]
#variations['polMagLower']=[
#        # Tried wordbreak \b and that doesnt really work searching for exact words.
#        #    we can instead make it select a full word by including the spaces
#        [f' {polmag000[0]:0.5f}',f' {polmag000[0]-polmag000[1]:0.5f}'],
#        [f' {polmag045[0]:0.5f}',f' {polmag045[0]-polmag045[1]:0.5f}'],
#        [f' {polmag090[0]:0.5f}',f' {polmag090[0]-polmag090[1]:0.5f}'],
#        [f' {polmag135[0]:0.5f}',f' {polmag135[0]-polmag135[1]:0.5f}']]
#variations['polMagUpper']=[
#        [f' {polmag000[0]:0.5f}',f' {polmag000[0]+polmag000[1]:0.5f}'],
#        [f' {polmag045[0]:0.5f}',f' {polmag045[0]+polmag045[1]:0.5f}'],
#        [f' {polmag090[0]:0.5f}',f' {polmag090[0]+polmag090[1]:0.5f}'],
#        [f' {polmag135[0]:0.5f}',f' {polmag135[0]+polmag135[1]:0.5f}']]
#variations['polOffset']=[
#        # offsets are 1.77, 2.85, 4.50, 3.43 degrees. There are syst and stat errors on them but ignore for now
#        # found the offsets in Jon Zarlings thesis which is taken from Alex's rho analysis.
#        # The regex syntax is complicated, it uses capture groups to search for lines that have some basic format
#        #    the lines that we are interested in start with amplitude and contains a 'polAngle polMag' string 
#        #    all polMag starts with 0.3 since they are always atleast 30% polarization and never greater than 40% 
#        #    the first (.*) captures any string between 'amplitude' and ' 0.0' and can be access by \1 
#        [f'amplitude(.*) 0.0 0.3(.*)',r'amplitude\1 1.77 0.3\2'],
#        [f'amplitude(.*) 45.0 0.3(.*)',r'amplitude\1 47.85 0.3\2'],
#        [f'amplitude(.*) 90.0 0.3(.*)',r'amplitude\1 94.50 0.3\2'],
#        [f'amplitude(.*) 135.0 0.3(.*)',r'amplitude\1 138.43 0.3\2']]

#max_pcwsbins=18
#for nbins in [16,15,14,13]:
#    maxm=f'{1.04+nbins*0.04:0.2f}'
#    variations[f'pcws{nbins}bins']=[]
#    for x in range(nbins,max_pcwsbins):
#        for e in ['Re','Im']:
#            for p in ['Neg','Pos']:
#                variations[f'pcws{nbins}bins'].append([f'.*pcwsBin_{x}{e}{p}[ \t].*\n', ''])
#                # [^ ] matches anything but a space
#                variations[f'pcws{nbins}bins'].append([
#                    fr'amplitude ([^ ]*) Piecewise ([^ ]*) [^ ]* [^ ]* (.*) \[pcwsBin_{nbins-1}Im{p}\].*',
#                    fr'amplitude \1 Piecewise \2 {maxm} {nbins} \3 [pcwsBin_{nbins-1}Im{p}]'])
#                variations[f'pcws{nbins}bins'].append([fr'Mpi0eta 1.04 ([^ ]*) ',fr'Mpi0eta 1.04 {maxm} '])
#                variations[f'pcws{nbins}bins'].append([fr'Mpi0eta_thrown 1.04 ([^\n]*)',fr'Mpi0eta_thrown 1.04 {maxm}'])

#for k,variation in variations.items():
#    ofiles=setVariation([])
#    for ofile in ofiles:
#        with open(ofile,"r") as out:
#            lines=out.readlines()
#        lines=''.join(lines)
#        with open(ofile,"w") as out:
#            for search,replace in variation:
#                lines=re.sub(search,replace,lines)
#            out.write(lines)
#
#    os.system(f"./reconfigureAndFit.py {k}")
#    os.system(f'mv 0*_{k} {ofolder}') 

################################################################
######  SYSTEMATICALLY VARY FORMAT 3: VARY ANCHOR EXPLICITY USING AN ARGUMENT IN RECONFIGEANDFIT.PY
################################################################
variations={}
#for x in [0,4,6,8,11]:
for x in [-1]:
    variations[f'anchorBin{x}']=x
for k,variation in variations.items():
    ofiles=setVariation([])
    os.system(f"./reconfigureAndFit.py {k} {variation}")
    os.system(f'mv 0*_{k} {ofolder}') 

####### DRAW ALL THE VARIATIONS #######
important_files=['checkQuality.py','overlayBins.C','overlayBins.py','run_overlayBins.py','checkFits.py']
[os.system(f'cp {important_file} {ofolder}') for important_file in important_files]
os.system(f'ln -snfr rootFiles {ofolder}')
#os.chdir(f'{ofolder}')
#os.system('./run_overlayBins.py') # this cannot actually be run now, since we have two different versions, ones for GPU another for nogpu (with more RAM)
#os.chdir('..')





