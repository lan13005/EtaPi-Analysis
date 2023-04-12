#!/usr/bin/python3

import os
import re
import numpy as np
from getYieldsFromConfig import getYield
import sys
import argparse
from reconfigureAndFit import reconfigureAndFit


parser = argparse.ArgumentParser()
parser.add_argument('alwaysOverwriteOutput', type=bool, nargs="?", default=False)
parser.add_argument('--tbins', type=int, nargs="+") # For example: ./runVariations.py --tbins 0 1 2
args = parser.parse_args()
alwaysOverwriteOutput=args.alwaysOverwriteOutput
tbinsChosen=args.tbins

ts=["010020","0200325","0325050","050075","075100"]
if tbinsChosen==None:
    tbinsChosen=[0,1,2,3,4]
ts=list(np.array(ts)[tbinsChosen])
#tbinsChosen=' '.join([str(i) for i in tbinsChosen])

veryNegNum=-999
veryPosNum=999

lowerMpi0eta=1.04
upperMpi0eta=1.72

class getNominal:
    def __init__(self,isgen):
        self.isgen=isgen
        self.nominal={
                'Mpi0eta':[lowerMpi0eta, upperMpi0eta],
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
                'Mpi0eta_thrown':[lowerMpi0eta, upperMpi0eta],
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
ofolder="shared_results/systematic_v29"
#ofolder="shared_results/trash"

if alwaysOverwriteOutput:
    pass
    #os.system(f'rm -r {ofolder}')
else:
    if os.path.exists(ofolder):
        option=input(f'{ofolder} exists. Do you want to overwrite yes/no or keep? (y/n/k): ')
        if option=='y':
            os.system(f'rm -r {ofolder}')
        elif option=='k':
            print('continuing...')
        else:
            print('exiting...')
print()
os.system(f'mkdir -p {ofolder}')
###########################

###############################################################
#####  SYSTEMATICALLY VARY FORMAT 1: THE EVENT SELECTIONS IN THE ROOTDataReaderFilter
###############################################################

variations={}
#variations['nominal']=[]
#variations['usL1']=[
#        ['unusedShowers',veryNegNum,1.5],
#        ['unusedEnergy',veryNegNum,veryPosNum] # unset the nominal unusedEnergy selection
#        ]
#variations['usL2']=[
#        ['unusedShowers',veryNegNum,2.5],
#        ['unusedEnergy',veryNegNum,veryPosNum] # unset the nominal unusedEnergy selection
#        ]
#variations['ueL1']=[['unusedEnergy',veryNegNum,0.25]]
#variations['ueL2']=[['unusedEnergy',veryNegNum,0.40]]
#variations['chiT']=[['chiSq',veryNegNum,10]]
#variations['chiL']=[['chiSq',veryNegNum,18]]
#variations['thetaBeamL']=[ # was set at 2.0
#            ['!photonTheta1',veryNegNum,2.1,10.3,11.9],
#            ['!photonTheta2',veryNegNum,2.1,10.3,11.9],
#            ['!photonTheta3',veryNegNum,2.1,10.3,11.9],
#            ['!photonTheta4',veryNegNum,2.1,10.3,11.9]]
#variations['thetaBeamT']=[ # was set at 3.0
#            ['!photonTheta1',veryNegNum,2.9,10.3,11.9],
#            ['!photonTheta2',veryNegNum,2.9,10.3,11.9],
#            ['!photonTheta3',veryNegNum,2.9,10.3,11.9],
#            ['!photonTheta4',veryNegNum,2.9,10.3,11.9]]
#variations['thetaTransL']=[ # looser selection originally just removed this selection
#            ['!photonTheta1',veryNegNum,2.5,10.4,11.7],
#            ['!photonTheta2',veryNegNum,2.5,10.4,11.7],
#            ['!photonTheta3',veryNegNum,2.5,10.4,11.7],
#            ['!photonTheta4',veryNegNum,2.5,10.4,11.7]]
#variations['thetaTransT']=[ # was set at 10.0, 12.5
#            ['!photonTheta1',veryNegNum,2.5,10.1,12.1],
#            ['!photonTheta2',veryNegNum,2.5,10.1,12.1],
#            ['!photonTheta3',veryNegNum,2.5,10.1,12.1],
#            ['!photonTheta4',veryNegNum,2.5,10.1,12.1]]
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
#variations['pZL']=[['proton_z',51,79]]
#variations['pZT']=[['proton_z',53,77]]
#variations['mmsqT1']=[['mmsq',-0.013,0.013]]
#variations['mmsqT2']=[['mmsq',-0.010,0.010]]

#for s,v in zip(['15','20','25','30','35','40'],[0.15,0.2,0.25,0.3,0.35,0.4]):
#for s,v in zip(['14','17'],[0.14,0.17]):
#    variations[f'ue{s}']=[['unusedEnergy',veryNegNum,v]]
#for s,v in zip(['11','12','14','16','18'],[11,12,14,16,18]):
#for s,v in zip(['18'],[18]):
#    variations[f'chi{s}']=[['chiSq',veryNegNum,v]]

### SCANS OF OTHER EVENT SELECTIONS ###
#for s,v in zip(['15','16','17','18','19','20'],[1.5,1.6,1.7,1.8,1.9,2.0]):
#for s,v in zip(['15','16','17','18'],[1.5,1.6,1.7,1.8]): # 1.9 and 2.0 change accidental subtracted yields by >10%
#    variations[f'Mpi0pGT{s}']=[['Mpi0p',v,veryPosNum],['pVH',veryNegNum,veryPosNum]]
##for s,v in zip(['00','02','04','06','08'],[0.0,0.2,0.4,0.6,0.8]):
#for s,v in zip(['00'],[0.0]):
#    variations[f'cosThetaLower{s}']=[['cosTheta_eta_hel',veryNegNum,v]]#,['pVH',veryNegNum,veryPosNum]]
##for s,v in zip(['00','02','04','06','08'],[0.0,-0.2,-0.4,-0.6,-0.8]):
#for s,v in zip(['00'],[0.0]):
#    variations[f'cosThetaUpper{s}']=[['cosTheta_eta_hel',v,veryPosNum]]#,['pVH',veryNegNum,veryPosNum]]
##for s,v in zip(['18','19','20','21','22','23','24','28','32'],[1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.8,3.2]):
##    variations[f'Mpi0pLT{s}']=[['Mpi0p',veryNegNum,v]]

print('variation: t-bin yields')
i=0
for k,vs in variations.items():
    ofiles=setVariation(vs)
    reconfigureAndFit(k,tbinsChosen=tbinsChosen)
    #os.system(f"./reconfigureAndFit.py {k} --tbins {tbinsChosen}")
    os.system(f'mv -f 0*_{k} {ofolder}') 
    

################################################################
######  SYSTEMATICALLY VARY FORMAT 2: SEARCH AND REPLACE VALUES
################################################################
## [polMag, uncertainty] of the phase1 dataset as calculated in 
## /d/grid17/ln16/dselector_v3/getPolMags/makePolValsV9/computeEtaPiPolMags.py 
variations={}
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

#max_pcwsbins=17 # this should be number of piecewise bins
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


#### Search and replace file locations
#variations['sbLaccN']=[[
#            'wUnusedShowers',
#            'wUnusedShowers_sbL_accN'
#            ]]
#variations['sbTaccN']=[[
#            'wUnusedShowers',
#            'wUnusedShowers_sbT_accN'
#            ]]
#variations['sbNaccL']=[[
#            'wUnusedShowers',
#            'wUnusedShowers_sbN_accL'
#            ]]
#variations['sbNaccT']=[[
#            'wUnusedShowers',
#            'wUnusedShowers_sbN_accT'
#            ]]

for k,variation in variations.items():
    ofiles=setVariation([])
    for ofile in ofiles:
        with open(ofile,"r") as out:
            lines=out.readlines()
        lines=''.join(lines)
        with open(ofile,"w") as out:
            for search,replace in variation:
                lines=re.sub(search,replace,lines)
            out.write(lines)

    reconfigureAndFit(k,tbinsChosen=tbinsChosen)
    #os.system(f"./reconfigureAndFit.py {k} --tbins {tbinsChosen}")
    os.system(f'mv -f 0*_{k} {ofolder}') 

################################################################
######  SYSTEMATICALLY VARY FORMAT 3: VARY ANCHOR EXPLICITY USING AN ARGUMENT IN RECONFIGEANDFIT.PY
################################################################
variations={}
#for x in [0,4,6,8,11]: # setting anchor to -1 means we do not change the anchor. basically is the nominal fit - just for testing
#    variations[f'anchorBin{x}']=x
for k,variation in variations.items():
    ofiles=setVariation([])
    reconfigureAndFit(k,anchor=variation,tbinsChosen=tbinsChosen)
    #os.system(f"./reconfigureAndFit.py {k} {variation} --tbins {tbinsChosen}")
    os.system(f'mv -f 0*_{k} {ofolder}') 


################################################################
######  SYSTEMATICALLY ADD WAVESETS: VARY ANCHOR EXPLICITY USING AN ARGUMENT IN RECONFIGEANDFIT.PY
################################################################
def loadValue(floc,search="bestMinimum"):
    with open(floc) as f:
        lines=[line.rstrip().lstrip() for line in f.readlines()]
        line=[float(line.split("\t")[1]) for line in lines if line.split("\t")[0]==search][0]
    return line

def incorporateWave(lines,waves):
    spect_map={"D": 2, "P": 1, "S": 0}
    def convertM(s):
        return (-1 if s[1]=='-' else 1)*int(s[0])
    for wave in waves:
        L=spect_map[wave[0]]
        m=convertM(wave[1:3])
        e="Positive" if wave[3]=='+' else "Negative"
        signs={
                "Positive": ['-1 -1', '+1 +1'],
                "Negative": ['-1 +1', '+1 -1']
                }
        reactions=['EtaPi0_000','EtaPi0_045','EtaPi0_090','EtaPi0_135']
        scalePols=["00","45","90","135"]
        for reaction,scalePol in zip(reactions,scalePols):
            for res1, res2 in zip(["","prime"],["","p"]):
                for part,sign in zip(["Im","Re"],signs[e]):
                    lines.append(f"amplitude {reaction}::{e}{part}::{res2}{wave} Zlm {L} {m} {sign} 0.0 0.3519\n")
                    lines.append(f"amplitude {reaction}::{e}{part}::{res2}{wave} BreitWigner [a2mass{res1}] [a2width{res1}] 2 2 3\n")
                    lines.append(f"initialize {reaction}::{e}{part}::{res2}{wave} cartesian 0.0 0.0\n")
                lines.append(f"scale {reaction}::{e}Im::{res2}{wave} [parScale{scalePol}]\n")
                lines.append(f"scale {reaction}::{e}Re::{res2}{wave} [parScale{scalePol}]\n")
                lines.append(f"constrain {reaction}::{e}Im::{res2}{wave} {reaction}::{e}Re::{res2}{wave}\n") # constrain Re and Im part
                for react in list(set(reactions)-set([reaction])):
                    lines.append(f"constrain {reaction}::{e}Im::{res2}{wave} {react}::{e}Re::{res2}{wave}\n") # constrain across polarizations Re
                    lines.append(f"constrain {reaction}::{e}Im::{res2}{wave} {react}::{e}Im::{res2}{wave}\n") # constrain across polarizations Im
                # amptools automatically initializes a waveset to 'cartesian 0 0' I think
    return lines

#allWaves=set(['D0++','D0+-','D1++','D1+-','D2++','D2+-','D1-+','D1--','D2-+','D2--']) # all potential waves
#baseWaves=set(['D2++','D1+-']) # waves in the minimal waveset
#usedWaves=set() # waves added so far
#while len(baseWaves)+len(usedWaves) < 10:
#    unused_waves=allWaves-baseWaves-usedWaves
#    variations={}
#    for wave in unused_waves:
#        ws=list(usedWaves)+[wave]
#        print(ws)
#        variations["waveset_"+'_'.join(ws)]=ws
#    waves=[] # I think dictionaries are unordered so when looping over the items later we might lose consistency if we fill waves beforehand
#    nlls=[]
#    print(f'1. variations: {variations}')
#    for k,variation in variations.items():
#        ofiles=setVariation([])
#        for ofile in ofiles:
#            with open(ofile,"r") as out:
#                lines=out.readlines()
#            with open(ofile,"w") as out:
#                lines=incorporateWave(lines,variation)
#                lines=''.join(lines)
#                out.write(lines)
#        reconfigureAndFit(k,tbinsChosen=tbinsChosen)
#        #os.system(f"./reconfigureAndFit.py {k} --tbins {tbinsChosen}")
#        waves.append(k.split("_")[-1]) # append only the newly added wave
#        nlls.append(loadValue(f'010020_{k}/etapi_result_0.fit')) # NEED TO UPDATE THIS LINE FOR OTHER T-BINS
#        os.system(f'mv -f 0*_{k} {ofolder}') 
#    bestWave=waves[np.argmin(nlls)]
#    usedWaves.add(bestWave)
#    print(f'2. NLLs: {list(zip(waves,nlls))}')
#    print(f'3. bestWave: {bestWave}')

################################################################
######  SYSTEMATICALLY ADD/SUBTRACT 1 ADDITIONAL WAVE
################################################################
variations={}
#waveMatch='D2\+\+' # will only rotate D2++ waves. Need to escape + signs
#waveMatch='p*D' # will match waves starting with 0+ "p" characters followed by a D
#waveMatch='pD' # will match waves starting with pD, or all the a2 primes
#variations['aD2mp_fixedPhase']=[['D2-+'],2] # "a" prefix means to add a wave whereas "m" prefix in the key means to minus
#variations['aD2mm_fixedPhase']=[['D2--'],2]
#variations['aD1mp_fixedPhase']=[['D1-+'],2]
#variations['aD2pm_fixedPhase']=[['D2+-'],2]
#variations['aD2mp_fixedMag']=[['D2-+'],3] 
#variations['aD2mm_fixedMag']=[['D2--'],3]
#variations['aD1mp_fixedMag']=[['D1-+'],3]
#variations['aD2pm_fixedMag']=[['D2+-'],3]
#variations['aD2mp_fixedMagPhase']=[['D2-+'],4] 
#variations['aD2mm_fixedMagPhase']=[['D2--'],4]
#variations['aD1mp_fixedMagPhase']=[['D1-+'],4]
#variations['aD2pm_fixedMagPhase']=[['D2+-'],4]
#variations['mD1mm']=['D1--'] # TMD waveset component
#variations['mD0pp']=['D0++'] # TMD waveset component
#variations['mD0pm']=['D0+-'] # TMD waveset component
#variations['mD1pp']=['D1++'] # TMD waveset component
#variations['mD1pm']=['D1+-'] # TMD waveset component
#variations['mD2pp']=['D2++'] # TMD waveset component

for k,variation in variations.items():
    ofiles=setVariation([])
    for ofile in ofiles:
        if k[0]=='a': # prefix with "a" if we wish to incorporate a wave
            with open(ofile,"r") as out:
                lines=out.readlines()
            with open(ofile,"w") as out:
                print(variation[0],variation[1])
                lines=incorporateWave(lines,variation[0])
                lines=''.join(lines)
                out.write(lines)
        elif k[0]=='m': # prefix with m if we wish to delete all rows of the file that has to do with a wave
            for v in variation[0]:
                os.system(f"sed -i '/{v}/d' {ofile}")
        else:
            print("Unexpected format for key when trying to +/- a wave")
            exit()
    reconfigureAndFit(k,tbinsChosen=tbinsChosen,method=variation[1],ignore_waves=variation[0],waveMatch=waveMatch) # ignore the wave we are adding
#    os.system(f"./reconfigureAndFit.py {k} --tbins {tbinsChosen} --method {variation[1]} --ignore_waves ")
    os.system(f'mv -f 0*_{k} {ofolder}') 



################################################################
######  SYSTEMATICALLY ZERO THE A2(1700) WAVES AND TRUNCATE THE MASS RANGE
################################################################
variations={}
#max_pcwsbins=17 # this should be number of piecewise bins
#for nbins in [17,16,15,14,13]:
#    maxm=f'{1.04+nbins*0.04:0.2f}'
#    variations[f'zero_pD_pcws{nbins}bins']=[['pD.*'],[]] # pD.* matches amplitude that starts with pD follow by any character any amount of times
#    for x in range(nbins,max_pcwsbins):
#        for e in ['Re','Im']:
#            for p in ['Neg','Pos']:
#                variations[f'zero_pD_pcws{nbins}bins'][1].append([f'.*pcwsBin_{x}{e}{p}[ \t].*\n', ''])
#                variations[f'zero_pD_pcws{nbins}bins'][1].append([
#                    fr'amplitude ([^ ]*) Piecewise ([^ ]*) [^ ]* [^ ]* (.*) \[pcwsBin_{nbins-1}Im{p}\].*',
#                    fr'amplitude \1 Piecewise \2 {maxm} {nbins} \3 [pcwsBin_{nbins-1}Im{p}]'])
#                variations[f'zero_pD_pcws{nbins}bins'][1].append([fr'Mpi0eta 1.04 ([^ ]*) ',fr'Mpi0eta 1.04 {maxm} '])
#                variations[f'zero_pD_pcws{nbins}bins'][1].append([fr'Mpi0eta_thrown 1.04 ([^\n]*)',fr'Mpi0eta_thrown 1.04 {maxm}'])
#for k in variations.keys():
#    variations[k][0]="|".join(variations[k][0])

def find_nth(haystack, needle, n):
    ''' Find nth occurence of the "needle" substr in the "haystack" string '''
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start

for k,variation in variations.items():
    ofiles=setVariation([])
    amp,vary = variation
    regexp_zeroer=re.compile(f'^initialize EtaPi0.*::.*::{amp}[-+]+.*')
    for ofile in ofiles:
        with open(ofile,"r") as out:
            lines=out.readlines()
        for i,line in enumerate(lines):
            if regexp_zeroer.search(line):
                line=line[:find_nth(line,' ',3)]
                line+=' 0.0 0.0 fixed\n'
                lines[i] = line
        lines=''.join(lines)
        for search,replace in vary:
            lines=re.sub(search,replace,lines)
        with open(ofile,"w") as out:
            out.write(lines)
    reconfigureAndFit(k,tbinsChosen=tbinsChosen)
    os.system(f'mv -f 0*_{k} {ofolder}') 




################################################################
######  SYSTEMATICALLY MODIFY THE RESONANCE PARAMETERS
################################################################
variations={}
# each map value is a list of lists. Each sub-list contains a search string and a value to be replaced at the specified index
#variations['a2massL']=[
#        ['parameter a2mass ', f'{0.0006*3:.5}', 5],  # search, value, idx. The value here is will be an integer multiple of the nominal
#        ]
#variations['a2widthL']=[
#        ['parameter a2width ',f'{0.0055*3:.5}', 5] 
#        ]
#variations['a2L']=[
#        ['parameter a2mass ', f'{0.0006*3:.5}', 5],  # search, value, idx
#        ['parameter a2width ', f'{0.0055*3:.5}', 5] 
#        ]
#variations['a2prmassL']=[
#        ['parameter a2massprime ', 'gaussian 1.698 0.04', 3],  # search, value, idx. The value here is will be an integer multiple of the nominal
#        ]
#variations['a2prwidthL']=[
#        ['parameter a2widthprime ', 'gaussian 0.265 0.060', 3] 
#        ]
#variations['a2prL']=[
#        ['parameter a2massprime ', 'gaussian 1.698 0.04', 3],  # search, value, idx
#        ['parameter a2widthprime ', 'gaussian 0.265 0.060', 3] 
#        ]


################################################################
######  SYSTEMATICALLY MODIFY THE RESONANCE PARAMETERS
################################################################
for k,variation in variations.items():
    ofiles=setVariation([])
    for ofile in ofiles:
        with open(ofile,"r") as out:
            lines=out.readlines()
        with open(ofile,"w") as out:
            for i in range(len(lines)):
                for v in variation:
                    search, value, idx = v
                    if lines[i].startswith(search):
                        line=lines[i].split(" ")
                        line[idx]=value
                        lines[i] = " ".join(line)+'\n'
            lines=''.join(lines)
            out.write(lines)
    reconfigureAndFit(k,tbinsChosen=tbinsChosen)
    #os.system(f"./reconfigureAndFit.py {k} --tbins {tbinsChosen}")
    os.system(f'mv -f 0*_{k} {ofolder}') 


################################################################
######  SYSTEMATICALLY MODIFY FCAL/BCAL EFFICIENCY
######  The purpose of this is to randomly subsample
################################################################
variations={}
variations['fcal99eff']=[-1,0.99] # first/second number is for bcal/fcal. -1 means we are not going to modify it
variations['fcal98eff']=[-1,0.98] # first/second number is for bcal/fcal. -1 means we are not going to modify it
variations['fcal97eff']=[-1,0.97] # first/second number is for bcal/fcal. -1 means we are not going to modify it
variations['fcal96eff']=[-1,0.96] # first/second number is for bcal/fcal. -1 means we are not going to modify it
variations['fcal95eff']=[-1,0.95] # first/second number is for bcal/fcal. -1 means we are not going to modify it
variations['bcal99eff']=[0.99,-1] # first/second number is for bcal/fcal. -1 means we are not going to modify it
variations['bcal98eff']=[0.98,-1] # first/second number is for bcal/fcal. -1 means we are not going to modify it
variations['bcal97eff']=[0.97,-1] # first/second number is for bcal/fcal. -1 means we are not going to modify it
variations['bcal96eff']=[0.96,-1] # first/second number is for bcal/fcal. -1 means we are not going to modify it
variations['bcal95eff']=[0.95,-1] # first/second number is for bcal/fcal. -1 means we are not going to modify it

def updateCfgWithRandomizer(floc,variation,seed='1992'):
    ### Replace folder
    folderSearch=r"rootFiles\(.*\)_selectGenTandM\/"
    folderReplace=r"rootFiles\1_selectGenTandM_nominal_wPhotonSyst\/"
    cmd=f"sed -i 's/{folderSearch}/{folderReplace}/' {floc}"
    os.system(cmd)
    ### Replace mass
    cmd=f"sed -i 's/104180/104172/g' {floc}"
    os.system(cmd)
    ### Replace file tags
    fileSearch="selected"
    fileReplace="selected_nominal_wPhotonSyst"
    cmd=f"sed -i 's/{fileSearch}/{fileReplace}/g' {floc}"
    os.system(cmd)
    ### Replace condition
    conditionSearch="Mpi0eta 1.04 1.72 pVH 0.5 999 unusedEnergy -999 0.01 chiSq -999 13.277 !photonTheta1 -999 2.5 !photonTheta1 10.3 11.9 !photonTheta2 -999 2.5 !photonTheta2 10.3 11.9 !photonTheta3 -999 2.5 !photonTheta3 10.3 11.9 !photonTheta4 -999 2.5 !photonTheta4 10.3 11.9 photonE1 0.1 999 photonE2 0.1 999 photonE3 0.1 999 photonE4 0.1 999 proton_momentum 0.3 999 proton_z 52 78 mmsq -0.05 0.05"
    conditionSearch.replace(' ','\ ')
    conditionReplace=f"{seed} "+conditionSearch
    if variation[0]!=-1:
        for i in range(4):
            conditionReplace+=f" ~photonSystem{i+1}_{variation[0]} -999 15"
    if variation[1]!=-1:
        for i in range(4):
            conditionReplace+=f" ~photonSystem{i+1}_{variation[1]} 15 999"
    conditionSearch.replace(' ','\ ')
    readerSearch='ROOTDataReaderFilter'
    readerReplace='ROOTDataReaderRandomFilter'
    for dtype in ['accmc']:
        cmd=f"sed -i '/^{dtype}/s/{conditionSearch}/{conditionReplace}/g' {floc}"
        os.system(cmd)
        ### Replace Reader except genmc
        cmd=f"sed -i '/^{dtype}/s/{readerSearch}/{readerReplace}/g' {floc}"
        os.system(cmd)

for k,variation in variations.items():
    ofiles=setVariation([])
    for ofile in ofiles:
        updateCfgWithRandomizer(ofile,variation)

    reconfigureAndFit(k,tbinsChosen=tbinsChosen)
    #os.system(f"./reconfigureAndFit.py {k} --tbins {tbinsChosen}")
    os.system(f'mv -f 0*_{k} {ofolder}') 



################################################################
######  SYSTEMATICALLY ALTER PIECEWISE BINNING
################################################################
variations={}
#variations['pcws40MeVbins']=1 # This number denotes how many times more piecewise bins to consider. Nominal is 40 MeV bins
#variations['pcws20MeVbins']=2 # This number denotes how many times more piecewise bins to consider. Nominal is 40 MeV bins
#variations['pcws13MeVbins']=3 # This number denotes how many times more piecewise bins to consider. Nominal is 40 MeV bins
#variations['pcws10MeVbins']=4 # This number denotes how many times more piecewise bins to consider. Nominal is 40 MeV bins

for k,variation in variations.items():
    ofiles=setVariation([])
    for ofile in ofiles:
        with open(ofile,"r") as out:
            lines=out.readlines()
        with open(ofile,"w") as out:
            ### Include additional lines to initialize new parameters according to searchString's search criteria ###
            searchStrings=[
                    ['parameter pcwsBin_',' ',1],  # Search for line starting with [0] and tack on additional tag at index [2] (after splitting on delimiter [1])
                    ['  pcwsBin_','\t',0], # use tab as the delimiter and dump the additional tag at index 0
                    ['parRange pcwsBin_',' ',1], 
            ]
            partSearcher=re.compile("Re|Im")
            deleteLines=[]
            for search, delim, loc in searchStrings:
                ##### DETERMINE WHICH LINES TO ANALYZE+MODIFY AND WHICH INDEX TO INSERT IT INTO ####
                for i,line in enumerate(lines):
                    if line.startswith(search):
                        splits=line.split("pcwsBin_")
                        partIdx=partSearcher.search(splits[1]).start()
                        ibin=splits[1][:partIdx]
                        additionalLines=[]
                        additionalLines.append(f'{splits[0]}pcwsBin_{variation*int(ibin)}{splits[1][partIdx:]}')
                        for j in range(1,variation):
                            additionalLines.append(f'{splits[0]}pcwsBin_{variation*int(ibin)+j}{splits[1][partIdx:]}')
                        lines[i]=''.join(additionalLines)
            ### Update Piecewise amplitude initialization to use the additional parameters ###
            for i,line in enumerate(lines):
                if line.startswith('amplitude') and 'Piecewise' in line:
                    line=line.split(' ')
                    ref=line[7]
                    pcwsBins=int(line[5])*variation
                    line[5]=str(pcwsBins) # scale the number of pcws bins
                    base=line[:9] # this is where the [pcwsBin_..] arguments starts 
                    for imass in range(pcwsBins):
                        base+=[f'[pcwsBin_{imass}Re{ref}] [pcwsBin_{imass}Im{ref}]']
                    lines[i]=' '.join(base)+'\n'
            lines=''.join(lines)
            out.write(lines)
    reconfigureAndFit(k,tbinsChosen=tbinsChosen)
    #os.system(f"./reconfigureAndFit.py {k} --tbins {tbinsChosen}")
    os.system(f'mv -f 0*_{k} {ofolder}') 

################################################################
####### DRAW ALL THE VARIATIONS #######
################################################################
important_files=['checkQuality.py','overlayBins.C','overlayBins.py','run_overlayBins.py','checkFits.py']
[os.system(f'cp {important_file} {ofolder}') for important_file in important_files]
os.system(f'ln -snfr rootFiles {ofolder}')
#os.chdir(f'{ofolder}')
#os.system('./run_overlayBins.py') # this cannot actually be run now, since we have two different versions, ones for GPU another for nogpu (with more RAM)
#os.chdir('..')





