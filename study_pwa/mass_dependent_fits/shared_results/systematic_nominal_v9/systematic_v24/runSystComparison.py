#!/usr/bin/python3

import os
import numpy as np

def cmd(call):
#    print(call)
    os.system(call)

minm=1.04
mwidth=0.04

tags=['nominal']
variations=[]
splits=[0]

# Which groups to calculate a maximum deviation in
#    ['event','counting'] will both first select barlow significant variations before calculating maximum deviation
#    all others will consider all variations and take the maximum deviation
classes=[] 

#### DEPRECATED: CAN TRY TO USED BOOTSTRAPPED ERRORS AS THE "NOMINAL STAT". THIS IS NOT THE WAY TO GO SINCE 
####    WE HAVE TO BE COMPARING MINUIT UNCERTAINTIES OF THE NOMINAL TO THE VARIATIONS WHEN DETERMINING BARLOW SIGNIFICANCE 
#bootstrap_folder='/d/grid17/ln16/dselector_v3/study_pwa/mass_dependent_fits/shared_results/systematic_nominal_v9/bootstrap/bootstrap_results'
bootstrap_folder=''

###########################
### Fiducial Selections
###########################
x=['thetaBeamL', 'thetaBeamT', 'thetaTransL', 'thetaTransT', 'eT1', 'eT2', 'pZL', 'pZT', 'mmsqT1', 'mmsqT2']
tags.extend(x)
variations.extend([
        r'$\theta_{\gamma,beam}>2.1^\circ$',r'$\theta_{\gamma,beam}>2.9^\circ$'
        ,r'$\theta_{\gamma,trans}$ cut [10.4,11.7]',r'$\theta_{\gamma,trans}$ cut $[10.1,12.1]^\circ$'
        ,r'$E_{\gamma}>0.12 GeV$',r'$E_{\gamma}>0.13 GeV$'
        ,r'$51<z_{proton}<79 cm$',r'$53<z_{proton}<77 cm$'
        ,r'$|MMsq|<0.013 GeV^2$',r'$|MMsq|<0.01 GeV^2$',
        ])
classes.extend(['event','event','event','event','event','event','event','event','event','event'])
print(len(tags),len(variations))

##########################
## Unused Shower
##########################
#x=['usL1','usL2']
#tags.extend(x)
#variations.extend(['<=1 unused showers', '<=2 unused showers'])
#classes.extend(['unusedEnergy','unusedEnergy']) # even though it is #unusedShowers it is still related to unusedEnergy
#splits.append(len(x))
#splits.append(len(x)+splits[-1])

###########################
### Unused Energy
###########################
##x=['15','20','25','30','35','40']
##xx=[0.15,0.2,0.25,0.3,0.35,0.4]
x=['14','15','17']
xx=[0.14,0.15,0.17]
for v,vv in zip(x,xx):
    tags.append(f'ue{v}')
    variations.append(fr'$E_{{unused}}<{vv} GeV$')
    classes.append('event')
print(len(tags),len(variations))

############################
#### ChiSq
############################
x=['11','12','14','16']#,'18']
xx=[11,12,14,16]#,18]
for v,vv in zip(x,xx):
    tags.append(f'chi{v}')
    variations.append(fr'$\chi^2<{vv}$')
    classes.append('event')
print(len(tags),len(variations))

###########################
### Baryon sensitive
###########################
x=[15,16,17,18]#,19,20]
xx=[1.5,1.6,1.7,1.8]#,1.9,2.0]
for v,vv in zip(x,xx):
    tags.append(f'Mpi0pGT{v}')
    variations.append(fr'$M(\pi p)>{vv}$')
    classes.append('event')
print(len(tags),len(variations))
splits.append(len(tags)-1) #### EVENT SELECTIONS COMPLETE

############################
#### Polarization Variations
############################
x=['polMagLower','polMagUpper','polOffset']
tags.extend(x)
variations.extend([r'$|P_{\gamma}|$ lower limit', r'$|P_{\gamma}|$ upper limit', r'$\phi_{\gamma}$ offset'])
classes.extend(['polarization','polarization','polarization'])
splits.append(len(tags)-1)
print(len(tags),len(variations))

###########################
### M(etapi) dependent
###########################
x=[16,15,14]#,13]
for i in x:
    tags.append(f'pcws{i}bins')
    variations.append(f'${minm}<M(\eta\pi)<{minm+mwidth*i:0.2f}$')
    classes.append('fitrange')
splits.append(len(tags)-1)

############################
### Mass sideband regions scan
############################
x=['sbLaccN','sbTaccN']
tags.extend(x)
variations.extend(['loose sideband','tight sideband'])
classes.extend(['counting','counting'])
print(len(tags),len(variations))

############################
### RF sideband regions scan
############################
x=['sbNaccL','sbNaccT']
tags.extend(x)
variations.extend(['+1 RF bunch','-1 RF bunch'])
classes.extend(['counting','counting'])
print(len(tags),len(variations))
splits.append(len(tags)-1)

############################
### anchor 
############################
#x=[0,6,8,11]
#for i in x:
#    tags.append(f'anchorBin{i}')
#    variations.append(f'anchor S-wave at {minm+mwidth*(i+0.5):0.2f}')
#    classes.append('anchor')
#splits.append(len(x)+splits[-1])
#print(len(tags),len(variations))

############################
### Pcws bin width
############################
x=['20']#, '13', '10']
for i in x:
    tags.append(f'pcws{i}MeVbins')
    variations.append(f'pcwsBinWidth {i}MeV')
    classes.append('pcws_width')
splits.append(len(tags)-1)
print(len(tags),len(variations))


###########################
##### Resonance parameters
###########################
x=['a2massL','a2widthL','a2L','a2prmassL','a2prwidthL','a2prL']
tags.extend(x)
variations.extend([
    '$a_2(1320)$ loose M',
    '$a_2(1320)$ loose $\Gamma$',
    '$a_2(1320)$ loose params',
    '$a_2(1700)$ loose M',
    '$a_2(1700)$ loose $\Gamma$',
    '$a_2(1700)$ loose params'
        ])
classes.extend([
    'bw_param',
    'bw_param',
    'bw_param',
    'bw_param',
    'bw_param',
    'bw_param',
    ])
splits.append(len(tags)-1)
print(len(tags),len(variations))

###########################
##### Additional Waveset systematics fixing particular contributions
###########################
#xs=['aD2mp','aD2mm','aD1mp','aD2pm']
#labels=['add $D_{-2}^{+}$', 'add $D_{-2}^{-}$', 'add $D_{-1}^{+}$', 'add $D_{2}^{-}$']
#for i,x,label in zip(range(len(xs)),xs,labels):
#    extension=[x,x+'_fixedPhase',x+'_fixedMag',x+"_fixedPhase_pD",x+"_fixedMag_pD"]
#    tags.extend(extension)
#    variations.extend([label,label+' $a_2$ $\phi$-fixed',label+' $a_2$ $R$-fixed',label+' $a_2\prime$ $\phi$-fixed',label+' $a_2\prime$  $R$-fixed'])
#    splits.append(len(extension) if i==0 else len(extension)+splits[-1])
#    classes.extend(['waveset']*len(extension))
#print(len(tags),len(variations))

###########################
##### Zeroing particular contributions and truncating M4g range
###########################
#minm=1.04
#mwidth=0.04
#x=[17,16,15,14,13]
#for i in x:
#    tags.append(f'zero_pD_pcws{i}bins')
#    variations.append(f'0 a2p ${minm}<M(\eta\pi)<{minm+mwidth*i:0.2f}$')
#    classes.append('zeroed_pD_Mpi0eta')
#splits.append(len(x)+splits[-1])
#print(len(tags),len(variations))


###########################
### Systematically change efficiency in MC
###########################
#x=[]
#for system in ['bcal','fcal']:
#    percentages=['95','96','97','98','99']
#    percentages.reverse()
#    for eff in percentages:
#        if int(eff)<97 and system=='fcal':
#            continue
#        x.append(f'{system}{eff}eff')
#        variations.append(f'{system.upper()} Eff. Reduce by {100-int(eff)}%')
#        classes.append('photonRecon')
#tags.extend(x)
#splits.append(len(x))

assert len(tags)==len(variations)+1

cmd(f"rm -rf SOURCE")
cmd(f"mkdir -p SOURCE")

ts=['010020','0200325','0325050','050075','075100']
allVariations=['nominal']+variations
with open(f'SOURCE/likelihoods.log','w') as likelihoodFile, open(f'SOURCE/variationsList.log','w') as variationFile:
    header= '# File is generated by runSystOverview.py     #\n'
    header+='# Contains a list of the variations performed #\n'  
    header+='# in order, matching the naming of fit files  #\n'
    header+='# VARIATION LIKELIHOOD                        #\n'
    variationFile.write(header)
    variationFile.write('\n'.join(allVariations))
    minima=[]
    for t in ts:
        cmd(f"mkdir -p SOURCE/{t}")
        for j,tag in enumerate(tags):
            if not os.path.exists(f'{t}_{tag}/etapi_result_0.fit'):
                print(f'{t}_{tag} did not contain a etapi_result_0.fit file! Might get an error from symlink but its ok!')
            fitdest=f'SOURCE/{t}/etapi_result_{j}.fit'
            cmd(f"ln -snfr {t}_{tag}/{t}_{tag}_0 SOURCE/{t}/{t}_{j}")
            cmd(f"ln -snfr {t}_{tag}/etapi_result_0.fit {fitdest}")
            cmd(f"ln -snfr {t}_{tag}/fit.log SOURCE/{t}/{t}_{j}/fit.log")
            cmd(f"ln -snfr {t}_{tag}/normInt* SOURCE/{t}")
            if os.path.exists(fitdest): # there could be broken symlinks
                with open(fitdest) as fit:
                    minimum = [float(line.split("\t")[-1]) for line in fit.readlines() if "bestMinimum" in line][0]
                    minima.append(minimum)
            else:
                print(f'{fitdest} link is missing or broken! Ignoring that file...')
                minima.append(np.inf) # -1 as a failure signal
    minima=np.array(minima).reshape(5,-1)
    for v in range(len(allVariations)):
        vv=len(allVariations)-1-v # additional -1 since python is 0 indexed
        line=f'{allVariations[vv]:<30}'
        for i in range(len(ts)):
            line+=f' {minima[i][vv]-minima[i][0]:<10.1f}'
        likelihoodFile.write(f'{line}\n')
        #likelihoodFile.write(f'{variation:<40} {minimum-minima[0]}\n')
    #likelihoodFile.write("\n")

keepOnlyConvergedFits='False'
#level=2
#cwd='/'.join(os.getcwd().split('/')[:-level])
#os.chdir('/'.join(['..']*level))
cwd=os.getcwd()
os.chdir('../..')
otag='syst'
variations=';'.join(variations)
classes=';'.join(classes)
splits=';'.join([str(s) for s in splits])
cmd(f"./compareSyst.py {cwd}/SOURCE '{bootstrap_folder}' {otag} '{variations}' '{classes}' '{splits}' {keepOnlyConvergedFits}")

os.chdir(cwd)






