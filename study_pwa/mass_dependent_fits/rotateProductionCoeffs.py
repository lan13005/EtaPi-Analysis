import re
import numpy as np

def simpleAmplitude(amplitude):
    '''
    Converts from EtaPi0_000::NegativeRe::D0+- to just D0+-
    '''
    k=amplitude.split('::')[-1]
    return k
    
def rotateProductionCoeffs(cfg,ignore_waves,waveMatch,method=0):
    '''
    "method" use cases:
    0. Do absolutely nothing
    1. Almost do nothing. We include parameters but fix them to zero. Simple test case. i.e. `parameter phase0 0 fixed`
        a. Free production coeff, phase fixed to 0
    2. Fix phase and free magnitude
        a. Rotated production coeff to real. phase set to the rotated phase and is fixed
    3. Fix magnitude and free phase
        a. Free production coeff, phase free

    ignore_waves: Ignore these waves when choosing which waves to rotate. Designed for waveset addition systematic - ignoring the added wave
        Default to some weird string that should never be matched
    '''
    if method==1:
        print("  rotateProductionCoeffs.py will fix new phase parameters to 0 as a test!")
    elif method==2:
        print("  rotateProductionCoeffs.py will fix phase and free magnitude!")
    elif method==3:
        print("  rotateProductionCoeffs.py will fix magnitude and free phase!")
    elif method==4:
        print("  rotateProductionCoeffs.py will fix the magnitude and phase!")
    else: # if 0, or anything else tbh
        print("  rotateProductionCoeffs.py will do nothing as requested...")
        return 0
        

    with open(cfg,'r') as f:
        lines=f.readlines()

    # TMD waveset is S+D waves. S is piecewise so production coeff is just (1,0). D-waves are the only ones that need rotating
    #waveMatch='D2\+\+' # will only rotate D2++ waves. Need to escape + signs
    #waveMatch='p*D' # will match waves starting with 0+ "p" characters followed by a D
    #waveMatch='pD' # will match waves starting with pD, or all the a2 primes
    REACTION='EtaPi0_\d{3}' # the reaction name. \d{3} will match any 3 digits

    # Search for "initialize" lines so we can rotate away the production coefficient
    #   Search for "amplitude" lines so we can incorporate the phase that was rotated away
    initialize={} # key:value = amplitude:line. amplitude of form :  EtaPi0_000::NegativeRe::D0+-
    zlm={} # key:value = amplitude:line
    ignore_waves=[w.replace('+','\+').replace('-','\-') for w in ignore_waves] # +/- needs to be escaped in the regex
    regexp_ignore=re.compile('|'.join(ignore_waves)) # join the words with an or statement
    regexp=re.compile(f'{REACTION}::(Negative|Positive)(Re|Im)::{waveMatch}.*')
    regexp_zlm=re.compile('^amplitude.*Zlm.*')
    regexp_init=re.compile('^initialize')
    for line in lines:
        if regexp.search(line):
            if regexp_zlm.search(line) and not regexp_ignore.search(line):
                amplitude=line.split(' ')[1]
                zlm[amplitude]=line
            if regexp_init.search(line) and not regexp_ignore.search(line):
                amplitude=line.split(' ')[1]
                initialize[amplitude]=line

    # Keep track of what to search and replace for initialize lines
    init_replace={} # key:value = search:replace
    amplitudeMap={} # key:value has form D2++:'phase0'. I.E. will include connect D2++ with parmaeter phase0
    additionalParameterLines='' # Keep track of the additional parameter lines to add
    iamp=0
    for amp,line in initialize.items():
        rePart,imPart=line.split(' ')[-2:]
        rePart,imPart=float(rePart),float(imPart)
        productionCoeff = complex(rePart,imPart)
        r=abs(productionCoeff)
        theta=np.angle(productionCoeff)
        #print(f'{simpleAmplitude(amp)} ReIm {rePart} {imPart} Polar {r} {theta}')
        if not (rePart==0 and imPart==0):
            # only modify the line if the wave is not newly initialized, (set to 0,0). Used for waveset systematic
            if method==1:
                init_replace[line]=line
            elif method==2:
                init_replace[line]=' '.join(line.rstrip().split(' ')[:-2])+f' {r} 0.0 real\n' # need to strip the newline and add it back
            elif method==3 or method==4:
                init_replace[line]=line.rstrip()+' fixed\n'
        simpAmp=simpleAmplitude(amp)
        if simpAmp not in amplitudeMap.keys():
            amplitudeMap[simpAmp]=f'phase{iamp}' # store the associated parameter name. zlm_replace will use it to match
            # Production Coeff = A+iB = Rexp(Theta), move exp(Theta) to Zlm amplitude, leave R as the production coeff fixed to be real
            if method==1 or method==4:
                additionalParameterLines+=f'parameter phase{iamp} 0.0'
                additionalParameterLines+=' fixed'
            elif method==2:
                additionalParameterLines+=f'parameter phase{iamp} {theta}'
                additionalParameterLines+=' fixed'
            elif method==3:
                additionalParameterLines+=f'parameter phase{iamp} 0.0'
            additionalParameterLines+='\n'
            iamp+=1

    #for k,v in init_replace.items():
    #    print(k,v)
        
    # Find initial parameter line and dump additionalParameterLines before it
    firstParamLine=min([i for i,line in enumerate(lines) if line.startswith('parameter')])
    lines.insert(firstParamLine,additionalParameterLines)

    # Keep track of what to search and replace for zlm amplitude lines
    zlm_replace={} # key:value = search:replace
    for amp,line in zlm.items():
        simpAmp=simpleAmplitude(amp)
        zlm_replace[line]=line.rstrip()+f' [{amplitudeMap[simpAmp]}]\n'

    lines=''.join(lines)
    for search,replace in init_replace.items():
        lines=lines.replace(search,replace.lstrip())
    for search,replace in zlm_replace.items():
        lines=lines.replace(search,replace.lstrip())

    with open(cfg,'w') as f:
        f.write(lines)

    return 0
    

