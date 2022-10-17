import uproot3 as uproot
import pandas as pd

def getYieldInFile(fileName,condition):
    ''' Function to load a root file '''
    tree=uproot.open(fileName)['kin']
    cols=condition.split(' ')[::3]
    cols=[x if x[0]!='!' else x[1:] for x in cols]+['weightASBS']
    cols=list(set(cols))
    #cols=['weightASBS','Mpi0eta','proton_z','proton_momentum','unusedEnergy','pVH','chiSq','photonTheta1','photonTheta2','photonTheta3','photonTheta4',
    #  'photonE1','photonE2','photonE3','photonE4','mmsq']
    df=tree.arrays(cols,outputtype=pd.DataFrame)

    vars1=condition.split(" ")[::3]
    mins1=[float(x) for x in condition.split(" ")[1::3]]
    maxs1=[float(x) for x in condition.split(" ")[2::3]]
    
    for var1,min1,max1 in zip(vars1,mins1,maxs1):
        if var1[0]=='!':
            var1=var1[1:]
            df=df[~((df[var1]>min1)&(df[var1]<max1))]
        else:
            df=df[((df[var1]>min1)&(df[var1]<max1))]

    return df['weightASBS'].sum()

def getYield(fname):
    with open(fname) as f:
        lines=f.readlines()
        lines=[l for l in lines if any([l.startswith(ftype) for ftype in ['data','bkgnd']])]
        fs=[l.split(' ')[3] for l in lines]
        conditions=[' '.join(l.split(' ')[4:]) for l in lines]
    
    totalYield=0
    for f,condition in zip(fs,conditions):
        a = 1 if 'data' else -1
        totalYield += a*getYieldInFile(f,condition)
        #print(f'current total {totalYield} after adding {f}')
    
    return totalYield
