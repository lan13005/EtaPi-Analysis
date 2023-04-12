#!/usr/bin/python3

import pandas as pd
import numpy as np

floc='./xsec_bootstrap_with_syst/results.csv'
df=pd.read_csv(floc)

def pretty(s):
    return f'{s:0.3f}'

new_df={}
new_df['t']=df['t']
new_df[r'\frac{d\sigma_{\text{pos}}}{dt}']=df['xsecpos'].apply(pretty)+'$\pm$'+ \
                                                df['xsecpos_staterr'].apply(pretty)+'$\pm$'+  \
                                                df['xsecpos_systerr'].apply(pretty) 
new_df[r'\frac{d\sigma_{\text{neg}}}{dt}']=df['xsecneg'].apply(pretty)+'$\pm$'+ \
                                                df['xsecneg_staterr'].apply(pretty)+'$\pm$'+ \
                                                df['xsecneg_systerr'].apply(pretty) 
new_df[r'\frac{d\sigma_{\text{total}}}{dt}']=df['xsec'].apply(pretty)+'$\pm$'+ \
                                                df['xsec_staterr'].apply(pretty)+'$\pm$'+ \
                                                df['xsec_systerr'].apply(pretty)
new_df=pd.DataFrame(new_df)
print(new_df.to_latex(escape=False, index=False))

tmin=df['t'].values[0]-df['tWidth'].values[0]/2
tmax=df['t'].values[-1]+df['tWidth'].values[-1]/2
print("-------------------------------")
print(f"Total cross section in {tmin:0.2f} < -t < {tmax:0.2f} GeV")
print("-------------------------------")
for base, label in zip(['xsecpos','xsecneg','xsec'],['Natural Production','Unnatural Production','Total']):
#for base, label in zip(['xsecpos'],['Natural Production']):
    nom=(df['tWidth']*df[base]).sum()
    err=np.sqrt(np.power(df['tWidth']*df[f'{base}_staterr'],2).sum())
    syst=np.sqrt(np.power(df['tWidth']*df[f'{base}_systerr'],2).sum())
    print(f'{label}: {nom:0.6f} +- {err:0.6f} +- {syst:0.6f}')
