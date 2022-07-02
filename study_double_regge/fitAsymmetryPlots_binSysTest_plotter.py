import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import mplhep as hep
mpl.rcParams.update(mpl.rcParamsDefault)
plt.style.use([hep.styles.ATLAS])
SMALL_SIZE = 24
MEDIUM_SIZE = 28
BIGGER_SIZE = 32

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

import pandas as pd
df=pd.read_csv("./results_bin_systematics.csv",delimiter=" ")

means=df.groupby("nbins")["par"].mean()
stds=df.groupby("nbins")["par"].std()
bins=df.nbins.unique()

print(means)
print(stds)
print(df[df.nbins==10])
print(df[df.nbins==150])
print(df[df.nbins==10]['par'].mean(), df[df.nbins==150]['par'].mean())
print(df[df.nbins==10]['par'].std(), df[df.nbins==150]['par'].std())

fig,ax=plt.subplots(1,1,figsize=(14,8))
ax.errorbar(bins,means,yerr=stds,fmt='.',c='black',linewidth=2)
ax.set_xlabel("Number of bins")
ax.set_ylabel("Correction")
ax.set_ylim(top=1.005)


axins = inset_axes(ax,width="50%",height="40%",loc=4,borderpad=3)
axins.errorbar(bins,means,yerr=stds,fmt='.',c='black',linewidth=2)
axins.set_ylim(bottom=0.985,top=1.005)
axins.set_xlabel("Number of bins")
axins.set_ylabel("Correction")
axins.set_title("Zoomed")

fig.savefig("results_bin_systematics.png")
