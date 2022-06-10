#!/usr/bin/python

import subprocess
import os

outFolder="results_v3_tbins5"

if os.path.isdir(outFolder):
    raise ValueError("Output folder already exists! Choose a different name!")

ts=["010020","0200325","0325050","050075","075100","010020"] # supposed to cycle back to the first element to reset everything
for i,t in enumerate(ts):
    searchStr='$t=.*'
    replaceStr='$t="'+t+'";'
    sedArgs=["sed","-i",'s@'+searchStr+'@'+replaceStr+'@g','divideData.pl']
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()

    if i<len(ts): 
        os.system("./divideData.pl")
        os.system("python runFits.py")
        os.system("mkdir -p "+outFolder+"/"+t)
        os.system("mv finalAmps "+outFolder+"/"+t)
