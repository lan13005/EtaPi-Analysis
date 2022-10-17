#!/usr/bin/python3

import os
import glob
import sys

if __name__ == "__main__":
    def help():
        print("Hide folders containing specified fit iterations so that drawOverview will not draw them")
        print("requires at least one argument: basefolder")
        print("requires at least one optional arg defined below")
        print(" optional args: --iterations - semicolon separated string containing the underscore separated iterations to keep for each t bin.")
        print(" optional args: --revert - True or False, to revert to the original subdirectory naming or not")
        exit()

    argc=len(sys.argv)
    if argc<2 or "-h" in sys.argv:
        help()
    floc=sys.argv[1]
    if not os.path.exists(floc):
        print(f"Specific baseFolder '{floc}' does not exist!")
        exit()

    revert='--revert' in sys.argv
    if '--iterations' in sys.argv:
        i=sys.argv.index('--iterations')
        iterations=sys.argv[i+1].split(";")
        iterations=[i.split("_") if "_" in i else [i] for i in iterations]
        numIters=len(iterations)
    else:
        numIters=0
    if '--iterations' not in sys.argv and not revert:
        help()

    ts=["010020","0200325","0325050","050075","075100"]
    if numIters!=len(ts) and not revert:
        print("Expects 5 sets of iterations for the 5 tbins")
        exit()

    cmds=[]
    for i,t in enumerate(ts):
        fs=glob.glob(floc+"/"+t+"/**/"+t+"*",recursive=True)
        iterationsInDir=[f.split("_")[1] for f in fs]
        if revert:
            cmds+=[f'mv {f} {"_".join(f.split("_")[:-1])}' for f,b in zip(fs,iterationsInDir) if "hidden" in f]
        else:
            cmds+=[f'mv {f} {f}_hidden' for f,b in zip(fs,iterationsInDir) if b not in iterations[i]]

    for cmd in cmds:
        print(cmd)
        os.system(cmd)
    
