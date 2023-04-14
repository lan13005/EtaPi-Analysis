#!/usr/bin/python3

import os
import sys
import numpy as np
import subprocess
import time
import random

def copyAndModifyCfg(cfg_input,sources):
    seed=random.randint(100000,999999)
    cfg_output=cfg_input.split("/")[-1].split(".")[0]+f"-seed{seed}.cfg"
    os.system(f"cp {cfg_input} {cfg_output}")
    for source in sources:
        # search for "source" and append (using $ to reach end of line) with the seed value. Sources {data,bkgnd,accmc}
        #   the seed is randomly sampled with randint, by default the RNG is seeded with the system time
        cmd=f"sed -i '/^{source}/ s/$/ {seed}/' {cfg_output}" 
        #print(cmd)
        os.system(cmd)
        cmd=f"sed -i '/^{source}/ s/rootFiles/\/geode2\/home\/u060\/malbrec\/BigRed200\/EtaPi0_lng_systematics\/etapi_samples_01_19_23/g' {cfg_output}"
        #print(cmd)
        os.system(cmd)
    cmd=f"sed -i '/^genmc/ s/rootFiles/\/geode2\/home\/u060\/malbrec\/BigRed200\/EtaPi0_lng_systematics\/etapi_samples_10_15_22/g' {cfg_output}"
    #print(cmd)
    os.system(cmd)
    return cfg_output

def main():
    ##########################################
    # Basic setup
    ##########################################
    nprocesses=9
    fitFileName="etapi_result.fit"
    niters=100
    workingDir=os.getcwd()

    variations=[
            ["./etapi_result_t010020_nominal.cfg", ["data","bkgnd","accmc"], "t010020"],
            ["./etapi_result_t0200325_nominal.cfg", ["data","bkgnd","accmc"], "t0200325"],
            ["./etapi_result_t0325050_nominal.cfg", ["data","bkgnd","accmc"], "t0325050"],
            ["./etapi_result_t050075_nominal.cfg", ["data","bkgnd","accmc"], "t050075"],
            ["./etapi_result_t075100_nominal.cfg", ["data","bkgnd","accmc"], "t075100"],
            ]
    for variation in variations:
        print("\nNew bootstrap variations")
        baseCfgFile, sources, otag = variation
        for i in range(niters):
            cfg_output=copyAndModifyCfg(baseCfgFile,sources)
            
            print("Starting fits")
            cmd="mpirun -np "+str(nprocesses)+" fitMPI -c "+cfg_output+" -m 1000000 -t 1.0 -x 1" 
            pipeCmd=' > fit.log'
            print(cmd+pipeCmd)
#            os.system(cmd+pipeCmd)
            
            # Move results to the desired folder
            ofolder=f"bootstrap_{otag}/bootstrap_{i}"
            os.system("mkdir -p "+ofolder)
            os.system("mv -f etapi_result*.fit "+ofolder+" 2>/dev/null")
            os.system("mv -f "+cfg_output+" "+ofolder+" 2>/dev/null")
            os.system("mv -f *.log "+ofolder+" 2>/dev/null")
            os.system("mv -f *.ni "+ofolder+" 2>/dev/null")
            #os.system("ln -snfr rootFiles "+ofolder+"/rootFiles")
        
if __name__ == "__main__":
    start_time=time.time()
    main()
    stop_time=time.time()
    print(f"Total elapsed time: {stop_time-start_time}s or {(stop_time-start_time)/60}m")
    
