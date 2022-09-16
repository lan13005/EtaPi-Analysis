#!/usr/bin/python3

import os
import re
import sys
import subprocess

'''
The purpose of this code is to reconfigure a .fit results file and seed the resulting cfg file with the best fit results 
  From here we can answer questions like how the modification of the data, i.e. by including a Delta cut, influences the amplitudes
'''

def replace(cfg,old,new):
    '''
    Replace the value of varName to varValue in the file called fileName. Depending on the value type we can include quotes or not
    '''
    sedArgs=["sed","-i","s@"+old+"@"+new+"@g",cfg]
    print(" ".join(sedArgs))
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()

def reconfigure(fname):
    '''
    This script creates a new configuration file based off the converged fit values in a .fit file
    The intention is to not change the actual values that are fitted for but you can modify the input data
        so that we can see the systemtatic changes things 
        1. Gauss constraint mean/mins
        2. Datasets
        3. The RootDataReaderFilter scheme
    '''

    with open(fname) as f:
        lines=f.readlines()
    
    cfgFirstLine=[i for i,v in enumerate(lines) if "## FIT CONFIGURATION ###" in v][0]
    fitResultLine=[i for i,v in enumerate(lines) if "+++ Parameter Values and Errors +++" in v][0]
    cfg=lines[cfgFirstLine+1:]
    fitresults=lines[fitResultLine:cfgFirstLine+1]
    fitresults=[re.sub(r" +","",line).rstrip().lstrip() for line in fitresults]
    # make sure the line is a parameterName and value that is separated by a tab
    fitresults=[line.split("\t") for line in fitresults if (len(line.split("\t"))==2)&(~line[0].isdigit())] 
    fitMap={k:v for k,v in fitresults}
    
    for i, line in enumerate(cfg):
        line=line.rstrip().lstrip()
        fields=line.split(" ")
        if fields[0]=="parameter": 
            # parameter lines have >=3 fields respectively
            # if the line has 4 fields then the 4th field is likely to be "fixed" so we dont need to change anything anyways
            if len(fields)==3:
                keyword, name, value = line.split() 
                cfg[i]=f"{keyword} {name} {fitMap[name]}\n"
            if len(fields)==6:
                keyword, name, value, constraint, constraintVar1, constraintVar2 = line.split()
                cfg[i]=f"{keyword} {name} {fitMap[name]} {constraint} {constraintVar1} {constraintVar2}\n"
        # initialize lines for amplitudes have 5 fields respectively. 
        # Sometimes it will have 6 where the 6th field is "fixed", so we wont have to change it anyways
        if len(fields)==5: 
            keyword, amp, coordinate, val1, val2 = line.split() 
            cfg[i]=f'{keyword} {amp} {coordinate}  {fitMap[amp+"_re"]} {fitMap[amp+"_im"]}\n'
    
    fout=fname.split(".")[0].split("/")[1]
    fout=fout+"_reconf.cfg"

    with open(fout,"w") as newcfg:
        newcfg.write("".join(cfg))

    return fout
    

if __name__ == "__main__":
    ts=["010020","0200325","0325050","050075","075100"]
    
    cfgs=[
            f"{t}/etapi_result.fit" 
             for t in ts
            ]
    
    nprocesses=9
    niters=1
    for t,cfgFile in zip(ts,cfgs):
        newCfgLoc=reconfigure(cfgFile)
        #replace(newCfgLoc,f"t{t}_m104180_vh",f"t{t}_m104180")

        ofolder=cfgFile.split("/")
        ofolder=ofolder[len(ofolder)-2:len(ofolder)-1][0]
        ofolder+="_reconfigured"
    
        print("Starting fits")
        cmd="mpirun -np "+str(nprocesses)+" fitMPI -c "+newCfgLoc+" -r "+str(niters)+" -m 1000000 -t 1.0 -x 0 -f 0.15" 
        pipeCmd=' > fit.log'
        print(cmd+pipeCmd)

        os.system(cmd+pipeCmd)

        print(f"moving output to {ofolder}")
        os.system(f"mkdir -p {ofolder}")
        os.system("mv -f etapi_result*.fit "+ofolder+" 2>/dev/null")
        os.system("mv -f "+newCfgLoc+" "+ofolder+" 2>/dev/null")
        os.system("mv -f *.log "+ofolder+" 2>/dev/null")
        os.system("mv -f *.ni "+ofolder+" 2>/dev/null")



