#!/usr/bin/python3

import os
import subprocess
import time
from multiprocessing import Pool

start = time.time()

##########################################
# Without compiling, running fitAsymmetryPlots over multiple process can interfere with each other
# the total run time actually increases due to the interference
##########################################
print("Deleting main and recompiling it")
subprocess.Popen("rm fitAsymmetryPlots", shell=True,  stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
compileMain=["g++","-o","fitAsymmetryPlots","fitAsymmetryPlots.C"]
try:
    rootFlags = subprocess.check_output(["root-config","--cflags","--glibs", "--libs"])
except:
    raise Exception("ROOT not loaded!")
rootFlags = rootFlags.decode(encoding="utf-8").rstrip().split(" ") #python3 requires decoding first, which python2 doesnt. But it doesn hurt
compileMain.extend(rootFlags)
print("\nStarting new compilation\n----------------------")
print(" ".join(compileMain))
out, err = subprocess.Popen(compileMain, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()
print(out)

if os.path.exists("logs"):
    os.system("rm -r logs")
os.system("mkdir -p logs")

def runCommand(cmd):
    os.system(cmd) # I could probably have use os.system as the command. 
    
#####################################
### Need to manually enter this, as parsing the file is complicated
# i is the ith element of eventSelects vector which contains a specific systematic you wish to run
js=range(17) # main systematics + default fits have 17 total variations
#js=[0,1]

cmds=["./fitAsymmetryPlots "+str(j)+" >> logs/output"+str(j)+".log 2>&1" for j in js]
print("\n********************************************\nExecuting all the following commands in 3 seconds!\n********************************************")
for cmd in cmds:
    print(cmd)
time.sleep(3)

with Pool(len(cmds)) as p:
    p.map(runCommand,cmds)

end = time.time()

print(f"\n***** FINISHED IN {end-start} seconds! *****")
