#!/usr/bin/python

import os

######################################
# The purpose of this script is to run the DSelector several times over different root files
#    that can also have different weight schemes. We can also run different DSelectors like
#    the one used for recon and thrown trees
######################################

def runSelector(inFileLoc, treeName, outFileName, choice, proof_Nthreads, cfiles):
    '''
    Run the DSelector over some files [inFileLoc] and output 3 root files: [outFileName]_{hist,tree,flat}.root
    The Weight branch will be filled following one of 3 schemes of your choice
        1. data
        2. bkgnd
        3. acc
    More information on the what these choices refer to can be found in the DSelector
    '''
    choiceToType={1:"data",2:"bkgnd",3:"acc"}
    csel, crun=cfiles
    os.system('sed -i "s/choice=[0-9];/choice='+str(choice)+';/g" '+csel)
    print("root -l -b -q '"+crun+"("+inFileLoc+","+'"'+treeName+'","'+outFileName+"_"+choiceToType[choice]+'",'+str(proof_Nthreads)+")'")
    os.system("root -l -b -q '"+crun+"("+inFileLoc+","+'"'+treeName+'","'+outFileName+"_"+choiceToType[choice]+'",'+str(proof_Nthreads)+")'")
    os.system("mv output_flat.root "+outFileName+"_"+choiceToType[choice]+"_flat.root")

def replaceTopology(topology):
    for afile in ["DSelector_etapi.C","DSelector_thrown.C"]:
        cmd="sed -i 's/topologyString=\".*\"/topologyString=\"%s\"/' %s" % (topology,afile)
        os.system(cmd)

def move(folder):
    os.system("mkdir -p %s" % folder)
    os.system("yes | mv *bkgndSample*.root %s" % folder)


proof_Nthreads=36
recon_cfiles=["DSelector_etapi.C", "runDSelector.C"]
thrown_cfiles=["DSelector_thrown.C", "runDSelector_thrown.C"]

reconTreeName="pi0eta__B4_M17_M7_Tree"
thrownTreeName="Thrown_Tree"

topologies=[
        "6#gammap[2#pi^{0},#eta]",
        "6#gammap[2#pi^{0},#eta]",
        "6#gammap[3#pi^{0}]",
        "6#gammap[2#pi^{0},#eta]",
        "3#gammap[#pi^{0}]",
        "4#gammap[2#pi^{0}]",
        "5#gammap[2#pi^{0}]",
        "5#gammap[2#pi^{0}]",
        "5#gammap[2#pi^{0}]"
        ]
ofolders=[
        "zDSelectedBkgndSamples/a2pi",
        "zDSelectedBkgndSamples/etap_to_etapipi",
        "zDSelectedBkgndSamples/eta_to_3pi",
        "zDSelectedBkgndSamples/f1_1285_to_etapipi",
        "zDSelectedBkgndSamples/omega_pi0g",
        "zDSelectedBkgndSamples/pi0pi0",
        "zDSelectedBkgndSamples/b1vps_2018_8",
        "zDSelectedBkgndSamples/b1vps_2018_1",
        "zDSelectedBkgndSamples/b1vps_2017_1"
        ]
ifolders=[
        '"/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/a2pi/tree_',
        '"/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/etap_to_etapipi/tree_',
        '"/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/eta_to_3pi/tree_',
        '"/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/f1_to_etapipi/tree_',
        '"/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/omega_pi0g/tree_',
        '"/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/pi0pi0/tree_',
        '"/d/grid17/ln16/rootFiles/omegapi_rad_massDepFit_2018_08_ver03.0_60M/root/merged/tree_',
        '"/d/grid17/ln16/rootFiles/omegapi_rad_massDepFit_2018_01_ver03.0/root/merged/tree_',
        '"/d/grid17/ln16/rootFiles/omegapi_rad_massDepFit_2017_01_ver03.0/root/merged/tree_'
        ]


assert len(ifolders)==len(ofolders) and len(ifolders)==len(topologies)


for ifolder,ofolder,topology in zip(ifolders,ofolders,topologies):
    replaceTopology(topology)
    runSelector(ifolder+'pi0eta*"',reconTreeName,"bkgndSample_recon",1,proof_Nthreads,recon_cfiles)
    runSelector(ifolder+'thrown*"',thrownTreeName,"bkgndSample_gen",1,proof_Nthreads,thrown_cfiles)
    move(ofolder)

replaceTopology("4#gammap[#pi^{0},#eta]")




