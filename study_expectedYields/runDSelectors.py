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
baseOutputFolder="zDSelectedBkgndSamples_noUE_t010100_m104172"


topologies=[]
ofolders=[]
ifolders=[]

##############################################
# MAIN STUDY: RUN OVER ALL BACKGROUND SIMULATIONS TO DETERMINE EXPECTED LEAKAGES
##############################################
topologies+=[
#        "6#gammap[2#pi^{0},#eta]",
#        "6#gammap[2#pi^{0},#eta]",
#        "6#gammap[3#pi^{0}]",
#        "6#gammap[2#pi^{0},#eta]",
        "3#gammap[#pi^{0}]",
        "3#gammap[#pi^{0}]",
        "3#gammap[#pi^{0}]",
#        "4#gammap[2#pi^{0}]",
#        "5#gammap[2#pi^{0}]",
        "5#gammap[2#pi^{0}]",
        "5#gammap[2#pi^{0}]",
        ]
ofolders+=[
#        "a2pi",
#        "etap_to_etapipi",
#        "eta_to_3pi",
#        "f1_1285_to_etapipi",
        "omega_pi0g",
        "omega_pi0g_2018_8_v1",
        "omega_pi0g_2018_8_v2",
#        "pi0pi0",
#        "b1vps_2018_8",
        "b1vps_2018_1",
        "b1vps_2017_1",
        ]
ifolders+=[
#        '"/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/a2pi/tree_',
#        '"/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/etap_to_etapipi/tree_',
#        '"/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/eta_to_3pi/tree_',
#        '"/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/f1_to_etapipi/tree_',
        '"/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/omega_pi0g/tree_',
        '"/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/omega_pi0g_2018_8/tree_',
        '"/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/omega_pi0g_2018_8_v2/tree_',
#        '"/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/pi0pi0/tree_',
#        '"/d/grid17/ln16/rootFiles/omegapi_rad_massDepFit_2018_08_ver03.0_60M/root/merged/tree_',
        '"/d/grid17/ln16/rootFiles/omegapi_rad_massDepFit_2018_01_ver03.0/root/merged/tree_',
        '"/d/grid17/ln16/rootFiles/omegapi_rad_massDepFit_2017_01_ver03.0/root/merged/tree_',
        ]


##############################################
# SECONDARY STUDY: RUN OVER OMEGA AND FLAT ETAPI MC TO UNDERSTAND THE EFFECTS OF LMAC SELECTION
# 1. we do not need to run over the thrown - turn it off
##############################################
#topologies=[
#        "3#gammap[#pi^{0}]",
#        "4#gammap[#pi^{0},#eta]",
#        ]
#ofolders=[
#        "omega_pi0g",
#        "flat_2017_1",
#        ]
#ifolders=[
#        '"/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/omega_pi0g/tree_',
#        '"/d/grid17/ln16/rootFiles/pi0eta/120921/2017_1_130M/merged/tree_'
#        ]


##############################################
# THIRD STUDY:  RUN OVER DATA
##############################################

#topologies+=[
#        "4#gammap[#pi^{0},#eta]",
#        "4#gammap[#pi^{0},#eta]",
#        "4#gammap[#pi^{0},#eta]"
#        ]
#ofolders+=[
#        "gluex_2017_1",
#        "gluex_2018_1",
#        "gluex_2018_8"
#        ]
#ifolders+=[
#        '"/d/grid17/ln16/dselector_v3/phase1_data_looseChiUE/D2017_1_loose_tree.root"',
#        '"/d/grid17/ln16/dselector_v3/phase1_data_looseChiUE/D2018_1_loose_tree.root"',
#        '"/d/grid17/ln16/dselector_v3/phase1_data_looseChiUE/D2018_8_loose_tree.root"'
#        ]

##############################################
# FOURTH STUDY:  RUN OVER FLAT MC TO HELP STUDY QFACTORS
##############################################
#topologies+=[
#        "4#gammap[#pi^{0},#eta]",
#        "4#gammap[#pi^{0},#eta]",
#        "4#gammap[#pi^{0},#eta]"
#        ]
#ofolders+=[
#        "flat_2017_1",
#        "flat_2018_1",
#        "flat_2018_8"
#        ]
#ifolders+=[
#        '"/d/grid17/ln16/rootFiles/pi0eta/120921/2017_1_130M/merged/tree_',
#        '"/d/grid17/ln16/rootFiles/pi0eta/120921/2018_1_400M/merged/tree_',
#        '"/d/grid17/ln16/rootFiles/pi0eta/120921/2018_8_260M_130M/merged/tree_'
#        ]


assert len(ifolders)==len(ofolders) and len(ifolders)==len(topologies)

######## THINGS TO KEEP IN MIND ###########
# 1. We do not want to select on the any region in Meta vs Mpi so we "choice" 3
# 2. DSelector_etapi should not throw out 0 weights so we have that line commented out
# This will give us access to the full mass region so we can see how, for instance, 
#   the b1 leakage looks like in Meta 
###########################################

for ifolder,ofolder,topology in zip(ifolders,ofolders,topologies):
    replaceTopology(topology)
    if ".root" in ifolder: # this should only be used for data! 
        runSelector(ifolder,reconTreeName,"bkgndSample_recon",3,proof_Nthreads,recon_cfiles)
    else:
        runSelector(ifolder+'pi0eta*"',reconTreeName,"bkgndSample_recon",3,proof_Nthreads,recon_cfiles)
        runSelector(ifolder+'thrown*"',thrownTreeName,"bkgndSample_gen",1,proof_Nthreads,thrown_cfiles)
    move(baseOutputFolder+"/"+ofolder)

replaceTopology("4#gammap[#pi^{0},#eta]")
