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

proof_Nthreads=48
recon_cfiles=["DSelector_etapi.C", "runDSelector.C"]
thrown_cfiles=["DSelector_thrown.C", "runDSelector_thrown.C"]

reconTreeName="pi0eta__B4_M17_M7_Tree"
thrownTreeName="Thrown_Tree"

#### PHASE 1 MONTE CARLO
#runSelector('"/d/grid17/ln16/rootFiles/pi0eta/120921/2017_1_130M/merged/tree_pi0eta*"',reconTreeName,"F2017_1_selected",3,proof_Nthreads,recon_cfiles)
#runSelector('"/d/grid17/ln16/rootFiles/pi0eta/120921/2018_1_400M/merged/tree_pi0eta*"',reconTreeName,"F2018_1_selected",3,proof_Nthreads,recon_cfiles)
#runSelector('"/d/grid17/ln16/rootFiles/pi0eta/120921/2018_8_260M_130M/merged/tree_pi0eta*"',reconTreeName,"F2018_8_selected",3,proof_Nthreads,recon_cfiles)

## FILTERED PHASE 1 DATA
#for i in [3]:
#    runSelector('"./phase1_data_looseChiUE/D2017_1_loose_tree.root"',reconTreeName,"D2017_1_selected",i,proof_Nthreads,recon_cfiles)
#    runSelector('"./phase1_data_looseChiUE/D2018_1_loose_tree.root"',reconTreeName,"D2018_1_selected",i,proof_Nthreads,recon_cfiles)
#    runSelector('"./phase1_data_looseChiUE/D2018_8_loose_tree.root"',reconTreeName,"D2018_8_selected",i,proof_Nthreads,recon_cfiles)

## UNFILTERED PHASE 1 DATA
for i in [3]:
    runSelector('"/d/grid17/ln16/rootFiles/pi0eta/RunPeriod-2017-01/analysis-ver52/tree_pi0eta__B4_M17_M7/merged/*"',reconTreeName,"D2017_1_selected",i,proof_Nthreads,recon_cfiles)
    runSelector('"/d/home/sdobbs/GlueX/gluex_data/RunPeriod-2018-01/analysis-ver02/tree_pi0eta__B4_M17_M7/merged/*"',reconTreeName,"D2018_1_selected",i,proof_Nthreads,recon_cfiles)
    runSelector('"/d/home/sdobbs/GlueX/gluex_data/RunPeriod-2018-08/analysis-ver02/tree_pi0eta__B4_M17_M7/merged/*"',reconTreeName,"D2018_8_selected",i,proof_Nthreads,recon_cfiles)

### Thrown trees generally need less threads to run over since there isn't much calculation being done
#runSelector('"/d/grid17/ln16/rootFiles/pi0eta/120921/2017_1_130M/merged/tree_thrown*"',thrownTreeName,"F2017_1_gen",1,proof_Nthreads,thrown_cfiles)
#runSelector('"/d/grid17/ln16/rootFiles/pi0eta/120921/2018_1_400M/merged/tree_thrown*"',thrownTreeName,"F2018_1_gen",1,proof_Nthreads,thrown_cfiles)
#runSelector('"/d/grid17/ln16/rootFiles/pi0eta/120921/2018_8_260M_130M/merged/tree_thrown*"',thrownTreeName,"F2018_8_gen",1,proof_Nthreads,thrown_cfiles)

#### KMATRIX MC
#reconTreeName="pi0eta__B4_M7_M17_Tree"
#runSelector('"/d/grid17/ln16/rootFiles/pi0eta/malte_kmatrix_10M_030322/root/trees_pt0/tree*"',reconTreeName,"kmatrix_selected_halved",1,proof_Nthreads,recon_cfiles)
#runSelector('"/d/grid17/ln16/rootFiles/pi0eta/malte_kmatrix_10M_030322/root/trees_pt0/tree*"',reconTreeName,"kmatrix_selected_halved",2,proof_Nthreads,recon_cfiles)
#runSelector('"/d/grid17/ln16/rootFiles/pi0eta/malte_kmatrix_10M_030322/root/thrown_pt0/tree*"',thrownTreeName,"kmatrix_gen_halved",1,proof_Nthreads,thrown_cfiles)














