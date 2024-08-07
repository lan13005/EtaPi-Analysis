#!/usr/bin/python

import os
import sys

major_version, minor_version, micro_version = sys.version_info[:3]
if major_version == 3 and minor_version > 9:
    print("rcdb (as of 6/26/24) requires python 3.9 or lower. Its a simple python import issue that can be fixed if you want to rebuild rcdb")
    print("Try using /usr/bin/python")
    print("EXITING...")
    exit(1)

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

proof_Nthreads=16
recon_cfiles=["DSelector_etapi.C", "runDSelector.C"]
thrown_cfiles=["DSelector_thrown.C", "runDSelector_thrown.C"]

reconTreeName="pi0eta__B4_M17_M7_Tree"
thrownTreeName="Thrown_Tree"

tag=""
#tag="_sbL_accN"
#tag="_sidebandStudy"
#tag="_wUnusedShowers_sbN_accT"
#tag="_bggen_2018_8"


### a0a2 bw study for NIFTy
# for i in [1,2,3]:
#     runSelector('"/w/halld-scshelf2101/lng/DATA/a0a2_2bw/root/trees/tree_pi0eta*"',reconTreeName,"a0a2_2bw",i,proof_Nthreads,recon_cfiles)
# runSelector('"/w/halld-scshelf2101/lng/DATA/a0a2_2bw/root/thrown/tree_thrown*"',thrownTreeName,"a0a2_2bw_gen",1,proof_Nthreads,thrown_cfiles)
# for i in [1,2,3]: # First MC sample. Turns out it is about 1:1 in size to the above a0a2 sample in the peak region
#     runSelector('"/w/halld-scshelf2101/lng/DATA/a0a2_flat/root/trees/tree_pi0eta*"',reconTreeName,"a0a2_flat",i,proof_Nthreads,recon_cfiles)
# runSelector('"/w/halld-scshelf2101/lng/DATA/a0a2_flat/root/thrown/tree_thrown*"',thrownTreeName,"a0a2_flat_gen",1,proof_Nthreads,thrown_cfiles)

# for i in [1,2,3]:
#     runSelector('"/w/halld-scshelf2101/lng/DATA/a0a2_bin_10_1M/trees/tree_pi0eta__B4_M17_M7_gen_amp/tree_pi0eta__B4_M17_M7_gen_amp_*.root"',reconTreeName,"a0a2_bin_10",i,proof_Nthreads,recon_cfiles)
# for i in [1]:
#     runSelector('"/w/halld-scshelf2101/lng/DATA/a0a2_bin_10_1M/thrown/tree_thrown_gen_amp_*"',thrownTreeName,"a0a2_bin10_gen",1,proof_Nthreads,thrown_cfiles)
# for i in [1,2,3]:
#     runSelector('"/w/halld-scshelf2101/lng/DATA/a0a2_flat_100M/trees/tree_pi0eta__B4_M17_M7_gen_amp/tree_pi0eta__B4_M17_M7_gen_amp_*"',reconTreeName,"flat",i,proof_Nthreads,recon_cfiles)
# for i in [1]:
#     runSelector('"/w/halld-scshelf2101/lng/DATA/a0a2_flat_100M/thrown/tree_thrown_gen_amp_*"',thrownTreeName,"flat_gen",1,proof_Nthreads,thrown_cfiles)


# FOR BORIS TEST 04/25/24
# runSelector('"/w/halld-scshelf2101/lng/DATA/tSlopeSampler/2017_130M/tree_builder/merged/tree_thrown_merged0.root"',thrownTreeName,"test_forBoris_2017_1_t010020_m080250",1,proof_Nthreads,thrown_cfiles)
# runSelector('"/w/halld-scshelf2101/lng/DATA/tSlopeSampler/2018_1_400M//tree_builder/merged/tree_thrown_merged0.root"',thrownTreeName,"test_forBoris_2018_1",1,proof_Nthreads,thrown_cfiles)
# runSelector('"/w/halld-scshelf2101/lng/DATA/tSlopeSampler/2018_8_260M//tree_builder/merged/tree_thrown_merged0.root"',thrownTreeName,"test_forBoris_2018_8",1,proof_Nthreads,thrown_cfiles)

# phase 1 data on ifarm
# runSelector('"/cache/halld/RunPeriod-2018-08/analysis/ver02/tree_pi0eta__B4_M17_M7/merged/*"',reconTreeName,"D2018_8_selected",3,proof_Nthreads,recon_cfiles)

### PHASE 1 MONTE CARLO
#runSelector('"/d/grid17/ln16/rootFiles/pi0eta/120921/2017_1_130M/merged/tree_pi0eta*"',reconTreeName,"F2017_1_selected"+tag,3,proof_Nthreads,recon_cfiles)
#runSelector('"/d/grid17/ln16/rootFiles/pi0eta/120921/2018_1_400M/merged/tree_pi0eta*"',reconTreeName,"F2018_1_selected"+tag,3,proof_Nthreads,recon_cfiles)
#runSelector('"/d/grid17/ln16/rootFiles/pi0eta/120921/2018_8_260M_130M/merged/tree_pi0eta*"',reconTreeName,"F2018_8_selected"+tag,3,proof_Nthreads,recon_cfiles)

## Another set of prefiltered datasets with loose sb and acc regions
#for i in [1,2,3]:
#    runSelector('"/w/halld-scshelf2101/lng/DATA/phase1_selected_v4/D2017_1_selected_wUnusedShowers_sbL_accN_acc_tree.root"',reconTreeName,"D2017_1_selected"+tag,i,proof_Nthreads,recon_cfiles)
#    runSelector('"/w/halld-scshelf2101/lng/DATA/phase1_selected_v4/D2018_1_selected_wUnusedShowers_sbL_accN_acc_tree.root"',reconTreeName,"D2018_1_selected"+tag,i,proof_Nthreads,recon_cfiles)
#    runSelector('"/w/halld-scshelf2101/lng/DATA/phase1_selected_v4/D2018_8_selected_wUnusedShowers_sbL_accN_acc_tree.root"',reconTreeName,"D2018_8_selected"+tag,i,proof_Nthreads,recon_cfiles)
#for i in [3]:
#    runSelector('"/w/halld-scshelf2101/lng/DATA/phase1_selected_v4/F2017_1_selected_wUnusedShowers_sbL_accN_acc_tree.root"',reconTreeName,"F2017_1_selected"+tag,i,proof_Nthreads,recon_cfiles)
#    runSelector('"/w/halld-scshelf2101/lng/DATA/phase1_selected_v4/F2018_1_selected_wUnusedShowers_sbL_accN_acc_tree.root"',reconTreeName,"F2018_1_selected"+tag,i,proof_Nthreads,recon_cfiles)
#    runSelector('"/w/halld-scshelf2101/lng/DATA/phase1_selected_v4/F2018_8_selected_wUnusedShowers_sbL_accN_acc_tree.root"',reconTreeName,"F2018_8_selected"+tag,i,proof_Nthreads,recon_cfiles)
### Thrown trees generally need less threads to run over since there isn't much calculation being done
#runSelector('"/w/halld-scshelf2101/lng/DATA/tSlopeSampler/2017_130M/tree_builder/merged/tree_thrown*"',thrownTreeName,"F2017_1_gen",1,proof_Nthreads,thrown_cfiles)
runSelector('"/w/halld-scshelf2101/lng/DATA/tSlopeSampler/2018_1_400M/tree_builder/merged/tree_thrown*"',thrownTreeName,"F2018_1_gen",1,proof_Nthreads,thrown_cfiles)
runSelector('"/w/halld-scshelf2101/lng/DATA/tSlopeSampler/2018_8_130M_260M/tree_builder/merged/tree_thrown*"',thrownTreeName,"F2018_8_gen",1,proof_Nthreads,thrown_cfiles)

# FILTERED LOOSE CHI UE PHASE 1 DATA. HAS NEAREST RF BUNCH SKIPPED
#for i in [3]:
#    runSelector('"./phase1_data_looseChiUE/D2017_1_loose_tree.root"',reconTreeName,"D2017_1_selected",i,proof_Nthreads,recon_cfiles)
#    runSelector('"./phase1_data_looseChiUE/D2018_1_loose_tree.root"',reconTreeName,"D2018_1_selected",i,proof_Nthreads,recon_cfiles)
#    runSelector('"./phase1_data_looseChiUE/D2018_8_loose_tree.root"',reconTreeName,"D2018_8_selected",i,proof_Nthreads,recon_cfiles)

### FOR THE DOUBLE REGGE STUDY
##for i in [3]:
##    runSelector('"./study_double_regge/rootFiles/D2017_1_selected_acc_tree.root"',reconTreeName,"D2017_1_selected",i,proof_Nthreads,recon_cfiles)
##    runSelector('"./study_double_regge/rootFiles/D2018_1_selected_acc_tree.root"',reconTreeName,"D2018_1_selected",i,proof_Nthreads,recon_cfiles)
##    runSelector('"./study_double_regge/rootFiles/D2018_8_selected_acc_tree.root"',reconTreeName,"D2018_8_selected",i,proof_Nthreads,recon_cfiles)
#

# UNFILTERED PHASE 1 DATA
#for i in [3]:
#    runSelector('"/d/grid17/ln16/rootFiles/pi0eta/RunPeriod-2017-01/analysis-ver52/tree_pi0eta__B4_M17_M7/merged/*"',reconTreeName,"D2017_1_selected",i,proof_Nthreads,recon_cfiles)
#    runSelector('"/d/home/sdobbs/GlueX/gluex_data/RunPeriod-2018-01/analysis-ver02/tree_pi0eta__B4_M17_M7/merged/*"',reconTreeName,"D2018_1_selected",i,proof_Nthreads,recon_cfiles)
#    runSelector('"/d/home/sdobbs/GlueX/gluex_data/RunPeriod-2018-08/analysis-ver02/tree_pi0eta__B4_M17_M7/merged/*"',reconTreeName,"D2018_8_selected",i,proof_Nthreads,recon_cfiles)

#### KMATRIX MC
#reconTreeName="pi0eta__B4_M7_M17_Tree"
#for i in [1,2,3]:
#    runSelector('"/d/grid17/ln16/rootFiles/pi0eta/malte_kmatrix_10M_030322/root/trees/tree*"',reconTreeName,"kmatrix_selected",i,proof_Nthreads,recon_cfiles)
#runSelector('"/d/grid17/ln16/rootFiles/pi0eta/malte_kmatrix_10M_030322/root/thrown/tree*"',thrownTreeName,"kmatrix_gen",1,proof_Nthreads,thrown_cfiles)

### PHASE 1 b1 MC
#runSelector('"/d/grid17/ln16/rootFiles/omegapi_rad_massDepFit_2017_01_ver03.0/root/merged/tree_pi0eta*"',
#        reconTreeName,"BOne2017_1_selected"+tag,3,proof_Nthreads,recon_cfiles)
#runSelector('"/d/grid17/ln16/rootFiles/omegapi_rad_massDepFit_2018_01_ver03.0/root/merged/tree_pi0eta*"',
#        reconTreeName,"BOne2018_1_selected"+tag,3,proof_Nthreads,recon_cfiles)
#runSelector('"/d/grid17/ln16/rootFiles/omegapi_rad_massDepFit_2018_08_ver03.0_60M/root/merged/tree_pi0eta*"',
#        reconTreeName,"BOne2018_8_selected"+tag,3,proof_Nthreads,recon_cfiles)
#runSelector('"/d/grid17/ln16/rootFiles/omegapi_rad_massDepFit_2017_01_ver03.0/root/merged/tree_thrown*"',
#        thrownTreeName,"BOne2017_1_gen"+tag,1,proof_Nthreads,thrown_cfiles)
#runSelector('"/d/grid17/ln16/rootFiles/omegapi_rad_massDepFit_2018_01_ver03.0/root/merged/tree_thrown*"',
#        thrownTreeName,"BOne2018_1_gen"+tag,1,proof_Nthreads,thrown_cfiles)
#runSelector('"/d/grid17/ln16/rootFiles/omegapi_rad_massDepFit_2018_08_ver03.0_60M/root/merged/tree_thrown*"',
#        thrownTreeName,"BOne2018_8_gen"+tag,1,proof_Nthreads,thrown_cfiles)

### f2 to pi0pi0 MC
#runSelector('"/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/pi0pi0/tree_pi0eta*"',
#        reconTreeName,"FTwo2017_1_selected"+tag,3,proof_Nthreads,recon_cfiles)
### etaprime to etapi0pi0 MC
#runSelector('"/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/etap_to_etapipi/tree_pi0eta*"',
#        reconTreeName,"etap6g2017_1_selected"+tag,3,proof_Nthreads,recon_cfiles)
### eta to 3pi0 MC
#runSelector('"/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/eta_to_3pi/tree_pi0eta*"',
#        reconTreeName,"eta6g2017_1_selected"+tag,3,proof_Nthreads,recon_cfiles)
### a2pi to etapi0pi0 MC
#runSelector('"/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/a2pi/tree_pi0eta*"',
#        reconTreeName,"a2pi6g2017_1_selected"+tag,3,proof_Nthreads,recon_cfiles)
### f1 to etapi0pi0 MC
#runSelector('"/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/f1_to_etapipi/tree_pi0eta*"',
#        reconTreeName,"f16g2017_1_selected"+tag,3,proof_Nthreads,recon_cfiles)

## Exotic review 2022
#    runSelector('"/d/grid17/ln16/rootFiles/pi0eta/exoticReview2022/ExoticReview2022_A2_EtaPi0_Spring2017_20221021063447pm/root/trees/tree_pi0eta*"',
#            reconTreeName,"ExoticReviewA2_2017_1_selected",3,proof_Nthreads,recon_cfiles)
#    runSelector('"/d/grid17/ln16/rootFiles/pi0eta/exoticReview2022/ExoticReview2022_PI1_EtaPi0_Spring2017_20221020103256am/root/trees/tree_pi0eta*"',
#            reconTreeName,"ExoticReview_2017_1_selected",3,proof_Nthreads,recon_cfiles)
#    runSelector('"/d/grid17/ln16/rootFiles/pi0eta/exoticReview2022/ExoticReview2022_A2_EtaPi0_Spring2018_20221021063729pm/root/trees/tree_pi0eta*"',
#            reconTreeName,"ExoticReviewA2_2018_1_selected",3,proof_Nthreads,recon_cfiles)
#    runSelector('"/d/grid17/ln16/rootFiles/pi0eta/exoticReview2022/ExoticReview2022_PI1_EtaPi0_Spring2018_20221020102347am/root/trees/tree_pi0eta*"',
#            reconTreeName,"ExoticReview_2018_1_selected",3,proof_Nthreads,recon_cfiles)
#    runSelector('"/d/grid17/ln16/rootFiles/pi0eta/exoticReview2022/ExoticReview2022_A2_EtaPi0_Fall2018_20221021064015pm/root/trees/tree_pi0eta*"',
#            reconTreeName,"ExoticReviewA2_2018_8_selected",3,proof_Nthreads,recon_cfiles)
#    runSelector('"/d/grid17/ln16/rootFiles/pi0eta/exoticReview2022/ExoticReview2022_PI1_EtaPi0_Fall2018_20221020095551am/root/trees/tree_pi0eta*"',
#            reconTreeName,"ExoticReview_2018_8_selected",3,proof_Nthreads,recon_cfiles)
#runSelector('"/d/grid17/ln16/rootFiles/pi0eta/exoticReview2022/ExoticReview2022_A2_EtaPi0_Spring2017_20221021063447pm/root/thrown/tree_thrown*"',thrownTreeName,
#        "ExoticReviewA2_2017_1_gen",1,proof_Nthreads,thrown_cfiles)
#runSelector('"/d/grid17/ln16/rootFiles/pi0eta/exoticReview2022/ExoticReview2022_PI1_EtaPi0_Spring2017_20221020103256am/root/thrown/tree_thrown*"',thrownTreeName,
#        "ExoticReview_2017_1_gen",1,proof_Nthreads,thrown_cfiles)
#runSelector('"/d/grid17/ln16/rootFiles/pi0eta/exoticReview2022/ExoticReview2022_A2_EtaPi0_Spring2018_20221021063729pm/root/thrown/tree_thrown*"',thrownTreeName,
#        "ExoticReviewA2_2018_1_gen",1,proof_Nthreads,thrown_cfiles)
#runSelector('"/d/grid17/ln16/rootFiles/pi0eta/exoticReview2022/ExoticReview2022_PI1_EtaPi0_Spring2018_20221020102347am/root/thrown/tree_thrown*"',thrownTreeName,
#        "ExoticReview_2018_1_gen",1,proof_Nthreads,thrown_cfiles)
#runSelector('"/d/grid17/ln16/rootFiles/pi0eta/exoticReview2022/ExoticReview2022_A2_EtaPi0_Fall2018_20221021064015pm/root/thrown/tree_thrown*"',thrownTreeName,
#        "ExoticReviewA2_2018_8_gen",1,proof_Nthreads,thrown_cfiles)
#runSelector('"/d/grid17/ln16/rootFiles/pi0eta/exoticReview2022/ExoticReview2022_PI1_EtaPi0_Fall2018_20221020095551am/root/thrown/tree_thrown*"',thrownTreeName,
#        "ExoticReview_2018_8_gen",1,proof_Nthreads,thrown_cfiles)

## BGGEN 2018_8
#runSelector('"/d/grid17/ln16/rootFiles/pi0eta/bggen_2018_8/symlinked_trees/tree_pi0eta__B4_M17_M7_*"',
#        reconTreeName,"BGGEN2018_8_selected",3,proof_Nthreads,recon_cfiles)













