#!/usr/bin/python

import os


####### MAKE SURE WE PULL DOWN THE LATEST FIRST #####
baseDir="/d/grid17/ln16/dselector_v3/"

ofolder="etapi0_results"
os.system("mkdir -p "+ofolder+"/a2_xsec")
os.system("mkdir -p "+ofolder+"/double_regge")
os.system("mkdir -p "+ofolder+"/event_selection")

a2_flocs=[
        ########
        "study_combinatorics/comparison.pdf", # Plots comparing kmatrix matched thrown to kmatrix reconstructed
        ########
        "study_pwa/jupyter_plotting/results/*pdf", # partial wave analysis plots
        ]

evtSel_flocs=[
        ########
        "study_expectedYields/zExpectedLeakage/montage/montage_Mpi0eta_AccWeight_sigRegion.png", # Expected leakage from background channels
        "study_expectedYields/zExpectedLeakage/montage/montage_Mpi0eta_weightASBS.png", # same as above but with sidebands subtracted
        "study_expectedYields/zExpectedLeakage/zB1_efficiency/efficiency_AccWeight_sigRegion.png", # Cross sections + expected leakage calculation for b1
        "study_expectedYields/zExpectedLeakage/zB1_efficiency/efficiency_weightASBS.png", # same as above but with sidebands subtracted
        ########
        "study_eventSelections/event_selections/*", # Plots for event selections - from jupyter notebook
        ########
        "study_lmac/results/lowMassAltComboSelect.png", # low mass alternative cut, comparing flat etapi MC to omega as etapi MC when applying lmac cut
        ]

dr_flocs=[
        ]

for floc in a2_flocs:
    os.system("cp -r "+baseDir+floc+" "+ofolder+"/a2_xsec")

for floc in dr_flocs:
    os.system("cp -r "+baseDir+floc+" "+ofolder+"/double_regge")

for floc in evtSel_flocs:
    os.system("cp -r "+baseDir+floc+" "+ofolder+"/event_selection")

