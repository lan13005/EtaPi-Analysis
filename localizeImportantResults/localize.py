#!/usr/bin/python

import os

baseDir="/d/grid17/ln16/dselector_v3/"

ofolder="etapi0_results"
os.system("mkdir -p "+ofolder)

flocs=[
        ########
        "./study_combinatorics/comparison.pdf", # Plots comparing kmatrix matched thrown to kmatrix reconstructed

        ########
        "./study_expectedYields/zExpectedLeakage/montage/montage_Mpi0eta_AccWeight_sigRegion.png", # Expected leakage from background channels
        "./study_expectedYields/zExpectedLeakage/montage/montage_Mpi0eta_weightASBS.png", # same as above but with sidebands subtracted
        "./study_expectedYields/zExpectedLeakage/zB1_efficiency/efficiency_AccWeight_sigRegion.png", # Cross sections + expected leakage calculation for b1
        "./study_expectedYields/zExpectedLeakage/zB1_efficiency/efficiency_weightASBS.png", # same as above but with sidebands subtracted

        ########
        "./study_eventSelections/event_selections/*", # Plots for event selections - from jupyter notebook

        ########
        "./study_lmac/lowMassAltComboSelect.png", # low mass alternative cut, comparing flat etapi MC to omega as etapi MC when applying lmac cut

        ########
        "./study_pwa/jupyter_plotting/results/*pdf", # partial wave analysis plots
        ]

for floc in flocs:
    os.system("cp -r "+baseDir+floc+" "+ofolder)

