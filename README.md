# EtaPi0 Analysis at GlueX

DSelector_etapi.C/h provide an implementation of a DSelector for the reaction 
\gamma p\rightarrow\pi^0\eta p \rightarrow 4\gamma p

To run the DSelector, runDSelector.C is used. Inside this macro, the physics analysis trees can be run over by Adding
them to a TChain. Once in a chain, PROOF can be used to run the DSelector multi-threaded for quicker runtimes. PROOF 
will output a root file of the reduced tree and a root file of the histograms made during the process. Once you are satisfied
with the program inputs, you can simply run it by

root -l -b -q runDSelector.C

Similarly, you can also create a DSelector for thrown trees and also another runDSelector for that. These are affixed with the thrown tag
To ease the running of these DSelectors, runDSelectors.py can be used. This is a simple python script that simply executes 
the runDSelector programs over some specified files

./runDSelectors.py

It is helpful to visualize the event selections to see what the effects are. DSelector_etapi makes various histograms
including most of the relevant variables that are used in the event selection process. It is also quite simple to draw
shaded boxes and functions to denote the selection regions. drawCuts.C is a macro that does this. It is another 
argument-less macro (therefore you need to update it to point to the correct files) and can be run simply as

root -l -b -q drawCuts.C

AMPTOOLS REQUIRES 4 TREES (data, bkgnd, accmc, genmc)
- data = selected trees of DATA where signal region has been selected, all weights = 1
- bkgnd = selected trees of DATA where sidebands have been selected, all weights = -weight
- accmc = selected trees of ACCEPTANCE MC where signal+sidebands have been selected, all weights = weight
- genmc = thrown trees created during simulation process
AMPTOOLS can also use polarization information to separate production mechanisms
If your data or bkgnd contains multiple polarizations (i.e. in this tutorial spring 2017 gluex data contains 4 different polarization orientations + amorphous runs)
then you can split the root tree into the different orientations using
Both these points are taken care of by split_flat_kinematics.C program which can split a single tree into separate parts.
For instance, (polarization, Mpi0eta, t). The program needs to be modified for your purposes 

root -l -b -q split_flat_kinematics.C
