Various channels can leak into the etapi->4g channel. For instance, b1->omegapi->5g where one photon is lost. The b1 
cross section has been measured at GlueX before. Simulations of b1->5g can be made to determine the efficiency 
of reconstructing a b1 as etapi. With the cross section and efficiency we have everything we need to determine an
expected leakage yield

This can also be done for other channels. Various simulations were made and some programs have been developed to 
run through all the trees, i.e. runSelectorOverBkgSamples.sh which will run the DSelectors and organize the results

Once all the trees exist, getEffs.py can be used to make plots of the expected leakage overlaid onto the phase 1 data 
NOTES:
1. Only some cross sections are measured at gluex. We are interested in upper limits so we can reuse cross sections if we knew the branching ratios of various processes
2. Most simulations that are currently used for the study is in the Fall 2018 configuration. Only b1 simulations include spring 2017/2018. This is the biggest leakage channel so a more detailed study is performed
3. DSelector_etapi.C keeps all Meta vs Mpi0. This allows us to see how, for instance, the b1 background shows up in Meta and Mpi 
4. getEffs.py program takes in cross section estimates and some recon+thrown trees of various MC samples to construct the efficiency. Knowing both of these makes it enough to extract an expected leakage yield
