Various channels can leak into the etapi->4g channel. For instance, b1->omegapi->5g where one photon is lost. The b1 
cross section has been measured at GlueX before. Simulations of b1->5g can be made to determine the efficiency 
of reconstructing a b1 as etapi. With the cross section and efficiency we have everything we need to determine an
expected leakage yield

This can also be done for other channels. Various simulations were made and some programs have been developed to 
run through all the trees, i.e. runSelectorOverBkgSamples.sh which will run the DSelectors and organize the results

Once all the trees exist, getEffs.py can be used to make plots of the expected leakage overlaid onto the phase 1 data 
