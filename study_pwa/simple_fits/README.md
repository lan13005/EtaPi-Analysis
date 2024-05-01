fit.C will perform a fit to the acceptance corrected mass distribution to extract the differential cross section of the a2(1320) so that it can be compared to the partial wave analysis results. Relativistic breit-wigner is implemented for use with RooFit and a convolution operation is performed to make a Voigtian. These distributions describe the mass distribution. The voigtian gaussian mass resolution is taken from another study I performed.

STEPS:
mkdir results

root -l
.L RooRelBreitWigner.C+
.q

root -l -b -q fit.C

. convertAndMerge.sh 