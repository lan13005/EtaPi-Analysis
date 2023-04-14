#ifndef DSELECTOR_HELPERFUNCS
#define DSELECTOR_HELPERFUNCS
// PUT THINGS HERE THAT ARE NEEDED BY THE DSELECTOR PROGRAM + OTHER PROGRAMS.
// FOR EXAMPLE FOR SIDEBAND SUBTRACTION WE NEED TO DEFINE THE REGIONS WHICH WILL BE USED TO CALCULATE THE SIDEBAND WEIGHTS BUT ALSO
// NEED TO THE REGIONS IN THE MAKEGRAPHS.C CODE TO DRAW THE LINES


//Variables for determining sideband subtraction regions.
double pi0Mean=0.135881;
double pi0Std=0.0076;
double pi0Sig=3;
double pi0Skip=1;

double etaMean=0.135881;
double etaStd=0.0076;
double etaSig=3;
double etaSkip=1;

double pi0SB=1;
double pi0SBweight= -1.5;
double etaSB=1;
double etaSBweight= -1.5;
#endif
