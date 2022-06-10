From the DSelector we have access to the truth information
We can determine if the beam photon is the true in-time beam photon using accidental subtraction
We can determine if the photon pairings that form the eta and pion candidate are the true ones using sideband subtraction

From here it is simple to compare the results. uproot is used to load the flat trees that already contain a lot of
kinematic variables, the accidental+sideband weights, and the truth information. We can simply overlay histograms
of the kinematic variables comparing the truth with the results from weighting

I have also done studies with various schemes for uniqueness tracking which is incompatible with the
current weighting schemes. At this point, accidental+sideband subtraction performs pretty well, so 
there is no need to consider tracking schemes
