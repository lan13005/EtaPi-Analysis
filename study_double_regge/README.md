Code to study the double regge process. The goal is to extract the beam asymmetry in various kinematic bins

fitAsymmetryPlots.C will loop through the provided data and perform selections (i.e. for event selection systematics). This selected data will then be used to fill the phi histograms which can then be fit to the "shifted-cosine" equation or the "yield asymmetry" equation. The specific variation, i.e. event selection, freeing phase offset, etc, can be set in the main program which then uses extractAsymmetries function to load the appropriate data. This program then uses constructAndFit function to fill the phi histograms and fit the distributions and saves all the final results.
 
Then fitAsymmetryPlots.ipynb is used to load all the data and organize the plot outputs.

### Versions
* Current version in folder: to ease the discussion, we will set chiSq<13.277 instead of 16. Will also remove the pi0 sidebands
* v0 - Final results before 03/16/23

### NOTES 
Include a selection on teta<1 OR tpi0<1 as we are only interested in fast pseudoscalars

MAINLY FOR EVENT SELECTION SYSTEMATICS
* No omega cut, not useful for DR studies
* No Delta cut, nor VH, nor t cuts
* 8.2 < E < 8.8 do not check systematics of coherent peak selection
* Select Double Regge region Mpi0eta [1.6, 3.0]
* ChiSq<40 or 16? Need to check again
* UE<10, basically freed
* The double regge region is pretty clean and it shows in the little # events in UE and chiSq that has larger than nominal selections
* Looser “standard” event selections
* Same photon E selection (since its very tight already)
* Same dEdx selection (weird to vary since its a functional form)
