
bootstrap_t010020 are made by me and I think results_t010020 are made by Malte
The files like "bootstrap_results/bootstrap_result_mean_std_t010020.log" do not have the uncertainties on the parameter values (like a2 mass, and bin 0 S-wave) only the fit fractions. This is the only thing that matters for the cross section measurement 

Use fit_bootstrap.py to run bootstrap fits with different seeds. It should result in a bunch of folders. 

You should use run_overlayBins.py program to extract the fit fractions. I forgot if this should be before or after reorganize.py.

Then you can use constructUncertainties.py to create a folder of files that looks a lot like the output of etapi_plotter but now contains the boostrapped statistical uncertainites. This format can be loaded easily by helper.py to load in the bootstrapped uncertainties 


02/08/23
Running bootstrap for the 5 different t-bins where the nominal fits
have been set to the best converged likelihood fit from malte_v3 (>250 randomized fits in 1.04-1.72 mass region)
Did not check for reconvergence, we think we can fix it using perturbation scans
