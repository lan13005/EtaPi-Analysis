The best fit files are taken from malte_v3_m104172
These are the best converged files, does not have to properly reconverge

etapi_hybrid_t010020_m104172/etapi_result_117.fit | nll: -1103254 : deltabestnll: 0 : minuit 0 : ematrix 3 : converged!
etapi_hybrid_t0200325_m104172/etapi_result_89.fit | NLL: -815841 : DeltaBestNLL: 3 : minuit 0 : ematrix 3 : converged!
etapi_hybrid_t0325050_m104172/etapi_result_172.fit | NLL: -535124 : DeltaBestNLL: 0 : minuit 0 : ematrix 3 : converged!
etapi_hybrid_t050075_m104172/etapi_result_181.fit | NLL: -393489 : DeltaBestNLL: 0 : minuit 0 : ematrix 3 : converged!
etapi_hybrid_t075100_m104172/etapi_result_289.fit | NLL: -267803 : DeltaBestNLL: 0 : minuit 0 : ematrix 3 : converged!

bootstrap_results contains fit results that malte ran for me with the configuration files that I made

systematics_v24 contains the first large scale systematics scan for this updated nominal fit results
A perturbation scan is performed so that we can atleast get converged fits for these systematic scans

bootstrap folder contains bootstrapped results that can be plotted and reoraganized into several csv files that contain a similar structure to etapi_plotter_output

cfg_original and cfg_wUnusedShowers are basically the same. The latter updates the root files to ones that include the # unused showers
