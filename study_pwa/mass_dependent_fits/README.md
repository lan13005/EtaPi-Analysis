setup_mass_dep_fits.py takes in a template config file and updates it to a functional form, where it will comment out 
unwanted lines and initiaze the complex amplitudes

Once a valid cfg file is made we can perform a mpi mass dependent fit using where -m is the maximum number of fit calls
and t is the tolerance. Our computer has 8 GPUs and we include an extra process as a leader node so np=num processes=9
mpirun -np 9 fitMPI -c etapi0_SD_TMD_allPol.cfg -m 80000 -t 0.1

A fit can fail for a number of reasons. We can reinitize the fit and refit in a while loop to search for convergence. Typically
we also care about running fits in multiple t-bins. To handle this fit.py program is made. At the end it will output a summary

./fit.py

summary2.py is made as a separate program that looks through some directories to check for statuses also. This is useful to check
how much times a fit has been run and to check the statuses. From here we can begin to diagnose any convergence issues

Once the fits have been performed we can use etapi_plotter program to draw the waves. The overlayBins set of programs will do this.
This set of programs will output a folder that compares various distributions with various partial wave components. From here you 
can quickly glance how the fits perform and how the partial waves look. The more important output of this program is a root file
that contains the shape of the partial wave. This can then be used to plot i.e. the piecewise s-wave and a breit-wigner D-wave 

./run_overlayBins.py

To get a cross section we need the integrated luminosity which we can obtain from plot_flux_ccdb.py program. getFluxs.py is
used to maintain the correct arguments to obtain the correct luminosity and to run the program


