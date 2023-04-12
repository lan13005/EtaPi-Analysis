./drawOverview.py systematic_nominal_v9/nominal -s systematic_nominal_v9/systematic_v24/SOURCE/systematic_syst/barlows.csv -b systematic_nominal_v9/bootstrap/bootstrap_results/
rm -rf systematic_nominal_v9/nominal/xsec_bootstrap_with_syst
mv systematic_nominal_v9/nominal/xsec systematic_nominal_v9/nominal/xsec_bootstrap_with_syst

#./drawOverview.py systematic_nominal_v9/nominal -b systematic_nominal_v9/bootstrap/bootstrap_results/
#rm -rf systematic_nominal_v9/nominal/xsec_bootstrap_no_syst
#mv systematic_nominal_v9/nominal/xsec systematic_nominal_v9/nominal/xsec_bootstrap_no_syst
#
#./drawOverview.py systematic_nominal_v9/nominal -s systematic_nominal_v9/systematic_v24/SOURCE/systematic_syst/barlows.csv
#rm -rf systematic_nominal_v9/nominal/xsec_minuit_with_syst
#mv systematic_nominal_v9/nominal/xsec systematic_nominal_v9/nominal/xsec_minuit_with_syst
#
#./drawOverview.py systematic_nominal_v9/nominal
#rm -rf systematic_nominal_v9/nominal/xsec_minuit_no_syst
#mv systematic_nominal_v9/nominal/xsec systematic_nominal_v9/nominal/xsec_minuit_no_syst
