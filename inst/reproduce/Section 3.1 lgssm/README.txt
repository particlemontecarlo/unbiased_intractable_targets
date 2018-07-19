This folder contains the code to reproduce the
experiments of Section 3.1 in
"Unbiased Markov chain Monte Carlo for intractable target distributions"
by Middleton, Deligiannidis, Doucet & Jacob, 2018

"lgssm_2params_100obs.R" includes the data used throughout this section,
though this can be equivalently simulated using "gen_data.R"

"exp_settings.R" includes the experiment settings used to obtain the plots.
It must be sourced before running any of the scripts in this folder.


Throughout we use the variable 'K' to denote the integer 'm' in the H_{k:m}
estimators.


####### Number of particles
The bulk of computation is performed in the following three scripts

1. "lg_ssm_different_N.R" - this is used to get meeting times on a grid of N (approx run time=2 days)
2. "lg_ssm_different_N_Hkm.R"  - this is used to unbiased estimators on a grid of N (approx run time=7 days)
3. "lg_ssm_different_N_serial.R" - this is used to perform the long runs of the serial algorithm (approx run time=2 days)

Each of these were run on a smaller server (~50 cores), so the run times are more indicative relative to each other.

Having run all these scripts, "lg_ssm_Hkm_process_results.R" produces
the plots of the outputs.

Running these as they are will be relatively expensive. All the coupled chains are saved. Then, having obtained them,
unbiased estimators of a certain function and plots are made in a single script (lg_ssm_Hkm_process_results.R).



####### Scaling with T
The final script "lg_ssm_scalingT.R" is used to produce the figures in
Section 3.1.2, examining how the distribution of meeting times change
with the number of observations. The script performs the computation and
produces the plots.

