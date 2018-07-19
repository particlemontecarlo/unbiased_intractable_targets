This folder contains the code to reproduce the 
experiments of Section 4.2 in 
"Unbiased Markov chain Monte Carlo for intractable target distributions"
by Middleton, Deligiannidis, Doucet & Jacob, 2018

File locations are specified in 'exp_settings.R' within a function.
This is called by each of the scripts in the folder.

Experiment settings are set to those used for the plots, though it may 
be advisable to play around with for example smaller values to begin with.
Variable 'K' corresponds to the integer m, when specifying estimator H_{k:m}

Computation is performed in two main scripts

1. ising.R - this is where unbiased estimators are obtained
2. ising_serial.R - this performs long runs of the serial algorithm to compare
inefficiencies. 

Data is already provided in the "data" folder though this can be created using
the "ising_gendata.R" script.  

Finally, having run both these, plots are created using ising_process_results.R
