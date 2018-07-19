This folder contains the code to reproduce the 
experiments of Section 3.2 in 
"Unbiased Markov chain Monte Carlo for intractable target distributions"
by Middleton, Deligiannidis, Doucet & Jacob, 2018

The data is in thaldata.csv
The file "binomialmodel.R" loads the data, and defines the model.
As well as the functions to run controlled Sequential Monte Carlo,
and particle marginal MH using controlled SMC. Finally it defines functions
to run coupled chains.

The file "binomialmodel.bpf.R" does the same but with bootstrap particle
filter instead of controlled SMC.

The script "comparecost.bpfvscsmc.R" runs both BPF and CSMC a number
of times to make sure that the costs are comparable, i.e. that the
numbers of particles for each are well chosen.


The script plot.datalikelihood.R creates Figure 4 of the manuscript.
It requires having run "csmc.2dgrid.R" before, which runs
controlled SMC on a 500x500 grid of parameter values.

The scripts "pmmh.csmc.binomial.*" run a long PMMH chain
using controlled SMC (CSMC) at each step, with two different choices
of proposal std deviation. Takes a few days to run.
The script pmmh.bpf.binomial.R does the same but with bootstrap
particle filter (BPF) instead of controlled SMC.

The files
run.odyssey.binomial.sh
run.odyssey.R
run.odyssey.bpf.R
load.odyssey.R
have to do with running scripts on the cluster. They cannot be used "as is"
but skimming through them will be informative to see how the files are produced.
For instance run.odyssey.R as it stands will create meeting times,
during a period of 23 hours, from coupled PMMH using controlled SMC.
Simple modifications of the script can be made to produce estimators,
for instance by replacing "m <- 0" by "m <- 10000".


plot.meetingtimes.R processes the output of runs, 
and creates plots for Figures 5 and 7 in the manuscript.

plot.traces.R processes the output of runs, 
and creates plots for Figure 6.

plot.efficiency.R processes the output of runs,
creates plots for Figure 8 and the numbers on the asymptotic inefficiency
reported in Section 3.2.4.
