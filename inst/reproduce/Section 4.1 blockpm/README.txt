This folder contains the code to reproduce the
experiments of Section 4.1 in
"Unbiased Markov chain Monte Carlo for intractable target distributions"
by Middleton, Deligiannidis, Doucet & Jacob, 2018

File locations are specified in 'exp_settings.R' within a function.
This is called by each of the scripts in the folder.

Experiment settings are set to those used for the plots, though it may
be advisable to play around with for example smaller values to begin with.
Variable 'K' corresponds to the integer m, when specifying estimator H_{k:m}

The package relies on "TruncatedNormal", developed by Zdravko Botev [2015].
It enables precise inversion of truncated Gaussian in the tails of the
distribution.

Computation is performed in these main scripts

1. election_vote_compare_meeting_times.R
2. election_vote_Hkm_estimators.R
3. election_vote_compare_serial.R


### 1. election_vote_compare_meeting_times.R
This vote gets meeting times for coupled PM and coupled block PM
The script performs the computation then produces the results, without
additional processing.

The script also estimates the variance of the log-likelihood used in the article
close to the mode.

### 2. election_vote_Hkm_estimators.R
This vote gets the unbiased estimators for both standard and block PM.
The script performs both the computation and produces the results at
the end.

### 3. election_vote_compare_serial.R
This performs a long run of serial block pseudo-marginal to estimate
the asymptotic variance used in the inefficiency.


To usefully use these scripts it is recommended to first gain familirity with
election_vote_compare_serial.R and then play around with the values  election_vote_Hkm_estimators.R


