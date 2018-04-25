# EFTofBOSS


Everything should be ready to run except for the files of the power spectra computed on the grid (very heavy). I will send one (LightConeDida, corresponding to a grid name LC0.61) through google drive.

For now, you should only run with a "simtype = LightConeDida". There are more row in the file "DataFrameCosmosims.csv" (which contains information on each simulation). You should also let the parameters "withBisp", "window" and "binning" as they are now.

If you have access to a lot of computer power, what you can do is write a code that runs Pierre's code (Power spectrum or bispectrum) on a grid of cosmological parameters.

Usuallt what I do is fix the ratio \omega_b/\omega_c to the true ratio (can be found in "DataFrameCosmosims.csv", or also in table 1 of arXiv:1706.02362) and run a grid for lnAs between 2.5 and 3.5 (100 points), Omega_m between 0.2 and 0.4 (50 points) and h between 0.63 and 0.77 (50 points). You can choose something else but be sure to document it.

Another idea would be to fix lnAs (that we don't constrain that well) and vary N_eff instead. If you have time that's worth exploring.

Concerning the MCMC itself, take your time to go through it. What it does is compute 4 chains in a row for minlength = 4000 steps, compute the Gelman Rubin criteria, and if it's not below a user-defined threshold (1 +\epsilon, where \epsilon = 0.06 here), it relaunches the sampler for 50 steps. This proceeds until convergence is reached.

For testing purposes, you might want to set \epsilon to very large value and minlength to a small one.