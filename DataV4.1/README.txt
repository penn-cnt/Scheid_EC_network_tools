Data folder readme:

Networks.mat  (11/20/19): All patients, networks generated from UEO clean data, gamma=0.5
Partitions.mat (11/20/19): Partitions generated from Networks.mat, gamma = 0.5, beta=0.01
Energy.mat (2/8/20): Energy with single-point control energy
EnergySOZ.mat (2/8/20): Energy with SOZ as control points, using average bandpower across entire preictal period as a reference, input is gaussian, A_s is NOT normalized!
EenrgySOZ_B1.mat (2/8/20): settings same as above, however B=1 at control points, 'relax=10^-5' elsewhere along diagonal, A_s was normalized, 1000 permutations 
EenrgySOZ_B1.mat (2/19/20): settings same as above, B=1 at control points, 'relax=0' elsewhere along diagonal, A_s was normalized, 100 permutations 
EnergySOZ_spread.mat (2/8/20): same as EnergySOZ_B1.mat, but used Gaussian Spread with parameter sigma, A_s was normalized, 1000 permutations


dataSets_clean.mat: 
- EECStart: earliest electrographic onset time on IEEG portal (sec)
- UEOstart: unequivical electrographic onset time on IEEG portal (sec)
- all_data: all data (cleaned?) from EEC to end across ALL channels
- data: all cleaned data from EEC to end across grid channels with toIgnore channels removed
- sozGrid: soz index on grid w/ toIgnore removed
- gridCoords: grid coorinates for each channel w/ toIgnore removed
- toIgnore: all channelnames to ignore
- cleanGridNames: names of grids with toIgnore removed

Metric_matrices.mat:


State_metrics.mat:
- Each metric column includes the average of metric values from all networks in each of three largest contiguous states.  
- The Zscored column includes the average of metric values from all networks in each of three largest contiguous states AFTER
all metric values have been zscored across the entire matrix. 

Energy.mat:
- x0: initial state, mean high gamma bandpower across all windows in each phase
- repMats: representative EC network for each seizure phase (NOT NORMALIZED). This is the network with the greatest average similarity to all other networks in the phase.  
- xf: final state, mean high gamma across all windows in preictal phase
- x0_z: x0 after zscore across xf and original x0
- xf_z: xf after zscore across xf and original x0
- t_traj: the vector of time horizons calculated
- rho: the vector of rho values used
NOTE: the following do not include i_ict for the 5 removed siezures in their computation:
- sxtrajErr: computation error for (i_driven node, i_traj, i_rho)
- sxNodeEnergy: node Energy arranged as (i_driven node, i_traj, i_rho)
- err_stats: compiled computational error across all i_ict subjects and states. 
- nullsets: indices of nodes randomly selected to be soz nodes
- SOZconfidence: The percentile of null soz energy wrt the real SOZ energy
- nodeEnergynull: optEnergy calculated for each null.

patients.mat: (not included, relevent version in DataV3.2), includes ALL channel names and ALL gridNames/coords (note ignored coords are included)
