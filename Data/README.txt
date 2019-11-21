Data folder readme:

Networks.mat: All patients, networks generated from UEO clean data, gamma=0.5
Partitions.mat: Partitions generated from Networks.mat, with beta=0.01. 
EnergySOZInf.mat: Energy with SOZ as control points, using average bandpower across entire preictal period as a reference
EenrgySOZInf_B1.mat: settings same as above, however B=1 at control points, zero else. 
EnergySOZInf_spread.mat: same as EnergySOZInf.mat... no idea. 
Partitions.mat: All patients, gamma = 0.5, beta=0.01

dataSets_clean.mat: 
- EECStart: earliest electrographic onset time on IEEG portal (sec)
- UEOstart: unequivical electrographic onset time on IEEG portal (sec)
- all_data: all data (cleaned?) from EEC to end across ALL channels
- data: all cleaned data from EEC to end across grid channels with toIgnore channels removed
- sozGrid: soz index on grid w/ toIgnore removed
- gridCoords: grid coorinates for each channel w/ toIgnore removed
- toIgnore: all channelnames to ignore
- cleanGridNames: names of grids with toIgnore removed


patients.mat: (not included, find in DataV3.2), includes ALL channel names and ALL gridNames/coords (note ignored coords are included)

