% Example EC Controllability Pipeline

addpath('helper_functions')

% load example grid recording data with artifact and strip channels removed
raw = load('HUP64_sz1_ictal.mat', 'd');

% process data
seizure = cleanRawSeizure(raw.d.data, raw.d.Fs);

% Get regularized partial correlation matrices across windowed data
timeseries_data = seizure; %(64 channels x 59255 samples)
windowSize = 500;
Networks = getNets('./', timeseries_data, windowSize, 0.5, [0, 0.1, .5]);

% Partition Network similarity matrix into desired number of communities
init_gamma = (0.8:0.05:1.2); % Set of initial spatial resolution parameters 
nTarget = 3;               % Target number of communities 
Qiter = 100;               % iterations of Louvain algorithm
Partitions = findSimComms(Networks.wSim_0_5, init_gamma, Qiter, nTarget);