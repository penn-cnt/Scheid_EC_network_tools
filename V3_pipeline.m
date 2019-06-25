%% EC Controllability- V3 Pipeline

Null='';

% 1) Load Data 
% run load_clean_data.m
load('Data/dataSets_clean.mat')
dataSets_clean= clean_EC_data;
load('Data/subjects.mat')

% 1.5) Generate Phase Randomized Null Data
% run genNullModel.m
load('Data/nullData.mat')
%%
% 2) Get Networks
% run get_networks.m
load(sprintf('Data/%sNetworks.mat', Null))

% 3) Find Partitions
% run dynamic_findSimComms.m
% run get_states.m
load(sprintf('Data/%sPartitions.mat', Null))

% 4) Compute Metrics
% run get_metrics.m
load(sprintf('Data/%sMetric_Matrices.mat', Null))
load(sprintf('Data/%sState_Metrics.mat', Null))

% 5) Compute Control Energy
% run get_opt_energy.mat
load('Data/Energy.mat')

% 6) Robustness
% assessModalThresh.m

%% Get coordinates of SOZ

for i_set=1:nSets
    
end

figure(4); clf; hold on;
plot(squeeze(sum(sum(icov>0,1)))./85^2)
plot(squeeze(sum(sum(icov<0,1)))./85^2)

figure(5); clf; hold on;
plot(squeeze(sum(sum(pcm>0,1)))./85^2)
plot(squeeze(sum(sum(pcm<0,1)))./85^2)




