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

%% Avg Distance from Node to SOZ nodes

% get distance 
for i_set=1:nIct

    if isempty(dataSets_clean(i_set).channels_soz)
        continue
    end
    % Get coordinates of soz nodes for subject
    pt_idx=strcmp(dataSets_clean(i_set).ID, {subjects.ID});
    sozIdx=contains(subjects(pt_idx).grids,dataSets_clean(i_set).channels_soz); 
    sozCoords=subjects(pt_idx).gridCoords(sozIdx,:);

    avgDist{i_set}=mean(pdist2(sozCoords, subjects(pt_idx).gridCoords))';
end


