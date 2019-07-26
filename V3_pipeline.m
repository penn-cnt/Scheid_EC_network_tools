%% EC Controllability- V3 Pipeline

cd('/Users/bscheid/Documents/LittLab/PROJECTS/p01_EC_controllability/v3/Code')
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
%run get_networks.m
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
avgDist=cell(39,1);
% get average distance from a given node to all SOZ nodes
% load('/Users/bscheid/Documents/LittLab/DATA/virtualResection_DataSet/dataSets_raw.mat')
for i_set=1:39

    if isempty(dataSets_clean(i_set).channels_soz)
        continue
    end
    % Get coordinates of soz nodes for subject
    pt_idx=strcmp(dataSets_clean(i_set).ID, {subjects.ID});
    r_id=strcmp(dataSets_clean(i_set).ID, {dataSets_raw.ID});
    r_blk=ismember([dataSets_raw.block],dataSets_clean(i_set).block);
    inds=find(r_id.*r_blk);
    
    nIgnore=~ismember(subjects(pt_idx).grids, dataSets_raw(inds(1)).toIgnore); 
    gridchans=subjects(pt_idx).grids(nIgnore); 
    sozIdx=ismember(subjects(pt_idx).grids,dataSets_clean(i_set).channels_soz);
    dataSets_clean(i_set).sozGrid=ismember(gridchans,dataSets_clean(i_set).channels_soz);
    dataSets_clean(i_set+39).sozGrid=ismember(gridchans,dataSets_clean(i_set).channels_soz);
     
    sozCoords=subjects(pt_idx).gridCoords(logical(sozIdx.*nIgnore),:);
    dataSets_clean(i_set).gridCoords=subjects(pt_idx).gridCoords(logical(nIgnore),:);
    dataSets_clean(i_set+39).gridCoords=subjects(pt_idx).gridCoords(logical(nIgnore),:);
    
    if size(sozCoords,1)==1
        avgDist{i_set}=pdist2(sozCoords, subjects(pt_idx).gridCoords(nIgnore,:))';
    else
        avgDist{i_set}=mean(pdist2(sozCoords, subjects(pt_idx).gridCoords(nIgnore,:)))';
    end
end

%save('Data/avgDist.mat', 'avgDist')

