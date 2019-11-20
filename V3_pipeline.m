%% EC Controllability- V3 Pipeline
dataFold='DataV3.2'
cd('/Users/bscheid/Documents/LittLab/PROJECTS/p01_EC_controllability/v3/Code')
Null='';
addpath('pipeline_scripts')
addpath('helper_functions')

% 1) Load Data 
% run load_clean_data.m
load(fullfile(dataFold,'dataSets_clean.mat'))
load(fullfile(dataFold,'subjects.mat'))

% 1.5) Generate Phase Randomized Null Data
% run genNullModel.m
%load(fullfile(dataFold,'nullData.mat'))
%%
% 2) Get Networks
%run get_networks.m
load(sprintf('%s/%sNetworks.mat', dataFold,Null))

% 3) Find Partitions
% run dynamic_findSimComms.m
% run get_states.m
for b=1:length(betas)
    ii = length(Network)*(b-1);
    for i_set=2:length(Network)
        d= Network(i_set);
        sim= Network(i_set).(sprintf('wSim_%s',betaStr{b}));
        n= findSimComms(sim, (0.8:.05:1.05), 100, 3);
        n.ID=d.ID; n.type=d.type; n.block=d.block; n.Fs=d.Fs;
        n.beta= betaStr{b}; 
        
        %Get assignments to three (or more) states
        nUnique= arrayfun(@(x)length(unique(n.consensusQcomms(:,x))),(1:length(n.gamma)));
        st= p.consensusQcomms(:,find(nUnique>=nStates,1,'first'))';

        % Check if median comms contains correct num of communities
        if isempty(st)
            nUnique=arrayfun(@(x)length(unique(n.medianQcomms(:,x))),(1:length(n.gamma)));
            st=n.medianQcomms(:,find(nUnique>=nStates,1,'first'))';
        end

        s=assignStates(st);
        Partitions(ii+i_set)= mergeStructs(n, s); 
    end
end






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

