% Initialize Project
cd '/Users/bscheid/Documents/LittLab/PROJECTS/p01_EC_controllability'

datafold='v4/DataV4.1';
addpath(genpath('~/Documents/CODE/'))
addpath('Code/helper_functions')
addpath('Code/pipeline_scripts')

if ~exist('dataSets_clean', 'var')
    load([datafold,'/dataSets_clean.mat'])
end
load(sprintf('v3/%s/subjects.mat', 'DataV3.2'));
load([datafold,'/Partitions.mat'])
load([datafold,'/Networks.mat'])
load([datafold,'/State_metrics.mat'])
load([datafold,'/Metric_matrices.mat'])
% load([datafold,'/Energy.mat'])
% load([datafold,'/EnergySOZ.mat'])

% SETTINGS



i_ict=find(strcmp({dataSets_clean.type},'ictal'));
i_preict=find(strcmp({dataSets_clean.type},'preictal'));
nSets=length(Partitions);
nIct=nSets/2;

rm=[31 36 33 35 29]; % Indices of seizures to remove in group analysis

cols=[[227,187,187]; [190,8,4]; [138,4,4];[140,42,195];[75,184,166];[242,224,43];[74,156,85];...
   [80,80,80]; [255,255,255]]/255;

wSim='wSim_0_01'; % default sim network. 

datafold='v4/DataV4.2'; % change datafold for saving data. 


