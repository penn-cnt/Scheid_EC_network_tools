%Function: load_data.m
%Project: EC Controllability- V3

% Load data
load('~/Documents/LittLab/DATA/PLoSCB_2015_RawDataset_Ictal_Preictal/datasets_all_clean.mat', ...
    'dataSets_clean')
load('Data/included.mat');
load('~/Documents/LittLab/DATA/PLoSCB_2015_RawDataset_Ictal_Preictal/clinical_metadata.mat')

[subject(:).ChannelsCoordCln]=subject.ChannelsCoord;
[subject(:).ChannelsCln]=subject.Channels;

% Remove channels
for i_set=(1:length(included))
    if ~isempty(included(i_set).rm_chan)
        toRm=true(1,size(dataSets_clean(i_set).data,1));
        toRm(included(i_set).rm_chan)=false;
        
        dataSets_clean(i_set).data(~toRm,:)=[];
        i_sub= strcmp({subject.ID},dataSets_clean(i_set).ID);
        subject(i_sub).ChannelsCln=subject(i_sub).Channels(toRm,:);
        subject(i_sub).ChannelsCoordCln= subject(i_sub).ChannelsCoord(toRm,:);
    end
end

% Remove Unwanted Siezures
dataSets_clean= dataSets_clean(logical([included.incl_sz]));

clear toRm i_sub i_set