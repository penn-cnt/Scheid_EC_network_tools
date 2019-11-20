% Generate Null model

% Create phase randomized models of each network, randomizing by time
% window
%load('Data/dataSets_clean.mat')
dataSets_clean= clean_EC_data;
load('Data/subjects.mat')

nullData=struct();

rng(20);
for i_set= 1:length(dataSets_clean)
    i_set
    nullData(i_set).ID=dataSets_clean(i_set).ID;
    nullData(i_set).type=dataSets_clean(i_set).type;
    nullData(i_set).block=dataSets_clean(i_set).block;
    Fs=dataSets_clean(i_set).Fs;
    nullData(i_set).Fs=Fs;
    randomized=zeros(size(dataSets_clean(i_set).data));
    
    for i=1:round(Fs):size(randomized,2)
        randomized(:,i:round(Fs)+i-1)= phaseRandomize(dataSets_clean(i_set).data(:,i:i+round(Fs)-1));
    end
    
    nullData(i_set).data=randomized;
end

save('Data/nullData.mat', 'nullData')