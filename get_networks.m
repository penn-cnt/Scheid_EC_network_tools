%% Get GLASSO Networks
% This script calls the getICov.R script located in EC_glasso to retrieve
% inverse covariance and partial correlation matrices across channels in
% time chunks of size Lwin. 

%load('Data/dataSets_clean.mat'); 
dataSet= dataSets_clean;
Null='';
gamma=.25
%load('Data/nullData.mat'); dataSet= nullData; 
%load('Data/subjects.mat')

%%

%Networks=struct();
for i_set=[1,2,7,9,10]
    fprintf('inds %d\n',i_set)
    %try

    data = dataSet(i_set).data;     % (channels x signal) matrix
    [N, l] = size(data);                   % number of channels
    Fs = round(dataSet(i_set).Fs);  % sample frequency
    Lwin = Fs;                            % num samples in a window
    TT= floor(l/(Lwin));                   % T: number of time windows
    ID=dataSet(i_set).ID;
    type=dataSet(i_set).type;
    block=dataSet(i_set).block;
    

    save('EC_glasso/tempSet.mat', 'data', 'TT', 'Lwin', 'gamma')

   
    % Call R script, this will take in the data set and return a network
    % and best rhos per time window
    cd EC_glasso
    system('Rscript getIcov.R')
    
    % This should return symmetricized and standardized icov (aka partial correlation)
    % networks, selected using extended bayesian information criteria (EBIC)
    load('glasso.mat')
    cd ../
   
    Networks(i_set).ID=ID;
    Networks(i_set).type=type;
    Networks(i_set).block=block;
    Networks(i_set).rhos=lambdas;
    Networks(i_set).pcm=nets; % Note, thess are actually PCMs
    Networks(i_set).icov=icovRho;
    
    Networks(i_set).config_pcm=[];
    Networks(i_set).config_icov=[];

    for t=1:TT
         pcm=nets(:,:,t);
         icov=icovRho(:,:,t);
         triu_i=find(triu(ones(size(pcm)),1)~=0);
         Networks(i_set).config_pcm(:,t)=pcm(triu_i);
         Networks(i_set).config_icov(:,t)=icov(triu_i);
    end
   
    Networks(i_set).sim=corrcoef(Networks(i_set).config_pcm).*~eye(TT);

end

% Weight similarity matrices with linear weighting 
for i_set=1:78
    N=Networks(i_set).sim;
    NN=zeros(length(N)); 
    for i =1:length(N)
        for j =1:length(N)
            NN(i,j)=N(i,j)*((1 - abs(i-j)/length(N))^2);    
        end
    end
    Networks(i_set).wSim=NN; 
end

eval([Null, 'Networks=Networks']);
save(sprintf('Data/%sNetworks.mat', Null), sprintf('%sNetworks', Null));

clear pcm triu_i t block type ID TT i_set networks N l Lwin lambdas Fs data
