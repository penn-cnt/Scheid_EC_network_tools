function [Networks]=getNets(path, data, Lwin, gamma, beta) 
% This script calls the getICov.R script located in EC_glasso to retrieve
% inverse covariance and partial correlation matrices across channels in

% data- Nxl data matrix, N channels, l-time samples
% Lwin- time window to create network over,
% gamma- EBIC parameter, 0.1-less regularization, 0.5-more regularization
% beta- tune the intensity of distance weighting for community detection

    Networks=struct();

    [~, l] = size(data);                   % number of channels
    TT= floor(l/(Lwin));                   % T: number of time windows
 
     save(fullfile(path, 'EC_glasso/tempSet.mat'), 'data', 'TT', 'Lwin', 'gamma')

    % Call R script, this will take in the data set and return a network
    % and best rhos per time window
    cd(fullfile(path, 'EC_glasso'))
    system('Rscript getIcov.R')
    cd ../
    
    % This should return symmetricized and standardized icov (aka partial correlation)
    % networks, selected using extended bayesian information criteria (EBIC)
    glasso=load(fullfile(path, 'EC_glasso/glasso.mat'));
   
    Networks.rhos=glasso.lambdas;
    Networks.pcm=glasso.nets; % Note, thess are actually PCMs
    Networks.icov=glasso.icovRho;
    
    Networks.config_pcm=[];
    Networks.config_icov=[];

    for t=1:TT
         pcm=glasso.nets(:,:,t);
         icov=glasso.icovRho(:,:,t);
         triu_i=find(triu(ones(size(pcm)),1)~=0);
         Networks.config_pcm(:,t)=pcm(triu_i);
         Networks.config_icov(:,t)=icov(triu_i);
    end
   
    Networks.sim=corrcoef(Networks.config_pcm).*~eye(TT);
    
    sim=Networks.sim;
    Wsim=zeros(length(sim)); 
    for i =1:TT
        for j =1:TT
            Wsim(i,j)=sim(i,j)*((1 - beta*abs(i-j)/length(sim))^2);    
        end
    end
    Networks.wSim=Wsim; 

end

