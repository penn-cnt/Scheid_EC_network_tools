function [p]= findSimComms(sim, gamma_init, Qiter, nTarget)
%
% Find a specified number of Similarity Communities using Louvaine-like modularity maximization
% and a Neuman-Girvan null model by tuning the structural resolution
% parameter, gamma. Algorithm overview: 
% 1) For each gamma value, the modularity maximization algorithm is run Qiter times
% 2) Consensus partitions over the Qiter partitions are computed for the initial set of gamma values
% 3) The range of gamma values is either extended or interpolated until nTarget communities are discovered. 
%
% INPUTS: 
% sim: symmetric adjacency matrix for finding communities
% gamma_init: vector of initial gammas to try
% Qiter: number of iterations of the Louvain modularity maximization
% algorithm.
% nTarget: number of communities to target
%
% OUTPUTS: 
%  Partitions struct, with the following fields: 
%   gamma: vector of gamma values used in modularity calculations
%   quantileQ: matrix containing 2.5, 25, 50, 75, and 9.75 percentile values of the modularity quality factor. 
%   quantileCommNumm: matrix containing percentiles of number of communities found
%   quantileCommSz: matrix contatining percentiles of average community size
%   medianQcomms: community partition for median Q value 
%   consensusQcomms: consensus partition over Qiter modularity maximization iterations
%   modInfInd: index of the inflection point of the modularity curve 
%   nTargInd: index where gamma produced a consensus parition of nTarget communities
    
maxIter= 5; % max # of iterations
ctr=1; 
    
gamma= gamma_init;

p.gamma= gamma;  
gamma_inds=1:length(gamma);  % Initial gamma Inds

N=size(sim,1);

%Initialize PartitionRho Struct
p.quantileQ= zeros(length(gamma),5);
p.quantileCommNum= zeros(length(gamma),5);
p.quantileCommSz= zeros(length(gamma),5);
p.medianQcomms= zeros(N, length(gamma));
p.consensusQcomms= zeros(N, length(gamma)); 

while ctr < maxIter  % At most complete 5 iterations
    fprintf('ctr %d\n', ctr)

        % Find communities for gamma at gamma_inds
        for i_gamma = gamma_inds
            fprintf('i_gamma %d\n',i_gamma)
            Qs=zeros(1,Qiter);
            Ss=zeros(N,Qiter);
            S_uniq=zeros(1,Qiter);
            S_avgSz=zeros(1,Qiter);

            k = full(sum(sim));

            for i_iter = 1:Qiter
                % Run Gen Louvain
                B = full(sim - gamma(i_gamma)*(k'*k)/sum(k));
                [S,Q]=genlouvain(B, 10000, 0);
                % add Partitions and modularity 
                Ss(:,i_iter)=S; 
                Qs(i_iter)=Q/sum(k);
                freq=histogram(S, length(unique(S))); 
                S_uniq(i_iter)=length(unique(S)); %num unique
                S_avgSz(i_iter)=mean(freq.Values); %avg comm size
            end

            % Add  modularity values and partition at max modularity
            i_median=find(Qs<=median(Qs), 1, 'last');
            p.medianQcomms(:,i_gamma)= Ss(:,i_median);                               % save median comm
            p.consensusQcomms(:,i_gamma)= consensus_similarity(Ss')';
            p.quantileQ(i_gamma,:)=quantile(Qs,[0.025 0.25 0.50 0.75 0.975]);
            p.quantileCommNum(i_gamma,:)=quantile(S_uniq,[0.025 0.25 0.50 0.75 0.975]);
            p.quantileCommSz(i_gamma,:)=quantile(S_avgSz,[0.025 0.25 0.50 0.75 0.975]);
        end

    % Dynamically find gamma for target comm size
    cNum=arrayfun(@(x)length(unique(p.consensusQcomms(:,x))),(1:length(p.gamma)));
    i_low= find(cNum<=nTarget, 1, 'last');
    i_high= find(cNum>=nTarget, 1, 'first');

    olgamma=gamma;

    if cNum(i_low) == cNum(i_high)
        ctr= maxIter;
        break; 
    end  % nTarget communities exists

    if isempty(i_low); gamma=(gamma(1)-.2:.05:gamma(1)-.05);            % search lower gamma
    elseif isempty(i_high); gamma=(gamma(end)+.05:.05:gamma(end)+.2);   % search higher gamma
    else; gamma=linspace(gamma(i_low), gamma(i_high), 5);               % Search gamma b/w values
    end

    % Assimilate new gamma values and round off small decimals
    gamma= sort(unique([olgamma, round(gamma, 6)]));
    p.gamma=gamma;

    gamma_inds=find(~ismember(gamma, olgamma));     % inds for new computations
    old_inds=find(ismember(gamma, olgamma));     %inds to move old comps to

    %Reassign old calculations
    MQ=zeros(N,length(gamma)); MQ(:,old_inds)=p.medianQcomms;
    CQ=zeros(N,length(gamma)); CQ(:,old_inds)=p.consensusQcomms;
    qQ=zeros(length(gamma),5); qQ(old_inds,:)=p.quantileQ;
    CN=zeros(length(gamma),5); CN(old_inds,:)=p.quantileCommNum;
    CZ=zeros(length(gamma),5); CZ(old_inds,:)=p.quantileCommSz;

    p.medianQcomms =MQ;
    p.consensusQcomms = CQ;
    p.quantileQ = qQ;
    p.quantileCommNum= CN;
    p.quantileCommSz= CZ;

    ctr=ctr+1;

end

% TODO: Put communities in order of newly explored gammas here instead of in the loop

% Find Q at inflection point
dt=max(min(diff(p.gamma)), 1e-4);
gi=(min(p.gamma):dt:max(p.gamma));
yi=interp1(p.gamma,p.quantileQ(:,3),gi);
[~,i_max]=max(diff(diff(yi)));
p.modInflection=yi(i_max+1);

infInd=find(p.quantileQ(:,3)<=p.modInflection, 1, 'first');
p.modInfInd=infInd;
p.nTargInd=find(cNum>=nTarget,1,'first');   

    
end