%% Find Similarity Communities
addpath(genpath('supportFcns'))
% Reorder at end... 

Null=''; 

betas={'0_01'};

load('Data/Networks.mat')
fprintf('Glasso gamma %0.02f \n', gamma)
    
Partitions=struct(); 

for s=1:length(betas)

maxIter= 5;           % max # of iterations
nTarget= 3;           % target Number of communities

gamma_init=(0.8:.05:1.05); % initial resolution parameter range
Qiter=100;            % number of mod. max iterations

n=length(Networks)*(s-1);

for i_set = 1:length(Networks)
    tic
    fprintf('inds %d\n',i_set)
    ctr=1; 
    
    %gamma=Partitions(i_set).gamma;
    
    gamma=gamma_init;
    % Set up Partitions.mat fields
    Partitions(n+i_set).ID=Networks(i_set).ID;    
    Partitions(n+i_set).type=Networks(i_set).type; 
    Partitions(n+i_set).block=Networks(i_set).block;
    Partitions(n+i_set).gamma=gamma;
    Partitions(n+i_set).simBeta=betas{s}; 
    
    gamma_inds=1:length(gamma);  % Initial gamma Inds
    if strcmp(betas{s}, '')
        sim= Networks(i_set).sim;
    else
        sim= Networks(i_set).(sprintf('wSim_%s',betas{s}));
    end
    
    N=size(sim,1);

   %Initialize PartitionRho Struct
   Partitions(n+i_set).quantileQ= zeros(length(gamma),5);
   Partitions(n+i_set).quantileCommNum= zeros(length(gamma),5);
   Partitions(n+i_set).quantileCommSz= zeros(length(gamma),5);
   Partitions(n+i_set).medianQcomms= zeros(N, length(gamma));
   Partitions(n+i_set).consensusQcomms= zeros(N, length(gamma)); 
    
    while ctr < maxIter  % At most complete 5 iterations
        fprintf('ctr %d\n', ctr)
    
        %if ctr~=1

            % Find communities for gamma at gamma_inds
            for i_gamma = gamma_inds
                fprintf('i_gamma %d\n',i_gamma)
                Qs=zeros(1,Qiter);
                Ss=zeros(N,Qiter);
                S_uniq=zeros(1,Qiter);
                S_avgSz=zeros(1,Qiter);

                A= sim;
                k = full(sum(A));

                for i_iter = 1:Qiter
                    % Run Gen Louvain
                    B = full(A - gamma(i_gamma)*(k'*k)/sum(k));
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
                Partitions(n+i_set).medianQcomms(:,i_gamma)= Ss(:,i_median);                               % save median comm
                Partitions(n+i_set).consensusQcomms(:,i_gamma)= consensus_similarity(Ss')';
                Partitions(n+i_set).quantileQ(i_gamma,:)=quantile(Qs,[0.025 0.25 0.50 0.75 0.975]);
                Partitions(n+i_set).quantileCommNum(i_gamma,:)=quantile(S_uniq,[0.025 0.25 0.50 0.75 0.975]);
                Partitions(n+i_set).quantileCommSz(i_gamma,:)=quantile(S_avgSz,[0.025 0.25 0.50 0.75 0.975]);
            end

        %end
        
        % Dynamically find gamma for target comm size
        p=Partitions(i_set);
        cNum=arrayfun(@(x)length(unique(p.consensusQcomms(:,x))),(1:length(p.gamma)))
        i_low= find(cNum<=nTarget, 1, 'last');
        i_high= find(cNum>=nTarget, 1, 'first');
        
        olgamma=gamma;

        if cNum(i_low) == cNum(i_high)
            ctr= maxIter;
            break; 
        end  % nTarget communities exists
        
        if isempty(i_low); gamma=(gamma(1)-.2:.05:gamma(1)-.05); % search lower gamma
        elseif isempty(i_high); gamma=(gamma(end)+.05:.05:gamma(end)+.2); % search higher gamma
        else; gamma=linspace(gamma(i_low), gamma(i_high), 5); %Search gamma b/w values
        end
        
        % Shift values over and round off small errors
        gamma= sort(unique([olgamma, round(gamma, 6)]));
        Partitions(i_set).gamma=gamma;
        
        gamma_inds=find(~ismember(gamma, olgamma)); % inds for new compuations
        old_inds=find(ismember(gamma, olgamma)); %inds to move old comps to
        
        %Reassign old calculations
        MQ=zeros(N,length(gamma)); MQ(:,old_inds)=Partitions(n+i_set).medianQcomms;
        CQ=zeros(N,length(gamma)); CQ(:,old_inds)=Partitions(n+i_set).consensusQcomms;
        qQ=zeros(length(gamma),5); qQ(old_inds,:)=Partitions(n+i_set).quantileQ;
        CN=zeros(length(gamma),5); CN(old_inds,:)=Partitions(n+i_set).quantileCommNum;
        CZ=zeros(length(gamma),5); CZ(old_inds,:)=Partitions(n+i_set).quantileCommSz;

        Partitions(n+i_set).medianQcomms =MQ;
        Partitions(n+i_set).consensusQcomms = CQ;
        Partitions(n+i_set).quantileQ = qQ;
        Partitions(n+i_set).quantileCommNum= CN;
        Partitions(n+i_set).quantileCommSz= CZ;
        
        ctr=ctr+1;
    
    end
    
    % TODO: Put communities in order here instead of later
    
    % Find Q at inflection point
    p=Partitions(n+i_set);
    dt=max(min(diff(p.gamma)), 1e-4);
    gi=(min(p.gamma):dt:max(p.gamma));
    yi=interp1(p.gamma,p.quantileQ(:,3),gi);
    [~,i_max]=max(diff(diff(yi)));
    Partitions(i_set).modInflection=yi(i_max+1);

    infInd=find(p.quantileQ(:,3)<=Partitions(i_set).modInflection, 1, 'first');
    Partitions(n+i_set).modInfInd=infInd;
    Partitions(n+i_set).nTargInd=find(cNum>=nTarget,1,'first');   
    toc
    
end %end i_set
end %end s

eval([Null, 'Partitions= Partitions;']);
save('Data/Partitions.mat', 'Partitions', 'gamma');
save(sprintf('Data/Robustness/Partitions_%s.mat',gammas{g}), sprintf('%sPartitions', Null))


%try; run ___.m; system('mail -s "Borel Finished" "brittany.h.scheid@gmail.com" <<< "Test Job Done"'); catch ME; system(sprintf('mail -s "Borel Error" "brittany.h.scheid@gmail.com" <<< "%s, %s"', ME.identifier, ME.message)); rethrow(ME);  end;   