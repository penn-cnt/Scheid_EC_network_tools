%% Compute Metrics (similarity, avgCtrl, etc) (get metric_matrices.mat)
addpath(genpath('supportFcns'))
Null='' % set string to 'Null' if computing on null model, else set empty

load(sprintf('Data/%sNetworks', Null))
load(sprintf('Data/%sPartitions', Null))
% 
% eval(['Networks=', Null, 'Networks;']);
% eval(['Partitions=', Null, 'Partitions;']);

nSets=length(Partitions);
metric_matrices= struct(); 
state_metrics=struct();

thresh=0.15;     % Threshold for transient/persistent mode selection
dt=1;        %

for i_set=1:nSets
    fprintf('inds %d\n',i_set)
    clear skewness kurtosis % clear to avoid func. ambiguity
    
    Net= Networks(i_set);
    
    % Initialize metric structs
    metric_matrices(i_set).ID= Net.ID; State_metrics(i_set).ID=Net.ID;    
    metric_matrices(i_set).type= Net.type; State_metrics(i_set).type=Net.type; 
    metric_matrices(i_set).block= Net.block; State_metrics(i_set).block=Net.block;
        
    config=  Net.config_pcm;
    pcm=     Net.pcm;
    [N, ~, T]= size(pcm);
   
   metrics={'globalCtrl', 'aveCtrl', 'modalCtrl', 'pModalCtrl', 'tModalCtrl',...
       'strength', 'strengthPos', 'strengthNeg', ...        % Network metrics %
       'skewness', 'kurtosis', 'degree', 'clustering3'}; %, ... % connection distribution %
       %'spreadEig', 'skeweig'};                             % eigenstructure distribution %
    
   % Initialize average, modal metrics
   [globalCtrl, aveCtrl, modalCtrl, pModalCtrl, tModalCtrl,...
       strength, strengthPos, strengthNeg, ...
       degree, clustering3]= deal(zeros(N, T));
   
   skewness=skewness(config);
   kurtosis=kurtosis(config);
    
   for t=1:T
       if mod(t,10)==0; fprintf('%s', '.'); end
       
       %%%%%%%%%%%%%%%%%%%%%%%%
       %%% Control  Metrics %%%
       %%%%%%%%%%%%%%%%%%%%%%%%
       
       % Normalize A
        A= pcm(:,:,t);
        normA= A./(1+svds(A,1))-eye(N);
       
        % Global Control
       for n=1:N
           B=eye(N).*1e-6; B(n,n)=1;                % set control points
           sys= ss(normA, B, [], []);
           gramian= gram(sys, 'c');
           globalCtrl(n,t)=min(eig(gramian)); 
       end
       
       % Avg Control
       sys=ss(normA,eye(N), [], []); 
       gramian=gram(sys, 'c');
       aveCtrl(:,t)= diag(gramian);
       
       %Modal Control
       modalCtrl(:,t)= modal_control(A, dt); %Note: norm happens in function
       % Note, find best threshold using "assessModalThresh.m", norm
       % happens in func:
      [pModalCtrl(:,t), tModalCtrl(:,t)]= modal_persist_trans(A, .15, dt);
    
       
       %%%%%%%%%%%%%%%%%%%%%%%%
       %%% Network Metrics %%%
       %%%%%%%%%%%%%%%%%%%%%%%%
       degree(:,t)=degrees_und(pcm(:,:,t));
       
       % Calculating strength as sum of pos and abs(neg)
       [nPos,nNeg,vpos,vneg]=strengths_und_sign(pcm(:,:,t));
       strength(:,t)=(nPos'-nNeg');
       strengthPos(:,t)=nPos;
       strengthNeg(:,t)=nNeg;
       
       % Calculating signed clustering using Perugini's signed extension
       [C_pos,C_neg,Ctot_pos,Ctot_neg] = clustering_coef_wu_sign(pcm(:,:,t),3);
       clustering3(:,t)=C_pos; 
       
       
       %%%%%%%%%%%%%%%%%%%%%%%%
       %%%  Eigenstructure  %%%
       %%%%%%%%%%%%%%%%%%%%%%%%
       %spreadEig(:,t)=eig(normA);
       
       % Calculating strength as sum of pos and abs(neg)
       [nPos,nNeg,vpos,vneg]=strengths_und_sign(pcm(:,:,t));
       strength(:,t)=(nPos'-nNeg');
       strengthPos(:,t)=nPos;
       strengthNeg(:,t)=nNeg;
       
       % Calculating signed clustering using Perugini's signed extension
       [C_pos,C_neg,Ctot_pos,Ctot_neg] = clustering_coef_wu_sign(pcm(:,:,t),3);
       clustering3(:,t)=C_pos; 


   end
    
    % Define function to compute state averages, compute averages from
    % longest run of each state type though. 
    st= Partitions(i_set).contigStates;
    stAvg= @(metric)cell2mat(arrayfun(@(x)mean(metric(:,st==x),2),...
        unique(st), 'UniformOutput', false)); 
    
    % Populate metric and state average structs
    for m=metrics
        eval(sprintf('metric_matrices(i_set).%s=%s;',m{1},m{1}));
        eval(sprintf('state_metrics(i_set).%s=stAvg(%s);',m{1},m{1}));
        % Get Zscore
        eval(sprintf('state_metrics(i_set).%sZ=stAvg((%s-mean(%s(:)))/std(%s(:)));',m{1},m{1},m{1},m{1}));
    end
   
end


eval([Null,'Metric_matrices= metric_matrices']);
eval([Null,'State_metrics= state_metrics']);
save(sprintf('Data/%sMetric_matrices.mat', Null), sprintf('%sMetric_matrices', Null))
save(sprintf('Data/%sState_metrics.mat', Null), sprintf('%sState_metrics', Null))

disp('done')
