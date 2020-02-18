%% Compute Metrics (similarity, avgCtrl, etc) (get metric_matrices.mat)

metric_matrices= struct(); 
State_metrics=struct();


thresh=0.15;     % Threshold for transient/persistent mode selection
dt=1;        %

% only compute for beta= 0.01
pinds=find(strcmp({Partitions.beta}, '0_01'))

for i_set=1:length(pinds)
    fprintf('inds %d\n',pinds(i_set))
    clear skewness kurtosis % clear to avoid func. ambiguity
    
    p= Partitions(pinds(i_set));
    % match ID
    netInd= strcmp(p.ID, {Network.ID}).*strcmp(p.type, {Network.type}); 
    Net=Network(find(netInd)); 
    
    % Initialize metric structs
    metric_matrices(i_set).ID= p.ID; State_metrics(i_set).ID=p.ID;    
    metric_matrices(i_set).type= p.type; State_metrics(i_set).type=p.type; 
    metric_matrices(i_set).block=p.block; State_metrics(i_set).block=p.block;
    metric_matrices(i_set).beta=p.beta; State_metrics(i_set).beta=p.beta;
        
    config=  Net.config_pcm;
    pcm=     Net.pcm;  
    [N, ~, T]=       size(pcm);
   
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
       
   end
   
    % Define function to compute state averages, compute averages from
    % longest run of each state type though. 
    st= Partitions(pinds(i_set)).contigStates;
    stAvg= @(metric)cell2mat(arrayfun(@(x)median(metric(:,st==x),2),...
        [1:3], 'UniformOutput', false)); 
    
    % Populate metric and state average structs
    for m=metrics
        eval(sprintf('metric_matrices(i_set).%s=%s;',m{1},m{1}));
        eval(sprintf('State_metrics(i_set).%s=stAvg(%s);',m{1},m{1}));
        % Get Zscore
        eval(sprintf('State_metrics(i_set).%sZ=stAvg((%s-mean(%s(:)))/std(%s(:)));',m{1},m{1},m{1},m{1}));
    end
   
end


Metric_matrices= metric_matrices;
% save(sprintf('Data/%sMetric_matrices.mat', Null), sprintf('%sMetric_matrices', Null))
% save(sprintf('Data/%sState_metrics.mat', Null), sprintf('%sState_metrics', Null))

disp('done')


% %% Correction for states
%    metrics={'globalCtrl', 'aveCtrl', 'modalCtrl', 'pModalCtrl', 'tModalCtrl',...
%        'strength', 'strengthPos', 'strengthNeg', ...        % Network metrics %
%        'skewness', 'kurtosis', 'degree', 'clustering3'}; 
% 
% State_metrics=struct();
% 
% for i_set=1:length(pinds)
%     st= Partitions(pinds(i_set)).contigStates;
%     stAvg= @(metric)cell2mat(arrayfun(@(x)median(metric(:,st==x),2),...
%     [1:3], 'UniformOutput', false)); 
% 
%     Net= Metric_matrices(i_set);
%     
%     % Initialize metric structs
%     State_metrics(i_set).ID=Net.ID;    
%     State_metrics(i_set).type=Net.type; 
%     State_metrics(i_set).block=Net.block;
%     
%     for m=metrics
%         met= Metric_matrices(i_set).(m{1});
%         State_metrics(i_set).(m{1})= stAvg(met);
%         State_metrics(i_set).([m{1},'Z'])= stAvg((met-mean(met(:)))/std(met(:))); 
%     end
% 
% end
% 
% disp('Done')
%     

