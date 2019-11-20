%% optEnergySOZ SOZ permutation test

% Get EnergySOZ
% load('Data/dataSets_clean')
% load('Data/WPartitionsUEO')
% load('Data/Networks')
load('Data/State_metrics')

i_ict=find(strcmp({dataSets_clean.type},'ictal'));
nSets=size(Networks,2);
i_soz=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 25 26 28];

%colors for plotting
cols=[[75,184,166];[255,168,231]; [36,67,152];[140,42,195];[121,29,38];[242,224,43];[74,156,85];...
   [80,80,80]; [255,255,255]]/255;

%Frequency Bands
alpha_theta=[5,15];
beta=[15,25];
low_g=[30,40];
high_g=[95,105];
broad_g=[30,105]; 


%% Part 1: Prepare EnergySOZ struct with initial and final values

for nWind=[Inf, 1, 5, 10]
EnergySOZ=struct();
%load('Data/EnergySOZ.mat')

% SETTINGS Select bandpower calculation range, nWind, info
freq=broad_g;
const=.5;
filename=sprintf('Data/EnergySOZ%d.mat', nWind);
info=sprintf(['x0: mean bandpower across first %d windows in phase, ',...
      'xf: zero state, ', ...
      'B: SOZ nodes=1, else %d'],nWind, const);

%line length function
llfun=@(x)sum(abs(diff(x,[],2)),2);

for i_set= i_ict
    Net= Networks(i_set);
    p=Partitions(i_set);
    i_set
    config= Net.config_pcm;
    [N,~,T]= size(Net.pcm);
   
    EnergySOZ(i_set).ID= Net.ID;   EnergySOZ(i_set).type= Net.type; 
    EnergySOZ(i_set).block= Net.block; EnergySOZ(i_set).x0=zeros(N,3);
    EnergySOZ(i_set).repMats=zeros(N,N,3);
    
    d=dataSets_clean(i_set);
    Fs=round(d.Fs);
    
    %Set up bandpower function
    bandfun=@(x)bandpower(x', Fs, high_g)';
    
    %Get preictal bandpower for state xf.
    %pi_d=dataSets_clean(i_set+nSets/2);
       
    stDiff=round(d.UEOStart-d.EECStart);
    pre_data= [dataSets_clean(i_set+nSets/2).data(:,1+2*(Fs*stDiff):end), d.data(:,1:Fs*stDiff)];
    EnergySOZ(i_set).xf=mean(MovingWinFeats(pre_data,Fs, 1, 1, bandfun),2);
    
    % Get ictal data from UEO
    data= d.data(:,1+Fs*stDiff:end);
   
    
    for s=1:3
        %Get representative network for state
        [~,m_ind]= max(sum(corrcoef(config(:,p.contigStates==s))));
        A_s= Net.pcm(:,:,p.contigStates==s); 
        EnergySOZ(i_set).repMats(:,:,s)= A_s(:,:,m_ind);
        
        %Get ECoG signal in longest contig. run of state
        st_inds=p.stateRuns(1,:)==s;
        [m_run, i_max]=max(p.stateRuns(2,st_inds));
        t2= sum(p.stateRuns(2,1:find(cumsum(st_inds)==i_max,1,'first')));
        t1= t2-p.stateRuns(2,find(cumsum(st_inds)==i_max,1,'first'))+1;

        signal=data(1:end,(Fs*(t1-1)+1:Fs*t2));
         
%       Compute average bandpower/1 sec window of state
        x0_t=MovingWinFeats(signal,Fs, 1, 1, bandfun);
        EnergySOZ(i_set).x0(:,s)=mean(x0_t(:,1:min(nWind,end)),2); %get first nWind of bandpower       
    end 
    
end
disp('done preparing energy struct')

% Part 2: Find ideal trajectory u* and optimal energy for SOZs
%load('Data/EnergySOZ.mat')

% EnergySOZ Parameters to Test
% power(10, linspace(-3, 2, 10)); 
t_traj= linspace(0.006, 0.15, 10); %power(10, linspace(-3, log10(5), 10)) %log10 distribution b/w 1-10
rho= power(10, linspace(-3, 2, 10)); % log10 distribution b/w 1-30

for i_set=i_ict
    i_set
    EnergySOZ(i_set).t_traj=t_traj;
    EnergySOZ(i_set).rho=rho;
    
    % Get SOZ nodes:
    if isempty(dataSets_clean(i_set).sozGrid); continue; end
    soz=dataSets_clean(i_set).sozGrid; 
    if sum(soz)==0; continue; end
    
    for s=1:3  % iterate through states
    
    A_s= EnergySOZ(i_set).repMats(:,:,s);
    x0= EnergySOZ(i_set).x0(:,s);
    xf= EnergySOZ(i_set).xf;
    [nodeEnergySOZ, trajErr]= deal(zeros(length(t_traj), length(rho)));
    compTime=zeros(length(t_traj), 1);
    
    for i_traj= 1:length(t_traj)
        tic
        T= t_traj(i_traj); 
        
        for i_rho=1:length(rho)
            r= rho(i_rho); 
             %Calculate energy per node
            try
                B=const*ones(length(x0),1); B(soz)=1; B=diag(B);
                [X_opt, U_opt, n_err] = optim_fun(A_s, T, B, x0, xf, r, eye(length(A_s)));
                nodeEnergySOZ(i_traj,i_rho)=sum(vecnorm(B*U_opt').^2);
                trajErr(i_traj,i_rho)=n_err;
            catch ME
                if strcmp(ME.identifier, 'MATLAB:svd:matrixWithNaNInf')
                    disp('err')
                    pause(0.2)
                    trajErr(:, i_traj,i_rho)=nan(length(x0),1);
                    nodeEnergySOZ(:,i_traj,i_rho)=nan(length(x0),1); 
                end
            end
            
        end % rho
        compTime(i_traj)=toc;
    end % t_traj
    EnergySOZ(i_set).(sprintf('s%dtrajErr',s))= trajErr;
    EnergySOZ(i_set).(sprintf('s%dNodeEnergySOZ',s))=nodeEnergySOZ;
    % Get idx of time where Err is min for each rho value
    [~,minErr]=min(trajErr); 
    EnergySOZ(i_set).T_min(:,s)=minErr';
    
    end % end states loop
end 

time=datetime(now,'ConvertFrom','datenum');
save(filename, 'EnergySOZ', 't_traj', 'rho', 'freq', 'info', 'time')
disp('Part 2 EnergySOZ Calc done')

% Part 3: Run SOZ permutation test using control parameters

%load('Data/EnergySOZ.mat')

percentRank = @(YourArray, TheProbes) reshape( mean( bsxfun(@le,...
    YourArray(:), TheProbes(:).') ) * 100, size(TheProbes) );

t_idx=6; %power(10, linspace(-3, log10(5), 10)) %log10 distribution b/w 1-10
rho_idx= 6; % log10 distribution b/w 1-30

nperms=500; % Number of permutations to select from

for i_set=i_ict
    fprintf('Finding null for i_set %d...\n', i_set)
    % Select the values of t and rho that were determined previously
    T=EnergySOZ(i_set).t_traj(t_idx);
    r=EnergySOZ(i_set).rho(rho_idx);

    % Get SOZ nodes:
    if isempty(dataSets_clean(i_set).sozGrid); continue; end
    soz=dataSets_clean(i_set).sozGrid; 
    if sum(soz)==0; continue; end
    
    % Get nullset
    nsoz=sum(soz); 
    nullset=find(~soz); 
    i_perm=zeros(nperms, nsoz);
    rng(5)
    for i=1:nperms
        i_perm(i,:)=nullset(randperm(sum(~soz),nsoz));
    end
    
    EnergySOZ(i_set).nullSets=i_perm; 
    EnergySOZ(i_set).SOZconfidence=zeros(1,3);
    EnergySOZ(i_set).nodeEnergynull=zeros(nperms,3);
       
    for s=1:3  % iterate through states
    
        A_s= EnergySOZ(i_set).repMats(:,:,s);
        x0= EnergySOZ(i_set).x0(:,s);
        xf= EnergySOZ(i_set).xf;
        sozEnergy=EnergySOZ(i_set).(sprintf('s%dNodeEnergySOZ',s))(t_idx, rho_idx); 
        nodeEnergynull= zeros(1,nperms);
     
            try
            for nt=1:nperms % Get null distribution of energy
                B=const*ones(length(x0),1); B(i_perm(nt,:))=1; B=diag(B);
                [X_opt, U_opt, n_err] = optim_fun(A_s, T, B, x0, xf, r, eye(length(A_s)));
                nodeEnergynull(nt)=sum(vecnorm(B*U_opt').^2);             
            end

        catch ME
            if strcmp(ME.identifier, 'MATLAB:svd:matrixWithNaNInf')
                disp('err')
                pause(0.2)
                nodeEnergynull(nt)=nan(length(x0),1); 
            end
            end %End try
            
    EnergySOZ(i_set).nodeEnergynull(:,s)=nodeEnergynull;
    
    % Get confidence interval
    EnergySOZ(i_set).SOZconfidence(s)=percentRank(nodeEnergynull, sozEnergy);
    
    end % end states loop
end 

time=datetime(now,'ConvertFrom','datenum');
save(filename, 'EnergySOZ', 'time', '-append')
%save(sprintf('Data/EnergySOZ_zero%d.mat',nWind), 'EnergySOZ', 't_traj', 'rho', 'freq')
disp('Part 3 EnergySOZ Calc done')


end

%% Visualize the SOZ signifcance level
% for each metric
nsoz=18; 

cb=97.5; %confidence bound (percent)

for j=[Inf]
    %load(sprintf('Data/EnergySOZ%d.mat', j));
    
figure(min(90,j)*17)
pcnts=[EnergySOZ.SOZconfidence];
sumSig=sum(pcnts>cb)+sum(pcnts<(100-cb));
%scatter(repmat([1:3],1,18)+normrnd(0,.05, [1,54]),  pcnts, 60, repmat(copper(18),3,1), 'filled')
scatter(repmat([1:3],1,nsoz)+normrnd(0,.05, [1,3*nsoz]),  pcnts, 25,'blue')
%title(sprintf('Energy SOZ Permutation Test (all windows) significance',j))
title(sprintf('Energy SOZ Permutation Test (%d windows) significance',j))
xlim([.5,3.5]); ylim([-5, 105])
xticklabels({'Phase1', 'Phase 2', 'Phase 3'})
xticks([1:3])
xtickangle(20)
hline(100-cb)
hline(cb)

figure(min(90,j)*17+1)
imagesc(reshape([EnergySOZ.SOZconfidence],3,nsoz)'<=(100-cb))
yticks([1:nsoz])
blk=cellfun(@num2str, {Metric_matrices(i_soz).block}', 'UniformOutput', false);
yticklabels(strcat({Metric_matrices(i_soz).ID}', {' '}, blk));
xticklabels({'Phase1', 'Phase 2', 'Phase 3'})
xticks([1:3])
title(sprintf('Energy SOZ Permutation Test (%d windows) significance',j))


end


%% Part 3: Add state EnergySOZ to State_Metrics

% Define optimal T and rho (use  functions below to work this out)
t_opt=6;
r_opt=6; 

for i_set=i_ict
    State_metrics(i_set).optEnergySOZ=[];
    State_metrics(i_set+length(i_ict)).optEnergySOZ=[];
    for s=1:3
     State_metrics(i_set).optEnergySOZ(:,s)=EnergySOZ(i_set).(sprintf('s%dNodeEnergySOZ',s))(:,t_opt, r_opt);
     State_metrics(i_set+length(i_ict)).optEnergySOZ(:,s)=EnergySOZ(i_set).(sprintf('s%dNodeEnergySOZ',s))(:,t_opt, r_opt);
    end
    State_metrics(i_set).optEnergySOZZ=(State_metrics(i_set).optEnergySOZ-mean(State_metrics(i_set).optEnergySOZ(:)))/...
        std(State_metrics(i_set).optEnergySOZ(:));
    State_metrics(i_set+length(i_ict)).optEnergySOZZ=(State_metrics(i_set).optEnergySOZ-mean(State_metrics(i_set).optEnergySOZ(:)))/...
        std(State_metrics(i_set).optEnergySOZ(:));
end

tOpt=t_traj(t_opt);
rOpt=rho(r_opt);
save('Data/State_metrics.mat', 'State_metrics', 'tOpt', 'rOpt');

disp('done')