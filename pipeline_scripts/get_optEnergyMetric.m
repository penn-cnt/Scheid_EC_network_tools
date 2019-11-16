%% optEnergyMetric for metric nodes permutation test

% Get EnergyMetric
load('Data/dataSets_clean')
load('Data/WPartitionsUEO')
load('Data/Networks')
load('Data/State_Metrics')

i_ict=find(strcmp({dataSets_clean.type},'ictal'));
nSets=size(Networks,2);

cols=[[75,184,166];[255,168,231]; [36,67,152];[140,42,195];[121,29,38];[242,224,43];[74,156,85];...
   [80,80,80]; [255,255,255]]/255;

%% Part 1: Prepare EnergyMetric struct with initial and final values

EnergyMetric=struct();
%load('Data/EnergyMetric.mat')

%Frequency Bands
alpha_theta=[5,15];
beta=[15,25];
low_g=[30,40];
high_g=[95,105];
broad_g=[30,105]; 

% Select bandpower calculation range
freq=broad_g;

%line length function
llfun=@(x)sum(abs(diff(x,[],2)),2);

for i_set= i_ict
    Net= Networks(i_set);
    p=Partitions(i_set);
    i_set
    config= Net.config_pcm;
    [N,~,T]= size(Net.pcm);
   
    EnergyMetric(i_set).ID= Net.ID;   EnergyMetric(i_set).type= Net.type; 
    EnergyMetric(i_set).block= Net.block; EnergyMetric(i_set).x0=zeros(N,3);
    EnergyMetric(i_set).repMats=zeros(N,N,3);
    
    d=dataSets_clean(i_set);
    Fs=round(d.Fs);
    
    %Set up bandpower function
    bandfun=@(x)bandpower(x', Fs, high_g)';
    
    %Get preictal bandpower for state xf.
    pi_d=dataSets_clean(i_set+nSets/2);
    stDiff=round(d.UEOStart-d.EECStart);
    
    pre_data= [dataSets_clean(i_set+nSets/2).data(:,1+2*(Fs*stDiff):end), d.data(:,1:Fs*stDiff)];

    xf=mean(MovingWinFeats(pre_data,Fs, 1, 1, bandfun),2);
    EnergyMetric(i_set).xf=mean(MovingWinFeats(pi_d.data,Fs, 1, 1, bandfun),2);
    
    % Get ictal data from UEO
    data= d.data(:,1+Fs*stDiff:end);
   
    
    for s=1:3
        %Get representative network for state
        [~,m_ind]= max(sum(corrcoef(config(:,p.contigStates==s))));
        A_s= Net.pcm(:,:,p.contigStates==s); 
        EnergyMetric(i_set).repMats(:,:,s)= A_s(:,:,m_ind);
        
        %Get ECoG signal in longest contig. run of state
        st_inds=p.stateRuns(1,:)==s;
        [m_run, i_max]=max(p.stateRuns(2,st_inds));
        t2= sum(p.stateRuns(2,1:find(cumsum(st_inds)==i_max,1,'first')));
        t1= t2-p.stateRuns(2,find(cumsum(st_inds)==i_max,1,'first'))+1;

        signal=data(1:end,(Fs*(t1-1)+1:Fs*t2));
         
%       Compute average bandpower/1 sec window of state
        x0_t=MovingWinFeats(signal,Fs, 1, 1, bandfun);
        EnergyMetric(i_set).x0(:,s)= mean(x0_t,2);         
    end 
    
end
disp('done preparing energy struct')

%% Part 2: Find ideal trajectory u* and optimal energy for nodes with high metric
%load('Data/EnergyMetric.mat')

% EnergyMetric Parameters to Test
%power(10, linspace(-3, 2, 10)); 
t_traj= linspace(0.006, 0.15, 10); %power(10, linspace(-3, log10(5), 10)) %log10 distribution b/w 1-10
rho= power(10, linspace(-3, 2, 10)); % log10 distribution b/w 1-30

metric='strength';   % metric to test energy values on
pcntN=.1;           % Percent nodes (if greater than .5, will take the highest value)

for i_set=i_ict
    i_set
    EnergyMetric(i_set).t_traj=t_traj;
    EnergyMetric(i_set).rho=rho;
    EnergyMetric(i_set).([metric, '_NOI'])=[];
    EnergyMetric(i_set).([metric, '_NullNodes'])=[];
    
    N_NOI=round(pcntN*size(State_metrics(i_set).(metric), 1)); 
    
    for s=1:3  % iterate through states
        
    % Get Metric nodes:
    % Node of Interest
    [~, i_srt]=sort(State_metrics(i_set).(metric)(:,s));
    i_NOI=i_srt(1:N_NOI); %get the smallest nodes N_NOI of metric
    EnergyMetric(i_set).([metric, '_NOI'])(:,s)=i_srt(1:N_NOI); 
    EnergyMetric(i_set).([metric, '_NullNodes'])(:,s)=i_srt(N_NOI+1:end); 
    
    A_s= EnergyMetric(i_set).repMats(:,:,s);
    x0= EnergyMetric(i_set).x0(:,s);
    xf= EnergyMetric(i_set).xf;
    [nodeEnergyMetric, trajErr]= deal(zeros(length(t_traj), length(rho)));
    compTime=zeros(length(t_traj), 1);
    
    for i_traj= 1:length(t_traj)
        tic
        T= t_traj(i_traj); 
        
        for i_rho=1:length(rho)
            r= rho(i_rho); 
             %Calculate energy per node
            try
                B=10e-5*ones(length(x0),1); B(i_NOI)=1; B=diag(B);
                [X_opt, U_opt, n_err] = optim_fun(A_s, T, B, x0, xf, r, eye(length(A_s)));
                nodeEnergyMetric(i_traj,i_rho)=sum(vecnorm(B*U_opt').^2);
                trajErr(i_traj,i_rho)=n_err;
             
            catch ME
                if strcmp(ME.identifier, 'MATLAB:svd:matrixWithNaNInf')
                    disp('err')
                    pause(0.2)
                    trajErr(:, i_traj,i_rho)=nan(length(x0),1);
                    nodeEnergyMetric(:,i_traj,i_rho)=nan(length(x0),1); 
                end
            end
            
        end % rho
        compTime(i_traj)=toc;
    end % t_traj
    EnergyMetric(i_set).(sprintf('%s_s%dtrajErr',metric, s))= trajErr;
    EnergyMetric(i_set).(sprintf('%s_s%dNodeEnergyMetric',metric,s))=nodeEnergyMetric;
    % Get idx of time where Err is min for each rho value
    [~,minErr]=min(trajErr); 
    EnergyMetric(i_set).(sprintf('%s_T_min',metric))(:,s)=minErr';
    
    end % end states loop
end 

save('Data/EnergyMetricStrength.mat', 'EnergyMetric', 't_traj', 'rho', 'freq')
disp('Part 2 EnergyMetric Calc done')

%% Part 3: Run permutation test using control parameters

%load('Data/EnergyMetric.mat')

percentRank = @(YourArray, TheProbes) reshape( mean( bsxfun(@le,...
    YourArray(:), TheProbes(:).') ) * 100, size(TheProbes) );

t_idx=6; %power(10, linspace(-3, log10(5), 10)) %log10 distribution b/w 1-10
rho_idx= 6; % log10 distribution b/w 1-30
metric='strength';
nperms=100; % Number of permutations to select from

for i_set=i_ict
    i_set
    % Select the values of t and rho that were determined previously
    T=EnergyMetric(i_set).t_traj(t_idx);
    r=EnergyMetric(i_set).rho(rho_idx);
    N=size(EnergyMetric(i_set).x0,1);

    % Get nullset:
    nullset=EnergyMetric(i_set).([metric, '_NullNodes']);
    n_null=size(nullset,1)
    i_perm=zeros(nperms, N-n_null);
    rng(5)
    % i_perm gets indices from 1:length(nullset).
    for i=1:nperms
        i_perm(i,:)=randperm(n_null,N-n_null);
    end
    
    EnergyMetric(i_set).nullSets=i_perm; 
    EnergyMetric(i_set).Metricconfidence=zeros(1,3);
    EnergyMetric(i_set).nodeEnergynull=zeros(nperms,3);
       
    for s=1:3  % iterate through states
    
        A_s= EnergyMetric(i_set).repMats(:,:,s);
        x0= EnergyMetric(i_set).x0(:,s);
        xf= EnergyMetric(i_set).xf;
        nullNodes=EnergyMetric(i_set).([metric, '_NullNodes'])(:,s);
        metricEnergy=EnergyMetric(i_set).(sprintf('%s_s%dNodeEnergyMetric',metric,s))(t_idx, rho_idx); 
        nodeEnergynull= zeros(1,nperms);
     
            try
            for nt=1:nperms % Get null distribution of energy
                B=10e-5*ones(length(x0),1); B(nullset(i_perm(nt,:),s))=1; B=diag(B);
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
            
    EnergyMetric(i_set).nodeEnergynull(:,s)=nodeEnergynull;
    
    % Get confidence interval
    EnergyMetric(i_set).Metricconfidence(s)=percentRank(nodeEnergynull, metricEnergy);
    
    end % end states loop
end 

save('Data/EnergyMetricStrength.mat', 'EnergyMetric', '-append')
disp('EnergyMetric Calc done')

%% Visualize the Metric Placement
% for each metric

figure(2); clf; 
figure(3); clf;

figure(2)
pcnts=[EnergyMetric.Metricconfidence];
sumSig=sum(pcnts>95)+sum(pcnts<5);
scatter(repmat([1:3],1,39)+normrnd(0,.05, [1,39*3]),  pcnts, 60, repmat(parula(39),3,1), 'filled')
title(sprintf('Energy, %d', sumSig))
xlim([.5,3.5])
ylim([-5, 105])
xticklabels({'Phase1', 'Phase 2', 'Phase 3'})
xticks([1:3])
xtickangle(20)
hline(5)
hline(95)

figure(3)
imagesc(reshape([EnergyMetric.Metricconfidence],3,39)'<=5)
yticks([1:39])
blk=cellfun(@num2str, {EnergyMetric([1:39]).block}', 'UniformOutput', false);
yticklabels(strcat({EnergyMetric([1:39]).ID}', {' '}, blk));
xticklabels({'Phase1', 'Phase 2', 'Phase 3'})
xticks([1:3])
yticks([1:39])


figure(2)
suptitle(sprintf('Preictal vs. phase, %s permutation test results', metric))


%% Part 3: Add state EnergyMetric to State_Metrics

% Define optimal T and rho (use  functions below to work this out)
t_opt=6;
r_opt=6; 

for i_set=i_ict
    State_metrics(i_set).optEnergyMetric=[];
    State_metrics(i_set+length(i_ict)).optEnergyMetric=[];
    for s=1:3
     State_metrics(i_set).optEnergyMetric(:,s)=EnergyMetric(i_set).(sprintf('s%dNodeEnergyMetric',s))(:,t_opt, r_opt);
     State_metrics(i_set+length(i_ict)).optEnergyMetric(:,s)=EnergyMetric(i_set).(sprintf('s%dNodeEnergyMetric',s))(:,t_opt, r_opt);
    end
    State_metrics(i_set).optEnergyMetricZ=(State_metrics(i_set).optEnergyMetric-mean(State_metrics(i_set).optEnergyMetric(:)))/...
        std(State_metrics(i_set).optEnergyMetric(:));
    State_metrics(i_set+length(i_ict)).optEnergyMetricZ=(State_metrics(i_set).optEnergyMetric-mean(State_metrics(i_set).optEnergyMetric(:)))/...
        std(State_metrics(i_set).optEnergyMetric(:));
end

tOpt=t_traj(t_opt);
rOpt=rho(r_opt);
save('Data/State_metrics.mat', 'State_metrics', 'tOpt', 'rOpt');

disp('done')
%% Visualization for zooming in 

figure(3); clf
i_set=38

for s=1:3

eng=EnergyMetric(i_set).(sprintf('%s_s%dNodeEnergyMetric',metric, s));
trajErr=EnergyMetric(i_set).(sprintf('%s_s%dtrajErr', metric, s));
    [~,minErr]=min(trajErr); 
    EnergyMetric(i_set).aveCtrl_T_min(:,s)=minErr';

st=1;
lim=10;
subset=trajErr(:, st:lim, st:lim);

figure(1)
[X,Y,Z] = ndgrid(1:size(subset,1), st:lim, st:lim);
pointsize = 30;
scatter3(X(:), Y(:), Z(:), pointsize, subset(:));
xlabel('nodes')
ylabel('T')
zlabel('rho')
set(gca, 'YTickLabels', round(t_traj(st:lim),2));
set(gca, 'ZTickLabels', round(rho(st:lim),2));
yticks((st:lim))
zticks((st:lim))

colorbar
colormap(autumn(64))
set(gca,'colorscale','log')

% try looking at variance b/w nodes?
figure(2)
imagesc(squeeze(mean(subset)))
set(gca,'colorscale','log')
set(gca,'YDir','normal')
xlabel('rho')
ylabel('T')
colormap(autumn(64))
set(gca, 'YTickLabels', round(t_traj(st:lim),2));
set(gca, 'XTickLabels', round(rho(st:lim),2));
colorbar

figure(3)
hold on
round(rho(st:lim),2)
% plot error by time
plot(squeeze(mean(subset)), 'color', cols(s,:))
xlabel('T')
ylabel('Error')
set(gca, 'XTickLabels', round(rho(st:lim),2));
xticks((st:lim)-st+1)
%legend(string(round(rho(st:lim),2)'))
% 
[~,idx]=min(squeeze(mean(subset)))
 pause(.3)
end

%% Find min and max of all error trajectories
% Get min, get max, choose value in the middle

engIdx=[EnergyMetric.T_min];

min(min(engIdx(1,:)))
max(max(engIdx(1,:)))

figure(2)
clf;
for i=1:10
    hold on
    histogram(t_traj(engIdx(i,:)),5)
    mean(t_traj(engIdx(i,:)))
end


%% Quantify error percentile for each patient and state

percentRank = @(YourArray, TheProbes) reshape( mean( bsxfun(@le,...
    YourArray(:), TheProbes(:).') ) * 100, size(TheProbes) );

mins=zeros(10,10,3);
maxes=zeros(10,10,3);


for t_opt= 1:10
    for r_opt= 1:10
        
    allRanks=zeros(length(i_ict),3);
    for i_set=i_ict
        for s=1:3
            % Average error across nodes
            errs=squeeze(mean(EnergyMetric(i_set).(sprintf('s%dtrajErr',s))));
           allRanks(i_set,s)=percentRank(errs(:,r_opt),errs(t_opt,r_opt));
        end
    end
    
    % Get min/max error across all siezures and states for param set
    mn=min(allRanks);
    mx=max(allRanks);
   
    
    mins(t_opt, r_opt,:)=mn;
    maxes(t_opt, r_opt,:)=mx;
   
    end  
end


figure(2)
imagesc(min(mins, [], 3))
caxis([0,100])
title('maximum percentile')
colorbar
figure(1)
clf
imagesc(max(maxes, [], 3))
caxis([0,100])
title('maximum error percentile across all siezures and phases')
xticklabels(EnergyMetric(1).t_traj)
yticklabels(strsplit(sprintf('%0.03f ',EnergyMetric(1).rho)))
colorbar


