% Get Energy
run('initProject.m')

%% Part 1: Prepare Energy struct with initial and final values

% For each patient, this section computes an x0 initial state for the three seizure
% phases as the average signal bandpower in "freq_band" over all windows in
% each phase, resulting in an N x 1 initial state vector. 
% xf is calculated as the average band power across all preictal windows.

% x0_z and xf_z are calculated by taking the zscore of all x0s and xf. 

Energy=struct();
%load('Data/Energy.mat')

%Frequency Bands
alpha_theta=[5,15];
beta=[15,25];
low_g=[30,40];
high_g=[95,105];

freq_band= high_g;
i_ict=find(strcmp({dataSets_clean.type},'ictal'));


info=['x0: mean bandpower across all windows in phase',...
      'xf: preict. bandpower', ...
      'B: single node, 1e-7 along diag'];

for i_set= i_ict % Note, should only iterate over ictal sets. 
    Net= Networks(i_set);
    p=Partitions(i_set);
    i_set
    config= Net.config_pcm;
    [N,~,T]= size(Net.pcm);
   
    Energy(i_set).ID= Net.ID;   Energy(i_set).type= Net.type; 
    Energy(i_set).block= Net.block; Energy(i_set).x0=zeros(N,3);
    Energy(i_set).repMats=zeros(N,N,3);
    
    d=dataSets_clean(i_set);
    Fs=round(d.Fs);
    
    %Set up bandpower function
    bandfun=@(x)bandpower(x', Fs, freq_band)';
    
    %Get preictal bandpower for state xf.
    pi_d=dataSets_clean(i_set+nSets/2);
    stDiff=round(d.UEOStart-d.EECStart);
    
    pre_data= [dataSets_clean(i_set+nSets/2).data(:,1+2*(Fs*stDiff):end), d.data(:,1:Fs*stDiff)];

    xf=mean(MovingWinFeats(pre_data,Fs, 1, 1, bandfun),2);
    Energy(i_set).xf=mean(MovingWinFeats(pi_d.data,Fs, 1, 1, bandfun),2);
    
    % Get ictal data from UEO
    data= d.data(:,1+Fs*stDiff:end);
   
    
    for s=1:3
        %Get representative network for state
        [~,m_ind]= max(sum(corrcoef(config(:,p.contigStates==s))));
        A_s= Net.pcm(:,:,p.contigStates==s); 
        Energy(i_set).repMats(:,:,s)= A_s(:,:,m_ind);
        
        %Get ECoG signal in longest contig. run of state
        t1= find(p.contigStates==s, 1, 'first');
        t2= find(p.contigStates==s, 1, 'last');
        signal=data(1:end,(Fs*(t1-1)+1:Fs*t2));
            
        %Uncomment to visualize and double check
%         figure(1); clf
%         subplot(121); hold on
%         plot(d.data'+[1:size(d.data,1)]*1000);
%         stem(d.Fs*(t1-1)+1, 50000, 'lineWidth', 2, 'color', 'red'); 
%         stem(d.Fs*t2, 50000, 'lineWidth', 2, 'color', 'red');
%         axis tight
%         subplot(122)
%         plot(signal'+[1:size(signal,1)]*1000);
%         axis tight; pause
%         figure(2)
%         imagesc(p.states) 
%         colormap(gca, cols([5,2,6],:));
%         set(gca, 'YTick', [], 'fontsize', 18)
   
        
%       Compute average bandpower/1 sec window of state
        x0_t=MovingWinFeats(signal,Fs, 1, 1, bandfun);
        Energy(i_set).x0(:,s)= mean(x0_t,2);
        
        % Uncomment to View x0 and xf over time windows:
%         figure(1)
%         imagesc(reshape(xf, 8, 8))
%         caxis([0,.5])
%         colorbar
% 
%         figure(2)
%         for i=1:size(x0,2)
%             imagesc(reshape(x0_t(:,i), 8, 8))
%             caxis([0,.5])
%             colorbar
%             pause(.1)
%         end
%         
    end
    
    % Calculate z-scored x0 and xf
    EZ= zscore([Energy(i_set).x0, Energy(i_set).xf], [],'all');
    Energy(i_set).x0_z=EZ(:,1:3); Energy(i_set).xf_z=EZ(:,4); 
    
end
disp('done')

%% Part 2: Find ideal trajectory u* and optimal energy

% Calculate control energy to drive state x0 to xf from each individual node
% with underlying connectivity as a represetative phase connectivity network A_s. 
% A_S gets stabilized as a continuous system

%load('Data/Energy.mat')

i_ict=find(strcmp({dataSets_clean.type},'ictal'));
relax= 10^-7; % value along non-driven diagonals of "relaxed" B matrix
%scale= 10^7; % value along non-driven diagonals of "input penalty" R matrix

% Energy Parameters to Test
% On first pass, run with large parameter set to find best parameters
tr= power(10, linspace(-2, log10(6), 10));
t_traj= round([tr(1:7),(1:5)],3); %log10 distribution b/w 0-1, then linear 1-6
rho= round([power(10, linspace(-2, 2, 10)), 110, 150],3); % log10 distribution b/w 1-30, rounded to 3 places

for i_set= i_ict(1:end)
    tic
    i_set
    Energy(i_set).t_traj=t_traj;
    Energy(i_set).rho=rho;
    xf= Energy(i_set).xf;
    
    for s=1:3  % iterate through states
    
    x0= Energy(i_set).x0(:,s);
    N=length(x0);
    A= Energy(i_set).repMats(:,:,s);  
    A_s= A./(1+svds(A,1))-eye(N);       % normalize representative matrix
    
    [nodeEnergy, trajErr]= deal(zeros(length(x0),...
        length(t_traj), length(rho)));

    xOpt=zeros(length(x0),length(t_traj), length(rho));

%     trajErr= Energy(i_set).(sprintf('s%dtrajErr',s));
%     nodeEnergy= Energy(i_set).(sprintf('s%dNodeEnergy',s));
%     xOpt = Energy(i_set).(sprintf('s%Xopt_dist',s));
    
    for i_traj= 1:length(t_traj) % try different time horizons
        i_traj
        T= t_traj(i_traj);
        for i_rho= 1:length(rho) % try different energy/distance tradeoffs
            r= rho(i_rho);
            try
                for n=1:length(x0) %Calculate energy per node
                    
                    B=relax*ones(length(x0),1); B(n)=1; B=diag(B); 
                    nodeEnergy(n,i_traj,i_rho)=sum(vecnorm(B*U_opt').^2);
                    [X_opt, U_opt, n_err] = optim_fun(A_s, T, B, x0, xf, r, eye(length(A_s)));
                    trajErr(n, i_traj,i_rho)=n_err;
                    xOpt(n, i_traj, i_rho)=mean((X_opt(end,1:N)'-xf).^2);
                    
                    % Uncoment to view output
%                     figure(12); imagesc(U_opt'); title('U'); 
%                     figure(22); subplot(1,2,1); imagesc(B); title('B');
%                                 subplot(1,2,2); imagesc(R); title('R');
%                     figure(32);
%                     subplot(1,6,[2:5]); imagesc([X_opt(:,1:N)'])% imagesc(X_opt(:,2:N)'); title('X');
%                     subplot(1,6,[2:5]); imagesc([X_opt(1,1:N)',X_opt(end,1:N)'])
%                     cl=caxis;
%                     subplot(1,6,1); imagesc(log(x0)); title('x0'); caxis(cl)
%                     subplot(1,6,6); imagesc(log(xf)); title('xf');  caxis(cl)

                end
                
            catch ME
                if strcmp(ME.identifier, 'MATLAB:svd:matrixWithNaNInf')
                    disp('err')
                    pause
                    trajErr(:, i_traj,i_rho)=nan(length(x0),1);
                    nodeEnergy(:,i_traj,i_rho)=nan(length(x0),1); 
                else; throw(ME)
                end              
            end % end try    
        end % rho       
    end % t_traj
    
    Energy(i_set).(sprintf('s%dtrajErr',s))= trajErr;
    Energy(i_set).(sprintf('s%dNodeEnergy',s))= nodeEnergy;
    Energy(i_set).(sprintf('s%Xopt_dist',s))= xOpt;
    
    
    end % end states loop
    toc
end 

save(fullfile(datafold, '/Energy_V4_1.mat'), 'Energy', 't_traj', 'rho', 'freq_band', 'relax', 'info', 'scale')
disp('Energy Calc done')

%% Part 2.5: Quantify error percentile for each patient and state

i_ict(ismember(i_ict,rm))=[];

percentRank = @(YourArray, TheProbes) reshape( mean( bsxfun(@le,...
    YourArray(:), TheProbes(:).') ) * 100, size(TheProbes) );

err_stats=zeros(length(t_traj),length(rho), length(i_ict)*3);
dist_err_stats=zeros(length(t_traj),length(rho), length(i_ict)*3);


allRanks=zeros(length(t_traj),length(rho),length(i_ict));
for i_set=1:length(i_ict)
    phaseRanks=zeros(length(t_traj),length(rho),3);
    for s=1:3
        % Average error across nodes
        errs=squeeze(mean(Energy(i_ict(i_set)).(sprintf('s%dtrajErr',s))));
        dist_errs=squeeze(mean(Energy(i_ict(i_set)).(sprintf('s%Xopt_dist',s))));
        
        % Error percentile out of entire "parameter square" for as single phase
       phaseRanks(:,:,s)=percentRank(errs,errs); %errs
       err_stats(:,:,(i_set-1)*3+s)=errs;
       
       dist_phaseRanks(:,:,s)=dist_errs;  %percentRank(errs,errs);
       dist_err_stats(:,:,(i_set-1)*3+s)=dist_errs;

    end
    allRanks(:,:,i_set)=max(phaseRanks,[],3); 
    dist_allRanks(:,:,i_set)=max(dist_phaseRanks,[],3);   % get max error at each cell across phases
end
    
max_percentiles=max(allRanks,[],3);
dist_max_percentiles=max(dist_allRanks,[],3);

figure(1)
clf
imagesc(max_percentiles')%percentRank(max_percentiles', max_percentiles'))
set(gca,'YDir','normal')
caxis([0,100])
title('maximum error percentile across all siezures and phases')
xticks([1:length(t_traj)])
xticklabels(strsplit(sprintf('%0.02f ',t_traj')))
yticklabels(strsplit(sprintf('%0.02f ',rho)))
yticks([1:length(rho)])
ylabel('\rho'); xlabel('\tau'); colorbar
xtickangle(45);

figure(2)
clf
imagesc(percentRank(dist_max_percentiles', dist_max_percentiles'))
set(gca,'YDir','normal')
caxis([0,100])
title('maximum distance error percentile across all siezures and phases')
xticks([1:length(t_traj)])
xticklabels(strsplit(sprintf('%0.02f ',t_traj')))
yticklabels(strsplit(sprintf('%0.02f ',rho)))
yticks([1:length(rho)])
ylabel('\rho'); xlabel('\tau'); colorbar
xtickangle(45);

[pct_err, i_err]=min(max_percentiles, [], 'all', 'linear')
[t_opt, r_opt]=ind2sub(size(max_percentiles), i_err)
mn_err= [mean(err_stats(t_opt, r_opt, :)),std(err_stats(t_opt, r_opt, :))] 

pctile_err= prctile(err_stats(t_opt, r_opt, :),[25 50 75]);

%save(fullfile(datafold, '/Energy.mat'), 'err_stats', '-append')
% saveas(gcf,'FigsV3.3/energy/energy_max_err_pcnt.png')
% saveas(gcf,'FigsV3.3/energy/energy_max_err_pcnt.fig')
% saveas(gcf,'FigsV3.3/energy/energy_max_distErr_pcnt.png')
% saveas(gcf,'FigsV3.3/energy/energy_max_distErr_pcnt.fig')

%% Part 3: Add state Energy to State_Metrics

% Define optimal T and rho (use  functions below to work this out)
t_opt=7;
r_opt=7; 
i_ict=find(strcmp({dataSets_clean.type},'ictal'));

for i_set=i_ict
    State_metrics(i_set).optEnergy=[];
    State_metrics(i_set+length(i_ict)).optEnergy=[];
    for s=1:3
     State_metrics(i_set).optEnergy(:,s)=Energy(i_set).(sprintf('s%dNodeEnergy',s))(:,t_opt, r_opt);
     State_metrics(i_set+length(i_ict)).optEnergy(:,s)=Energy(i_set).(sprintf('s%dNodeEnergy',s))(:,t_opt, r_opt);
    end
    State_metrics(i_set).optEnergyZ= zscore(State_metrics(i_set).optEnergy, [], 'all');
    State_metrics(i_set+length(i_ict)).optEnergyZ= zscore(State_metrics(i_set).optEnergy, [], 'all');
end

tOpt=t_traj(t_opt);
rOpt=rho(r_opt);
%save('Data/State_metrics.mat', 'State_metrics', 'tOpt', 'rOpt', '-append');

disp('done')
%% Visualization for zooming in 

for i_set=i_ict
figure(3); clf
traj_lim=[1:10];
rho_lim=[1:10];
i_set

for s=1:3

eng=Energy(i_set).(sprintf('s%dNodeEnergy',s));
trajErr=Energy(i_set).(sprintf('s%dtrajErr', s));
    [~,minErr]=min(squeeze(mean(trajErr))); 
    Energy(i_set).T_min(:,s)=minErr'

subset=trajErr(:, traj_lim, rho_lim);

figure(1)
[X,Y,Z] = ndgrid(1:size(subset,1), traj_lim, rho_lim);
pointsize = 30;
scatter3(X(:), Y(:), Z(:), pointsize, subset(:));
xlabel('nodes')
ylabel('T')
zlabel('rho')
set(gca, 'YTickLabels', round(t_traj(traj_lim),2));
set(gca, 'ZTickLabels', round(rho(rho_lim),2));
yticks((traj_lim))
zticks((rho_lim))

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
yticks(traj_lim); xticks(rho_lim);
set(gca, 'YTickLabels', round(t_traj(traj_lim),3));
set(gca, 'XTickLabels', round(rho(rho_lim),2));
colorbar

figure(3)
% plot error by time
semilogy(squeeze(mean(subset)), 'color', cols(s,:))
hold on
xlabel('T')
ylabel('Error')
set(gca, 'XTickLabels', round(t_traj(traj_lim),2));
xticks((rho_lim)-rho_lim(1)+1)
% 
[~,idx]=min(squeeze(mean(subset)), [], 'all', 'linear')
 %pause
end
pause
end
%% Find min and max of all error trajectories
% Get min, get max, choose value in the middle

i_ict(ismember(i_ict,rm))=[];
engIdx=[Energy(i_ict).T_min];

min(min(engIdx(1,:)))
max(max(engIdx(1,:)))

figure(2)
clf;
for i=1:10
    hold on
    histogram(t_traj(engIdx(i,:)),5)
    mean(t_traj(engIdx(i,:)))
end





