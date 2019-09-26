%% optEnergySOZ SOZ permutation test



% Get EnergySOZ
load('Data/dataSets_clean')
load('Data/WPartitionsUEO')
load('Data/Networks')
load('Data/State_Metrics')

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

for nWind=[Inf]
EnergySOZ=struct();
%load('Data/EnergySOZ.mat')

% SETTINGS Select bandpower calculation range, nWind, info
freq=broad_g;
sigma=1.6332; % input control energy spread
filename=sprintf('Data/EnergySOZ%d_spread.mat', nWind);
info=sprintf(['x0: mean bandpower across first %d windows in phase, ',...
      'xf: preict state, ', ...
      'B: spread with sigma'],nWind, sigma);

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
    B= getSpreadControl(dataSets_clean(i_set).gridCoords, soz);
    %B=const*ones(length(x0),1); B(soz)=1; B=diag(B);
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
                B= getSpreadControl(dataSets_clean(i_set).gridCoords, i_perm(nt,:));
                %B=const*ones(length(x0),1); B(i_perm(nt,:))=1; B=diag(B);
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
hold on
imagesc(reshape([EnergySOZ.SOZconfidence],3,nsoz)');
plot([1:3],(reshape([EnergySOZ.SOZconfidence],3,nsoz)'<=(100-cb)).*repmat([1:18]',1,3), '*', 'color', 'red')
yticks([1:nsoz])
blk=cellfun(@num2str, {EnergySOZ(i_soz).block}', 'UniformOutput', false);
yticklabels(strcat({EnergySOZ(i_soz).ID}', {' '}, blk));
xticklabels({'Phase1', 'Phase 2', 'Phase 3'})
xticks([1:3])
title(sprintf('Energy SOZ Permutation Test (%d windows) significance',j))
axis tight
ylim([.5, 18.5])
set(gca, 'YDir', 'reverse')


end


%% Part 3: Add state EnergySOZ to State_Metrics

% Define optimal T and rho (use  functions below to work this out)
t_opt=6;
r_opt=6; 

for i_set=i_soz
    State_metrics(i_set).optEnergySOZ=[];
    State_metrics(i_set+length(i_ict)).optEnergySOZ=[];
    for s=1:3
     State_metrics(i_set).optEnergySOZ(:,s)=EnergySOZ(i_set).(sprintf('s%dNodeEnergySOZ',s))(t_opt, r_opt);
     State_metrics(i_set+length(i_ict)).optEnergySOZ(:,s)=EnergySOZ(i_set).(sprintf('s%dNodeEnergySOZ',s))(t_opt, r_opt);
    end
    State_metrics(i_set).optEnergySOZZ=(State_metrics(i_set).optEnergySOZ-mean(State_metrics(i_set).optEnergySOZ(:)))/...
        std(State_metrics(i_set).optEnergySOZ(:));
    State_metrics(i_set+length(i_ict)).optEnergySOZZ=(State_metrics(i_set).optEnergySOZ-mean(State_metrics(i_set).optEnergySOZ(:)))/...
        std(State_metrics(i_set).optEnergySOZ(:));
end

tOpt=t_traj(t_opt);
rOpt=rho(r_opt);
time=datetime(now,'ConvertFrom','datenum');
save('Data/State_metrics.mat', 'State_metrics', 'tOpt', 'rOpt', 'time');

disp('done')
%% Spatial visual of node energies

clf;
cmap=colormap; 
for i_set=3 %a
    figure(1)
    clf; hold on
    if isempty(dataSets_clean(i_set).sozGrid)
        continue
    end
    
    e=EnergySOZ(i_set);
    N=size(e.x0,1);
    
    % Calculate average energy values
    avgEng=nan(N,3);
    for node=1:N
       [I,~] = ind2sub(size(e.nullSets),find(e.nullSets==node))
        if ~isempty(I)
            %avgEng(node,:)=e.nodeEnergynull(I(randi(length(I),1)),:);
            avgEng(node,:)=mean(e.nodeEnergynull(I,:));
        end
    end
    
    soz=dataSets_clean(i_set).sozGrid; 
    avgEng(soz,1)=e.s1NodeEnergySOZ(t_opt, r_opt);
    avgEng(soz,2)=e.s2NodeEnergySOZ(t_opt, r_opt); 
    avgEng(soz,3)=e.s3NodeEnergySOZ(t_opt, r_opt); 
    
        st=avgEng
        c_met=1+reshape(round(63*rescale(st(:))), size(st));
        clims=[min(st(:)), max(st(:))];
        X=dataSets_clean(i_set).gridCoords(:,1);
        Y=dataSets_clean(i_set).gridCoords(:,2);
        for s=3:-1:1
            subplot(1,3,s)
            colormap parula
            hold on
            scatter(X(dataSets_clean(i_set).sozGrid), ...
                 Y(dataSets_clean(i_set).sozGrid),120, 'r', 'filled')
            scatter(X,Y, 80, st(:,s), 'filled')           
            xticklabels([]); yticklabels([])
            xlim([.5,8.5]); ylim([.5,8.5]);
            xlabel(sprintf('phase %d', s))
            %caxis(clims)
            %set(gca, 'colorscale', 'log')
            %ylabel(metrics{m})
            colorbar
            %if s==3; colorbar; end
      
        
    end
    %suptitle(sprintf('Spatial Metrics for %s, block %d', dataSets_clean(i_set).ID, dataSets_clean(i_set).block))
   % pause
end

%% Plot trajectories 
clear energy
for i_set=3 %i_ict
    i_set
    % Select the values of t and rho that were determined previously
    T=EnergySOZ(i_set).t_traj(t_idx);
    r=EnergySOZ(i_set).rho(rho_idx);

    for s=1:3
        A_s= EnergySOZ(i_set).repMats(:,:,s);
        x0= EnergySOZ(i_set).x0(:,s);
        xf= EnergySOZ(i_set).xf;
        B=10e-5*ones(length(x0),1); B(soz)=1; B=diag(B);
        %T=1;
        [X_opt, U_opt, n_err] = optim_fun(A_s, T, B, x0, xf, r, eye(length(A_s)));
        nodeEnergySOZ(i_traj,i_rho)=sum(vecnorm(B*U_opt').^2);
        energy(i_set).X_opt(:,:,s)=X_opt(:,end/2:end)';
        energy(i_set).U_opt(:,:,s)=U_opt';
        energy(i_set).n_err(:,:,s)=n_err;
    end     
end
%% Get distance from final location of x opt. 
i_set=3;
plot(mean(energy(i_set).U_opt(:,:,3)))
xtraj=zeros(size(energy(i_set).X_opt(:,:,3),1),3);
for j=1:size(energy(i_set).X_opt(:,:,3),2)
    xtraj(j,1)=norm(energy(i_set).X_opt(:,j,1)-energy(i_set).X_opt(:,end,1));
    xtraj(j,2)=norm(energy(i_set).X_opt(:,j,2)-energy(i_set).X_opt(:,end,2));
    xtraj(j,3)=norm(energy(i_set).X_opt(:,j,3)-energy(i_set).X_opt(:,end,3));   
end

clf
plot([1:size(energy(i_set).U_opt,2)],xtraj)

%% Get optimal control energy over time

i_set=3;
plot(mean(energy(i_set).U_opt(:,:,3)))
xtraj=zeros(size(energy(i_set).X_opt,1),3);
for j=1:size(energy(i_set).X_opt,2)
    xtraj(j,1)=norm(B*energy(i_set).U_opt(:,j,1));
    xtraj(j,2)=norm(B*energy(i_set).U_opt(:,j,2));
    xtraj(j,3)=norm(B*energy(i_set).U_opt(:,j,3));
end

clf


%% Plot Results
clf; hold on;
plot([1:size(energy(i_set).U_opt,2)]/1000, xtraj(:,1), 'color', cols(1,:))
plot([1:size(energy(i_set).U_opt,2)]/1000, xtraj(:,2), 'color', cols(2,:))
plot([1:size(energy(i_set).U_opt,2)]/1000, xtraj(:,3), 'color', cols(3,:))
xlabel('T')
ylabel('||X(t)-X_T|| (a.u.)o')
legend('Phase 1', 'Phase 2', 'Phase 3')


%% Look at correlation of mean control metric of nodes in regions at top of distribution.
colormap('winter')
figure(304)
clf;
figure(7)
clf; 
metrics={'aveCtrl', 'pModalCtrl', 'tModalCtrl', 'strength', 'strengthNeg'};
[SOZMetp,SOZMetcorr]=deal(zeros(18, 3)); 

for m=1:length(metrics)
    figure(304)
    for i_set=i_soz
        met=State_metrics(i_set).(metrics{m});
        for s=1:3
            subplot(1,3,s)
            hold on
            % Get node index of lowest controll energy sets, sorted low to high
            energyNull=EnergySOZ(i_set).nodeEnergynull(:,s);
            [nullE, i_eNull]=sort(energyNull);
            nodeSets=EnergySOZ(i_set).nullSets(i_eNull, :); 

            % Get metric values for nodes
            setMet=mean(reshape(met(nodeSets,1), size(nodeSets)),2); 
            [rr, pp]=corr(nullE, setMet); 
%             mdl = polyfit(nullE, setMet, 1);
%             f=polyval(mdl, nullE); 
%             scatter(nullE, setMet)
%             plot(nullE,f); 
%             title(sprintf('r^2=%0.3f, p=%0.3f', rr^2, pp))
            SOZMetcorr(i_set,s)=rr; 
            SOZMetp(i_set,s)=pp;
        end
        %suptitle(sprintf('Energy and %s correlation, HUP %s, %d', metrics{m}, EnergySOZ(i_set).ID, EnergySOZ(i_set).block))
        %pause
    end
    
    figure(7)
    subplot(1,length(metrics), m)
    hold on;
    imagesc(SOZMetcorr(i_soz,:))
    plot([1:3], (SOZMetp(i_soz,:)<=0.05).*repmat([1:length(i_soz)]',1, 3), '*', 'color', 'black');
    title( metrics{m})
    axis tight
    blk=cellfun(@num2str, {EnergySOZ(i_soz).block}', 'UniformOutput', false);
    yticklabels(strcat({EnergySOZ(i_soz).ID}', {' '}, blk));
    yticks([1:18])
    set(gca,'TickLabelInterpreter','none')
    xticks([1:3])
    xticklabels({'Phase 1', 'Phase 2', 'Phase 3'})
    set(gca,'Ydir', 'reverse')
    ylim([.5,18.5])
end

%% Visualization for zooming in 

figure(3); clf
i_set=38

for s=1:3

eng=EnergySOZ(i_set).(sprintf('s%dNodeEnergySOZ',s));
trajErr=EnergySOZ(i_set).(sprintf('s%dtrajErr', s));
    [~,minErr]=min(squeeze(mean(trajErr))); 
    EnergySOZ(i_set).T_min(:,s)=minErr';

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

engIdx=[EnergySOZ.T_min];

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
        
    allRanks=zeros(length(a),3);
    for i_set=i_soz
        for s=1:3
            % Average error across nodes
            errs=EnergySOZ(i_set).(sprintf('s%dtrajErr',s));
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
xticklabels(EnergySOZ(1).t_traj)
yticklabels(strsplit(sprintf('%0.03f ',EnergySOZ(1).rho)))
colorbar


