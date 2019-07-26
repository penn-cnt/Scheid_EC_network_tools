% Get Energy
load('Data/dataSets_clean')
load('Data/Partitions')
load('Data/Networks')
load('Data/State_Metrics')

i_ict=find(strcmp({dataSets_clean.type},'ictal'));
nSets=size(Networks,2);

cols=[[75,184,166];[255,168,231]; [36,67,152];[140,42,195];[121,29,38];[242,224,43];[74,156,85];...
   [80,80,80]; [255,255,255]]/255;

%% Part 1: Prepare Energy struct with initial and final values

%Energy=struct();
%load('Data/Energy.mat')

%Frequency Bands
alpha_theta=[5,15];
beta=[15,25];
low_g=[30,40];
high_g=[95,105];

freq=high_g;

%line length function
llfun=@(x)sum(abs(diff(x,[],2)),2);

for i_set= i_ict
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
    bandfun=@(x)bandpower(x', Fs, high_g)';
    
    %Get preictal bandpower for state xf.
    pi_d=dataSets_clean(i_set+nSets/2);
    xf=mean(MovingWinFeats(pi_d.data,Fs, 1, 1, bandfun),2);
    Energy(i_set).xf=mean(MovingWinFeats(pi_d.data,Fs, 1, 1, bandfun),2);
   
    
    for s=1:3
        %Get representative network for state
        [~,m_ind]= max(sum(corrcoef(config(:,p.states==s))));
        A_s= Net.pcm(:,:,p.states(m_ind)); 
        Energy(i_set).repMats(:,:,s)= A_s;
        
        %Get ECoG signal in longest contig. run of state
        st_inds=p.stateRuns(1,:)==s;
        [m_run, i_max]=max(p.stateRuns(2,st_inds));
        t2= sum(p.stateRuns(2,1:find(cumsum(st_inds)==i_max,1,'first')));
        t1= t2-p.stateRuns(2,find(cumsum(st_inds)==i_max,1,'first'))+1;
        signal=d.data(1:end,(Fs*(t1-1)+1:Fs*t2));
        
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
    
end
disp('done')

%% Part 2: Find ideal trajectory u* and optimal energy
%load('Data/Energy.mat')

% Energy Parameters to Test
%power(10, linspace(-3, 2, 10)); 
t_traj= linspace(0.006, 0.15, 10); %power(10, linspace(-3, log10(5), 10)) %log10 distribution b/w 1-10
rho= power(10, linspace(-3, 2, 10)); % log10 distribution b/w 1-30

for i_set=differs(1:end/2) %i_ict(29:end)
    i_set
    Energy(i_set).t_traj=t_traj;
    Energy(i_set).rho=rho;
    
    for s=1:3  % iterate through states
    
    A_s= Energy(i_set).repMats(:,:,s);
    x0= Energy(i_set).x0(:,s);
    xf= Energy(i_set).xf;
    [nodeEnergy, trajErr]= deal(zeros(length(x0),...
        length(t_traj), length(rho)));
    compTime=zeros(length(t_traj), 1);
    
    for i_traj= 1:length(t_traj)
        tic
        T= t_traj(i_traj); 
        
        for i_rho=1:length(rho)
            r= rho(i_rho); 
             %Calculate energy per node
            try
                for n=1:length(x0)
                    B=10e-5*ones(length(x0),length(x0)); B(n,n)=1;
                    [X_opt, U_opt, n_err] = optim_fun(A_s, T, B, x0, xf, r, eye(length(A_s)));
                    nodeEnergy(n,i_traj,i_rho)=sum(vecnorm(B*U_opt').^2);
                    trajErr(n, i_traj,i_rho)=n_err;
                end
                
            catch ME
                if strcmp(ME.identifier, 'MATLAB:svd:matrixWithNaNInf')
                    disp('err')
                    pause(0.2)
                    trajErr(:, i_traj,i_rho)=nan(length(x0),1);
                    nodeEnergy(:,i_traj,i_rho)=nan(length(x0),1); 
                end
            end
            
        end % rho
        compTime(i_traj)=toc;
    end % t_traj
    Energy(i_set).(sprintf('s%dtrajErr',s))= trajErr;
    Energy(i_set).(sprintf('s%dNodeEnergy',s))=nodeEnergy;
    Energy(i_set).compuTime=compTime
    % Get idx of time where Err is min
    [~,minErr]=min(squeeze(mean(trajErr))); 
    Energy(i_set).T_min(:,s)=minErr';
    
    end % end states loop
end 

save('Data/Energy.mat', 'Energy', 't_traj', 'rho', 'freq')
disp('Energy Calc done')


%% Part 3: Add state Energy to State_Metrics

% Define optimal T and rho (use  functions below to work this out)
t_opt=6;
r_opt=6; 

for i_set=i_ict
%     State_metrics(i_set).optEnergy=[];
%     State_metrics(i_set+length(i_ict)).optEnergy=[];
%     eval(sprintf('State_metrics(i_set).%sZ=stAvg((%s-mean(%s(:)))/std(%s(:)));',m{1},m{1},m{1},m{1}));
%     for s=1:3
%      State_metrics(i_set).optEnergy(:,s)=Energy(i_set).(sprintf('s%dNodeEnergy',s))(:,t_opt, r_opt);
%      State_metrics(i_set+length(i_ict)).optEnergy(:,s)=Energy(i_set).(sprintf('s%dNodeEnergy',s))(:,t_opt, r_opt);
%     end
    State_metrics(i_set).optEnergyZ=(State_metrics(i_set).optEnergy-mean(State_metrics(i_set).optEnergy(:)))/...
        std(State_metrics(i_set).optEnergy(:));
    State_metrics(i_set+length(i_ict)).optEnergyZ=(State_metrics(i_set).optEnergy-mean(State_metrics(i_set).optEnergy(:)))/...
        std(State_metrics(i_set).optEnergy(:));
end

tOpt=t_traj(t_opt);
rOpt=rho(r_opt);
save('Data/State_metrics.mat', 'State_metrics', 'tOpt', 'rOpt');

disp('done')
%% Visualization for zooming in 

figure(3); clf
i_set=38

for s=1:3

eng=Energy(i_set).(sprintf('s%dNodeEnergy',s));
trajErr=Energy(i_set).(sprintf('s%dtrajErr', s));
    [~,minErr]=min(squeeze(mean(trajErr))); 
    Energy(i_set).T_min(:,s)=minErr';

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

engIdx=[Energy.T_min];

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
            errs=squeeze(mean(Energy(i_set).(sprintf('s%dtrajErr',s))));
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
xticklabels(Energy(1).t_traj)
yticklabels(strsplit(sprintf('%0.03f ',Energy(1).rho)))
colorbar


