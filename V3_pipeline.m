%% EC Controllability- V3 Pipeline

Null='';

% 1) Load Data 
% run load_clean_data.m
load('Data/dataSets_clean.mat')
load('Data/subjects.mat')

% 1.5) Generate Phase Randomized Null Data
% run genNullModel.m
load('Data/nullData.mat')
%%
% 2) Get Networks
% run get_networks.m
load(sprintf('Data/%sNetworks.mat', Null))

% 3) Find Partitions
% run dynamic_findSimComms.m
% run get_states.m
load(sprintf('Data/%sPartitions.mat', Null))

% 4) Compute Metrics
% run get_metrics.m
load(sprintf('Data/%sMetric_Matrices.mat', Null))
load(sprintf('Data/%sState_Metrics.mat', Null))

% 5) Compute Control Energy
%load('Data/energy.mat')

% 6) Robustness
% assessModalThresh.m

%% Find best energy values by changing all three parameters

energy=struct();

%Prep for energy analysis. Get x0, xT states and rep. EC network

%line length function
llfun=@(x)sum(abs(diff(x,[],2)),2);

for i_set= i_ict
    
    Net= Networks(i_set);
    p=Partitions(i_set);
    
    config= Net.config_pcm;
    [N,~,T]= size(Net.pcm);
   
    energy(i_set).ID= Net.ID;   
    energy(i_set).type= Net.type; 
    energy(i_set).block= Net.block; 
    energy(i_set).x0=zeros(N,3);
    energy(i_set).repMats=zeros(N,N,3);
    
%     'Get soz nodes'
%     subj=subject(strcmp({subject.ID},partitions(i_set).ID));
%     soz=ismember(string(subj.Channels(1:end-1)), string(dataSets_clean(i_set).channels_soz));
    
    pi_d=dataSets_clean(i_set+nSets/2);
    energy(i_set).xf=mean(MovingWinFeats(pi_d.data,pi_d.Fs, 1, 1, llfun),2);
    %subplot(2,2,1)
    %imagesc(xf);
    
%     figure(1)
%     imagesc(p.states) 
%     colormap(gca, cols([5,2,6],:));
%     set(gca, 'YTick', [], 'fontsize', 18)
    
    for s=1:3
        % Get representative network for state
        [~,m_ind]= max(sum(corrcoef(config(:,p.states==s))));
        A_s= Net.pcm(:,:,p.states(m_ind)); 
        energy(i_set).repMats(:,:,s)= A_s;
        
        % Get ECoG signal in longest contig. run of state
        st_inds=p.stateRuns(1,:)==s;
        [m_run, i_max]=max(p.stateRuns(2,st_inds));
        t2= sum(p.stateRuns(2,1:find(cumsum(st_inds)==i_max,1,'first')));
        t1= t2-p.stateRuns(2,find(cumsum(st_inds)==i_max,1,'first'))+1;
        d=dataSets_clean(i_set);
        signal=d.data(1:end,(d.Fs*(t1-1)+1:d.Fs*t2));
        
%         figure(2)
%         clf
%         subplot(121)
%         hold on
%         plot(d.data'+[1:size(d.data,1)]*1000);
%         stem(d.Fs*(t1-1)+1, 50000, 'lineWidth', 2, 'color', 'red'); 
%         stem(d.Fs*t2, 50000, 'lineWidth', 2, 'color', 'red');
%         axis tight
%         subplot(122)
%         plot(signal'+[1:size(signal,1)]*1000);
%         axis tight
%         pause
        
        %Compute average linelength/1 sec window of state
        x0=mean(MovingWinFeats(signal,d.Fs, 1, 1, llfun),2);
        energy(i_set).x0(:,s)= x0;
    end 
    
end
disp('done')


%%
% Test Energy Parameters 
t_traj= power(10, (0:6)*.2-.5) %log10 distribution b/w 1-10
rho= power(10, (0:6)*.2-.5) % log10 distribution b/w 1-30

s=2;

for i_set=1:i_ict
    
    A_s= energy(i_set).repMats(:,:,s);
    x0= energy(i_set).x0(:,s);
    xf= energy(i_set).xf;
    [nodeEnergy, trajErr]= deal(zeros(length(x0),...
        length(t_traj), length(rho)));
    
    for i_traj= 1:length(t_traj)
        tic
        T= t_traj(i_traj); 
        
        for i_rho=1:length(rho)
            i_rho
            r= rho(i_rho); 
             %Calculate energy per node
            for n=1:length(x0)
                B=10e-5*ones(length(x0),length(x0)); B(n,n)=1;
                [X_opt, U_opt, n_err] = optim_fun(A_s, T, B, x0, xf, r, eye(length(A_s)));
                nodeEnergy(n,i_traj,i_rho)=sum(vecnorm(B*U_opt').^2);
                trajErr(n, i_traj,i_rho)=n_err;
            end
        end
        toc
    end
    energy(i_set).s1trajErr= trajErr;
    energy(i_set).s1NodeEnergy=nodeEnergy;
end

eng=nodeEnergy;

[X,Y,Z] = ndgrid(1:size(eng,1), 1:size(eng,2), 1:size(eng,3));
pointsize = 30;
scatter3(X(:), Y(:), Z(:), pointsize, trajErr(:));
xlabel('nodes')
ylabel('T')
zlabel('rho')
set(gca,'colorscale','log')

figure(2)
for i=1:5
    i
%imagesc(trajErr(:,:,i))
plot(mean(trajErr))
ylabel('n')
xlabel('T')
colorbar
pause(.5)
end
%%

figure(4); clf; hold on;
plot(squeeze(sum(sum(icov>0,1)))./85^2)
plot(squeeze(sum(sum(icov<0,1)))./85^2)

figure(5); clf; hold on;
plot(squeeze(sum(sum(pcm>0,1)))./85^2)
plot(squeeze(sum(sum(pcm<0,1)))./85^2)




