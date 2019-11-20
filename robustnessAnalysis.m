%% Display Robustness Outcomes
%% Sparsity tuning affect on density and controllability

%Density vs. rho for each patient
clear d;
cd('/Users/bscheid/Documents/LittLab/PROJECTS/p01_EC_controllability/v3/Code')
addpath('helper_functions')
addpath('pipeline_scripts')

betas= [0, .01, .05, .1 ];
gammas= [ 0, .1, .25, .5, .75];
npts=12;

betaStr=replace(string(betas), '.', '_');
gammaStr=replace(string(gammas), '.', '_');


for g=1:length(gammas)
    for i_pt=1:length(Network)
        load(sprintf('Data/Robustness/Network_%s.mat',gammaStr{g}))
        pcm= Network(i_pt).pcm;
        [N,~,T]=size(pcm);
        % Get density, ignore zeros. 
        d(g,i_pt)=(sum(pcm(:)~=0))./(numel(pcm)-N*T);
        avg_rho(g,i_pt)=mean(Network(i_pt).rhos);
        std_rho(g,i_pt)=std(Network(i_pt).rhos);
    end
end

%%

figure(3); clf; hold on
plot(gammas, d(:,1:npts))
xlabel('gamma'); ylabel('density')
xticks(gammas);
title('Density plot')

figure(2); clf; hold on
for z=1:npts
    errorbar(gammas, avg_rho(:,z), std_rho(:,z));
end
xlabel('gamma'); ylabel('rho values')
xticks(gammas);
title('avg. rho')

% How does sparsity affect average controllability values?

%% Contiguity Tuning (beta)

load(sprintf('Data/Robustness/Partitions_0_25'))
figure(1); clf
figure(2); clf
 
for i_pt=1:npts

p= Partitions(i_pt);
inds= find(strcmp(p.ID, {Partitions.ID}));
rns=cell(length(betas),1);

for b= 1:(2*length(betas))
    for i=1:size(Partitions(inds(b)).medianQcomms,2)
        s=assignStates(Partitions(inds(b)).medianQcomms(:,i));
        %rns{b,1}(i,:)=[Partitions(inds(b)).gamma(i), size(s.stateRuns,2)];  % # of temporal clusters
        rns{b,1}(i,:)=[Partitions(inds(b)).gamma(i), max(s.contigStates)]; % total # of communities  
    end
    figure(mod(b,2)+1); hold on;
    subplot(2,ceil(npts/2),i_pt)
    hold on; 
    plot(rns{b,1}(:,1), rns{b,1}(:,2), '--o', 'linewidth', .5)
     
end

% plot(p.gamma, rns)
figure(1)
xlabel('gamma'); ylabel('#comms'); title(p.ID)
ylim([0,10])

figure(2); xlabel('gamma'); ylabel('#comms'); title(p.ID)
ylim([0,10])

end

figure(1); suptitle('preictal')
figure(2); suptitle('ictal')

legend(betaStr, "Interpreter", "none")
%% View network for different beta

figure(2)
subplot(2,3,1)
imagesc(Network(i_pt).sim)
for i=1:length(betas)
    subplot(2,3,i+1)
    imagesc(Network(1).(['wSim_', betaStr{i}]))
end

%% State arrangement with different beta

load(sprintf('Data/Robustness/Partitions_0_1'))
npts= length(unique({Partitions.ID}));
figure(1); clf
figure(2); clf
 
for i_pt=1:npts

p= Partitions(i_pt);
ict= find(strcmp(p.ID, {Partitions.ID}) .* strcmp({'ictal'}, {Partitions.type}));
preict= find(strcmp(p.ID, {Partitions.ID}) .* strcmp({'preictal'}, {Partitions.type}));

for b= 1:length(betas)
    figure(1); hold on;
    subplot(length(betas),npts, npts*(b-1)+i_pt); hold on; 
    imagesc(Partitions(preict(b)).states) 
    axis tight; xticklabels([]); yticklabels([])
    caxis([0,4])
    if i_pt==1, ylabel(betas(b)); end
    if b==1, title(p.ID); end
    
    figure(2); hold on;
    subplot(length(betas),npts, npts*(b-1)+i_pt); hold on; 
    imagesc(Partitions(ict(b)).states)   
    axis tight; xticklabels([]); yticklabels([])
    caxis([0,4])
    if i_pt==1, ylabel(betas(b)); end
    if b==1, title(p.ID); end
    
end

end

figure(1); suptitle('preictal')
figure(2); suptitle('ictal')


