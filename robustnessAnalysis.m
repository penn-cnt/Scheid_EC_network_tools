%% Display Robustness Outcomes

cd('/Users/bscheid/Documents/LittLab/PROJECTS/p01_EC_controllability/v3/Code')
addpath('helper_functions')
addpath('pipeline_scripts')
addpath(genpath('~/Documents/CODE'))

betas= [0, .01, .05, .1 ];
gammas= [ 0, .1, .25, .5, .75];
npts=14;
nStates=3;

betaStr=replace(string(betas), '.', '_');
gammaStr=replace(string(gammas), '.', '_');

cols=[[227,187,187]; [190,8,4]; [138,4,4];[140,42,195];[75,184,166];[242,224,43];[74,156,85];...
   [80,80,80]; [255,255,255]]/255;

%% How does sparsity affect density, rho values, distribution of values?

d=[];
d_mn=[]; d_std=[];
sign_mn=[]; sign_std=[];

figure(3); clf;
set(groot,'defaultAxesColorOrder','default')

for g=1:length(gammas)
        load(sprintf('Data/Robustness/Network_%s.mat',gammaStr{g}))
    for i_pt= npts+(1:npts)
        pcm= Network(i_pt).pcm;
        [N,~,T]=size(pcm);
        % Get total density, ignore zeros. 
        %d(g,i_pt)=(sum(pcm(:)~=0))./(numel(pcm)-N*T);
        
        d_mn(g,i_pt)=mean(squeeze(sum(pcm~=0, [1,2]))./(N^2-N));
        d_std(g,i_pt)=std(squeeze(sum(pcm~=0, [1,2]))./(N^2-N));
        
        avg_rho(g,i_pt)=mean(Network(i_pt).rhos);
        std_rho(g,i_pt)=std(Network(i_pt).rhos);
        
        sign_mn(g,i_pt)=mean(squeeze(sum(pcm<0, [1,2]))./(N^2-N));
        sign_std(g,i_pt)=std(squeeze(sum(pcm<0, [1,2]))./(N^2-N));
   
    end
end

%%
cmap=brewermap(66,'PuBuGn');      
set(groot,'defaultAxesColorOrder',cmap(20:6:66,:))

% Ictal Density
figure(10); clf; hold on
for z=1:npts
errorbar(gammas, d_mn(:,z), d_std(:,z),'linewidth', 1.5)
end
xlabel('\gamma', 'fontsize', 15); ylabel('Average network density')
xticks(gammas); xlim([gammas(1)-.05, gammas(end)+.05]);

% Rho values
figure(11); clf; hold on
for z=1:npts
    errorbar(gammas, avg_rho(:,z), std_rho(:,z), 'linewidth', 1.5);
end
xlabel('\gamma'); ylabel('Average \rho value')
xticks(gammas); xlim([gammas(1)-.05, gammas(end)+.05]);

% 
% % Sign
% figure(4); clf; hold on
% for z=1:npts
% errorbar(gammas, sign_mn(:,z), sign_std(:,z),'linewidth', 1.5)
% end
% xlabel('\gamma'); ylabel('density')
% xticks(gammas); xlim([gammas(1)-.05, gammas(end)+.05]);
% title('Edge Sign')

% How does sparsity affect average controllability values?

%% Contiguity Tuning 1 (beta vs. gamma (?))

load(sprintf('Data/Robustness/Partitions_0_25'))
figure(1); clf
figure(2); clf
 
for i_pt=1:npts

p= Partitions(i_pt);
inds= find(strcmp(p.ID, {Partitions.ID}));
rns=cell(length(betas),1);

for b= 1:(2*length(betas))
    for i=1:size(Partitions(inds(b)).medianQcomms,2)
        %s=assignStates(Partitions(inds(b)).medianQcomms(:,i));
        rns{b,1}(i,:)=[p(inds(b)).gamma(i), size(p(inds(b)).stateRuns,2)];  % # of temporal clusters
        %rns{b,1}(i,:)=[Partitions(inds(b)).gamma(i), max(p(inds(b)).contigStates)]; % total # of communities  
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

%% Contiguity Tuning 2 (beta)

% How do number of temoral demarkations change with increasing beta for all
% patients (ictal vs. interical) 

%load(sprintf('Data/Robustness/Partitions_0_25.mat'))

contigs_ict=zeros(npts,length(betas));
contigs_preict=zeros(npts,length(betas));
 
for i_pt=1:npts

p= Partitions(i_pt);
inds= find(strcmp(p.ID, {Partitions.ID}));

for b= 1:2:length(betas)*2
    contigs_ict(i_pt,(b+1)/2)=size(Partitions(inds(b)).stateRuns,2);  % # of temporal clusters
    contigs_preict(i_pt, (b+1)/2)=size(Partitions(inds(b+1)).stateRuns,2);
    
    % longest contig state/ total state length
    pct_ict(i_pt,(b+1)/2)=min(arrayfun(@(s)sum(Partitions(inds(b)).contigStates==s)/sum(Partitions(inds(b)).states==s), [1:3]))
    pct_preict(i_pt, (b+1)/2)=min(arrayfun(@(s)sum(Partitions(inds(b+1)).contigStates==s)/sum(Partitions(inds(b+1)).states==s), [1:3]))
end

end

% cmap=parula;
% %set(groot,'defaultAxesColorOrder',summer)
% colororder(cmap)

cmap=brewermap(66,'PuBuGn');      
set(groot,'defaultAxesColorOrder',cmap(20:6:66,:))
         
% Total number of temporal groups
figure(80); clf; hold on
subplot(1,2,1)
plot(betas, contigs_preict+normrnd(0,.9, size(contigs_preict)), 'o-', 'linewidth', 1,...
        'markeredge', 'black')
ylim([0, max(contigs_preict(:))+5])
xlabel('\beta', 'fontsize',  15)
xlim([betas(1)-.005, betas(end)+.005])
title('Preictal')
ylabel('# temporally contig. groups')

subplot(1,2,2)
plot(betas, contigs_ict+normrnd(0,.05, size(contigs_ict)), 'o-', 'linewidth', 1, ...
    'markeredge', 'black')
ylim([0, max(contigs_ict(:))])
xlabel('\beta', 'fontsize',  15)
xlim([betas(1)-.005, betas(end)+.005])
title('Ictal','fontweight','normal')
ylim([0,10])

% Total number of temporal groups
figure(81); clf; hold on
subplot(1,2,1)
plot(betas,  pct_preict+normrnd(0,.005, size(pct_preict)), 'o-', 'linewidth', 1,...
        'markeredge', 'black')
ylim([0, max( pct_preict(:))])
xlabel('beta')
title('Preictal','fontweight','normal')
ylabel('# temporal groups')
ylim([0,1.01])

subplot(1,2,2)
plot(betas, pct_ict+normrnd(0,.001, size(pct_ict)), 'o-', 'linewidth', 1,...
        'markeredge', 'black')
xlabel('beta')
title('Ictal')
ylabel('# temporal groups')
% Contiguous length/total length
ylim([0,1.01])

%% View network for different beta

figure(2)
subplot(2,3,1)
imagesc(Network(i_pt).sim)
for i=1:length(betas)
    subplot(2,3,i+1)
    imagesc(Network(1).(['wSim_', betaStr{i}]))
end

%% State arrangement with different beta

load(sprintf('Data/Robustness/Partitions_0_5'))
npts= length(unique({Partitions.ID}));
figure(1); clf
figure(2); clf
 
for i_pt= 1:npts

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

% Show only a single instance

figure(70); clf; hold on;
i_pt=npts
for b=1:length(betas)
    subplot(length(betas),2,1+2*(b-1))
    imagesc(Partitions(preict(b)).states)
    colormap(gca, cols(1:3,:)*1.1)
    axis tight; xticklabels([]); yticklabels([])
    ylabel(betas(b), 'fontsize', 14);
    
    subplot(length(betas),2,2+2*(b-1))
    imagesc(Partitions(ict(b)).states)
    colormap(gca, cols(1:3,:)*1.1)
    axis tight; xticklabels([]); yticklabels([])

end

%% Metric Comparison with different Betas

load(sprintf('Data/Robustness/Partitions_0_5'))
load(sprintf('Data/Robustness/State_metrics_0_5'))
npts= length(unique({Partitions.ID}));
metrics={'aveCtrl', 'modalCtrl', 'tModalCtrl', 'pModalCtrl'}; 

figure(1); clf
figure(2); clf

alpha=0.05;
ph=1;

[p_diff, rho_diff, p_diff_ict, rho_diff_ict]= deal(zeros(length(betas)-1, nStates, length(metrics)));
ict_met=[]; preict_met=[];
 
ctr=1; 
for m=1:length(metrics)
    for b= 1:length(betas)        
        i_ict= (1:npts)+(2*npts)*(b-1);
        i_preict= i_ict+npts; 
        
        tmp_ict=ict_met;
        tmp_preict=preict_met;
        

       ict_metZ= cell2mat(cellfun(@mean,...
            {State_metrics(i_ict).([metrics{m},'Z'])}, 'UniformOutput', false)');
       preict_metZ= cell2mat(cellfun(@mean,...
            {State_metrics(i_preict).([metrics{m},'Z'])}, 'UniformOutput', false)');
       
        
       ict_met= cell2mat(cellfun(@mean,...
            {State_metrics(i_ict).(metrics{m})}, 'UniformOutput', false)');
       preict_met= cell2mat(cellfun(@mean,...
            {State_metrics(i_preict).(metrics{m})}, 'UniformOutput', false)');
       
        % Perform Wilcoxon signed rank test
       if b>1
           p_diff(b-1,:,m)= arrayfun(@(x)signrank(tmp_ict(:,x), ict_met(:,x)), [1:3]) 
           p_diff_preict(b-1,:,m)= arrayfun(@(x)signrank(tmp_preict(:,x), preict_met(:,x)), [1:3]) 
           
           [~,~,sss]=arrayfun(@(x)signrank(tmp_ict(:,x), ict_met(:,x)), [1:3],'UniformOutput', false)
           sss=cell2mat(sss);
           [~,~,sss_preict]= arrayfun(@(x)signrank(tmp_preict(:,x), preict_met(:,x)), [1:3]) 
           
           rho_diff(b-1,:,m)= [sss.signedrank];
           rho_diff_preict(b-1,:,m)= [sss_preict.signedrank];
       end

        figure(1); hold on;
        subplot(length(betas),length(metrics), ctr); hold on;    
        h= boxplot(ict_metZ, {'Phase 1', 'Phase 2', 'Phase 3'},'Colors', cols(1:3,:));
        set(h,{'linew'},{2})
        if m==1; ylabel(sprintf('beta: %0.2f', betas(b))); end
        if b==1; title(sprintf('%s', metrics{m})); end
        
        [~,~,cc]=friedman(ict_met,1, 'off'); [ee]=multcompare(cc, 'display', 'off');
        sigstar(num2cell(ee(ee(:,6)<=alpha,[1:2]),2))
        
        figure(2); hold on;
        subplot(length(betas),length(metrics), ctr); hold on;    
        h= boxplot(preict_metZ, {'Phase 1', 'Phase 2', 'Phase 3'},'Colors', cols(1:3,:));
        set(h,{'linew'},{2})
        if m==1; ylabel(sprintf('beta: %0.2f', betas(b))); end
        if b==1; title(sprintf('%s', metrics{m})); end
        
        [~,~,cc]=friedman(preict_met,1, 'off'); [ee]=multcompare(cc,'display', 'off');
        sigstar(num2cell(ee(ee(:,6)<=alpha,[1:2]),2))
        
        ctr=ctr+1;

    end    
end


figure(1); suptitle('ictal')
figure(2); suptitle('preictal')

% significance of correlation between b=0 and b=0.01 values for all metrics
% and phases

squeeze(p_diff(1,:,:)) % Wilcox signed rank test p-value
squeeze(p_diff_preict(1,:,:))

squeeze(rho_diff(1,:,:)) % Wilcox signed rank test W statistic
squeeze(rho_diff_preict(1,:,:))

%% Metric change with beta for individuals
figure(90); clf; hold on
metric= 'aveCtrl';
s=2;

beta_metric=zeros(npts, length(betas));
beta_metric_preict=zeros(npts, length(betas));

for b=1:length(betas)
    i_ict= (1:npts)+(2*npts)*(b-1);
    i_preict= i_ict+npts; 
    
    ict_metZ= cell2mat(cellfun(@mean, {State_metrics(i_ict).([metric,'Z'])},...
        'UniformOutput', false)');
    pre_metZ= cell2mat(cellfun(@mean, {State_metrics(i_preict).([metric,'Z'])},...
        'UniformOutput', false)');
    
    beta_metric(:,b)=ict_metZ(:,s);
    beta_metric_preict(:,b)=pre_metZ(:,s)
end

subplot(1,2,1)
plot(betas,beta_metric_preict, 'o-', 'linewidth', 1.5, 'markeredge', 'black')
xlabel('\beta', 'fontsize',  15)
title('Preictal','fontweight','normal'); ylabel(sprintf('%s (z-scored)',metric))
xlim([betas(1)-.005, betas(end)+.005])
ylim([min(beta_metric_preict(:))-.05, max(beta_metric_preict(:)+.05)])

subplot(1,2,2)
plot(betas,beta_metric, 'o-', 'linewidth', 1.5, 'markeredge', 'black')
xlabel('\beta', 'fontsize',  15)
title('Ictal','fontweight','normal')
xlim([betas(1)-.005, betas(end)+.005])
ylim([min(beta_metric(:))-.05, max(beta_metric(:)+.05)])

%% Single Metric beta change

load(sprintf('Data/Robustness/Partitions_0_5'))
load(sprintf('Data/Robustness/State_metrics_0_5'))
npts= length(unique({Partitions.ID}));
metrics={'aveCtrl', 'modalCtrl', 'tModalCtrl', 'pModalCtrl'}; 

figure(1); clf
figure(2); clf

alpha=0.05;
 
ctr=1; 
for b= 1:length(betas)
    for m=1:length(metrics)
        
        i_ict= (1:npts)+(2*npts)*(b-1);
        i_preict= i_ict+npts; 

       ict_metZ= cell2mat(cellfun(@mean,...
            {State_metrics(i_ict).([metrics{m},'Z'])}, 'UniformOutput', false)');
       preict_metZ= cell2mat(cellfun(@mean,...
            {State_metrics(i_preict).([metrics{m},'Z'])}, 'UniformOutput', false)');
        
       ict_met= cell2mat(cellfun(@mean,...
            {State_metrics(i_ict).(metrics{m})}, 'UniformOutput', false)');
       preict_met= cell2mat(cellfun(@mean,...
            {State_metrics(i_preict).(metrics{m})}, 'UniformOutput', false)');

        figure(1); hold on;
        subplot(length(betas),length(metrics), ctr); hold on;    
        h= boxplot(ict_metZ, {'Phase 1', 'Phase 2', 'Phase 3'},'Colors', cols(1:3,:));
        set(h,{'linew'},{2})
        if m==1; ylabel(sprintf('beta: %0.2f', betas(b))); end
        if b==1; title(sprintf('%s', metrics{m})); end
        
        [~,~,cc]=friedman(ict_met,1, 'off'); [ee]=multcompare(cc, 'display', 'off');
        sigstar(num2cell(ee(ee(:,6)<=alpha,[1:2]),2))
        
        figure(2); hold on;
        subplot(length(betas),length(metrics), ctr); hold on;    
        h= boxplot(preict_metZ, {'Phase 1', 'Phase 2', 'Phase 3'},'Colors', cols(1:3,:));
        set(h,{'linew'},{2})
        if m==1; ylabel(sprintf('beta: %0.2f', betas(b))); end
        if b==1; title(sprintf('%s', metrics{m})); end
        
        [~,~,cc]=friedman(preict_met,1, 'off'); [ee]=multcompare(cc,'display', 'off');
        sigstar(num2cell(ee(ee(:,6)<=alpha,[1:2]),2))
        
        ctr=ctr+1;

    end
end


figure(1); suptitle('ictal')
figure(2); suptitle('preictal')

%% Metric comparison across gamma values (beta=0.01)

metrics={'aveCtrl', 'modalCtrl', 'tModalCtrl', 'pModalCtrl'}; 

figure(1); clf
figure(2); clf

alpha=0.05;
ph=1;
b=betas(2); % Select beta value

[p_diff, rho_diff, p_diff_ict, rho_diff_ict]= deal(zeros(length(betas)-1, nStates, length(metrics)));
[ict_met, preict_met, ict_metZ, preict_metZ]= deal(zeros(npts, 3, length(metrics)));
 
ctr=1; 
i_ict= (1:npts); i_preict= i_ict+npts; 
s=3

for g= 1:length(gammas)    
    load(sprintf('Data/Robustness/State_metrics_%s.mat',gammaStr{g}))
     
    tmp_ict=ict_met;
    tmp_preict=preict_met;
    
    for m=1:length(metrics)

        ict_metZ(:,:,m)= cell2mat(cellfun(@mean, {State_metrics(i_ict).([metrics{m},'Z'])}, 'UniformOutput', false)');
        preict_metZ(:,:,m)= cell2mat(cellfun(@mean, {State_metrics(i_preict).([metrics{m},'Z'])}, 'UniformOutput', false)');

        ict_met(:,:,m)= cell2mat(cellfun(@mean,{State_metrics(i_ict).(metrics{m})}, 'UniformOutput', false)');
        preict_met(:,:,m)= cell2mat(cellfun(@mean,{State_metrics(i_preict).(metrics{m})}, 'UniformOutput', false)');

%         figure(1); hold on;
%         subplot(length(gammas),length(metrics), ctr); hold on;    
%         h= boxplot(ict_metZ(:,:,m), {'Phase 1', 'Phase 2', 'Phase 3'},'Colors', cols(1:3,:));
%         set(h,{'linew'},{2})
%         if m==1; ylabel(sprintf('gamma: %0.2f', gammas(g))); end
%         if g==1; title(sprintf('%s', metrics{m})); end
% 
%         [~,~,cc]=friedman(ict_met(:,:,m),1, 'off'); [ee]=multcompare(cc, 'display', 'off');
%         sigstar(num2cell(ee(ee(:,6)<=alpha,[1:2]),2))
% 
%         figure(2); hold on;
%         subplot(length(gammas),length(metrics), ctr); hold on;    
%         h= boxplot(preict_metZ(:,:,m), {'Phase 1', 'Phase 2', 'Phase 3'},'Colors', cols(1:3,:));
%         set(h,{'linew'},{2})
%         if m==1; ylabel(sprintf('gamma: %0.2f', gammas(g))); end
%         if g==1; title(sprintf('%s', metrics{m})); end
% 
%         [~,~,cc]=friedman(preict_met(:,:,m),1, 'off'); [ee]=multcompare(cc,'display', 'off');
%         sigstar(num2cell(ee(ee(:,6)<=alpha,[1:2]),2))
% 
%         ctr=ctr+1;
%         
%        % Perform Wilcoxon signed rank test
       if g > 1
           p_diff(g-1,:,m)= arrayfun(@(x)signrank(tmp_ict(:,x,m), ict_met(:,x,m)), [1:3]) 
           p_diff_preict(g-1,:,m)= arrayfun(@(x)signrank(tmp_preict(:,x), preict_met(:,x)), [1:3]) 
           
           [~,~,sss]=arrayfun(@(x)signrank(tmp_ict(:,x,m), ict_met(:,x,m)), [1:3],'UniformOutput', false)
           sss=cell2mat(sss);
           [~,~,sss_preict]= arrayfun(@(x)signrank(tmp_preict(:,x,m), preict_met(:,x,m)), [1:3]) 
           
           rho_diff(g-1,:,m)= [sss.signedrank];
           rho_diff_preict(g-1,:,m)= [sss_preict.signedrank];
       end          
    end  
    
        modal_metric(:,g)= ict_metZ(:,s,1);
        modal_metric_preict(:,g)= preict_metZ(:,s,1);
    
end


figure(1); suptitle('ictal')
figure(2); suptitle('preictal')

% significance of correlation between b=0 and b=0.01 values for all metrics
% and phases

squeeze(p_diff(3,:,:)) % Wilcox signed rank test p-value
squeeze(p_diff_preict(1,:,:))

squeeze(rho_diff(3,:,:)) % Wilcox signed rank test W statistic
squeeze(rho_diff_preict(1,:,:))

figure(32)
subplot(1,2,1)
plot(gammas,modal_metric_preict, 'o-', 'linewidth', 1.5, 'markeredge', 'black')
xlabel('\gamma', 'fontsize',  15)
title('Preictal','fontweight','normal'); ylabel(sprintf('%s (z-scored)',metrics{1}))
xlim([gammas(1)-.005, gammas(end)+.005])
%ylim([min(beta_metric_preict(:))-.05, max(beta_metric_preict(:)+.05)])

subplot(1,2,2)
plot(gammas,modal_metric, 'o-', 'linewidth', 1.5, 'markeredge', 'black')
xlabel('\gamma', 'fontsize',  15)
title('Ictal','fontweight','normal')
xlim([gammas(1)-.005, gammas(end)+.005])
% ylim([min(beta_metric(:))-.05, max(beta_metric(:)+.05)])

%% Persistent Transient Modal Controllability 

load('Data/Robustness/modalThresholds.mat')

% figure(62); clf
% figure(63); clf
% figure(2); clf
% figure(3); clf

alpha=0.017;
ph=1;
b=betas(2); % Select beta value

%[p_diff, rho_diff, p_diff_ict, rho_diff_ict]= deal(zeros(length(betas)-1, nStates, length(metrics)));
[ict_met, preict_met, ict_metZ, preict_metZ]= deal(zeros(npts, 3, length(thresholds)));

thresholds=[.1, .15, .2]; 
ctr=1; 

metrics= {'pModalCtrl', 'tModalCtrl'};
  
for m= 1:2
    for th=1:length(thresholds)
    
        % Get indices
        th_inds= (round([modalStateRobustness.thresh],2) == thresholds(th));
        i_ict= find(th_inds.*strcmp({modalStateRobustness.type}, 'ictal'));
        i_preict=find(th_inds.*strcmp({modalStateRobustness.type}, 'preictal'));
        
%         tmp_ict=ict_met;
%         tmp_preict=preict_met;
        
       ict_metZ(:,:,th)= cell2mat(cellfun(@mean, {modalStateRobustness(i_ict).([metrics{m},'Z'])}, 'UniformOutput', false)');
       preict_metZ(:,:,th)= cell2mat(cellfun(@mean,{modalStateRobustness(i_preict).([metrics{m},'Z'])}, 'UniformOutput', false)');
        
       ict_met(:,:,th)= cell2mat(cellfun(@mean, {modalStateRobustness(i_ict).(metrics{m})}, 'UniformOutput', false)');
       preict_met(:,:,th)= cell2mat(cellfun(@mean, {modalStateRobustness(i_preict).(metrics{m})}, 'UniformOutput', false)');

        % Show boxplots
        figure(62); hold on;
        subplot(length(metrics), length(thresholds), ctr); hold on;    
        h= boxplot(squeeze(ict_metZ(:,:,th)), {'Phase 1', 'Phase 2', 'Phase 3'},'Colors', cols(1:3,:));
        set(h,{'linew'},{2})
        if m==1; title(sprintf('thresh: %0.2f',thresholds(th))); end
        if th==1; ylabel(sprintf('%s', metrics{m})); end
        
        [~,statt,cc]=friedman(squeeze(ict_met(:,:,th)),1, 'off'); [ee]=multcompare(cc, 'display', 'off');
        
        sigstar(num2cell(ee(ee(:,6)<=alpha,[1:2]),2))
        ict_pvals(m,:,th)=ee(:,6);
        
        figure(63); hold on;
        subplot(length(metrics),length(thresholds), ctr); hold on;    
        h= boxplot(squeeze(preict_metZ(:,:,th)), {'Phase 1', 'Phase 2', 'Phase 3'},'Colors', cols(1:3,:));
        set(h,{'linew'},{2})
        
        if m==1; title(sprintf('thresh: %0.2f',thresholds(th))); end
        if th==1; ylabel(sprintf('%s', metrics{m})); end
        
        [~,~,cc]=friedman(squeeze(preict_met(:,:,th)),1, 'off'); [ee]=multcompare(cc,'display', 'off');
        sigstar(num2cell(ee(ee(:,6)<=alpha,[1:2]),2))
        
         ctr=ctr+1;
%               
%         % Perform wilcoxon signed rank test
%         if th>1
%            p_diff(th-1,:,m)= arrayfun(@(x)signrank(tmp_ict(:,x), ict_met(:,x)), [1:3]) 
%            p_diff_preict(th-1,:,m)= arrayfun(@(x)signrank(tmp_preict(:,x), preict_met(:,x)), [1:3]) 
%            
%            [~,~,sss]=arrayfun(@(x)signrank(tmp_ict(:,x), ict_met(:,x)), [1:3],'UniformOutput', false)
%            sss=cell2mat(sss);
%            [~,~,sss_preict]= arrayfun(@(x)signrank(tmp_preict(:,x), preict_met(:,x)), [1:3]) 
%            
%            rho_diff(th-1,:,m)= [sss.signedrank];
%            rho_diff_preict(th-1,:,m)= [sss_preict.signedrank];
%         end

    
    end
    
    % Display actual values
%     figure(m)
%     for s=1:3  
%         subplot(3,2,2*(s-1)+1)
%         plot(thresholds,squeeze(preict_met(:,s,:)), 'o-', 'linewidth', 1.5, 'markeredge', 'black')
%         xlabel('threshold', 'fontsize',  15)
%         title('Preictal','fontweight','normal'); ylabel(sprintf('%s (z-scored)',metrics{m}))
% 
%         subplot(3,2,2*(s-1)+2)
%         plot(thresholds,squeeze(ict_met(:,s,:)), 'o-', 'linewidth', 1.5, 'markeredge', 'black')
%         xlabel('threshold', 'fontsize',  15)
%         title('Ictal','fontweight','normal')
%     end
    
end


 figure(62); suptitle('ictal')
 figure(63); suptitle('preictal')
% 
% squeeze(p_diff(1,:,:)) % Wilcox signed rank test p-value
% squeeze(p_diff_preict(1,:,:))
% 
% squeeze(rho_diff(1,:,:)) % Wilcox signed rank test W statistic
% squeeze(rho_diff_preict(1,:,:))

%% Get all final figures

% gamma figures

% figure(10) % beta contig
% set(gcf, 'Position', [184   663   290   236])
% saveas(gcf, 'FigsV3.3/robustness/gamma_density.png')
% saveas(gcf, 'FigsV3.3/robustness/gamma_density.fig')
% 
% figure(11) % beta contig
% set(gcf, 'Position', [184   663   290   236])
% saveas(gcf, 'FigsV3.3/robustness/gamma_rho.png')
% saveas(gcf, 'FigsV3.3/robustness/gamma_rho.fig')

% figure(33) % beta contig
% set(gcf, 'Position', [184   634   582   265])
% saveas(gcf, 'FigsV3.3/robustness/gamma_metric_s3.png')
% saveas(gcf, 'FigsV3.3/robustness/gamma_metric_s3.fig')


% beta figures 

% figure(80) % beta contig
% set(gcf, 'Position', [1119 690 511 264])
% saveas(gcf, 'FigsV3.3/robustness/beta_contig.png')
% saveas(gcf, 'FigsV3.3/robustness/beta_contig.fig')

% figure(70) % beta contig
% set(gcf, 'Position', [1300 494 346 461])
% saveas(gcf, 'FigsV3.3/robustness/beta_states.png')
% saveas(gcf, 'FigsV3.3/robustness/beta_states.fig')

% figure(90) % beta contig
% set(gcf, 'Position', [1119 690 511 264])
% saveas(gcf, 'FigsV3.3/robustness/beta_metric.png')
% saveas(gcf, 'FigsV3.3/robustness/beta_metric.fig')

% transient figures

% figure(62) % modal threshold
% set(gcf, 'Position', [411   390   654   368])
% saveas(gcf, 'FigsV3.3/robustness/modal_thresh.png')
% saveas(gcf, 'FigsV3.3/robustness/modal_thresh.fig')











