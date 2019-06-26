%V3_analysis

Null='';
load('Data/dataSets_clean.mat');
load('Data/subjects.mat')
load(sprintf('Data/%sNetworks.mat', Null)); 
load(sprintf('Data/%sPartitions.mat', Null));
load(sprintf('Data/%sMetric_matrices.mat', Null)); 
load(sprintf('Data/%sState_metrics.mat', Null))

eval(['Networks=', Null, 'Networks'])
eval(['Partitions=', Null, 'Partitions'])
eval(['Metric_matrices=', Null, 'Metric_matrices'])
eval(['State_metrics=', Null, 'State_metrics'])

if strcmp(Null,'Null')
    load('Data/nullData.mat')
    dataSets_clean=nullData;
end

addpath(genpath('~/Documents/CODE/'))

cols=[[75,184,166];[255,168,231]; [36,67,152];[140,42,195];[121,29,38];[242,224,43];[74,156,85];...
   [80,80,80]; [255,255,255]]/255;

i_ict=find(strcmp({Partitions.type},'ictal'));
i_preict=find(strcmp({Partitions.type},'preictal'));

nSets=length(Partitions);
nIct=nSets/2;

                %% States Analysis %%
                
% Display Median State Times vs.Lengths of Largest 3 states
% note, must run above section first z(tran)=m

meds_preict=cell2mat({Partitions(i_preict).stateMedian}');
meds_ict=cell2mat({Partitions(i_ict).stateMedian}');
lens_preict=cell2mat({Partitions(i_preict).runLen}');
lens_ict=cell2mat({Partitions(i_ict).runLen}');
labels(1:nIct)={'State 1'}; labels(nIct+1:2*nIct)={'State 2'}; labels(2*nIct+1:3*nIct)={'State 3'};
  
figure(1); clf; hold on
scatterhist(meds_ict(:),lens_ict(:),'Group',labels','Kernel','on','Location','NorthEast',...
    'Direction','out','Color',cols(1:3,:),'LineStyle',{'-','-','-'},...
    'LineWidth',[2,2,2],'Marker','ddd','MarkerSize',[5,5,5]);
ylabel(''); xlabel('');
set(gca,'xaxisLocation','bottom', 'yaxislocation', 'left', 'fontsize', 18)

figure(2); clf; hold on
scatterhist(meds_preict(:),lens_preict(:),'Group',labels','Kernel','on','Location','NorthEast',...
    'Direction','out','Color',cols(1:3,:),'LineStyle',{'-','-','-'},...
    'LineWidth',[2,2,2],'Marker','ddd','MarkerSize',[5,5,5]);
ylabel(''); xlabel('');
set(gca,'xaxisLocation','bottom', 'yaxislocation', 'left','fontsize', 18)

%% Show table of mean and std statistics for ictal and preictal states
mn_len=mean(lens_ict); std_len=std(lens_ict);
mn_med=mean(meds_ict);std_med=std(meds_ict);
mn_lenpi=mean(lens_preict);std_lenpi=std(lens_preict);
mn_medpi=mean(meds_preict);std_medpi=std(meds_preict);

stateTable=table([mn_len', std_len'], [mn_med', std_med'], [mn_lenpi', std_lenpi'], [mn_medpi', std_medpi'],...
    'VariableNames', {'ictal_longest_run', 'ictal_med_value','preictal_longest_run',...
    'preictal_med_value'})
writetable(stateTable,'Data/stateTable.csv')

clear mn_len mn_med mn_lenpi mn_medpi labels x_ict y_ict x_preict y_preict meds_preict lens_preict ...
    meds_ict lens_ict

            %% Metric Trends (individual) %%

% Select metrics to plot
n=2; m=4;
all=false; % show both preictal and ictal
metrics={'aveCtrl', 'modalCtrl', 'tModalCtrl','pModalCtrl'}; 
%metrics={'degree', 'aveCtrl', 'modalCtrl'}; 

for i_set=i_ict
    
    figure(3); clf; 
    p= Partitions(i_set);
    
    if isempty(p.states)
        continue
    end
    
    st= p.states; mm= Metric_matrices(i_set);
    
    if all
        st=[Partitions(i_set+39).states+3, st];
        mm_pre= Metric_matrices(i_set+39); 
    end
    
    [N, T]=size(mm.degree);
    suptitle(sprintf('Network Metrics for %s, %s %s', p.ID, p.type, p.block))
 
    % Display Metrics
    for i=1:length(metrics)
        if all
            met=[mm_pre.(metrics{i}), mm.(metrics{i})];
        else 
            met= mm.(metrics{i});
        end
        
        subplot(n,m,i); hold on
        imagesc(met);
        stem((diff(st)~=0)*N,'Marker', 'none', 'lineWidth', 2, 'color', 'red')
        yyaxis right
        plot(mean(met), 'color', 'cyan');
        title(metrics{i}); axis tight
    end

    % Display signal
    subplot(n,m,(n*m)-2); hold on
    Fs= dataSets_clean(i_set).Fs;
    data= dataSets_clean(i_set).data'+(1:1000:1000*(N));
    plot((0:length(data)-1)/Fs,data);
    %stem([(diff(st)~=0)]*1000*N,'Marker', 'none', 'lineWidth', 2, 'color', 'red')
    title('Signal'); axis tight

    % Config matrix
    subplot(n,m,(n*m)-1); hold on
    imagesc(Networks(i_set).config_pcm)
    set(gca,'colorscale','log')
    stem((diff(st)~=0)*(N*(N-2))/2,'Marker', 'none', 'lineWidth', 2, 'color', 'red')
    title('Icov'); axis tight
    
    % Similarity Matrix
    subplot(n,m,(n*m)); hold on
    imagesc(Networks(i_set).sim)
    set(gca,'colorscale','log')
    stem((diff(st)~=0)*T,'Marker', 'none', 'lineWidth', 2, 'color', 'red')
    title('Sim'); axis tight 

    figure(4)
    %set(gca, 'Position', get(gca, 'Position')+[0,.05,0,0])
    imagesc(st) 
    colormap(gca, cols(1:3,:));
    set(gca, 'YTick', [], 'fontsize', 18)

    pause
end

%% Is there a difference between states?

ctype='bonferroni';
display='off';

analysis=struct();
metrics={'aveCtrl', 'modalCtrl', 'tModalCtrl','pModalCtrl', 'strength'}; 

for i_set=1:nSets
    
    s=State_metrics(i_set);
    if isempty(s.aveCtrl)
        continue
    end
    
    for i=1:length(metrics)
        % Perform Friedman's analysis & post hoc (non parametric rank test)
        eval(sprintf('[~,~,stats_%s]=friedman(s.%s(:,1:3),1,display);', metrics{i}, metrics{i}));
        eval(sprintf('[c_%s,m_%s]= multcompare(stats_%s, ''ctype'', ctype, ''display'', display);', ...
            metrics{i}, metrics{i}, metrics{i}));
        
        % Save P value, difference, and global average
        eval(sprintf('analysis(i_set).P%s=c_%s(:,6);',metrics{i}, metrics{i}));
        eval(sprintf('diffs%s(i_set,:)=c_%s(:,4);', metrics{i}, metrics{i}));
        eval(sprintf('glob%s(i_set,:)=mean(s.%s(:,1:3));', metrics{i}, metrics{i})); 
    end
end

display='off';
metrics=[metrics, 'optEnergy'];
% Get group level averages
for i=1:length(metrics)
    
    % ictal 
    eval(sprintf('[~,~,stats_g%s]=friedman(glob%s(i_ict,1:3),1,display)', metrics{i},  metrics{i}));
    eval(['[c_i_ict_glob',metrics{i},',m_ict_g',metrics{i},...
        ']= multcompare(stats_g',metrics{i},', ''ctype'', ctype, ''display'', display);'])
    %preicatl
    eval(sprintf('[~,~,stats_g%s]=friedman(glob%s(i_preict,1:3),1,display)', metrics{i},  metrics{i}));
    eval(['[c_i_preict_glob',metrics{i},',m_ict_g',metrics{i},...
        ']= multcompare(stats_g',metrics{i},', ''ctype'', ctype, ''display'', display);'])
end
    
disp('done')
%% Display results of Friedman's test individually
ctr=1;

for type=[i_ict', i_preict']
    figure(ctr)
    alpha=0.01;
    nMeas=3;   % number of measures
    nSamp=length([analysis(type).PmodalCtrl]);
    
    ma=reshape([analysis(type).PmodalCtrl],3, nSamp)'<alpha;
    sa=reshape([analysis(type).PaveCtrl],3, nSamp)'<alpha;
    st=reshape([analysis(type).Pstrength],3, nSamp)'<alpha; 

    dma=diffsmodalCtrl(type,:).*ma;
    dsa=diffsaveCtrl(type,:).*sa;
    ddst=diffsstrength(type,:).*st;

    clf; hold on;
    colormap(gca,cols([9,5,2,6],:));
    imagesc([sa,ma,st].*repmat([1:3],nSamp,nMeas))
    diffs=[dsa,dma,ddst];

    plot((diffs>0).*repmat((1:nMeas*3),nSamp,1),repmat((1:nSamp)',1,nMeas*3),'v','color', 'cyan')
    plot((diffs<0).*repmat((1:nMeas*3),nSamp,1),repmat((1:nSamp)',1,nMeas*3),'^','color', 'white')

    yticks([1:nSamp])
    yticklabels({Metric_matrices(type).ID});
    xticks([2:3:nMeas*3-1])
    xticklabels({'avg', 'modal','strength'})
    stem(3*[1:nMeas]+.5,ones(1,nMeas)*length(sa),...
        'Marker', 'none', 'lineWidth', 0.5, 'color', 'black')
    title({'Significant effect of seizure state on metric'})
    %'1-2, 1-3, 2-3'
    axis tight
    yyaxis right
    ylim([1,nSamp])
    yticks((1:nSamp))
    %yticklabels({numSOZ(sz_pairs(:,1))'./numChannels(sz_pairs(:,1))'*100});
    %yticklabels(sprintf('%d/%d\n', [numSOZ(i_ict)',numChannels(i_ict)']'))
    ylabel('#SOZ nodes/ # Channels')
    xlim([0.5, Inf])
    ctr= ctr+1;
end
%% Group Level Trends Ictal preictal, and null model

ctr=1;

metrics={'aveCtrl', 'modalCtrl', 'strength'};
for type={'i_ict', 'i_preict'} %, 'i_null'
    figure(ctr); clf;
    for i=1:nMeas
        subplot(1,3,i);
        hold on
        h=boxplot(eval(['glob',metrics{i},'(', type{1},',:)']),{'State 1', 'State 2', 'State 3'}, 'Colors', cols(1:3,:));
        set(h,{'linew'},{2})
        ylabel(metrics{i})
        title(metrics{i})
        cs=eval(['c_',type{1},'_glob',metrics{i}]);
        sigstar(num2cell(cs(cs(:,6)<=0.05,[1:2]),2));
    end
    suptitle([type{1}(3:end), 'al'])
    ctr= ctr+1;
end


%% Community Curves %%

for i_set=1:nSets
    
    p=Partitions(i_set);

if strcmp(p.type, 'preictal')==1
    color=cols(1,:);
else
    color= cols(2,:);
end

figure(1); clf; hold on
plot(p.gamma, p.quantileQ(:,3), 'color', color)
scatter(p.gamma(p.modInfInd), p.quantileQ(p.modInfInd,3))


title({'Louvain maximized modularity and inflection point',sprintf( '%s, %s phase #%s', p.ID, p.type, p.block)})
xlabel('gamma (\gamma)')
ylabel('modularity (Q)')

figure(2); clf; hold on
plot(p.gamma, p.quantileCommNum(:,3), 'color', color)
infInd1= p.modInfInd; %find(p.quantileCommNum(:,3)>=3,1,'first');
scatter(p.gamma(infInd1), p.quantileCommNum(infInd1,3),...
    80, cols(1), 'HandleVisibility','off')
text(p.gamma(infInd1)-1e-2, p.quantileCommNum(infInd1,3),...
    sprintf('%d', p.quantileCommNum(infInd1,3)))

try
infInd2=p.modInfInd; %find(p.quantileCommNum(:,3)>=4,1,'first');
scatter(p.gamma(infInd2), p.quantileCommNum(infInd2,4),...
    80, cols(1), 'HandleVisibility','off', 'color', 'red')
text(p.gamma(infInd2)-1e-2, p.quantileCommNum(infInd2,4),...
    sprintf('%d', p.quantileCommNum(infInd2,4)))
catch ME
    disp(ME)
end

title({'NumComs',sprintf( '%s, %s phase #%s', p.ID, p.type, p.block)})
xlabel('gamma (\gamma)')
ylabel('numComms')

%saveas(gcf,sprintf('Figures/Mods/%s_%s.png',p.ID ,p.block))

i_set
ncomms=unique(p.quantileCommNum(infInd2,3))
pause

end

%% Connection Density, distance correlation

[Density, RPosNeg_icov, RPosNeg_pcm, ...
    Corr_icov, Corr_pcm]=deal(zeros(nSets,2));

f_mnstd=@(x)[mean(x), std(x)]; 

for i_set=1:nSets
    config_icov= Networks(i_set).config_icov;
    config_pcm= Networks(i_set).config_pcm;
    pcm= Networks(i_set).pcm;
    
    % Compute matrix density
    density= sum(config_icov~=0)./size(config_icov,1);
    Density(i_set,:)=f_mnstd(density);

    % Compute ratio of pos to neg connections
    rPosNeg_icov= sum(config_icov>0)./sum(config_icov<0);
    rPosNeg_pcm= sum(config_pcm>0)./sum(config_pcm<0);
    
    RPosNeg_icov(i_set,:)=f_mnstd(rPosNeg_icov);
    RPosNeg_pcm(i_set,:)=f_mnstd(rPosNeg_pcm);
    
    % Compute distance correlation
    subj= subject(strcmp({subject.ID}, Networks(i_set).ID));
    triu_i=find(triu(ones(size(pcm,1)),1)~=0);
    dist=pdist2(subj.gridCoords(1:end,:),subj.gridCoords(:,1:end));
    dist=dist(triu_i); 
    try
    [corr_icov,pcorr_icov]=corr(dist(~isnan(dist)), config_icov(~isnan(dist),:));
    [corr_pcm,pcorr_pcm]=corr(dist(~isnan(dist)), config_pcm(~isnan(dist),:));
    catch ME
        disp(ME)
    end
    
    %Fisher transform to calculate mean, then transform back
    Corr_icov(i_set,:)=[tanh(mean(atanh(corr_icov))), sum(pcorr_icov<(0.05/nSets))];
    Corr_pcm(i_set,:)=[tanh(mean(atanh(corr_pcm))), sum(pcorr_pcm<(0.05/nSets))];
    
end

NetworkTable=table({Networks.ID}', {Networks.type}', {Networks.block}', Density,...
    RPosNeg_icov, RPosNeg_pcm, Corr_icov, Corr_pcm ,...
    'VariableNames', {'ID', 'type', 'block', 'Density', 'RposNeg_icov',...
    'RposNeg_pcm', 'corricov', 'corrpcm'});
writetable(NetworkTable,'Data/networkTable.csv')





