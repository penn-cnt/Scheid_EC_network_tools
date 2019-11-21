%V3_analysis

cd '/Users/bscheid/Documents/LittLab/PROJECTS/p01_EC_controllability/v3/Code'

Null='';
datafold='Data';
load(sprintf('%s/dataSets_clean.mat', datafold));
load(sprintf('%s/subjects.mat', 'DataV3.2'));
load(sprintf('%s/%sNetworks.mat',datafold, Null)); 
load(sprintf('%s/%sPartitions.mat', datafold, Null));
load(sprintf('%s/%sMetric_matrices.mat', datafold, Null)); 
load(sprintf('%s/%sState_metrics.mat', datafold, Null))
load('%s/avgDist.mat');

eval(['Networks=', Null, 'Networks'])
eval(['Partitions=', Null, 'Partitions'])
eval(['Metric_matrices=', Null, 'Metric_matrices'])
eval(['State_metrics=', Null, 'State_metrics'])

% if strcmp(Null,'Null')
%     load('Data/nullData.mat')
%     dataSets_clean=nullData;
% end

addpath(genpath('~/Documents/CODE/'))

% cols=[[75,184,166];[255,168,231]; [36,67,152];[140,42,195];[121,29,38];[242,224,43];[74,156,85];...
%    [80,80,80]; [255,255,255]]/255;
cols=[[227,187,187]; [190,8,4]; [138,4,4];[140,42,195];[75,184,166];[242,224,43];[74,156,85];...
   [80,80,80]; [255,255,255]]/255;

i_soz=true(1,39); i_soz(16)=false; i_soz= find(i_soz);

i_ict=find(strcmp({Partitions.type},'ictal'));
i_preict=find(strcmp({Partitions.type},'preictal'));
wSim='wSim_0_01'; % default sim network. 

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
title('ictal')
set(gca,'xaxisLocation','bottom', 'yaxislocation', 'left', 'fontsize', 18)

figure(2); clf; hold on
scatterhist(meds_preict(:),lens_preict(:),'Group',labels','Kernel','on','Location','NorthEast',...
    'Direction','out','Color',cols(1:3,:),'LineStyle',{'-','-','-'},...
    'LineWidth',[2,2,2],'Marker','ddd','MarkerSize',[5,5,5]);
ylabel(''); xlabel('');
title('preictal')
set(gca,'xaxisLocation','bottom', 'yaxislocation', 'left','fontsize', 18)

%% Show table of mean and std statistics for ictal and preictal states
mn_len=mean(lens_ict); std_len=std(lens_ict);
mn_med=mean(meds_ict);std_med=std(meds_ict);
mn_lenpi=mean(lens_preict);std_lenpi=std(lens_preict);
mn_medpi=mean(meds_preict);std_medpi=std(meds_preict);

avg_sz_length=mean(cellfun(@(x)size(x,3), {Networks(i_ict).pcm}))
std_sz_length=std(cellfun(@(x)size(x,3), {Networks(i_ict).pcm}))

avg_sz_N=mean(cellfun(@(x)size(x,2), {Networks(i_ict).pcm}))
std_sz_N=std(cellfun(@(x)size(x,2), {Networks(i_ict).pcm}))

% Calculate % coverage of top 3 contiguous communities
ictal_contig_coverage=arrayfun(@(x)sum(Partitions(x).runLen), [1:nIct]); 
preictal_contig_coverage=arrayfun(@(x)sum(Partitions(x).runLen), [nIct+1:nSets]);

[mean(ictal_contig_coverage), std(ictal_contig_coverage)]
[mean(preictal_contig_coverage), std(preictal_contig_coverage)]

stateTable=table([mn_len', std_len'], [mn_med', std_med'], [mn_lenpi', std_lenpi'], [mn_medpi', std_medpi'],...
    'VariableNames', {'ictal_longest_run', 'ictal_med_value','preictal_longest_run','preictal_med_value'}, ...
    'RowNames', {'st1', 'st2', 'st3'})
sum(table2array(stateTable(:,:)))

writetable(stateTable,'Data/stateTable.csv')

clear mn_len mn_med mn_lenpi mn_medpi labels x_ict y_ict x_preict y_preict meds_preict lens_preict ...
    meds_ict lens_ict



            %% Metric Trends (individual) %%

% Select metrics to plot
n=2; m=4;
all=false; % show both preictal and ictal
metrics={'aveCtrl', 'modalCtrl', 'tModalCtrl'}; 
%metrics={'degree', 'aveCtrl', 'modalCtrl'}; 

for i_set=1:nSets
    i_set
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
    suptitle(sprintf('Network Metrics for %s, %s %d', p.ID, p.type, p.block))
 
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
    Fs=round( dataSets_clean(i_set).Fs);
    if strcmp(dataSets_clean(i_set).type, 'ictal')
        stDiff=round(dataSets_clean(i_set).UEOStart-dataSets_clean(i_set).EECStart);
        data= dataSets_clean(i_set).data(:,1+Fs*stDiff:end)';
    else
        stDiff=round(dataSets_clean(i_set-39).UEOStart-dataSets_clean(i_set-39).EECStart);
        data= [dataSets_clean(i_set).data(:,1+2*(Fs*stDiff):end), dataSets_clean(i_set-39).data(:,1:Fs*stDiff)]';
    end
    sozGrid=find(dataSets_clean(i_set).sozGrid);
    plot((0:length(data)-1)/Fs,data+(1:1000:1000*(N)));
    if ~isempty(sozGrid)
        plot((0:length(data)-1)/Fs,data(:, sozGrid)+1000*(sozGrid-1)', 'c');
    end
    stem([(diff(st)~=0)]*1000*N,'Marker', 'none', 'lineWidth', 2, 'color', 'red')
    title('Signal'); axis tight

    % Config matrix
    subplot(n,m,(n*m)-1); hold on
    imagesc(Networks(i_set).config_pcm)
    set(gca,'colorscale','log')
    stem((diff(st)~=0)*(N*(N-2))/2,'Marker', 'none', 'lineWidth', 2, 'color', 'red')
    title('regular. PCM'); axis tight
    
    % Similarity Matrix
    subplot(n,m,(n*m)); hold on
    imagesc(Networks(i_set).(wSim))
    set(gca,'colorscale','log')
    stem((diff(st)~=0)*T,'Marker', 'none', 'lineWidth', 2, 'color', 'red')
    title('Sim'); axis tight 

    figure(4)
    %set(gca, 'Position', get(gca, 'Position')+[0,.05,0,0])
    imagesc(st) 
    colormap(gca, cols(unique(st),:));
    set(gca, 'YTick', [], 'fontsize', 18)

    pause
end

%% Show SOZ in data
for i_set=1:nSets %i_ict
    i_set
    st= Partitions(i_set).states;
        Fs=round(dataSets_clean(i_set).Fs);
        d=dataSets_clean(i_set);
    if strcmp(dataSets_clean(i_set).type, 'ictal')
        stDiff=round(dataSets_clean(i_set).UEOStart-dataSets_clean(i_set).EECStart);
        data= dataSets_clean(i_set).data(:,1+Fs*stDiff:end)';
    else
        stDiff=round(dataSets_clean(i_set-39).UEOStart-dataSets_clean(i_set-39).EECStart);
        data= [dataSets_clean(i_set).data(:,1+2*(Fs*stDiff):end), dataSets_clean(i_set-39).data(:,1:Fs*stDiff)]';
    end
    N=size(data,2);
    figure(9); clf; hold on
    sozGrid=find(dataSets_clean(i_set).sozGrid);
    plot((0:length(data)-1)/Fs,data+(1:1000:1000*(N)));
    if ~isempty(sozGrid)
        sum(sozGrid)
        plot((0:length(data)-1)/Fs,data(:, sozGrid)+1000*(sozGrid-1)', 'c');
    end
    stem([(diff(st)~=0)]*1000*N,'Marker', 'none', 'lineWidth', 2, 'color', 'red')
    title(sprintf('HUP %s, %d', d.ID, d.block)); axis tight
    
    figure(10); 
    imagesc(Partitions(i_set).contigStates); axis off
    pause
    
    
end

%% Display Kurtosis
for i_set=i_ict
    figure(1)
    clf;
    p= Partitions(i_set);
    st= p.states; mm= Metric_matrices(i_set);
    hold on
    plot(mm.kurtosis)
    stem((diff(st)~=0)*max(mm.kurtosis),'Marker', 'none', 'lineWidth', 2, 'color', 'red')
    ylim([0,20])
   
     figure(2)
    %set(gca, 'Position', get(gca, 'Position')+[0,.05,0,0])
    imagesc(st) 
    colormap(gca, cols(1:3,:));
    set(gca, 'YTick', [], 'fontsize', 18)
    
        figure(3)
    clf;
    p= Partitions(i_set+39);
    st= p.states; mm= Metric_matrices(i_set+39);
    hold on
    plot(mm.kurtosis)
    stem((diff(st)~=0)*max(mm.kurtosis),'Marker', 'none', 'lineWidth', .5, 'color', 'red')
    ylim([0,20])
   
     figure(5)
    %set(gca, 'Position', get(gca, 'Position')+[0,.05,0,0])
    imagesc(st) 
    colormap(gca, cols(1:3,:));
    set(gca, 'YTick', [], 'fontsize', 18)
    
    pause
end

%% Is there a difference between states? (creates glob,c_i_preict_glob, i_ict_glob)

ctype='bonferroni';
display= 'off';

analysis=struct();
diffs=struct();
glob=struct(); c_i_preict_glob=struct(); c_i_ict_glob=struct();
metrics={'strength', 'aveCtrl', 'modalCtrl', 'tModalCtrl','pModalCtrl'}; %'optEnergySOZ'}; %'strength', 'clustering3', 'optEnergy', 'kurtosis', 'skewness'};

for i_set=1:nSets
    
    s= State_metrics(i_set);
    if isempty(s.aveCtrl)
        continue
    end
    
    for i=1:length(metrics)
        
        % Skip paatients without values for a specific metric
        if isempty(s.(metrics{i}))
            eval(sprintf('glob.%s(i_set,:)=nan;', metrics{i})); 
            continue
        end
        % Take care of metrics without nodal values
        if size(s.(metrics{i}),1)==1
            i_set
            eval(sprintf('glob.%s(i_set,:)=s.%s(1:3);', metrics{i}, metrics{i})); 
            continue
        end
            
        % Perform Friedman's analysis & post hoc (non parametric rank test)
        eval(sprintf('[~,~,stats_%s]=friedman(s.%s(:,1:3),1,display);', metrics{i}, metrics{i}));
        eval(sprintf('[c_%s,m_%s]= multcompare(stats_%s, ''ctype'', ctype, ''display'', display);', ...
            metrics{i}, metrics{i}, metrics{i}));
        
        % Save P value, difference, and global average
        eval(sprintf('analysis(i_set).%s=c_%s(:,6);',metrics{i}, metrics{i}));
        eval(sprintf('diffs%s(i_set,:)=c_%s(:,4);', metrics{i}, metrics{i}));
        % Save average of ZSCORED data
        eval(sprintf('glob.%s(i_set,:)=mean(s.%sZ(:,1:3));', metrics{i}, metrics{i})); 
        
    end
end

display='off';
%metrics={'optEnergy'}
% Get group level averages
for i=1:length(metrics)
    if contains(metrics{i}, 'SOZ')
        %SOZ
        % ictal 
        eval(sprintf('[~,table1,stats_g%s]=friedman(glob.%s(i_soz,1:3),1,display);', metrics{i},  metrics{i}))
        eval(['[c_i_ict_glob.',metrics{i},',m_ict_g',metrics{i},...
            ']= multcompare(stats_g',metrics{i},', ''ctype'', ctype, ''display'', display);'])
        disp(table1)
        %preictal
        eval(sprintf('[~,~,stats_g%s]=friedman(glob.%s(i_soz+39,1:3),1,display);', metrics{i},  metrics{i}));
        eval(['[c_i_preict_glob.',metrics{i},',m_ict_g',metrics{i},...
            ']= multcompare(stats_g',metrics{i},', ''ctype'', ctype, ''display'', display);'])
        continue 
        pause
    end
    
    % ictal 
    eval(sprintf('[~,table1,stats_g%s]=friedman(glob.%s(i_ict,1:3),1,display);', metrics{i},  metrics{i}))
    eval(['[c_i_ict_glob.',metrics{i},',m_ict_g',metrics{i},...
        ']= multcompare(stats_g',metrics{i},', ''ctype'', ctype, ''display'', display);'])
    disp(table1)
    %preictal
    eval(sprintf('[~,~,stats_g%s]=friedman(glob.%s(i_preict,1:3),1,display);', metrics{i},  metrics{i}));
    eval(['[c_i_preict_glob.',metrics{i},',m_ict_g',metrics{i},...
        ']= multcompare(stats_g',metrics{i},', ''ctype'', ctype, ''display'', display);'])
    

end
    
disp('done')
%% Display results of Friedman's test individually

ctr=1;
alpha=.05; 

metrics= {'aveCtrl', 'modalCtrl', 'strength'} % optEnergy'


for type=[i_ict', i_preict']
    figure(ctr)
    nMeas=length(metrics);   % number of measures
    nSamp=length([analysis(type).modalCtrl]);
    
    % Get all metrics and 
    met=[];
    diffs=[];
    for m=1:length(metrics)
        ma=reshape([analysis(type).(metrics{m})],3, nSamp)'<alpha;
        dma=diffsmodalCtrl(type,:).*ma;
        met=[met, ma];
        diffs=[diffs, dma]; 
    end

    clf; hold on;
    colormap(gca,cols([9,1,2,3],:));
    imagesc(met.*repmat([1:3],nSamp,nMeas))

    plot((diffs>0).*repmat((1:nMeas*3),nSamp,1),repmat((1:nSamp)',1,nMeas*3),'^','color', 'cyan')
    plot((diffs<0).*repmat((1:nMeas*3),nSamp,1),repmat((1:nSamp)',1,nMeas*3),'v','color', 'white')

    yticks([1:nSamp])
    %yticklabels({Metric_matrices(type).block});
    blk=cellfun(@num2str, {Metric_matrices(type).block}', 'UniformOutput', false);
    yticklabels(strcat({Metric_matrices(type).ID}', {' '}, blk));
    xticks([2:3:nMeas*3-1])
    xticklabels(metrics)
    set(gca,'TickLabelInterpreter','none')
    stem(3*[1:nMeas]+.5,ones(1,nMeas)*size(met,1),...
        'Marker', 'none', 'lineWidth', 0.5, 'color', 'black')
    title({'Significant effect of seizure state on metric'})
    if ctr==1
        suptitle('ictal')
    else; suptitle('preictal')
    end
    %'1-2, 1-3, 2-3'
%     axis tight
%     yyaxis right
%     ylim([1,nSamp])
%     yticks((1:nSamp))
%     yticklabels({numSOZ(sz_pairs(:,1))'./numChannels(sz_pairs(:,1))'*100});
%     yticklabels(sprintf('%d/%d\n', [numSOZ(i_ict)',numChannels(i_ict)']'))
%     ylabel('#SOZ nodes/ # Channels')
    xlim([0.5, Inf])
    ctr= ctr+1
end
disp('done')
%% Group Level Trends Ictal preictal, and null model

% 
% Boxplot visualization of group level trends using Friedman's test. 

fig_ctr=6;
alpha=0.05;

metrics= {'aveCtrl', 'modalCtrl', 'tModalCtrl', 'pModalCtrl'};
%metrics={'optEnergySOZ'}
for type={'i_ict', 'i_preict'} %, 'i_null'
    figure(fig_ctr); clf;
    for i=1:length(metrics)
        subplot(1,length(metrics),i);
        hold on
        outs=sort(eval(['glob.',metrics{i}, '(isoutlier(glob.',metrics{i},'))']));
        sorted=sort(eval(['glob.', metrics{i}, '(:)']));
        
        % Set datalimits if outliers are > 2 std away
        maxLim=outs(find(zscore(outs)>2, 1, 'first')); 
        minLim=outs(find(zscore(outs)<-2, 1, 'last'));
        
        if isempty(minLim); minLim= -Inf;
        else; minLim=sorted(find(sorted>minLim, 1, 'first')); end
        if isempty(maxLim); maxLim= Inf; 
        else; maxLim=sorted(find(sorted<maxLim, 1, 'last')); end
        
        h=boxplot(eval(['glob.',metrics{i},'(', type{1},',:)']), ...
            {'Phase 1', 'Phase 2', 'Phase 3'}, 'Colors', cols(1:3,:), ...
            'Symbol', 'x', ...
            'dataLim', [minLim, maxLim]);
        set(h,{'linew'},{2})
        ylabel([metrics{i}, ' zscore'])
        title(metrics{i})
        cs=eval(['c_',type{1},'_glob.',metrics{i}]);
        sigstar(num2cell(cs(cs(:,6)<=alpha,[1:2]),2));
        %ylim([min(eval(['glob.',metrics{i},'(:)']))*1.1, max(eval(['glob.',metrics{i},'(:)'])*1.2)])
    end
    suptitle([type{1}(3:end), 'al'])
    fig_ctr= fig_ctr+2;
end

%% Fit a Linear Mixed model to metrics

% Get table with (ID, Phase, Metric)
metrics= {'aveCtrl', 'modalCtrl','optEnergy', 'strength'};
for m= metrics
    varnames= {m{1},'Phase', 'ID'};
    tbl=table([],[],[],'VariableNames', varnames); 
    for i=i_ict
        tbl= [tbl; table(mean(State_metrics(i).(m{1})(:,1:3))', [1:3]', repmat(string(Metric_matrices(i).ID),3,1), ...
            'VariableNames', varnames)];
    end

    tbl.Phase = categorical(tbl.Phase); 
    tbl.ID = categorical(tbl.ID);
    
lme = fitlme(tbl, sprintf('%s ~ 1 + Phase + (1|ID)', m{1}))

end


%% Is there a correlation between metrics?
clf
metric1='optEnergysoz'; %aveCtrl', 'modalCtrl', 'optEnergy', 'tModalCtrl','pModalCtrl', 'strength'};
metric2='strength';

sozCorr=Inf(39, 3);
sozCorrP=Inf(39,3);
figure(1)
for i_set=i_ict
    i_set
    met1=State_metrics(i_set).(metric1); %avgDist{i_set}
    met2=State_metrics(i_set).(metric2);
%     if isempty(avgDist{i_set})
%         disp('skipping')
%         continue
%     end
    
    soz=dataSets_clean(i_set).sozGrid;
    
for s=1:3
    subplot(1,3,s)
    hold on
    % Get correlation between soz nodes and others per state
    [cor, pval]=corr(met1(:,s), met2(:,s));
    sozCorr(i_set, s)=cor;
    sozCorrP(i_set, s)=pval;
    
    scatter(met1(:,s), met2(:,s));
    scatter(met1(soz,s),met2(soz,s), 'r');
    mdl = polyfit(met1(:,s), met2(:,s), 1);
    f=polyval(mdl, met1(:,s)); 
    plot(met1(:,s),f); 
 
    title(sprintf('pval=%0.3f', pval))
    %title(sprintf('pval=%0.3f, R^2=%0.3f', mdl.Coefficients.pValue(2), mdl.Rsquared.Adjusted))
    xlabel('optEnergy'); ylabel('modalCtrl')
end

suptitle(sprintf('Avg. Ctrl. SOZ distance corr, %s %d', State_metrics(i_set).ID, State_metrics(i_set).block))
%pause
end


figure(2)
clf
hold on
imagesc(sozCorr(i_ict, :).^2)
plot([1:3], (sozCorrP(i_ict,:)<=0.05).*repmat([1:length(i_ict)]',1, 3), '*', 'color', 'black');
blk=cellfun(@num2str, {Metric_matrices(i_ict).block}', 'UniformOutput', false);
yticklabels(strcat({Metric_matrices(i_ict).ID}', {' '}, blk));
yticks([1:39])
set(gca,'TickLabelInterpreter','none')
xticks([1:3])
xticklabels({'Phase 1', 'Phase 2', 'Phase 3'})
title(sprintf('Corr values between %s and %s', metric1, metric2))
% figure(3)
% imagesc((sozCorrP<0.05).*(sozCorr<0))

%% View Nodes and Values spatially  %%
figure(1)
clf;
metrics={'optEnergySOZ'};
cmap=colormap; 
for i_set=i_soz %1:39
    figure(1)
    clf; hold on
    if isempty(dataSets_clean(i_set).sozGrid)
        continue
    end
    for m=1:length(metrics)
        %st=State_metrics(i_set).strengthPos+State_metrics(i_set).strengthNeg;
        st=State_metrics(i_set).(metrics{m});
        %st=abs(State_metrics(i_set).(metrics{m})-State_metrics(i_set+39).(metrics{m})(:,1));
        c_met=1+reshape(round(63*rescale(st(:))), size(st));
        clims=[min(st(:)), max(st(:))];
        X=dataSets_clean(i_set).gridCoords(:,1);
        Y=dataSets_clean(i_set).gridCoords(:,2);
        
        soz=dataSets_clean(i_set).sozGrid; 
        ctrl=diag(getSpreadControl([X,Y],soz)); 
        
        for s=3:-1:1
            subplot(length(metrics),3,((m-1)*3+s))
            colormap parula
            hold on
           
            % Show control gradient
            scatter(X, Y, 500*ctrl, 'b', 'filled')
            
            % Show metric
            %scatter(X,Y,c_met(:,s)+80, st(:,s), 'filled')
            
              % Show control metric value
%             scatter(X(EnergyMetric(i_set).aveCtrl_NOI(:,s)), ...
%                  Y(EnergyMetric(i_set).aveCtrl_NOI(:,s)),100, 'r')

              % Show soz
              scatter(X(soz), Y(soz),100, 'r')
              

              
              
            xticklabels([]); yticklabels([])
            xlim([.5,8.5]); ylim([.5,8.5]);
            xlabel(sprintf('phase %d', s))
            caxis(clims)
            ylabel(metrics{m})
            %if s==3; colorbar; end
        end
        
    end
    suptitle(sprintf('Spatial Metrics for %s, block %d', dataSets_clean(i_set).ID, dataSets_clean(i_set).block))
    pause
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
    subj= subjects(strcmp({subjects.ID}, Networks(i_set).ID));
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

%%
figure(5)
subplot(121)
histogram(Density(1:39,1),10)
title('Average ictal PCM density')
subplot(122)
histogram(Density(40:end,1),10)
title('Average preictal PCM density')

figure(6)
subplot(121)
histogram(RPosNeg_pcm(1:39,1),10)
title('ictal')
subplot(122)
histogram(RPosNeg_pcm(40:end,1),10)
title('preictal')
suptitle('average positive/negative edge ratio')


%% Permutation test Method 1: See if SOZ nodes are relevant

percentRank = @(YourArray, TheProbes) reshape( mean( bsxfun(@le,...
    YourArray(:), TheProbes(:).') ) * 100, size(TheProbes) );

metrics={'strength','aveCtrl', 'modalCtrl', 'tModalCtrl','pModalCtrl'}
nperms=2000;  %number of 
 
sozRanks=struct(); 
disp('Calculating sozRanks...')

for i_set=i_ict
    i_set
    % Grab SOZ nodes from pre ict and ictal
    if isempty(dataSets_clean(i_set).sozGrid); continue; end
    soz=dataSets_clean(i_set).sozGrid; 
    if sum(soz)==0; continue; end
    
    % Get nullset
    nsoz=sum(soz); 
    nullset=find(~soz); 
    i_perm=zeros(nperms, nsoz);
    for i=1:nperms
        i_perm(i,:)=nullset(randperm(sum(~soz),nsoz));
    end

    d=dataSets_clean(i_set); 
    st= Partitions(i_set).states;
    %figure(1); clf; hold on
        Fs=round( dataSets_clean(i_set).Fs);
    if strcmp(dataSets_clean(i_set).type, 'ictal')
        stDiff=round(dataSets_clean(i_set).UEOStart-dataSets_clean(i_set).EECStart);
        data= dataSets_clean(i_set).data(:,1+Fs*stDiff:end)';
    else
        stDiff=round(dataSets_clean(i_set-39).UEOStart-dataSets_clean(i_set-39).EECStart);
        data= [dataSets_clean(i_set).data(:,1+2*(Fs*stDiff):end), dataSets_clean(i_set-39).data(:,1:Fs*stDiff)]';
    end
    N=size(data,2);
    
%     figure(1); clf; hold on
      sozGrid=find(dataSets_clean(i_set).sozGrid);
%     plot((0:length(data)-1)/Fs,data+(1:1000:1000*(N)));
%     if ~isempty(sozGrid)
%         plot((0:length(data)-1)/Fs,data(:, sozGrid)+1000*(sozGrid-1)', 'c');
%     end
%     stem([(diff(st)~=0)]*1000*N,'Marker', 'none', 'lineWidth', 2, 'color', 'red')
%     title('Signal'); axis tight

    for m=metrics
%         figure(3); clf;
        for s=1:3
        % Difference between ictal and preictal phase
            %diffs=State_metrics(i_set).(m{1})(:,s)-State_metrics(i_set+39).(m{1})(:,1);
            diffs=State_metrics(i_set).(m{1})(:,s); 
            soz_diffs= diffs(soz);
            dist_diffs=diffs(i_perm)';

            % Get absolute value of sum
            soz_diffs=mean(abs(soz_diffs),1); 
            dist_diffs=mean(abs(dist_diffs),1); 
                        
%             subplot(1,3,s)
%             hold on
%             histogram(dist_diffs)
%             stem(quantile(dist_diffs, [.05,.25,.5,.75,.95]), ones(1,5)*100, 'linewidth', 2)
%             stem(soz_diffs, 100, 'linewidth', 2)
              sozRanks(i_set).(m{1})(s)=percentRank(dist_diffs, soz_diffs);
        end   
%         suptitle(sprintf('%s, HUP %s, siezure %d', m{1}, d.ID,d.block))
        %pause
        
    end

end
disp('done with sozRanks calculation')

%% Display permutation test result
a=[1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    25    26    28];

% for each metric
metrics={'aveCtrl', 'modalCtrl', 'tModalCtrl','pModalCtrl', 'strength'}
figure(2); clf; 
figure(3); clf;
for m=1:length(metrics)
    figure(2)
    subplot(1,length(metrics), m)
    pcnts=[sozRanks.(metrics{m})];
    sumSig=sum(pcnts>95)+sum(pcnts<5);
    scatter(repmat([1:3],1,18)+normrnd(0,.05, [1,54]),  pcnts, 60, repmat(parula(18),3,1), 'filled')
    title(sprintf('%s', metrics{m}))
    xlim([.5,3.5])
    ylim([-5, 105])
    xticklabels({'Phase1', 'Phase 2', 'Phase 3'})
    xtickangle(20)
    hline(2.5)
    hline(97.5)
    
    figure(3)
    subplot(1,length(metrics), m)
    imagesc(double(reshape([sozRanks.(metrics{m})],3,18)'>=97.5)...
        +2*double(reshape([sozRanks.(metrics{m})],3,18)'<=2.5))
    yticks([1:18])
    blk=cellfun(@num2str, {Metric_matrices(a).block}', 'UniformOutput', false);
    yticklabels(strcat({Metric_matrices(a).ID}', {' '}, blk));
end

figure(2)
suptitle('Preictal vs. phase, SOZ permutation test results')



%% Permutation test Method 2: See which nodes are on the edges of the distribution

percentRank = @(YourArray, TheProbes) reshape( mean( bsxfun(@le,...
    YourArray(:), TheProbes(:).') ) * 100, size(TheProbes) );

metrics={'aveCtrl', 'modalCtrl', 'tModalCtrl','pModalCtrl', 'strength', 'clustering3'};
nperms=2000;  %number of 
 
kRanks=struct(); 
disp('Calculating distributions...')

for i_set=i_ict
    i_set

    d=dataSets_clean(i_set);
    st= Partitions(i_set).states;
    Fs=round(d.Fs);
    N=size(d.data,1); 
    
    % Get nullset
    k=round(N*.1); 
    i_perm=zeros(nperms, k);
    rng(4)
    for i=1:nperms
        i_perm(i,:)=randperm(N,k);
    end
    
    [all_freq, all_edges]=histcounts(i_perm(:),[0:64]+.5);
%     figure(1); clf; hold on
%     if strcmp(dataSets_clean(i_set).type, 'ictal')
%         stDiff=round(dataSets_clean(i_set).UEOStart-dataSets_clean(i_set).EECStart);
%         data= dataSets_clean(i_set).data(:,1+Fs*stDiff:end)';
%     else
%         stDiff=round(dataSets_clean(i_set-39).UEOStart-dataSets_clean(i_set-39).EECStart);
%         data= [dataSets_clean(i_set).data(:,1+2*(Fs*stDiff):end), dataSets_clean(i_set-39).data(:,1:Fs*stDiff)]';
%     end
% 
%     figure(1); clf; hold on
%     plot((0:length(data)-1)/Fs,data+(1:1000:1000*(N)));
%     if ~isempty(sozGrid)
%         plot((0:length(data)-1)/Fs,data(:, sozGrid)+1000*(sozGrid-1)', 'c');
%     end
%     stem([(diff(st)~=0)]*1000*N,'Marker', 'none', 'lineWidth', 2, 'color', 'red')
%     title('Signal'); axis tight

    for m=metrics
        for s=1:3
        % Difference between ictal and preictal phase
            diffs=State_metrics(i_set).(m{1})(:,s)-State_metrics(i_set+39).(m{1})(:,1);
            dist_diffs=diffs(i_perm)';

            % Get absolute value of sum
            dist_diffs=mean(abs(dist_diffs),1);
            
            % Get ordering/rank
            [~,rnk]=sort(dist_diffs); 
            kRanks(i_set).(m{1})(s,:)= rnk/nperms; 
            figure(1)
            clf
            hold on
            % Which nodes appear the most in the bottom percentile?
            sml=(rnk/nperms)<0.025; lg=(rnk/nperms)>0.975; 
            [sml_freq]=histcounts(i_perm(sml,:),[0:64]+.5);
            
            bar(sml_freq./all_freq)
            [lg_freq]=histcounts(i_perm(lg,:),[0:64]+.5);
            
            b2=bar(lg_freq./all_freq);
            b2.FaceAlpha = 0.5;
            
        end   
    end
    
end
disp('done with sozRanks calculation')

%% Generate Figures
%close all

% PCM figure
figure(1); set(gcf,'Name', 'pcm')
imagesc(Networks(1).pcm(:,:,6))
colormap('winter'); colorbar
xlabel('nodes'); ylabel('nodes');

% methods metric
figure(3); clf; set(gcf,'Name', 'methodsMetric')
subplot(10,1,[1:8])
imagesc(Metric_matrices(27).aveCtrl)
xlabel('nodes'); ylabel('nodes');
colormap('winter'); colorbar
subplot(10,1,10)
imagesc(Partitions(27).states) 
colormap(gca, cols(1:3,:));
xticklabels(''); yticklabels('');




