%% V3_analysis  data

run('Code/initProject.m')

% i_soz=true(1,39); i_soz(16)=false; i_soz= find(i_soz);


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

tst_id=table({datestr(now())}, "StateStats", 'VariableNames', {'Date', 'test'});
writetable(tst_id, 'Data/resultsTable.xlsx', 'Range', 'A1:B2');
writetable(stateTable,'Data/resultsTable.xlsx', 'Range', 'A3:H6')

clear mn_len mn_med mn_lenpi mn_medpi labels x_ict y_ict x_preict y_preict meds_preict lens_preict ...
    meds_ict lens_ict

%% Friedmans- Is there a difference between states? (creates glob,c_i_preict_glob, i_ict_glob)

i_ict=find(strcmp({Partitions.type},'ictal'));
i_preict=find(strcmp({Partitions.type},'preictal'));

ctype='bonferroni';
display= 'off'; %DON'T TOUCH!

groupOn= false;

analysis=struct();
diffs=struct();
rnks=struct();
glob=struct(); c_i_preict_glob=struct(); c_i_ict_glob=struct();
metrics={'aveCtrl', 'modalCtrl', 'tModalCtrl','pModalCtrl', 'strength', 'optEnergy'}; %'strength', 'clustering3', 'optEnergy', 'kurtosis', 'skewness'};

lstID=State_metrics(1).ID; ctr=1; 
for i_set=1:nSets
    
    s= State_metrics(i_set);
    if isempty(s.aveCtrl)
        continue
    end
    
    for i=1:length(metrics)
        
        % Skip patients without values for a specific metric
        if isempty(s.(metrics{i}))
            glob.(metrics{i})(i_set,:)=nan; 
            continue
        end
        % Take care of metrics without nodal values
        if size(s.(metrics{i}),1)==1
            i_set
            glob.(metrics{i})(i_set,:)=s.(metrics{i})(1:3);
            continue
        end
            
        % Perform Friedman's analysis & post hoc (non parametric rank test)
        [~,~,stats.(metrics{i})]=friedman(s.(metrics{i})(:,1:3),1,display);
        multcomp.(metrics{i})= multcompareBS(stats.(metrics{i}), 'ctype', ctype, 'display', display);
        %multcomp.(metrics{i})= postHocSignedRank(s.(metrics{i})(:,1:3), 'ttest');
        
        % Save P value, difference, and global network average
        analysis(i_set).(metrics{i})=multcomp.(metrics{i})(:,6);
        diffs.(metrics{i})(i_set,:)=multcomp.(metrics{i})(:,4);
        
        [~, rnks.(metrics{i})(i_set,:)]=ismember(stats.(metrics{i}).meanranks,sort(stats.(metrics{i}).meanranks));
        % Save average of ZSCORED data
        eval(sprintf('glob.%s(i_set,:)=mean(s.%sZ(:,1:3));', metrics{i}, metrics{i}));
        
    end
end

% Random removal of 5 seizures from study026 to balance data set. 
display='off';
rm=[31    36    33    35    29];
i_ict(ismember(i_ict,rm))=[];
i_preict(ismember(i_preict,rm+39))=[];
subinds=[i_ict,i_preict];


if groupOn
    % Get global avg grouped by sz type
    cls=strcat({dataSets_clean([i_ict, i_preict]).type}', {dataSets_clean([i_ict, i_preict]).ID}', ...
        num2str([dataSets_clean([i_ict, i_preict]).sz_subtype]'));
    [uns, bb]=unique(string(cls),'rows');

    grp_glob=struct();
    ctr=1;
    for u=1:length(uns)
        inds=strcmp(cls, uns(u,:));
        for i=1:length(metrics)
            grp_glob.(metrics{i})(u,:)=mean(glob.(metrics{i})(subinds(inds),1:3),1);
        end
    end

    n_gp=length(grp_glob.(metrics{1}))/2;
    i_ict=[1:n_gp]; i_preict=[1:n_gp]+n_gp;
end

stats_ict=zeros(length(metrics),4);
stats_preict=zeros(length(metrics),4);
xl_ptr=1; if groupOn; xl_ptr= 5*length(metrics)+1; end

%metrics={'optEnergy'}
% Get group level averages
for i=1:length(metrics)
    if contains(metrics{i}, 'SOZ')
        %SOZ
        % ictal 
        [~,table1,stats_g.(metrics{i})]=friedman(glob.(metrics{i})(i_soz,1:3),1,display);
        c_i_ict_glob.(metrics{i})= multcompareBS(stats_g.(metrics{i}),'ctype', ctype, 'display', display);
        %c_i_ict_glob.(metrics{i})= postHocSignedRank(glob.(metrics{i})(i_soz,1:3));
       
        %preictal
        [~,~,stats_g.(metrics{i})]=friedman(glob.(metrics{i})(i_soz+39,1:3),1,display);
        c_i_preict_glob.metrics{i}= multcompareBS(stats_g.metrics{i},'ctype', ctype, 'display', display);
       % c_i_ict_glob.(metrics{i})= postHocSignedRank(glob.(metrics{i})(i_soz+39,1:3));
       
        continue 
        pause
    end
        
    if groupOn
        % Group Independent
        
        % ictal 
        [~,tbl_ict,stats_gict]=friedman(grp_glob.(metrics{i})([1:n_gp],1:3),1,display);
        c_i_ict_glob.(metrics{i})= multcompareBS(stats_gict,'ctype', ctype, 'display', display);
        %c_i_ict_glob.(metrics{i})= postHocSignedRank(grp_glob.(metrics{i})([1:n_gp],1:3));
        
        %preictal
        [~,tbl_preict,stats_gpreict]=friedman(grp_glob.(metrics{i})([1:n_gp]+n_gp,1:3),1,display);
        c_i_preict_glob.(metrics{i})= multcompareBS(stats_gpreict,'ctype', ctype, 'display', display);
        %c_i_preict_glob.(metrics{i})= postHocSignedRank(grp_glob.(metrics{i})([1:n_gp]+n_gp,1:3));
    else
        % ictal 
        [~,tbl_ict,stats_gict]=friedman(glob.(metrics{i})(i_ict,1:3),1,display);
        c_i_ict_glob.(metrics{i})=multcompareBS(stats_gict,'ctype', ctype, 'display', display);  
%       c_i_ict_glob.(metrics{i})= postHocSignedRank(glob.(metrics{i})(i_ict,1:3), 'ttest');
                
        
        %preictal
        [~,tbl_preict,stats_gpreict]=friedman(glob.(metrics{i})(i_preict,1:3),1,display);
        c_i_preict_glob.(metrics{i})= multcompareBS(stats_gpreict,'ctype', ctype, 'display', display);       
%       c_i_preict_glob.(metrics{i})= postHocSignedRank(glob.(metrics{i})(i_preict,1:3), 'ttest');
        
    end
    
    % Write results to table
    writecell({['ictal ', metrics{i}, ' groupOn: ', string(groupOn)]},...
        [datafold, '/resultsTable.xlsx'], 'Sheet', 'Grp Metrics', 'Range', sprintf('L%d', xl_ptr));
    writecell({['preictal ', metrics{i},' groupOn: ', string(groupOn)]},...
        [datafold, '/resultsTable.xlsx'], 'Sheet', 'Grp Metrics', 'Range', sprintf('S%d', xl_ptr)); xl_ptr=xl_ptr+1; 
    writematrix(c_i_ict_glob.(metrics{i}),[datafold, '/resultsTable.xlsx'],'Sheet','Grp Metrics','Range', sprintf('L%d:R%d', xl_ptr, xl_ptr+3));   
    writematrix(c_i_preict_glob.(metrics{i}),[datafold, '/resultsTable.xlsx'],'Sheet','Grp Metrics','Range', sprintf('S%d:Y%d', xl_ptr, xl_ptr+3)); 
    xl_ptr=xl_ptr+4; 
    
    stats_ict(i,:)= [tbl_ict{2,3},tbl_ict{2,5},stats_gict.n, tbl_ict{2,6}];
    stats_preict(i,:)= [tbl_preict{2,3},tbl_preict{2,5},stats_gpreict.n, tbl_preict{2,6}];

end
 
% Show table and write to results document

statTbl_ict= array2table(stats_ict, 'VariableNames', {'df','chi2', 'n', 'pval'})
statTbl_preict= array2table(stats_preict, 'VariableNames', {'df','chi2', 'n', 'pval'})

if groupOn
    writecell({datestr(now()), 'GroupLevel_metrics-grouped by subject'}, [datafold, '/resultsTable.xlsx'], 'Sheet', 'Grp Metrics', 'Range', 'F1:G2');
    writetable([table(metrics', 'VariableNames', {'Ictal_Metrics'}),statTbl_ict],[datafold, '/resultsTable.xlsx'], ...
        'Sheet', 'Grp Metrics', 'Range', sprintf('F3:L%d',height(statTbl_ict)+10))
    writetable([table(metrics', 'VariableNames', {'Preictal_Metrics'}), statTbl_preict],[datafold, '/resultsTable.xlsx'], ...
        'Sheet', 'Grp Metrics', 'Range', sprintf('F%d:L%d',height(statTbl_ict)+4,4+2*height(statTbl_ict)))
else
    writecell({datestr(now()), 'GroupLevel_metrics- sz indep.'}, [datafold, '/resultsTable.xlsx'], 'Sheet', 'Grp Metrics', 'Range', 'A1:B2');
    writetable([table(metrics', 'VariableNames', {'Ictal_Metrics'}), statTbl_ict],[datafold, '/resultsTable.xlsx'], ...
        'Sheet', 'Grp Metrics', 'Range', sprintf('A3:E%d',height(statTbl_ict)+10))
    writetable([table(metrics', 'VariableNames', {'Preictal_Metrics'}), statTbl_preict],[datafold, '/resultsTable.xlsx'], ...
        'Sheet', 'Grp Metrics', 'Range', sprintf('A%d:E%d',height(statTbl_ict)+4,4+2*height(statTbl_ict)))
end

    
disp('done')

%% Display results of Friedman's test individually

ctr=0;
alpha=.017; 

metrics= {'aveCtrl', 'modalCtrl', 'tModalCtrl', 'pModalCtrl', 'strength'} % 'optEnergy'}

for type=[i_ict', i_preict']
    figure(ctr+1)
    nMeas=length(metrics);   % number of measures
    nSamp=length([analysis(type).modalCtrl]);
    
    % Get all metrics and 
    met=[];

    figure(ctr+1); clf;
    for m=1:length(metrics)
        ma=reshape([analysis(type).(metrics{m})],3, nSamp)'<alpha;
        
        drnk=rnks.(metrics{m})(type,:);
        
        i_pp=[find(ma(:,1));find(ma(:,2));find(ma(:,3))]
        pp=[drnk(ma(:,1),[2,1]); drnk(ma(:,2),[3,1]); drnk(ma(:,3),[3,2])]
        ph=[repmat([2,1],sum(ma(:,1)),1);repmat([3,1],sum(ma(:,2)),1);repmat([3,2],sum(ma(:,3)),1)]
        
        plot(3.25*(m-1)+ph', ((2*i_pp)+(pp*.5))', 'black')
                
         figure(ctr+1)
         hold on
         plot(3.25*(m-1)+[(1:3),1], [drnk,drnk(:,1)]*.5+2*[1:nSamp]', 'color', 'black', 'lineWidth', .5)
         X=repmat(3.25*(m-1)+[(1:3),1],nSamp,1)'; Y=([drnk,drnk(:,1)]*.5+2*[1:nSamp]')';
         scatter(X(:), Y(:), 60, repmat(cols([(1:3),1],:),floor(numel(X)/4),1), 'filled')
        
    end
    
    % Plot lines
     figure(ctr+1)
     blk=cellfun(@num2str, {Metric_matrices(type).block}', 'UniformOutput', false);
     xticks((2:3.25:nMeas*3.5-1)); xticklabels(metrics)
     yticks((1:nSamp)*2+1); yticklabels(strcat({Metric_matrices(type).ID}', {' '}, blk));
     xlim([0,nMeas*3.5-.5])

     ctr= ctr+1
    if ctr==1, suptitle('ictal'), else; suptitle('preictal'); end

     

end
disp('done')

%% Group Level Trends Ictal/Preictal (must run Friedman's block first)
% 
% Boxplot visualization of group level trends using Friedman's test. 

fig_ctr=6;
alpha=0.016;

metrics= {'aveCtrl', 'modalCtrl', 'tModalCtrl', 'pModalCtrl'};

%metrics={'optEnergySOZ'}
for type={'i_ict', 'i_preict'} %, 'i_null'
    figure(fig_ctr); clf;
    for i=1:length(metrics)
        if groupOn
            glb= grp_glob.(metrics{i});
        else
            glb= glob.(metrics{i});
        end
        subplot(1,length(metrics),i);
        hold on
        outs= sort(glb(isoutlier(glb)));
        sorted=sort(glb(:));
        
        % Set datalimits if outliers are > 2 std away
        maxLim=outs(find(zscore(outs)>2, 1, 'first')); 
        minLim=outs(find(zscore(outs)<-2, 1, 'last'));
        
        if isempty(minLim); minLim= -Inf;
        else; minLim=sorted(find(sorted>minLim, 1, 'first')); end
        if isempty(maxLim); maxLim= Inf; 
        else; maxLim=sorted(find(sorted<maxLim, 1, 'last')); end
        
        h=boxplot(eval(['glb(', type{1},',:)']), ...
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
    fig_ctr= fig_ctr+2;
    suptitle(sprintf('%sal, alpha: %0.3f grouped: %s', type{1}(3:end),alpha, string(groupOn)))
end

% set(gcf, 'Position', [560   706   387   242])

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

%% Individual Energy Analysis

i_set= 2; 

p=Partitions(i_set);
d=dataSets_clean(i_set);
stDiff=round(d.UEOStart-d.EECStart);
data= d.data(:,1+Fs*stDiff:end);
Fs=round(d.Fs);

xf=Energy(i_set).xf;
bandfun=@(x)bandpower(x', Fs, freq_band)';

for s=1:3
    t1= find(p.contigStates==s, 1, 'first');
    t2= find(p.contigStates==s, 1, 'last');
    signal=data(1:end,(Fs*(t1-1)+1:Fs*t2));
    
    x0=Energy(i_set).x0(:,s);
    
    x0_t=MovingWinFeats(signal,Fs, 1, 1, bandfun);
    
% Look at Energy Trajectory- Distance over time
figure(1); clf
subplot(121); hold on
plot(d.data'+[1:size(d.data,1)]*1000);
vline(d.Fs*(t1-1)+1)
vline(d.Fs*t2)
axis tight
subplot(122)
plot(signal'+[1:size(signal,1)]*1000);
axis tight; 
figure(2)
imagesc(p.states) 
colormap(gca, cols([5,2,6],:));
set(gca, 'YTick', [], 'fontsize', 18)

figure(3)
imagesc(reshape(xf, 8, 8))
title('avg final bandpower (preictal)'), caxis([0,.5])
colorbar

figure(4)
imagesc(reshape(x0, 8, 8))
title(sprintf('avg bandpower ictal phase %d (%s %d)', s, d.ID, d.block))
caxis([0,.5])
colorbar
pause
end


    figure(1); imagesc(U_opt'); title('U'); 
    figure(2);
    subplot(1,6,[2:5]); imagesc([X_opt(:,1:N)'])% imagesc(X_opt(:,2:N)'); title('X');
    subplot(1,6,[2:5]); imagesc([X_opt(1,1:N)',X_opt(end,1:N)'])
    cl=caxis;
    subplot(1,6,1); imagesc(log(x0)); title('x0'); caxis(cl)
    subplot(1,6,6); imagesc(log(xf)); title('xf');  caxis(cl)

%% View Nodes and Values spatially  %%
figure(5)
clf;
metrics={'optEnergy'};
cmap=colormap; 
for i_set=1:34
    i_set
    figure(5)
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
            %scatter(X, Y, 500*ctrl, 'b', 'filled')
            
            % Show metric
            scatter(X,Y,c_met(:,s)+80, st(:,s), 'filled')
            
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

%% Permutation test Method 1: See if SOZ nodes are relevant

percentRank = @(YourArray, TheProbes) reshape( mean( bsxfun(@le,...
    YourArray(:), TheProbes(:).') ) * 100, size(TheProbes) );

metrics={'aveCtrl','modalCtrl', 'tModalCtrl','pModalCtrl', 'optEnergy'}
nperms=1000;  %number of 
 
sozRanks=struct(); 
disp('Calculating sozRanks...')

for i_set=i_ict
    i_set
    % Grab SOZ nodes from ictal phase
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
        % figure(1); clf;
        for s=1:3
        % Difference between ictal and preictal phase
            %diffs=State_metrics(i_set).(m{1})(:,s)-State_metrics(i_set+39).(m{1})(:,1);
            
            alldiffs=State_metrics(i_set).(m{1})(:,s); 
            soz_diffs= alldiffs(soz);
            dist_diffs=alldiffs(i_perm)';

            % Get absolute value of sum
            soz_diffs=mean(abs(soz_diffs),1); 
            dist_diffs=mean(abs(dist_diffs),1); 
            sozRanks(i_set).(m{1})(s)=percentRank(dist_diffs, soz_diffs);
                        
%             subplot(1,3,s)
%             hold on
%             histogram(dist_diffs)
%             stem(quantile(dist_diffs, [.05,.25,.5,.75,.95, .983]), ones(1,6)*100, 'linewidth', 2)
%             stem(soz_diffs, 100, 'linewidth', 2)
              
        end   
         %suptitle(sprintf('%s, HUP %s, siezure %d', m{1}, d.ID,d.block))
        %pause
        
    end

end
disp('done with sozRanks calculation')

%% Display permutation test result

ptSOZ=cellfun(@sum, {dataSets_clean(1:end/2).sozGrid})>0;

alpha=5; %.2/3*100/2;

% for each metric
metrics={'optEnergy'}
figure(2); clf; 
figure(3); clf;
for m=1:length(metrics)
        
    figure(2)
    subplot(1,length(metrics), m)
    pcnts=[sozRanks.(metrics{m})];
    sumSig=sum(pcnts>95)+sum(pcnts<5);
    scatter(repmat([1:3],1,length(pcnts)/3)+normrnd(0,.05, [1,numel(pcnts)]),  pcnts, 60, repmat(parula(length(pcnts)/3),3,1), 'filled')
    title(sprintf('%s', metrics{m}))
    xlim([.5,3.5])
    ylim([-5, 105])
    xticklabels({'Phase1', 'Phase 2', 'Phase 3'})
    xtickangle(20)
    hline(2.5)
    hline(97.5)
    
    figure(3)
    subplot(1,length(metrics), m)
    imagesc(double(reshape([sozRanks.(metrics{m})],3,length(pcnts)/3)'>=(100-alpha))...
        +2*double(reshape([sozRanks.(metrics{m})],3,length(pcnts)/3)'<=alpha))
    yticks([1:length(ptSOZ)])
    blk=cellfun(@num2str, {Metric_matrices(ptSOZ).block}', 'UniformOutput', false);
    yticklabels(strcat({Metric_matrices(ptSOZ).ID}', {' '}, blk));
    title(metrics{m})
    
    sigs=double(reshape([sozRanks.(metrics{m})],3,length(pcnts)/3)'>=(100-alpha))+double(reshape([sozRanks.(metrics{m})],3,length(pcnts)/3)'<=alpha)
    ids=find(ptSOZ)'; in=ids(logical(sigs(:,1))); out=ids(~logical(sigs(:,1)))
    
    [p,~,stats]=ranksum(cellfun(@sum,{dataSets_clean(in).sozGrid}), cellfun(@sum,{dataSets_clean(out).sozGrid}))


end

figure(2)
suptitle('Preictal vs. phase, SOZ permutation test results')

figure(3)
suptitle(sprintf('Green: SOZ metric >= %0.3f perms, Yellow: SOZ metric <=%0.3f perms', 1-alpha/100, alpha/100))

%% Permutation test Method 2: See which nodes are on the edges of the distribution

percentRank = @(YourArray, TheProbes) reshape( mean( bsxfun(@le,...
    YourArray(:), TheProbes(:).') ) * 100, size(TheProbes) );

metrics={'pModalCtrl'}%, %'aveCtrl', 'modalCtrl', 'tModalCtrl','pModalCtrl', 'strength', 'clustering3'};
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
            alldiffs=State_metrics(i_set).(m{1})(:,s)-State_metrics(i_set+39).(m{1})(:,1);
            dist_diffs=alldiffs(i_perm)';

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

    pause
end
disp('done with sozRanks calculation')

%% Strength percentages

pp=cell2mat({State_metrics.PosRatio}')
boxplot(pp, {'Phase 1', 'Phase 2', 'Phase 3'}, 'Colors', cols(1:3,:), 'Symbol', 'x')



%
%% NOT USED: Fit a Linear Mixed model to metrics

% Get table with (ID, Phase, Metric)
metrics= {'aveCtrl'};
%, 'modalCtrl', 'strength'};
for m= metrics
    varnames= {m{1},'Phase', 'ID'};
    tbl=table([],[],[],'VariableNames', varnames); 
    for i=i_ict
        tbl= [tbl; table(mean(State_metrics(i).(m{1})(:,1:3))', [1:3]',...
            repmat(string([Metric_matrices(i).ID+"_"+dataSets_clean(i).sz_subtype]),3,1), ...
            'VariableNames', varnames)];
    end

    tbl.Phase = categorical(tbl.Phase); 
    tbl.ID = categorical(tbl.ID);
    
lme = fitlme(tbl, sprintf('%s ~ Phase + (1|ID)', m{1}))
lme2 = fitlme(tbl, sprintf('%s ~ Phase + (1|ID) + (Phase-1|ID)', m{1}))
lme3 = fitlme(tbl, sprintf('%s ~ Phase + (Phase|ID)', m{1}))

end

%% NOT USED: Connection Density, distance correlation

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

%% Generate Figures
%close all

% PCM figure
figure(1); set(gcf,'Name', 'pcm')
imagesc(Networks(1).pcm(:,:,6))
colormap('winter'); colorbar
xlabel('nodes'); ylabel('nodes');

% Example similarity matrix
A= [[rand(5,5)*.6, rand(5,5)];[rand(5,5), rand(5,5)*.4]]; A=A+A'; A=A.*~eye(10)+eye(10);
imagesc(A)
colormap(bone)
axis off

% Example ec matrix
A= [[rand(10,10), rand(10,10)];[rand(10,10), rand(10,10)]]; A=A+A'; A=A.*~eye(20)+eye(20);
imagesc(A)
colormap(gray)
axis off

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

% Slow/fast mode impulse response
x=linspace(0,1); 
rgb_array =...
[[0.2235    0.6118    0.6275];
    [0.8235    0.9765    0.9843];...
    [0.4588    0.8196    0.8353];...
    [0.0863    0.4863    0.5059];...
    [0.0196    0.1961    0.2039];...
    [0.2941    0.4078    0.6902];...
    [0.8431    0.8824    0.9882];...
    [0.5137    0.6157    0.8627];...
    [0.1412    0.2627    0.5608];...
    [0.0392    0.0941    0.2235];...
    [0.2824    0.8000    0.3176];...
    [0.8275    0.9922    0.8392];...
    [0.4980    0.9098    0.5255];...
    [0.1059    0.6549    0.1412];...
    [0.0235    0.2627    0.0392];...
];

figure(39)
clf; hold on
plot(exp(-9*x),'color', rgb_array(9,:))
plot(exp(-6*x),'color', rgb_array(8,:))
plot(exp(-2*x),'color', rgb_array(7,:))
plot(exp(-1.2*x),'color', rgb_array(3,:))
plot(exp(-1*x),'color', rgb_array(4,:))
colormap(rgb_array([9,6,8,7,2,3,1,4,5]',:))

% figure(11)
% set(gcf, 'Position', [944   184   705   771])
% saveas(gcf, 'FigsV3.3/controllability/subject_level_metrics_ictal.png')
% saveas(gcf, 'FigsV3.3/controllability/subject_level_metrics_ictal.fig')
% 
% figure(12)
% set(gcf, 'Position', [944   184   705   771])
% saveas(gcf, 'FigsV3.3/controllability/subject_level_metrics_preictal.png')
% saveas(gcf, 'FigsV3.3/controllability/subject_level_metrics_preictal.fig')




