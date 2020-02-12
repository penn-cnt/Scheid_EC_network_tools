% V3_viewing

% These sections are for scripts used to view data before/during final
% analysis

% Note, run the first sectino in V3_analysis first to load all metrics. 

nSets=length(Partitions);
nIct=nSets/2;
            %% Metric Trends (individual visualizations) %%

% Select metrics to plot
n=2; m=4;
all=false; % show both preictal and ictal
metrics={'strength', 'aveCtrl', 'modalCtrl', 'tModalCtrl'}; 
%metrics={'degree', 'aveCtrl', 'modalCtrl'}; 

for i_set=i_preict %1:nSets
    i_set
    figure(3); clf; 
    p= Partitions(i_set);
    
    if isempty(p.states)
        continue
    end
    
    st= p.contigStates; mm= Metric_matrices(i_set);
    
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

%% Show SOZ channels in data

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

%% Zoomed State Data View

i_set= 9; 

p=Partitions(i_set);
d=dataSets_clean(i_set);
stDiff=round(d.UEOStart-d.EECStart);
data= d.data(:,1+Fs*stDiff:end);
Fs=round(d.Fs);

xf=Energy(i_set).xf;

for s=1:3
    t1= find(p.contigStates==s, 1, 'first');
    t2= find(p.contigStates==s, 1, 'last');
    signal=data(1:end,(Fs*(t1-1)+1:Fs*t2));
    
    x0=Energy(i_set).x0(:,s);
    
    % Look at Energy Trajectory- Distance over time
    figure(1); clf
    subplot(121); hold on
    plot(d.data'+[1:size(d.data,1)]*1000);
    vline(d.Fs*(t1-1)+1); vline(d.Fs*t2); axis tight
    subplot(122); plot(signal'+[1:size(signal,1)]*1000);
    axis tight; 
    
    figure(2)
    imagesc(p.states) 
    colormap(gca, cols([5,2,6],:)); set(gca, 'YTick', [], 'fontsize', 18)

    figure(3)
    imagesc(reshape(xf, 8, 8)); colorbar
    title('avg final bandpower (preictal)'), caxis([0,.5])
    
    figure(4)
    imagesc(reshape(x0, 8, 8)); colorbar; caxis([0,.5])
    title(sprintf('avg bandpower ictal phase %d (%s %d)', s, d.ID, d.block))
      
    pause
end
