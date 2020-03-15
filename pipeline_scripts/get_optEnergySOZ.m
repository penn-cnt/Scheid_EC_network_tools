%% optEnergySOZ SOZ permutation test

run('initProject.m')
% Note, get_optEnergy.m must be completed first. This will ensure consistency
% between xf, x0 and representative connectivity matrices.

%Frequency Bands
alpha_theta=[5,15];
beta=[15,25];
low_g=[30,40];
high_g=[95,105];
broad_g=[30,105]; 


%% Part 1: Prepare EnergySOZ struct with initial and final values

fldnames=fieldnames(Energy);
EnergySOZ= rmfield(Energy, fldnames(~ismember(fldnames,...
    {'ID', 'type', 'block', 'x0','repMats', 'xf', 'x0_z', 'xf_z','t_traj', 'rho'})));

sigma=1.6332; % input control energy spread
info= sprintf(['x0: mean bandpower',...
      'xf: preict state, ', ...
      'B: spread with sigma %0.04f'],sigma);
  
%% Part 2: Find ideal trajectory u* and optimal energy for SOZs

relax=1;
scale=10^5; 

for i_set=i_ict
    disp(i_set)
    
    % Get SOZ nodes:
    soz=dataSets_clean(i_set).sozGrid; 
    if sum(soz)==0; continue; end
    
    for s=1:3  % iterate through states
    
    x0= EnergySOZ(i_set).x0(:,s);
    xf= EnergySOZ(i_set).xf;
    
    A= EnergySOZ(i_set).repMats(:,:,s);  
    A_s= A./(1+svds(A,1))-eye(length(x0)); 
    
    %B= getSpreadControl(dataSets_clean(i_set).gridCoords, soz, sigma);
    B=relax*ones(length(x0),1); B(soz)=1; B=diag(B);
    R=scale*ones(length(x0),1); R(soz)=1; R=diag(R); 
    [nodeEnergySOZ, trajErr]= deal(zeros(length(t_traj), length(rho)));
    xOpt= zeros(length(t_traj), length(rho),length(soz));
    
    for i_traj=1:length(t_traj)
        tic
        T= t_traj(i_traj); 
        
        for i_rho= 1:length(rho)
            r= rho(i_rho); 
             %Calculate energy per node
            try   
                [X_opt, U_opt, n_err] = optim_fun_input_constrained(A_s, T, B, x0, xf, r, R, eye(length(A_s)));
                nodeEnergySOZ(i_traj,i_rho)=sum(vecnorm(B*U_opt').^2);
                trajErr(i_traj,i_rho)=n_err;
                xOpt(i_traj, i_rho,:)=X_opt(end,1:length(soz))';
            catch ME
                disp('ERROR!')
                if strcmp(ME.identifier, 'MATLAB:svd:matrixWithNaNInf')
                    disp('err')
                    pause(0.2)
                    trajErr(:, i_traj,i_rho)=nan(length(x0),1);
                    nodeEnergySOZ(:,i_traj,i_rho)=nan(length(x0),1);
                else 
                    throw(ME)
                end
            end
            
        end % rho
    end % t_traj
    
    EnergySOZ(i_set).(sprintf('s%dtrajErr',s))= trajErr;
    EnergySOZ(i_set).(sprintf('s%dNodeEnergySOZ',s))=nodeEnergySOZ;
    EnergySOZ(i_set).(sprintf('s%dX_optf',s))=xOpt;
    
    end % end states loop
end 


save('Data/EnergySOZ_V4_1.mat', 'EnergySOZ', 't_traj', 'rho', 'info', 'relax', 'scale')
disp('Part 2 EnergySOZ Calc done')

%% Part 3: Run SOZ permutation test using control parameters

% Only run for seizures included in the analysis
empties=find(cellfun(@isempty,{EnergySOZ.s1trajErr}));
i_ict(ismember(i_ict,[rm,empties]))=[];

percentRank = @(YourArray, TheProbes) reshape( mean( bsxfun(@le,...
    YourArray(:), TheProbes(:).') ) * 100, size(TheProbes) );

t_idx= 7; %power(10, linspace(-3, log10(5), 10)) %log10 distribution b/w 1-10
rho_idx= 7; % log10 distribution b/w 1-30

nperms= 100; % Number of permutations to select from

for i_set=i_ict
    fprintf('Finding null for i_set %d...\n', i_set)
    % Select the values of t and rho that were determined previously
    T= EnergySOZ(i_set).t_traj(t_idx);
    r= EnergySOZ(i_set).rho(rho_idx);

    % Get SOZ nodes:
    if isempty(dataSets_clean(i_set).sozGrid); continue; end
    soz= dataSets_clean(i_set).sozGrid; 
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
    
        %A_s= EnergySOZ(i_set).repMats(:,:,s);
        
        x0= EnergySOZ(i_set).x0(:,s);
        xf= EnergySOZ(i_set).xf;
        
        A= Energy(i_set).repMats(:,:,s);  
        A_s= A./(1+svds(A,1))-eye(length(x0));  
        
        sozEnergy=EnergySOZ(i_set).(sprintf('s%dNodeEnergySOZ',s))(t_idx, rho_idx); 
        nodeEnergynull= zeros(1,nperms);
     
            try
            for nt=1:nperms % Get null distribution of energy
                %B= getSpreadControl(dataSets_clean(i_set).gridCoords, i_perm(nt,:), sigma);
                B=relax*ones(length(x0),1); B(i_perm(nt,:))=1; B=diag(B);
                R=scale*ones(length(x0),1); R(i_perm(nt,:))=1; R=diag(R); 
                [X_opt, U_opt, n_err] = optim_fun_input_constrained(A_s, T, B, x0, xf, r, R, eye(length(A_s)));
                nodeEnergynull(nt)=sum(vecnorm(B*U_opt').^2);             
            end

        catch ME
            if strcmp(ME.identifier, 'MATLAB:svd:matrixWithNaNInf')
                disp('err')
                pause(0.2)
                nodeEnergynull(nt)=nan(length(x0),1); 
            end
            throw(ME)
            end %End try
            
    EnergySOZ(i_set).nodeEnergynull(:,s)=nodeEnergynull;
    
    % Get confidence interval
    EnergySOZ(i_set).SOZconfidence(s)=percentRank(nodeEnergynull, sozEnergy);
    
    end % end states loop
end 

%save('Data/EnergySOZ_V4_1.mat', 'EnergySOZ', 't_idx', 'rho_idx', '-append')
%disp('Part 3 EnergySOZ Calc done')

%% Part 3.5: Quantify error percentiles for each patient and state

fctr=3;

empties=find(cellfun(@isempty,{EnergySOZ.s1trajErr}));
i_ict(ismember(i_ict,[rm,empties]))=[];

percentRank = @(YourArray, TheProbes) reshape( mean( bsxfun(@le,...
    YourArray(:), TheProbes(:).') ) * 100, size(TheProbes) );

err_stats_soz= zeros(length(t_traj),length(rho), length(i_ict)*3);
final_diff_err= zeros(length(t_traj),length(rho), length(i_ict)*3);

allRanks=zeros(length(t_traj),length(rho),length(i_ict));
for i_set=2; %1:length(i_ict)
    phaseRanks=zeros(length(t_traj),length(rho),3);
    for s=1:3
        % Average computational error across nodes
        errs=EnergySOZ(i_ict(i_set)).(sprintf('s%dtrajErr',s));
        
        % mean squared error of xf and X_optf 
        xopt=EnergySOZ(i_ict(i_set)).(sprintf('s%dX_optf',s));
        diff= shiftdim(xopt,2)-repmat(EnergySOZ(i_ict(i_set)).xf, 1, length(t_traj));
        
      %  imagesc([squeeze(xopt(1,1,:)), EnergySOZ(i_ict(i_set)).xf])
        
        
        % calcualtioni Error percentile out of entire "parameter square" for as single phase
       phaseRanks(:,:,s)= errs; %percentRank(errs,errs); %Use if taking max of percentiles
       err_stats_soz(:,:,(i_set-1)*3+s)= phaseRanks(:,:,s); 
       final_diff_err(:,:,(i_set-1)*3+s)= squeeze(mean(diff.^2));
       
    end
    allRanks(:,:,i_set)=max(phaseRanks,[],3); % get max error at each cell across phases
end
    
max_percentiles=max(allRanks,[],3);
max_finaldist= max(final_diff_err,[], 3);

figure(1+fctr); clf
imagesc(percentRank(max_percentiles',max_percentiles'))
set(gca,'YDir','normal')
caxis([0,100])
title('maximum calculation error percentile across all siezures and phases')
xticks([1:length(t_traj)])
xticklabels(strsplit(sprintf('%0.02f ',t_traj')))
yticklabels(strsplit(sprintf('%0.02f ',rho)))
yticks([1:length(rho)])
ylabel('\rho'); xlabel('\tau'); colorbar
xtickangle(45); 

figure(2+fctr)
clf
imagesc(percentRank(max_finaldist',max_finaldist'))
set(gca,'YDir','normal')
caxis([0,100])
title('maximum final distance error percentile across all siezures and phases')
xticks([1:length(t_traj)])
xticklabels(strsplit(sprintf('%0.02f ',t_traj')))
yticklabels(strsplit(sprintf('%0.02f ',rho)))
yticks([1:length(rho)])
ylabel('\rho'); xlabel('\tau'); colorbar
xtickangle(45);

[pct_err, i_err]=min(max_percentiles, [], 'all', 'linear');
[t_opt, r_opt]=ind2sub(size(max_percentiles), i_err);
mn_err= [mean(err_stats_soz(t_opt, r_opt, :)),std(err_stats_soz(t_opt, r_opt, :))];

% saveas(gcf,'FigsV3.3/energy/energySOZ_B1_max_distErr_pcnt.png')
% saveas(gcf,'FigsV3.3/energy/energySOZ_B1_max_distErr_pcnt.fig')
% saveas(gcf,'FigsV3.3/energy/energySOZ_B1_max_err_pcnt.png')
% saveas(gcf,'FigsV3.3/energy/energySOZ_B1_max_err_pcnt.fig')

%save(fullfile(datafold, '/Energy.mat'), 'err_stats', '-append')

%% Visualize the SOZ signifcance level
% for each metric
nsoz=22; 

cb=97.5; %confidence bound (percent)
    
figure(4); clf
pcnts=[EnergySOZ.SOZconfidence];
%scatter(repmat([1:3],1,18)+normrnd(0,.05, [1,54]),  pcnts, 60, repmat(copper(18),3,1), 'filled')
%scatter(repmat([1:3],1,nsoz)+normrnd(0,.05, [1,3*nsoz]),  pcnts, 25,'grey')
scatter(repmat([1:3],1,nsoz)+repelem([1-nsoz:2:nsoz]/1e2, 3), pcnts, 40, 'filled','MarkerEdgeColor','black','MarkerFaceColor',[1,1,1]*.8)

%title(sprintf('Energy SOZ Permutation Test (all windows) significance',j))
title('Energy SOZ Permutation Test significance')
xlim([.5,3.5]); ylim([-5, 105])
xticklabels({'Phase1', 'Phase 2', 'Phase 3'})
xticks([1:3])
xtickangle(20)
hline(100-cb,'black:')

%hline(cb,'black:')

figure(5)
hold on

i_soz=~cellfun(@isempty,{EnergySOZ.SOZconfidence});

imagesc(reshape([EnergySOZ.SOZconfidence],3,nsoz)');
plot([1:3],(reshape([EnergySOZ.SOZconfidence],3,nsoz)'<=(100-cb)).*repmat([1:nsoz]',1,3), '*', 'color', 'red')
yticks([1:nsoz])
blk=cellfun(@num2str, {EnergySOZ(i_soz).block}', 'UniformOutput', false);
yticklabels(strcat({EnergySOZ(i_soz).ID}', {' '}, blk));
xticklabels({'Phase1', 'Phase 2', 'Phase 3'})
xticks([1:3])
title('Energy SOZ Permutation Test significance')
axis tight
ylim([.5, nsoz+.5])
set(gca, 'YDir', 'reverse')


% U-test of significance
nsozs=cellfun(@sum, {dataSets_clean(i_ict).sozGrid});
sigs= reshape([EnergySOZ.SOZconfidence],3,nsoz)'<=(100-cb)
[pp, hh, stats]=ranksum(nsozs(sigs(:,1)), nsozs(~sigs(:,1)))

%% Part 5: Add state EnergySOZ to State_Metrics

i_soz=find(~cellfun(@isempty, {EnergySOZ.s1trajErr}));
% Define optimal T and rho (use  functions below to work this out)
t_opt=7;
r_opt=7; 

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
save('Data/State_metrics.mat', 'State_metrics', 'tOpt', 'rOpt', '-append');

disp('done')

%% Spatial visual of node energies

r_opt=7; t_opt=7;
clf;
cmap=colormap; 
for i_set=8 %a
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
            colormap(brewermap([],'*PuBuGn'))
            colorbar
            %if s==3; colorbar; end
      
        
    end
    %suptitle(sprintf('Spatial Metrics for %s, block %d', dataSets_clean(i_set).ID, dataSets_clean(i_set).block))
   % pause
end

%export_fig('~/Desktop/try.png', '-transparent', '-r1000')

%% 1. Plot trajectories- Get energy 
clear energy

t_idx=7;
rho_idx=7;

i_set=8 % i_ict
    i_set
    % Select the values of t and rho that were determined previously
    T=EnergySOZ(i_set).t_traj(t_idx);
   % T=2
    r=EnergySOZ(i_set).rho(rho_idx);
    soz=dataSets_clean(i_set).sozGrid;
    xf= EnergySOZ(i_set).xf;

    for s=1:3
    
        A= EnergySOZ(i_set).repMats(:,:,s);  
        A_s= A./(1+svds(A,1))-eye(size(A,1));       % normalize representative matrix
        x0= EnergySOZ(i_set).x0(:,s);
        
        B=relax*ones(length(x0),1); B(soz)=1; B=diag(B);
        R=scale*ones(length(x0),1); R(soz)=1; R=diag(R); 
       % R=scale*ones(length(x0),1); R(i_perm(nt,:))=1; R=diag(R); 
        
        %B= getSpreadControl(dataSets_clean(i_set).gridCoords, soz, sigma);

        [X_opt, U_opt, n_err] = optim_fun_input_constrained(A_s, T, B, x0, xf, r, R, eye(length(A_s)));
%         nodeEnergySOZ(i_traj,i_rho)=sum(vecnorm(B*U_opt').^2);
        energy(i_set).X_opt(:,:,s)=X_opt;
        energy(i_set).U_opt(:,:,s)=U_opt';
        energy(i_set).n_err(:,:,s)=n_err;
    end     


figure(1); clf; imagesc(U_opt'); title('U'); 
figure(2); clf; imagesc(X_opt(:,1:length(soz))'); title('xopt');
figure(3); clf; imagesc(X_opt(:,end/2+1:end)'); title('p_star');

figure(4); clf;
%subplot(1,6,[2:5]); imagesc([X_opt(:,1:length(soz))'])% imagesc(X_opt(:,2:N)'); title('X');
subplot(1,6,[2:5]); imagesc(log([X_opt(1,1:length(soz))',X_opt(end,1:length(soz))'])); 
cl=caxis;
subplot(1,6,1); imagesc(log(x0)); title('x0'); caxis(cl)
subplot(1,6,6); imagesc(log(xf)); title('xf'); caxis(cl)

% figure(5); clf;
% subplot(2,1,1); imagesc(energy(i_set).X_opt(:,:,1))
% subplot(2,1,2); imagesc(energy(i_set).U_opt(:,:,1))
%% 2a Get distance from final location of x opt. 

plot(mean(energy(i_set).U_opt(:,:,3)))
xtraj=zeros(size(energy(i_set).X_opt(:,:,3),1),3);

xtraj(:,1)=sum(sqrt((energy(i_set).X_opt(:,1:length(soz),1)'-xf).^2));
xtraj(:,2)=sum(sqrt((energy(i_set).X_opt(:,1:length(soz),2)'-xf).^2));
xtraj(:,3)=sum(sqrt((energy(i_set).X_opt(:,1:length(soz),3)'-xf).^2));   

figure(5); clf
plot(xtraj)
%% 2b Get optimal control energy over time

plot(mean(energy(i_set).U_opt(:,:,3)))
xtraj=zeros(size(energy(i_set).X_opt,1),3);
for j=1:size(energy(i_set).U_opt,2)
    xtraj(j,1)=norm(B*energy(i_set).U_opt(:,j,1));
    xtraj(j,2)=norm(B*energy(i_set).U_opt(:,j,2));
    xtraj(j,3)=norm(B*energy(i_set).U_opt(:,j,3));
end

clf
%% 3 Plot Results
figure(2); clf; 
 semilogy([1:size(xtraj,1)]/1000, xtraj(:,1), 'color', cols(1,:))
 hold on;
 semilogy([1:size(xtraj,1)]/1000, xtraj(:,2), 'color', cols(2,:))
 semilogy([1:size(xtraj,1)]/1000, xtraj(:,3), 'color', cols(3,:))
 xlabel('\tau') 

ylabel('||X(t)-X_\tau|| (a.u.)')
%legend('Phase 1', 'Phase 2', 'Phase 3')
xlim([-.01,t_traj(t_idx)+.01])
box off


figure(1); clf; hold on
plot([1:size(xtraj,1)]/1000, xtraj(:,1), 'color', cols(1,:))
plot([1:size(xtraj,1)]/1000, xtraj(:,2), 'color', cols(2,:))
plot([1:size(xtraj,1)]/1000, xtraj(:,3), 'color', cols(3,:))
xlabel('\tau')
ylabel('||X(t)-X_T|| (a.u.)')
legend('Phase 1', 'Phase 2', 'Phase 3')
xlim([0,t_traj(t_idx)])
box off

%% Look at correlation of mean control metric of nodes in regions at top of distribution.
colormap('winter')
figure(304)
clf;
figure(7)
clf; 
%metrics={'aveCtrl', 'pModalCtrl', 'tModalCtrl', 'strength', 'strengthNeg'};
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
i_set=1

for s=1:3

eng=EnergySOZ(i_set).(sprintf('s%dNodeEnergySOZ',s));
trajErr=EnergySOZ(i_set).(sprintf('s%dtrajErr', s));
    [~,minErr]=min(squeeze(mean(trajErr))); 
    EnergySOZ(i_set).T_min(:,s)=minErr';

st=1;
lim=10;
subset=trajErr(st:lim, st:lim);

figure(1)
[X,Y] = ndgrid(st:lim, st:lim);
pointsize = 30;
scatter(X(:), Y(:), pointsize, log(subset(:)));
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


