run('initProjectBorel.m')

disp('running Energy SOZ B1')

fldnames=fieldnames(Energy);
EnergySOZ= rmfield(Energy, fldnames(~ismember(fldnames,...
    {'ID', 'type', 'block', 'x0','repMats', 'xf', 'x0_z', 'xf_z','t_traj', 'rho'})));

%sigma=1.6332; % input control energy spread
relax= 1e-5;
info= ['x0: mean bandpower',...
      'xf: preict state, ', ...
      'B: 1 at SOZ points, 1e-5 else'];
  
%% Part 2: Find ideal trajectory u* and optimal energy for SOZs

tr= power(10, linspace(-2, log10(6), 10));
t_traj= round([tr(1:7),(1:6)],3)

for i_set=i_ict
    disp(i_set)
    EnergySOZ(i_set).t_traj=t_traj;
    % Get SOZ nodes:
%     if ~any(dataSets_clean(i_set).sozGrid); continue; end
    soz=dataSets_clean(i_set).sozGrid; 
    if sum(soz)==0; continue; end
    
    for s=1:3  % iterate through states
    
   % Normalize A_s
    
    A= Energy(i_set).repMats(:,:,s);  
    A_s= A./(1+svds(A,1))-eye(length(soz));   
    
    %B= getSpreadControl(dataSets_clean(i_set).gridCoords, soz, sigma);
    x0= EnergySOZ(i_set).x0(:,s);
    xf= EnergySOZ(i_set).xf;
    B=relax*ones(length(x0),1); B(soz)=1; B=diag(B);
    [nodeEnergySOZ, trajErr]= deal(zeros(length(t_traj), length(rho)));
    xOpt= zeros(length(t_traj), length(rho),length(soz));
    
    for i_traj=1:length(t_traj)
        tic
        T= t_traj(i_traj); 
        
        for i_rho= 1:length(rho)
            r= rho(i_rho); 
             %Calculate energy per node
            try   
                [X_opt, U_opt, n_err] = optim_fun(A_s, T, B, x0, xf, r, eye(length(A_s)));
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


save('Data/EnergySOZ_B1.mat', 'EnergySOZ', 't_traj', 'rho', 'info', 'relax')
disp('Part 2 EnergySOZ Calc done')

%%
% % Only run for seizures included in the analysis
% empties=find(cellfun(@isempty,{EnergySOZ.s1trajErr}));
% i_ict(ismember(i_ict,[rm,empties]))=[];
% 
% percentRank = @(YourArray, TheProbes) reshape( mean( bsxfun(@le,...
%     YourArray(:), TheProbes(:).') ) * 100, size(TheProbes) );
% 
% t_idx= 7; %power(10, linspace(-3, log10(5), 10)) %log10 distribution b/w 1-10
% rho_idx= 7; % log10 distribution b/w 1-30
% 
% nperms= 1000; % Number of permutations to select from
% 
% for i_set=i_ict
%     fprintf('Finding null for i_set %d...\n', i_set)
%     % Select the values of t and rho that were determined previously
%     T= EnergySOZ(i_set).t_traj(t_idx);
%     r= EnergySOZ(i_set).rho(rho_idx);
% 
%     % Get SOZ nodes:
%     if isempty(dataSets_clean(i_set).sozGrid); continue; end
%     soz= dataSets_clean(i_set).sozGrid; 
%     if sum(soz)==0; continue; end
%     
%     % Get nullset
%     nsoz=sum(soz); 
%     nullset=find(~soz); 
%     i_perm=zeros(nperms, nsoz);
%     rng(5)
%     for i=1:nperms
%         i_perm(i,:)=nullset(randperm(sum(~soz),nsoz));
%     end
%     
%     EnergySOZ(i_set).nullSets=i_perm; 
%     EnergySOZ(i_set).SOZconfidence=zeros(1,3);
%     EnergySOZ(i_set).nodeEnergynull=zeros(nperms,3);
%        
%     for s=1:3  % iterate through states
%     
%         A_s= EnergySOZ(i_set).repMats(:,:,s);
%         x0= EnergySOZ(i_set).x0(:,s);
%         xf= EnergySOZ(i_set).xf;
%         sozEnergy=EnergySOZ(i_set).(sprintf('s%dNodeEnergySOZ',s))(t_idx, rho_idx); 
%         nodeEnergynull= zeros(1,nperms);
%      
%             try
%             for nt=1:nperms % Get null distribution of energy
%                 B= getSpreadControl(dataSets_clean(i_set).gridCoords, i_perm(nt,:), sigma);
%                 %B=const*ones(length(x0),1); B(i_perm(nt,:))=1; B=diag(B);
%                 [X_opt, U_opt, n_err] = optim_fun(A_s, T, B, x0, xf, r, eye(length(A_s)));
%                 nodeEnergynull(nt)=sum(vecnorm(B*U_opt').^2);             
%             end
% 
%         catch ME
%             if strcmp(ME.identifier, 'MATLAB:svd:matrixWithNaNInf')
%                 disp('err')
%                 pause(0.2)
%                 nodeEnergynull(nt)=nan(length(x0),1); 
%             end
%             end %End try
%             
%     EnergySOZ(i_set).nodeEnergynull(:,s)=nodeEnergynull;
%     
%     % Get confidence interval
%     EnergySOZ(i_set).SOZconfidence(s)=percentRank(nodeEnergynull, sozEnergy);
%     
%     end % end states loop
% end 
% 
% save('Data/EnergySOZ', 'EnergySOZ', 't_traj', 'rho', 'freq')
% disp('Part 3 EnergySOZ Calc done')