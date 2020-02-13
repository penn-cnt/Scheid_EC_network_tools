% For Borel optEnergy
run('initProject.m')

%% Part 1: Prepare Energy struct with initial and final values

% For each patient, this section computes an x0 initial state for the three seizure
% phases as the average signal bandpower in "freq_band" over all windows in
% each phase, resulting in an N x 1 initial state vector. 
% xf is calculated as the average band power across all preictal windows.

% x0_z and xf_z are calculated by taking the zscore of all x0s and xf. 

Energy=struct();
%load('Data/Energy.mat')

%Frequency Bands
alpha_theta=[5,15];
beta=[15,25];
low_g=[30,40];
high_g=[95,105];

freq_band= high_g;


info=['x0: mean bandpower across all windows in phase',...
      'xf: preict. bandpower', ...
      'B: single node, 1e-5 along diag'];

for i_set= i_ict % Note, should only iterate over ictal sets. 
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
    bandfun=@(x)bandpower(x', Fs, freq_band)';
    
    %Get preictal bandpower for state xf.
    pi_d=dataSets_clean(i_set+nSets/2);
    stDiff=round(d.UEOStart-d.EECStart);
    
    pre_data= [dataSets_clean(i_set+nSets/2).data(:,1+2*(Fs*stDiff):end), d.data(:,1:Fs*stDiff)];

    xf=mean(MovingWinFeats(pre_data,Fs, 1, 1, bandfun),2);
    Energy(i_set).xf=mean(MovingWinFeats(pi_d.data,Fs, 1, 1, bandfun),2);
    
    % Get ictal data from UEO
    data= d.data(:,1+Fs*stDiff:end);
   
    
    for s=1:3
        %Get representative network for state
        [~,m_ind]= max(sum(corrcoef(config(:,p.contigStates==s))));
        A_s= Net.pcm(:,:,p.contigStates==s); 
        Energy(i_set).repMats(:,:,s)= A_s(:,:,m_ind);
        
        %Get ECoG signal in longest contig. run of state
        t1= find(p.contigStates==s, 1, 'first');
        t2= find(p.contigStates==s, 1, 'last');
        signal=data(1:end,(Fs*(t1-1)+1:Fs*t2));
            
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
    
    % Calculate z-scored x0 and xf
    EZ= zscore([Energy(i_set).x0, Energy(i_set).xf], [],'all');
    Energy(i_set).x0_z=EZ(:,1:3); Energy(i_set).xf_z=EZ(:,4); 
    
end
disp('done')

%% Part 2: Find ideal trajectory u* and optimal energy

% Calculate control energy to drive state x0 to xf from each individual node
% with underlying connectivity as a represetative phase connectivity network A_s. 
% A_S gets stabilized as a continuous system

%load('Data/Energy.mat')

i_ict=find(strcmp({dataSets_clean.type},'ictal'));
relax= 1e-5; % value along non-driven diagonals of "relaxed" B matrix

% Energy Parameters to Test
% On first pass, run with large parameter set to find best parameters
tr= power(10, linspace(-2, log10(6), 10));
t_traj= round([tr(1:7),(1:6)],3) %log10 distribution b/w 0-1, then linear 1-6
rho= round([power(10, linspace(-2, 2, 10)), 110, 150],3); 

for i_set=i_ict
    tic
    i_set
    Energy(i_set).t_traj=t_traj;
    Energy(i_set).rho=rho;
    xf= Energy(i_set).xf;
    
    for s=1:3  % iterate through states
    
    x0= Energy(i_set).x0(:,s);
    N=length(x0);
    A= Energy(i_set).repMats(:,:,s);  
    A_s= A./(1+svds(A,1))-eye(N);       % normalize representative matrix
    
    trajErr= Energy(i_set).(sprintf('s%dtrajErr',s));
    nodeEnergy= Energy(i_set).(sprintf('s%dNodeEnergy',s));
    
%     [nodeEnergy, trajErr]= deal(zeros(length(x0),...
%         length(t_traj), length(rho)));
    xOpt=zeros(length(x0),length(t_traj), length(rho));
    
    for i_traj= [1:7] %1:length(t_traj) % try different time horizons
        tic
        T= t_traj(i_traj)
        for i_rho= 1:length(rho) % try different energy/distance tradeoffs
            r= rho(i_rho);
            try
                for n=1:length(x0) %Calculate energy per node
                    
                    B=relax*ones(length(x0),1); B(n)=1; B=diag(B); 
                    [X_opt, U_opt, n_err] = optim_fun(A_s, T, B, x0, xf, r, eye(length(A_s)));
                    nodeEnergy(n,i_traj,i_rho)=sum(vecnorm(B*U_opt').^2);
                    trajErr(n, i_traj,i_rho)=n_err;
                    xOpt(n, i_traj, i_rho)=mean((X_opt(end,1:N)'-xf).^2);
                    
                    % Uncoment to view output
%                     figure(1); imagesc(U_opt'); title('U'); 
%                     figure(2); imagesc(B); title('B');
%                     figure(3);
%                     subplot(1,6,[2:5]); imagesc([X_opt(:,1:N)'])% imagesc(X_opt(:,2:N)'); title('X');
%                     subplot(1,6,[2:5]); imagesc([X_opt(1,1:N)',X_opt(end,1:N)'])
%                     cl=caxis;
%                     subplot(1,6,1); imagesc(log(x0)); title('x0'); caxis(cl)
%                     subplot(1,6,6); imagesc(log(xf)); title('xf');  caxis(cl)

                end
                
            catch ME
                if strcmp(ME.identifier, 'MATLAB:svd:matrixWithNaNInf')
                    disp('err')
                    pause
                    trajErr(:, i_traj,i_rho)=nan(length(x0),1);
                    nodeEnergy(:,i_traj,i_rho)=nan(length(x0),1); 
                else; throw(ME)
                end
                
            end
            
        end % rho
        toc
    end % t_traj
    
    Energy(i_set).(sprintf('s%dtrajErr',s))= trajErr;
    Energy(i_set).(sprintf('s%dNodeEnergy',s))= nodeEnergy;
    Energy(i_set).(sprintf('s%Xopt_dist',s))= xOpt;
    
    end % end states loop
    toc
end 

save('Data/Energy2.mat', 'Energy', 't_traj', 'rho', 'freq_band', 'relax', 'info')
disp('Energy Calc done')