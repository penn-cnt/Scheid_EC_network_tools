% Get Energy
load('Data/dataSets_clean')
load('Data/Partitions')
loda('Data/Networks')

dataSets_clean= clean_EC_data;
i_ict=find(strcmp({dataSets_clean.type},'ictal'));

%% Part 1: Prepare Energy struct with initial and final values

Energy=struct();

% Frequency Bands
alpha_theta=[5,15];
beta=[15,25];
low_g=[30,40];
high_g=[95,105];

%line length function
llfun=@(x)sum(abs(diff(x,[],2)),2);

for i_set= i_ict
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
    
    % Set up bandpower function
    bandfun=@(x)bandpower(x', Fs, high_g)';
    
    % Get preictal bandpower for state xf.
    pi_d=dataSets_clean(i_set+nSets/2);
    xf=mean(MovingWinFeats(pi_d.data,Fs, 1, 1, bandfun),2);
    Energy(i_set).xf=mean(MovingWinFeats(pi_d.data,Fs, 1, 1, bandfun),2);
   
    
    for s=1:3
        % Get representative network for state
        [~,m_ind]= max(sum(corrcoef(config(:,p.states==s))));
        A_s= Net.pcm(:,:,p.states(m_ind)); 
        Energy(i_set).repMats(:,:,s)= A_s;
        
        % Get ECoG signal in longest contig. run of state
        st_inds=p.stateRuns(1,:)==s;
        [m_run, i_max]=max(p.stateRuns(2,st_inds));
        t2= sum(p.stateRuns(2,1:find(cumsum(st_inds)==i_max,1,'first')));
        t1= t2-p.stateRuns(2,find(cumsum(st_inds)==i_max,1,'first'))+1;
        signal=d.data(1:end,(Fs*(t1-1)+1:Fs*t2));
        
        % Uncomment to visualize and double check
%         figure(1); clf
%         subplot(121); hold on
%         plot(d.data'+[1:size(d.data,1)]*1000);
%         stem(d.Fs*(t1-1)+1, 50000, 'lineWidth', 2, 'color', 'red'); 
%         stem(d.Fs*t2, 50000, 'lineWidth', 2, 'color', 'red');
%         axis tight
%         subplot(122)
%         plot(signal'+[1:size(signal,1)]*1000);
%         axis tight; pause
    %     figure(2)
    %     imagesc(p.states) 
    %     colormap(gca, cols([5,2,6],:));
    %     set(gca, 'YTick', [], 'fontsize', 18)
        
        %Compute average bandpower/1 sec window of state
        x0_t=MovingWinFeats(signal,Fs, 1, 1, bandfun);
        Energy(i_set).x0(:,s)= mean(x0_t,2);
        
%         % Uncomment to View x0 and xf over time windows:
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
    
end
disp('done')


%% Part 2: Find ideal parameters

% Energy Parameters to Test
t_traj= power(10, linspace(-3, 1, 10)) %log10 distribution b/w 1-10
rho= power(10, linspace(-3, 1, 10)) % log10 distribution b/w 1-30

for i_set=2;% 1:i_ict
    
    for s=1:3  % iterate through states
    
    A_s= Energy(i_set).repMats(:,:,s);
    x0= Energy(i_set).x0(:,s);
    xf= Energy(i_set).xf;
    [nodeEnergy, trajErr]= deal(zeros(length(x0),...
        length(t_traj), length(rho)));
    compTime=zeros(length(t_traj), 1);
    
    for i_traj= 1:length(t_traj)
        tic
        T= t_traj(i_traj); 
        
        for i_rho=1:length(rho)
            r= rho(i_rho); 
             %Calculate energy per node
            try
                for n=1:length(x0)
                    B=10e-5*ones(length(x0),length(x0)); B(n,n)=1;
                    [X_opt, U_opt, n_err] = optim_fun(A_s, T, B, x0, xf, r, eye(length(A_s)));
                    nodeEnergy(n,i_traj,i_rho)=sum(vecnorm(B*U_opt').^2);
                    trajErr(n, i_traj,i_rho)=n_err;
                end
                
            catch ME
                if strcmp(ME.identifier, 'MATLAB:svd:matrixWithNaNInf')
                    disp('err')
                    trajErr(:, i_traj,i_rho)=nan(length(x0),1);
                    nodeEnergy(:,i_traj,i_rho)=nan(length(x0),1); 
                else; throw ME
                end
            end
            
        end % rho
        compTime(i_traj)=toc
    end % t_traj
    
    Energy(i_set).(sprintf('s%dtrajErr',s))= trajErr;
    Energy(i_set).(sprintf('s%dNodeEnergy',s))=nodeEnergy;
    Energy(i_set).compuTime=compTime;
    
    end % end states loop
end 

%% Visualization for zooming in 
close all
for s=1:3

eng=Energy(i_set).(sprintf('s%dNodeEnergy',s));
trajErr=Energy(i_set).(sprintf('s%dtrajErr', s));

st=2;
lim=9;
subset=trajErr(:, st:lim, st:lim);

figure(1)
[X,Y,Z] = ndgrid(1:size(x0), st:lim, st:lim);
pointsize = 30;
scatter3(X(:), Y(:), Z(:), pointsize, subset(:));
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
plot(squeeze(mean(subset)), 'color', cols(s,:))
xlabel('T')
ylabel('Error')
set(gca, 'XTickLabels', round(rho(st:lim),2));
xticks((st:lim)-st+1)
%legend(string(round(rho(st:lim),2)'))
% 
 pause(.3)
end
%%
figure(2)
for i=1:5
    i
imagesc(log10(trajErr(:,:,i)))
%plot(mean(trajErr))
ylabel('n')
xlabel('T')
colormap(autumn(64))
colorbar
pause(.5)
end
