%% Robustness Measures

% SETTINGS %
load('Data/dataSets_clean.mat');
addpath('helper_functions'); 
addpath('pipeline_scripts')

%subset=[1,2,7,9,10,14,16,18,20,21,23,25,27,37] % seizures included
all=[1:78];
%subset=all(~ismember(all, subset))

betas= [0, .01, .05, .1];
gammas= [ 0, .1, .25, .5, .75];

betaStr=replace(string(betas), '.', '_');
gammaStr=replace(string(gammas), '.', '_');
nSets=28; 

%%  Run Single Trial 
% i_set= 5;
% d= dataSets_clean(i_set); 
% network= getNets(pwd, d.data, d.Fs, gamma, betas); 
% sim = network.sim; 
% b=.05
% Wsim= zeros(length(sim)); 
% for i =1:length(sim)
%     for j =1:length(sim)
%         Wsim(i,j)=sim(i,j)*((1 - b*abs(i-j)/length(sim)));    
%     end
% end
% 
% partitions= findSimComms(Wsim, (0.8:.05:1.05), 100, 3);
% s= assignStates(partitions, 3);
% partitions= mergeStructs(partitions, s);
% 
% figure(2)
% imagesc(partitions.states)

%% Generate Networks for different Gammas and Betas
for gamma = gammas
    
     ctr=25;
     load(sprintf('Data/Robustness/Network_%s.mat',...
            replace(string(gamma), '.', '_')), 'Network')
        
     for i_set=[subset, subset+39]
        % Get different Gammas
        d=dataSets_clean(i_set);
        n=getNets(pwd, d.data, d.Fs, gamma, betas);
        n.ID=d.ID; n.type=d.type; n.block=d.block; n.Fs=d.Fs;
        if ismember('wSim_1', fieldnames(Network))
            n.wSim_1=[]; n.wSim_0_25=[]; n.wSim_0_75=[]; n.wSim_0_5=[];
        end
        if ~ismember('wSim_0_5', fieldnames(Network))
            if ismember('wSim_0_5', fieldnames(n)), n=rmfield(n,'wSim_0_5'), end
        end
        Network(ctr)=orderfields(n, Network(1));
        ctr=ctr+1;  
     end % end i_set
         
    save(sprintf('Data/Robustness/Network_%s.mat',...
        replace(string(gamma), '.', '_')), 'Network')
end

% 
% % Reorder at end
% for gamma = gammas
%     
%      load(sprintf('Data/Robustness/Network_%s.mat',...
%             replace(string(gamma), '.', '_')), 'Network')
%      
%         T = struct2table(Network); % convert the struct array to a table
%         sortedT = sortrows(T, 'ID'); sortedT = sortrows(sortedT, 'type');
%         Network = table2struct(sortedT);
%         
%             save(sprintf('Data/Robustness/Network_%s.mat',...
%         replace(string(gamma), '.', '_')), 'Network')
% end

%% Add/update Wsims

for g=1:length(gammas)
    fprintf('g=%d \t', g)

load(sprintf('Data/Robustness/Network_%s.mat',gammaStr{g}))

    for i_set=1:length(Network)
        sim= Network(i_set).sim;
        for b=betas
            Wsim= zeros(length(sim)); 
            for i =1:length(sim)
                for j =1:length(sim)
                    Wsim(i,j)=sim(i,j)*((1 - b*abs(i-j)/length(sim))^2);    
                end
            end
            Network(i_set).(sprintf('wSim_%s', replace(string(b), '.', '_')))= Wsim;
        end

    %save(sprintf('Data/Robustness/Network_%s.mat',gammaStr{g}), 'Network')
    end
end

disp('done updating WSims')

%% Community Detection on Networks 

nStates=3;

for g=1:length(gammas)

fprintf('Glasso gamma %s \n', gammaStr{g})
load(sprintf('Data/Robustness/Network_%s.mat',gammaStr{g}))

% Partitions=struct(); 
% clear Partitions
load(sprintf('Data/Robustness/Partitions_%s.mat',gammaStr{g}))

ctr=101

for b=1:length(betas)% Note: beta=2 has been done
    ii = length(Network)*(b-1);
    for i_set=[2, 12, 16,26]%1:length(Network)
        d= Network(i_set);
        sim= Network(i_set).(sprintf('wSim_%s',betaStr{b}));
        
        n= findSimComms(sim, (0.8:.05:1.05), 100, 3);
        n.ID=d.ID; n.type=d.type; n.block=d.block; n.Fs=d.Fs;
        n.beta= betaStr{b}; 
        
        %Get assignments to three (or more) states
        nUnique= arrayfun(@(x)length(unique(n.consensusQcomms(:,x))),(1:length(n.gamma)));
        st= n.consensusQcomms(:,find(nUnique==nStates,1,'first'))';
        %st= p.consensusQcomms(:,find(nUnique>=nStates,1,'first'))';

        % Check if median comms contains correct num of communities
        if isempty(st)
            nUnique=arrayfun(@(x)length(unique(n.medianQcomms(:,x))),(1:length(n.gamma)));
            st=n.medianQcomms(:,find(nUnique>=nStates,1,'first'))';
        end

        s=assignStates(st);
        Partitions(ctr)= mergeStructs(n, s); 
        ctr=ctr+1
    end
end

    %save(sprintf('Data/Robustness/Partitions_%s.mat',gammaStr{g}), 'Partitions')
end

%% State (re)Assignment
for g=1:length(gammas)

fprintf('Glasso gamma %s \n', gammaStr{g})
load(sprintf('Data/Robustness/Partitions_%s.mat',gammaStr{g}))

for b=1:length(betas)
    ii = nSets*(b-1);
    for i_set=1:nSets
        p=Partitions(ii+i_set); 

        
        %Get assignments to three (or more) states
        nUnique= arrayfun(@(x)length(unique(p.consensusQcomms(:,x))),(1:length(p.gamma)));
        st= p.consensusQcomms(:,find(nUnique==nStates,1,'first'))';
        %st= p.consensusQcomms(:,find(nUnique>=nStates,1,'first'))';

        % Check if median comms contains correct num of communities
        if isempty(st)
            nUnique=arrayfun(@(x)length(unique(p.medianQcomms(:,x))),(1:length(p.gamma)));
            st=p.medianQcomms(:,find(nUnique==nStates,1,'first'))';
        end
        disp(nUnique(find(nUnique==nStates,1,'first')));

         s=assignStates(st);
         Partitions(ii+i_set)=mergeStructs(s,Partitions(ii+i_set)); 
    end
end

%     save(sprintf('Data/Robustness/Partitions_%s.mat',gammaStr{g}), 'Partitions')
end

%% Reorder Partitions
% Reorder at end
for gamma = gammas
    
     load(sprintf('Data/Robustness/Partitions_%s.mat',...
            replace(string(gamma), '.', '_')), 'Partitions')
     
        T = struct2table(Partitions); % convert the struct array to a table
        sortedT = sortrows(T, 'ID'); sortedT = sortrows(sortedT, 'type'); sortedT= sortrows(sortedT,'beta'); 
        Partitions = table2struct(sortedT);
        
            save(sprintf('Data/Robustness/Partitions_%s.mat',...
        replace(string(gamma), '.', '_')), 'Partitions')
end

%% Calculate All Metrics 

% Get metrics
for g=2:length(gammas)

    fprintf('Glasso gamma %s \n', gammaStr{g})
    load(sprintf('Data/Robustness/Partitions_%s.mat',gammaStr{g}))
    load(sprintf('Data/Robustness/Network_%s.mat',gammaStr{g}))
    
    % Note, metrics.m has been changed for use with robustness measures
    run get_metrics.m
    
    save(sprintf('Data/Robustness/Metric_matrices_%s.mat',gammaStr{g}), 'Metric_matrices')
    save(sprintf('Data/Robustness/State_metrics_%s.mat',gammaStr{g}), 'State_metrics')

end

%% Transient/Persistent Metric Calculation

thresholds=[.05:.05:.25]; 

modalMetRobustness=struct();
modalStateRobustness=struct();

g=4

% Get metrics
%for g=1:length(gammas)
    fprintf('Glasso gamma %s \n', gammaStr{g})
    load(sprintf('Data/Robustness/Partitions_%s.mat',gammaStr{g}))
    load(sprintf('Data/Robustness/Network_%s.mat',gammaStr{g}))
    
    for th=1:length(thresholds)
    % Note, metrics.m has been changed for use with robustness measures
    [mm, sm]= getMets(Network, Partitions, '0_01', {'pModalCtrl', 'tModalCtrl'}, thresholds(th));  
    
    [mm(1:28).thresh]=deal(thresholds(th));
    [sm(1:28).thresh]=deal(thresholds(th));
    
    [mm(1:28).gamma]=deal(gammas(g));
    [sm(1:28).gamma]=deal(gammas(g));
    
        if isempty(fieldnames(modalMetRobustness)) 
            modalMetRobustness= mm; 
            modalStateRobustness= sm; 
        else 
        modalMetRobustness=[modalMetRobustness, mm]; 
        modalStateRobustness=[modalStateRobustness, sm]; 
        end
       
    end
    
%end

%save('Data/Robustness/modalThresholds.mat', 'modalMetRobustness', 'modalStateRobustness', 'thresholds'); 
