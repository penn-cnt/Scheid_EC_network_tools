%% Robustness Measures

% SETTINGS %
load('Data/dataSets_clean.mat');

%subset=[1,2,7,9,10,14,16,20,21,23,27,37]
subset=[1,7,10,14,20,21,37];
all=[1:78];
subset=all(~ismember(all, subset))

betas= [0, .01, .05, .1, .5];
gammas= [ 0, .1, .25, .5, .75];

betaStr=replace(string(betas), '.', '_');
gammaStr=replace(string(gammas), '.', '_');

%%  Run Single Trial: 
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
    
     ctr=6;
%      load(sprintf('Data/Robustness/Network_%s.mat',...
%             replace(string(gamma), '.', '_')), 'Network')

     for i_set=[[21,37], subset+39]
        % Get different Gammas
        d=dataSets_clean(i_set);
        n=getNets(pwd, d.data, d.Fs, gamma, betas);
        n.ID=d.ID; n.type=d.type; n.block=d.block; n.Fs=d.Fs;
        Network(ctr)=n; 
        ctr=ctr+1;  
     end % end i_set
         
    save(sprintf('Data/Robustness/Network_%s.mat',...
        replace(string(gamma), '.', '_')), 'Network')
end

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

    save(sprintf('Data/Robustness/Network_%s.mat',gammaStr{g}), 'Network')
    end
end

disp('done updating WSims')
%% Community Detection on Networks 

for g=1:length(gammas)

fprintf('Glasso gamma %s \n', gammaStr{g})
load(sprintf('Data/Robustness/Network_%s.mat',gammaStr{g}))

%Partitions=struct(); 
%clear Partitions
load(sprintf('Data/Robustness/Partitions_%s.mat',gammaStr{g}))

for b=1:length(betas)
    ii = length(Network)*(b-1);
    for i_set=2:length(Network)
        d= Network(i_set);
        sim= Network(i_set).(sprintf('wSim_%s',betaStr{b}));
        n= findSimComms(sim, (0.8:.05:1.05), 100, 3);
        n.ID=d.ID; n.type=d.type; n.block=d.block; n.Fs=d.Fs;
        n.beta= betaStr{b}; 
        
        %Get assignments to three (or more) states
        nUnique= arrayfun(@(x)length(unique(n.consensusQcomms(:,x))),(1:length(n.gamma)));
        st= p.consensusQcomms(:,find(nUnique>=nStates,1,'first'))';

        % Check if median comms contains correct num of communities
        if isempty(st)
            nUnique=arrayfun(@(x)length(unique(n.medianQcomms(:,x))),(1:length(n.gamma)));
            st=n.medianQcomms(:,find(nUnique>=nStates,1,'first'))';
        end

        s=assignStates(st);
        Partitions(ii+i_set)= mergeStructs(n, s); 
    end
end

    save(sprintf('Data/Robustness/Partitions_%s.mat',gammaStr{g}), 'Partitions')
end
