% get_states.m

Null=''; % set to 'Null' to perform on null model, and '' otherwise

% load(sprintf('Data/%sPartitions.mat', Null))
% eval(['Partitions=', Null, 'Partitions;']) 
% nSets=length(Partitions);

for i_set=1:nSets
    p=Partitions(i_set);
    
    %Get assignments to three (or more) states
    nUnique= arrayfun(@(x)length(unique(p.consensusQcomms(:,x))),(1:length(p.gamma)));
    st= p.consensusQcomms(:,find(nUnique>=nStates,1,'first'))';
    
    % Check if median comms contains correct num of communities
    if isempty(st)
        nUnique=arrayfun(@(x)length(unique(p.medianQcomms(:,x))),(1:length(p.gamma)));
        st=p.medianQcomms(:,find(nUnique>=nStates,1,'first'))';
    end
    
    if isempty(st)
        continue
    end
    
    l=length(st);
    u=unique(st);
    [runstates, runLength]=RunLength(st(:)');
    stateRuns=[runstates;runLength];
    
    % Get total size, max run length, and median of each state
    stateLens=arrayfun(@(x)sum(st==x), u)./l*100;
    stateMedian=(arrayfun(@(x)round(median(find(st==x))), u)-.5)./l;
    runLen=arrayfun(@(x)max(runLength(runstates==x)), u);
    
    % Select top 3 states
    [~,idx]=sort(stateLens,'descend'); %idx of top states 
    topStates=idx(1:nStates); 
   
    % Get temporal order transformation for top3 states.
    tran=zeros(1,length(u));
    f=arrayfun(@(x)find(ismember(stateRuns',[x,runLen(x)],'rows'),1),topStates);
    [~,fin]=sort(f);
    tran(topStates(fin))=[1:nStates]; % Transformed ordering of states 
    if length(u)>nStates; tran(tran==0)=u(nStates+1:end); end
    [~, inv]=sort(tran);  % inverted state transformation 

    Partitions(i_set).states=tran(st(:));
    Partitions(i_set).stateLens=stateLens(inv(1:nStates));
    Partitions(i_set).stateMedian=stateMedian(inv(1:nStates))*100;
    Partitions(i_set).runLen=runLen(inv(1:nStates))./l*100;
    stateRuns(1,:)= tran(stateRuns(1,:));
    Partitions(i_set).stateRuns=stateRuns;

end

%save('DataBackup/Partitions.mat', 'Partitions', 'nStates'); 

eval([Null,'Partitions=Partitions;'])
save(sprintf('Data/%sWPartitions.mat', Null), sprintf('%sPartitions', Null),...
    'nStates')
disp('done')
clear tran u f idx l p st nUnique top3 stateRuns runstates runLength 