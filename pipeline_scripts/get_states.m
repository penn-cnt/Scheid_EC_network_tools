% get_states.m

Null=''; % set to 'Null' to perform on null model, and '' otherwise
nStates=3;
load(sprintf('Data/%sPartitions.mat', Null))
% eval(['Partitions=', Null, 'Partitions;']) 

for i_set=1:length(Partitions)
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
        fprintf('skipping %d \n',i_set)
        continue
    end
    
    l=length(st);
    u=unique(st);
    [runstates, runLength]=RunLength(st(:)');
    stateRuns=[runstates;runLength];
    
    % Get total size, max run length, and median of each state
    stateLens=arrayfun(@(x)sum(st==x), u)./l*100;
    stateMedian=(arrayfun(@(x)round(median(find(st==x))), u)-.5)./l;
    [runLen, mxRunIdx]=arrayfun(@(x)max(runLength(runstates==x)), u);
    
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
    
    % get contiguous states (maximum runs)
    contigSt=zeros(1,length(st));
    bnds= [[1,cumsum(stateRuns(2,1:end-1))+1]', cumsum(stateRuns(2,:))']; % run bounds
    for x=u
        bnd=bnds(stateRuns(1,:)==x,:); 
        mxbnd=bnd(mxRunIdx(x),:);
        contigSt(mxbnd(1):mxbnd(2))=tran(x);
    end

    Partitions(i_set).states=tran(st(:));    
    Partitions(i_set).contigStates=contigSt; 
    Partitions(i_set).stateLens=stateLens(inv(1:nStates));
    Partitions(i_set).stateMedian=stateMedian(inv(1:nStates))*100;
    Partitions(i_set).runLen=runLen(inv(1:nStates))./l*100;
    stateRuns(1,:)= tran(stateRuns(1,:));
    Partitions(i_set).stateRuns=stateRuns;

end

disp('Done')
% eval([Null,'Partitions=Partitions;'])
save('Data/Partitions.mat', sprintf('%sPartitions', Null),'nStates', '-append')

% disp('done')
% clear tran u f idx l p st nUnique top3 stateRuns runstates runLength 