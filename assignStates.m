function [s]= assignStates(st)
% st: the original partition assignment
% assignStates will ensure that states are chronological. 
   
    if isempty(st)
        error('nCommunities < nStates')
    end
    
    l=length(st);
    u=unique(st);
    nStates=length(u);
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
    for x=u'
        bnd=bnds(stateRuns(1,:)==x,:); 
        mxbnd=bnd(mxRunIdx(x),:);
        contigSt(mxbnd(1):mxbnd(2))=tran(x);
    end

    s.states=tran(st(:));    
    s.contigStates=contigSt; 
    s.stateLens=stateLens(inv(1:nStates));
    s.stateMedian=stateMedian(inv(1:nStates))*100;
    s.runLen=runLen(inv(1:nStates))./l*100;
    stateRuns(1,:)= tran(stateRuns(1,:));
    s.stateRuns=stateRuns;

end