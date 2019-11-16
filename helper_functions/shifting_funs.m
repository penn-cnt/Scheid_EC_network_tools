% Shift to use UEO instead of EEO

%% Ensure that all data lengths have the correct values. 

for i_set=1:78
    oldL=size(Networks(i_set).rhos,1)
    newL=size(Networks(i_set).pcm,3)
    
    dd=oldL-newL
    
    Fs=round(dataSets_clean(i_set).Fs);
    
    % is dataset long enough?
    if length(dataSets_clean(i_set).data)/Fs ~= newL
        plot(dataSets_clean(i_set).data'+[1:size(dataSets_clean(i_set).data,1)]*1000)
        if strcmp(Networks(i_set).type, 'ictal')
            dataSets_clean(i_set).data=dataSets_clean(i_set).data(1:(end-Fs*dd))
        else 
            dataSets_clean(i_set).data=dataSets_clean(i_set).data((1+Fs*dd):end)
        end
    end
            
   % does block end match?
   if floor(dataSets_clean(i_set).blockEnd-dataSets_clean(i_set).EECStart) ~= newL
        if strcmp(Networks(i_set).type, 'ictal')
            dataSets_clean(i_set).blockEnd=dataSets_clean(i_set).blockEnd-dd
        else 
            dataSets_clean(i_set).EECStart=dataSets_clean(i_set).EECStart+dd
        end
   end
    
end
   
%% Compute UEOStart Networks       

for i_set=1:39
    
    stDiff=round(dataSets_clean(i_set).UEOStart-dataSets_clean(i_set).EECStart)
    Fs=round(dataSets_clean(i_set).Fs);
    
    sz=size(Networks(i_set).pcm,3)
    sub=floor(dataSets_clean(i_set).blockEnd-dataSets_clean(i_set).EECStart-stDiff)
    
    figure(1) 
    clf
    hold on
    plot(dataSets_clean(i_set).data'+[1:size(dataSets_clean(i_set).data,1)]*1000)
    stem(Fs*stDiff, 1000*size(dataSets_clean(i_set).data,1), 'linewidth', 2)
    
    if size(Networks(i_set).pcm,3)<=floor(dataSets_clean(i_set).blockEnd-dataSets_clean(i_set).EECStart-stDiff)
        fprintf('%d shifted to UEO \n', i_set)
        pause
        continue
    end
        
    pause
    
    if stDiff>0
        % shift stDiff networks to preictal
        Networks(i_set+39).pcm=cat(3,Networks(i_set+39).pcm(:,:,2*stDiff+1:end), Networks(i_set).pcm(:,:,1:stDiff));      
        Networks(i_set).pcm=Networks(i_set).pcm(:,:,stDiff+1:end);
        
        Networks(i_set+39).icov=cat(3,Networks(i_set+39).icov(:,:,2*stDiff+1:end), Networks(i_set).icov(:,:,1:stDiff));      
        Networks(i_set).icov=Networks(i_set).icov(:,:,stDiff+1:end);
        
        Networks(i_set+39).config_pcm=cat(2,Networks(i_set+39).config_pcm(:,2*stDiff+1:end), Networks(i_set).config_pcm(:,1:stDiff));      
        Networks(i_set).config_pcm=Networks(i_set).config_pcm(:,stDiff+1:end);
        
        Networks(i_set+39).config_icov=cat(2,Networks(i_set+39).config_icov(:,2*stDiff+1:end), Networks(i_set).config_icov(:,1:stDiff));      
        Networks(i_set).config_icov=Networks(i_set).config_icov(:,stDiff+1:end); 
    end
end
%save('NetworksUEO.mat', 'mat')