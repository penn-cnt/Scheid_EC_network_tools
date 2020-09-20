function cleanSz = cleanRawSeizure(rawSz, fs)
% Pre-process raw seizure data with sampling rate fs. Enforces a common
% average reference, as well as a 4th order butterworth highpass, lowpass
% and bandpass filters. 

% Filter 
fc_line= [58, 61]; % cutoff frequency, Hz
fc_low = 120; 
fc_high = 1; 
order=4; 

data_0ref= rawSz-mean(rawSz);  
[b_stop,a_stop]= butter(order/2, 2*fc_line/fs, 'stop');  %Note order param is div 2
[b_low,a_low] = butter(order,(2*fc_low)/fs,'low');
[b_high,a_high] = butter(order,(2*fc_high)/fs,'high'); 

 cleanSz= zeros(size(rawSz));

 %filter all channels sequentially
 for i_chan = 1:size(rawSz, 1)
      cleanSz(i_chan,:) = filtfilt(b_stop,a_stop, data_0ref(i_chan, :)); 
      cleanSz(i_chan,:) = filtfilt(b_high,a_high, cleanSz(i_chan,:));
      cleanSz(i_chan, :)= filtfilt(b_low,a_low, cleanSz(i_chan,:)); 
 end