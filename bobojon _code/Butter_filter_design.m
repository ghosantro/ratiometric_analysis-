function[b,a] = Butter_filter_design(pass_freq, stop_freq,frameperiod)
%This function filtered the signal trace to remove noise built in
%butterworth filter is used 
% optimize this filter 

nyquist_frequency = 1000/frameperiod/2; % half of the sampling frequency 

Wp = pass_freq/nyquist_frequency;
Ws = stop_freq/nyquist_frequency;
[n,Wn] = buttord(Wp,Ws,3,60);

[b,a] = butter(n,Wn,'low');
% freqz(b,a,512,1000/frameperiod); 
% title('n=5 Butterworth Lowpass Filter')


