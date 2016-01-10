function filtered_sig = ButterFilter(trace,b,a)
%This function filtered the signal trace to remove noise built in
%butterworth filter is used 


filtered_sig = filtfilt(b, a, trace);
%filtered_sig = smooth(filtered_sig);
%figure; hold on; plot(trace,'r');plot(filtered_sig);hold off
%hold on;  plot(trace,'b');hold off;  