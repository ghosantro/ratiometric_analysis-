function [upstroke_time_points,rapid_depolar_end] = Get_upstroke(trace,frameperiod, plot_index)
%% This function get the beginning and end points of the upstroke 

gradient = (trace(2:end)-trace(1:end-1))./frameperiod;
index_depolar = find (gradient>0.0035); % Capture the fast depolarisation phase 
subs_index = index_depolar(2:end)-index_depolar(1:end-1);
discontinious_index = find(subs_index>1)+3; %distinguish depolarisations (to get rid of littel humps set the discontinous index to 10 instead of 1)  
phase_0_start = [index_depolar(1);index_depolar(discontinious_index)];
%This needs to be further filtered by a big rise in aplitude 
rapid_rise_end = [index_depolar(discontinious_index-1);index_depolar(end)];

rise_amplitude = trace(rapid_rise_end)-trace(phase_0_start); 
significant_rise_index = find(rise_amplitude>0.35 );   %IF THE APD MAP HAS LOTS OF HOLES CHANGE THE THRESHOLD HERE
phase_0_start = phase_0_start(significant_rise_index);
rapid_rise_end = rapid_rise_end(significant_rise_index);

% if the two activation points are found to be too close to in each other
% remove the later one and this must be paired up with one false rapid rise
% end being found, therefore remove the first one from the rapid rise
% vector 
remove_index = []; 
for i = 1: length(phase_0_start)-1
    if phase_0_start(i+1)-phase_0_start(i)<30/frameperiod 
       remove_index = [remove_index,i+1]; 
    end 
end 
if isempty(remove_index)~=1
    phase_0_start(remove_index) = []; 
    rapid_rise_end(remove_index-1) = []; 
end 

upstroke_time_points = phase_0_start; 
rapid_depolar_end = rapid_rise_end; 
   
if plot_index ==1
   figure; 
   hold on; 
   plot(trace, 'b'); 
   plot(upstroke_time_points, trace(upstroke_time_points),'*'); 
   plot(rapid_depolar_end, trace(rapid_depolar_end),'o'); 
   hold off; 
end 

