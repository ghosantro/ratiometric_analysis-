function [upstroke_time_points,rapid_depolar_end] = Get_rough_upstroke(trace,frameperiod, plot_index,row_trace)
%% This function get roughly the upstroke for normalized trace 

% Normalize the trace first 
trace_original = trace; 
row_trace_original = row_trace; 
trace= trace-min(trace); 
trace = trace./max(trace);
row_trace= row_trace-min(row_trace); 
row_trace = row_trace./max(row_trace);


gradient = (trace(2:end)-trace(1:end-1))./frameperiod; % get first deravative 


index_depolar = find (gradient>0.0045); % Capture the fast depolarisation phase 
subs_index = index_depolar(2:end)-index_depolar(1:end-1);
discontinious_index = find(subs_index>1)+1; %distinguish depolarisations (to get rid of littel humps set the discontinous index to 10 instead of 1)index=10 now   
if isempty(index_depolar)~=1
    phase_0_start = [index_depolar(1),index_depolar(discontinious_index)];

    %This needs to be further filtered by a big rise in aplitude 
    rapid_rise_end = [index_depolar(discontinious_index-1),index_depolar(end)];

    rise_amplitude = trace(rapid_rise_end)-trace(phase_0_start); 
    if plot_index ==1
        figure; 
       hist(rise_amplitude) 
    end 
    significant_rise_index = find(rise_amplitude>0.35 );   %IF THE APD MAP HAS LOTS OF HOLES CHANGE THE THRESHOLD HERE changed from 0.25 to 0.20; 0.15 works well
    phase_0_start = phase_0_start(significant_rise_index);
    rapid_rise_end = rapid_rise_end(significant_rise_index);

    % if the two activation points are found to be too close to in each other
    % remove the later one and this must be paired up with one false rapid rise
    % end being found, therefore remove the first one from the rapid rise
    % vector 
    % if  pacing_fequency ~=0 % this means the pacing frequency is known 
    %    remove_cutoff = 1000/ pacing_fequency/1.96*0.6;
    % else 
    %     remove_cutoff = 80; 
    % end 
    %    
    % remove_index = []; 
    % for i = 1: length(phase_0_start)-1
    %     if phase_0_start(i+1)-phase_0_start(i)<remove_cutoff 
    %        remove_index = [remove_index,i+1]; 
    %     end 
    % end 
    % if isempty(remove_index)~=1
    %     phase_0_start(remove_index) = []; 
    %     rapid_rise_end(remove_index) = []; 
    % end 

    upstroke_time_points = phase_0_start; 
    rapid_depolar_end = rapid_rise_end; 
   
    if plot_index ==1
       figure; 
       hold on; 
       plot(1:frameperiod:frameperiod*length(row_trace_original),row_trace_original, 'b'); 
       plot(1:frameperiod:frameperiod*length(trace_original),trace_original,'r');
       plot(upstroke_time_points*frameperiod, trace_original(upstroke_time_points),'*'); 
       plot(rapid_depolar_end*frameperiod, trace_original(rapid_depolar_end),'o'); 
       xlabel('time (ms)') ; 
       ylabel('normalized Vm'); 
       legend('row trace','filtered trace','onset upstroke','upstroke end'); 
       hold off; 
       pause
    end 
else 
    upstroke_time_points = []; 
    rapid_depolar_end = []; 
end 
