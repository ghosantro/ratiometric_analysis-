function [APD,peak_value] = Find_all_APD(trace,upstroke_point, rapid_depolar_end, baseline,frameperiod,poly_coef, plot_index)
%% This function calculate APD with 95 percentile as the peak value 

%  This function takes in: 
%  Trace for one action potential (fitted prefered) 
%  Upstroke time point for this action potential 
%  End of uprstoke point for this action potential 
%  Baseline (holding potential value) 
%  plot_index: 1 = plot, 0 = no plot 
%% Code 
% find the 'peak' value 
peak_value = quantile(trace,0.95);  

[APD50_time_point,threshold_50] = Get_APD_time_point (50, upstroke_point,rapid_depolar_end,peak_value, baseline, frameperiod, poly_coef);
[APD80_time_point,threshold_80] = Get_APD_time_point (80, upstroke_point,rapid_depolar_end,peak_value, baseline, frameperiod, poly_coef); 

%More APD can be added 

% Check with a plot 
if plot_index ==1 
    figure;
    hold on; 
   plot(0:frameperiod:frameperiod*(length(trace)-1),trace,'bo'); 
   plot(APD50_time_point, polyval(poly_coef, APD50_time_point),'r*'); 
   plot(APD80_time_point, polyval(poly_coef, APD80_time_point),'ro');
   plot(0:frameperiod:frameperiod*(length(trace)-1),ones(1,length(trace)).*peak_value,'k');
   plot(0:frameperiod:frameperiod*(length(trace)-1),ones(1,length(trace)).*threshold_80,'k');
   plot(0:frameperiod:frameperiod*(length(trace)-1),ones(1,length(trace)).*threshold_50,'k');
   xlabel('time (ms)'); 
   ylabel('vol sig'); 
   legend('trace','APD50','APD80','peakvalue','threshold 80','threshold 50'); 
   hold off; 
end 
   

APD = struct('APD80',APD80_time_point,'APD50',APD50_time_point); 


