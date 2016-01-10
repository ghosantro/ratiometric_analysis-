function [APD,peak_value] = Find_all_CaTD(trace,row_trace, upstroke_point, rapid_depolar_end, baseline,frameperiod,poly_coef, pacing_interval, plot_index)
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

[APD50_time_point,threshold_50] = Get_CaTD_time_point (50, upstroke_point,rapid_depolar_end,peak_value, baseline, frameperiod, poly_coef);
[APD80_time_point,threshold_80] = Get_CaTD_time_point (80, upstroke_point,rapid_depolar_end,peak_value, baseline, frameperiod, poly_coef);
[APD90_time_point,threshold_90] = Get_CaTD_time_point (90, upstroke_point,rapid_depolar_end,peak_value, baseline, frameperiod, poly_coef);
not_from_poly_80=0;
if isempty(APD80_time_point)==1
    if rapid_depolar_end-upstroke_point>0
        
        APD80_time_index= find(trace(rapid_depolar_end-upstroke_point:end)<threshold_80);
        if isempty(APD80_time_index)==0
            APD80_time_point = min(APD80_time_index)+ (rapid_depolar_end-upstroke_point); 
            not_from_poly_80=1;
        else
            APD80_time_point= pacing_interval; 
        end 
    else 
        
         APD80_time_point = 0; %find(trace(1:end),threshold_80,'first'); 
    end       
   
else 
    not_from_poly_80=2; 
end

not_from_poly_50=0;
if isempty(APD50_time_point)==1   
   if rapid_depolar_end-upstroke_point>0
       APD50_time_index= find(trace(rapid_depolar_end-upstroke_point:end)<threshold_50);
       if isempty(APD50_time_index)==0       
          APD50_time_point = min(APD50_time_index)+ (rapid_depolar_end-upstroke_point); 
          not_from_poly_50=1;
       else 
          APD50_time_point=pacing_interval; 
       end 
    else 
        APD50_time_point = 0;%find(trace(1:end),threshold_50,'first'); 
    end 
   
   
else 
    not_from_poly_50=2; 
end 


not_from_poly_90=0;
if isempty(APD90_time_point)==1
    if rapid_depolar_end-upstroke_point>0
        
        APD90_time_index= find(trace(rapid_depolar_end-upstroke_point:end)<threshold_90);
        if isempty(APD90_time_index)==0
            APD90_time_point = min(APD90_time_index)+ (rapid_depolar_end-upstroke_point); 
            not_from_poly_90=1;
        else
            APD90_time_point= pacing_interval; 
        end 
    else 
        
         APD90_time_point = 0; %find(trace(1:end),threshold_80,'first'); 
    end       
   
else 
    not_from_poly_90=2; 
end


%More APD can be added 

% Check with a plot 
if plot_index ==1 
    figure;
    hold on; 
   plot(0:frameperiod:frameperiod*(length(trace)-1),trace,'bo');
   if not_from_poly_50==2
      plot(APD50_time_point, polyval(poly_coef, APD50_time_point),'r*'); 
   elseif not_from_poly_50==1
      plot(APD50_time_point,trace(floor(APD80_time_point/frameperiod)),'r*');
   else 
      plot(APD50_time_point, threshold_50,'r*'); 
   end 
   if not_from_poly_80==2
      plot(APD80_time_point, polyval(poly_coef, APD80_time_point),'ro');
   elseif not_from_poly_80==1
      plot(APD80_time_point,trace(floor(APD80_time_point/frameperiod)),'ro');
   else 
      plot(APD80_time_point, threshold_80,'ro');
   end 
   
   if not_from_poly_90==2
      plot(APD90_time_point, polyval(poly_coef, APD90_time_point),'ro');
   elseif not_from_poly_90==1
      plot(APD90_time_point,trace(floor(APD90_time_point/frameperiod)),'ro');
   else 
      plot(APD90_time_point, threshold_90,'ro');
   end 
   
   plot(0:frameperiod:frameperiod*(length(trace)-1),ones(1,length(trace)).*peak_value,'k');
   plot(0:frameperiod:frameperiod*(length(trace)-1),ones(1,length(trace)).*threshold_80,'k');
   plot(0:frameperiod:frameperiod*(length(trace)-1),ones(1,length(trace)).*threshold_50,'k');
   plot(0:frameperiod:frameperiod*(length(trace)-1),ones(1,length(trace)).*threshold_90,'k');
   plot(0:frameperiod:frameperiod*(length(row_trace)-1),row_trace,'k*');
   xlabel('time (ms)'); 
   ylabel('vol sig'); 
   legend('trace','APD50','APD80','APD90','peakvalue','threshold 80','threshold 50','threshold 90'); 
   hold off; 
end 
   

APD = struct('APD80',APD80_time_point,'APD50',APD50_time_point,'APD90',APD90_time_point); 


