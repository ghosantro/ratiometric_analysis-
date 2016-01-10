function [fitted_Vm,p_repolar,baseline_intersection_index] = fitt_AP_curve(trace,row_trace,upstroke_point,rapid_rise_end, baseline,frameperiod, plot_index,signal_swing) 
%% This function fit curve to the AP repolarisation phase  

% input: trace of ONE AP (this is the trace with high frequency noise filtered out) 
%        row trace 
%        index of first upstroke point (after second correction) 
%        end of rapid depolarisation phase point (after correction) 
%        baseline value of this AP 
%        frame period 
%        1 = want a plot 
%        0 = don't want a plot 
%        signal swing for this AP (estimates how good the AP signal is)

% Output: fitted trace    

%         coefficient of the fitting 
%         further corrected upstroke time point index 

%% Code 

% we want to find the first point it's back from repolarisation_threshold 
% Find the intersection of the baseline with the repolarisation trace 
peak_value = max(trace);
 
if rapid_rise_end-upstroke_point<1
   intersection_point = find(trace(1:end)<baseline+abs(peak_value-baseline)*0.05,1,'first'); 
   %disp('warning: this upstroke point is not well identified'); 
else 
      
    intersection_point = rapid_rise_end-upstroke_point+find(trace(rapid_rise_end-upstroke_point+1:end)<baseline+abs(peak_value-baseline)*0.05,1,'first');
end 

if isempty(intersection_point)==1
    intersection_point = rapid_rise_end-upstroke_point+find(trace(rapid_rise_end-upstroke_point+1:end)<baseline+abs(peak_value-baseline)*0.15,1,'first');
    %disp('difficult to find intersection point')
end 

% for the repolarisation phase a cubic curve is fitted to the data
T_repolar = (rapid_rise_end-upstroke_point:intersection_point).*frameperiod;
if signal_swing >=0.05
    if rapid_rise_end-upstroke_point>=1
        [p_repolar, s_repolar] = polyfit(T_repolar, row_trace(rapid_rise_end-upstroke_point:intersection_point),4);
    else
        T_repolar =(1:intersection_point).*frameperiod; 
        [p_repolar, s_repolar] = polyfit(T_repolar, row_trace(1:intersection_point),4);
    end 
else 
    if rapid_rise_end-upstroke_point>=1
        [p_repolar, s_repolar] = polyfit(T_repolar, trace(rapid_rise_end-upstroke_point:intersection_point),4); 
    else 
        T_repolar =(1:intersection_point).*frameperiod; 
        [p_repolar, s_repolar] = polyfit(T_repolar, trace(1:intersection_point),4); 
    end
end 

%look for hte intersection of this polynomial with the baseline 
cubic_coef = p_repolar;
cubic_coef(end) = p_repolar(end)-baseline; 
root_vector = roots(cubic_coef);
%set the complex roots to zero
root_vector = sort(root_vector);
for i = 1 :length(root_vector) 
    if isreal(root_vector(i))==0
        root_vector(i)=0; 
    end 
end 
%root_vector
baseline_intersection_time = root_vector(find(root_vector >(rapid_rise_end-upstroke_point)*frameperiod,1,'first'));
if isempty(baseline_intersection_time)~=1&&baseline_intersection_time>frameperiod &&baseline_intersection_time <length(trace)*frameperiod 
    baseline_intersection_index = round(baseline_intersection_time/frameperiod);
else 
    baseline_intersection_index = intersection_point; 
end 
fit_index = [rapid_rise_end, upstroke_point+baseline_intersection_index];
% evaluate the voltage after curve fitting 
fitted_Vm = row_trace; 
if rapid_rise_end-upstroke_point>=1
    T_fit = (rapid_rise_end-upstroke_point:baseline_intersection_index).*frameperiod;
    fitted_Vm (rapid_rise_end-upstroke_point:baseline_intersection_index) = polyval(p_repolar, T_fit);
else 
    T_fit = (1:baseline_intersection_index).*frameperiod; 
    fitted_Vm (1:baseline_intersection_index) = polyval(p_repolar, T_fit);
end 


fitted_Vm (baseline_intersection_index:end) = baseline ;

if plot_index ==1
    figure; 
    hold on; 
    plot(0:frameperiod:frameperiod*(length(trace)-1),trace,'r');
    plot(0:frameperiod:frameperiod*(length(trace)-1),row_trace,'k'); 
    plot(0:frameperiod:frameperiod*(length(trace)-1),fitted_Vm,'b');
    plot(0:frameperiod:frameperiod*(length(trace)-1),ones(1,length(trace)).*baseline,'k--')
    plot(round(intersection_point)*frameperiod, trace(round(intersection_point)),'ro'); 
    plot(baseline_intersection_index*frameperiod, fitted_Vm(baseline_intersection_index),'r*'); 
      
    xlabel('time(ms)'); 
    ylabel('F/F0'); 
    legend('filtered trace','row trace','fitted trace','baseline');%'intersection point','poly intersection'); 
    
    hold off; 
end  



