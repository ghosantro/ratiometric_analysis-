function [fitted_Vm,p_depolar,further_corrected_upstroke] = fitt_AP_upstroke(trace,row_trace,upstroke_point,rapid_rise_end, baseline,frameperiod, plot_index,signal_swing) 
%% This function fit a straight line to the AP upstroke phase and further correct the onset of upstroke time  

% input: trace of ONE AP (this is the trace with high frequency noise filtered out) 
%        row trace 
%        index of first upstroke point (after first correction) 
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
% get the fitting points 

start_fitting_index = round((rapid_rise_end-upstroke_point)/4*2);% choose an inter median point for the starting point use 2 thrids way up
end_fitting_index = round((rapid_rise_end-upstroke_point)/5*4);
T_straight = (start_fitting_index:end_fitting_index).*frameperiod; %This is the coordinate on which the fitting will be done (the reference here is the old upstroke point: set as 1)

% fit a straight line for the upstroke phase 
if signal_swing >=0.05
    
    [p_depolar, s_depolar] = polyfit(T_straight, row_trace(start_fitting_index:end_fitting_index),1);
    
    if p_depolar(1) > 0 
    
    % find the intersection point of the stratight line with the baseline for
    % to further correct the onset of upstroke 

        fitted_Vm = row_trace; 
        linear_coef = p_depolar;
        linear_coef(end) = p_depolar(end)-baseline; 
        root_vector_depolar = roots(linear_coef);
        if root_vector_depolar > 0 % i.e. there's some correction done 
           further_corrected_upstroke = round(root_vector_depolar/frameperiod+upstroke_point); 
           % correct the upstroke part of the trace 
           if further_corrected_upstroke-upstroke_point>0
              T_fit_linear = (further_corrected_upstroke-upstroke_point:rapid_rise_end-upstroke_point).*frameperiod; 
              fitted_Vm(further_corrected_upstroke-upstroke_point:rapid_rise_end-upstroke_point) = polyval(p_depolar,T_fit_linear);
           else 
              T_fit_linear = (1:rapid_rise_end-upstroke_point).*frameperiod; 
              fitted_Vm(1:rapid_rise_end-upstroke_point) = polyval(p_depolar,T_fit_linear);
           end
        else % there's no further correction on the upstroke point we reject the fitting 
            further_corrected_upstroke =upstroke_point; 
        end


        if plot_index ==1
            figure; 
            hold on; 
            plot(0:frameperiod:frameperiod*(length(trace)-1),trace,'r');
            plot(0:frameperiod:frameperiod*(length(trace)-1),row_trace,'k'); 
            plot(0:frameperiod:frameperiod*(length(trace)-1),fitted_Vm,'b');
            plot(0:frameperiod:frameperiod*(length(trace)-1),ones(1,length(trace)).*baseline,'k--')
             plot(start_fitting_index*frameperiod, row_trace(start_fitting_index),'ro'); 
            plot(end_fitting_index*frameperiod, row_trace(end_fitting_index),'r+');
            %plot(intersection_point*frameperiod, trace(intersection_point),'ro'); 
            %plot(baseline_intersection_index*frameperiod, fitted_Vm(baseline_intersection_index),'r*');
           
            if root_vector_depolar > 0
                if further_corrected_upstroke-upstroke_point==0
                   plot(floor((further_corrected_upstroke-upstroke_point+1))*frameperiod, row_trace(floor(further_corrected_upstroke-upstroke_point+1)),'r*');
                else
                    plot(floor((further_corrected_upstroke-upstroke_point))*frameperiod, row_trace(floor(further_corrected_upstroke-upstroke_point)),'r*');
                end 
            end

            xlabel('time(ms)'); 
            ylabel('F/F0'); 
            legend('filtered trace','row trace','fitted trace','baseline','start fitting point','end fitting point','corrected upstroke point');%'intersection point','poly intersection'); 

            hold off; 
        end  
       
    else 
        fitted_Vm = trace; 
        p_depolar = zeros(1,2); 
        further_corrected_upstroke = upstroke_point; 
    end 
else 
    [p_depolar, s_depolar] = polyfit(T_straight, trace(start_fitting_index:end_fitting_index),1);
       
    if p_depolar(1) > 0 
    
    % find the intersection point of the stratight line with the baseline for
    % to further correct the onset of upstroke 

        fitted_Vm = row_trace; 
        linear_coef = p_depolar;
        linear_coef(end) = p_depolar(end)-baseline; 
        root_vector_depolar = roots(linear_coef);
        if root_vector_depolar > 1 % i.e. there's some correction done 
           further_corrected_upstroke = round(root_vector_depolar/frameperiod+upstroke_point); 
           % correct the upstroke part of the trace 
           if further_corrected_upstroke-upstroke_point>0
               T_fit_linear = (further_corrected_upstroke-upstroke_point:rapid_rise_end-upstroke_point).*frameperiod; 
               fitted_Vm(further_corrected_upstroke-upstroke_point:rapid_rise_end-upstroke_point) = polyval(p_depolar,T_fit_linear);
           else
               T_fit_linear = (1:rapid_rise_end-upstroke_point).*frameperiod; 
               fitted_Vm(1:rapid_rise_end-upstroke_point) = polyval(p_depolar,T_fit_linear);
           end
        else % there's no further correction on the upstroke point we reject the fitting 
            further_corrected_upstroke =upstroke_point; 
        end

        
        if plot_index ==1
            figure; 
            hold on; 
            plot(0:frameperiod:frameperiod*(length(trace)-1),trace,'r');
            plot(0:frameperiod:frameperiod*(length(trace)-1),row_trace,'k'); 
            plot(0:frameperiod:frameperiod*(length(trace)-1),fitted_Vm,'b');
            plot(0:frameperiod:frameperiod*(length(trace)-1),ones(1,length(trace)).*baseline,'k--')
            %plot(intersection_point*frameperiod, trace(intersection_point),'ro'); 
            %plot(baseline_intersection_index*frameperiod, fitted_Vm(baseline_intersection_index),'r*');
           
            plot(start_fitting_index*frameperiod, row_trace(start_fitting_index),'ro'); 
            plot(end_fitting_index*frameperiod, row_trace(end_fitting_index),'r+');
            if root_vector_depolar > 1
                plot((further_corrected_upstroke-upstroke_point)*frameperiod, trace(further_corrected_upstroke-upstroke_point),'r*');
            end

            xlabel('time(ms)'); 
            ylabel('F/F0'); 
            legend('filtered trace','row trace','fitted trace','baseline','start fitting index','end of fitting index','corrected upstroke_point');%'intersection point','poly intersection'); 

            hold off; 
        end  
       
    else 
        fitted_Vm = trace; 
        p_depolar = zeros(1,2); 
        further_corrected_upstroke = upstroke_point; 
    end 
end 
    