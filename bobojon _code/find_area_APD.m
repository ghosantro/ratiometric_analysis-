function [area_under_AP,area_ratio]= find_area_APD(AP_trace,baseline,frame_period,baseline_intresection,APD80, peak_value,plot_index )
%% This function evaluates area under the action potential curve and hence finds APD50 and APD90 according to the under-curve area 
% This function use trapezium rule to estimate the area under curve
% This function takes:fitted trace of one action potential(i.e. from the further corrected upstroke point to the next upstroke point)
%                     
% 
%% code 

accumulate_area_vector = []; 
if APD80 ~=0  %only look for area under the curve for well fitted curve 
    if floor(APD80/frame_period)>2 % only find the area if AP is long enough to use trapezium rule 
       % disp('should get area')
        for i = 2:baseline_intresection
            area = (min(AP_trace(i),peak_value)-baseline+min(AP_trace(i-1),peak_value)-baseline)*frame_period/2; %trapezium 
            if isempty(accumulate_area_vector)~=1
               accumulate_area_vector = [accumulate_area_vector, area+accumulate_area_vector(end)]; 
            else 
               accumulate_area_vector= [accumulate_area_vector,area];
            end 
        end 
    else 
        accumulate_area_vector =0; 
        disp('this AP is too short to find area')
        area_ratio=0; 
    end 

    
    if isempty(accumulate_area_vector)~=1
        area_under_AP = accumulate_area_vector(end); 
    else 
        area_under_AP =0; 
        disp('no area found ')
        area_ratio=0;
    end
     
    square_area = (peak_value-baseline)*(baseline_intresection-1)*frame_period;
    area_ratio = area_under_AP/square_area; 
    if area_ratio >1
        disp('there is something wrong in area under the curve') 
    end 
    if plot_index ==1 
        T1 = 0:frame_period:(baseline_intresection-1)*frame_period;
        
        %T1 = 0.1:frame_period:baseline_intresection*frame_period;% for the test
        figure; 
        hold on; 
        
        plot(frame_period.*[1, baseline_intresection], [baseline, baseline], 'r');
        plot(frame_period.*[1, baseline_intresection], [peak_value, peak_value], 'r');
        plot(frame_period.*[1, 1], [baseline, peak_value], 'r');
        plot(frame_period.*[baseline_intresection, baseline_intresection], [baseline, peak_value], 'r');
        plot(T1,AP_trace(1:baseline_intresection),'b'); 
%         plot(T1,ones(1,baseline_intresection)*baseline ,'r');
%         plot(T1,ones(1,baseline_intresection)*peak_value ,'r');
        
        hold off; 
    end  
        
else 
    area_under_AP = 0; 
    area_ratio = 0; 
end 

      