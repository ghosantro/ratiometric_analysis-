function [alternans_index,alternans_score] = test_alternans(peak_value, upstroke_time_point, APD80, AP_area, known_pacing)
%% This function looks for alternans in the trace 
%  This function takes: fitted action potential trace 
%                       upstroke time point (further corrected)  
%                       APD80
%                       area under curve for each action potential 
%                       known_pacing index: 1 = with pacing 0 = unkown pacing            

%  This function ouputs: alterans_index 1 = has alternansn, 0 = no alternans 
%% Code 

% find inter upstroke interval time outlier and exclude it if pacing
% frequency is known 


if known_pacing ==1 
    inter_upstroke_time= upstroke_time_point(2:end)-upstroke_time_point(1:end-1); 
    mean_upstroke_interval = mean(inter_upstroke_time);
    std_upstroke_interval = std(inter_upstroke_time);

    if std_upstroke_interval<5
        outlier_index = find(inter_upstroke_time>mean_upstroke_interval+8*std_upstroke_interval);
    elseif std_upstroke_interval>=5&& std_upstroke_interval<25
        outlier_index = find(inter_upstroke_time>mean_upstroke_interval+2.5*std_upstroke_interval);
    else
        outlier_index = find(inter_upstroke_time>mean_upstroke_interval+1*std_upstroke_interval);
    end 
        
   
    upstroke_time_point(outlier_index) = []; 
    peak_value(outlier_index) = []; 
    APD80(outlier_index) = []; 
    AP_area(outlier_index) = []; 
end 

% remove the beats with no APD 
no_APD80 = find(APD80==0);

if isempty(find(no_APD80==1))==0&&isempty(APD80)~=1
    
    no_APD80(1)=[]; 
    upstroke_time_point(1) = []; 
    APD80(1) = []; 
    peak_value(1)= [];
    AP_area(1) = [];
    no_APD80= no_APD80-1; % as the first one is removed
end 
no_APD80= no_APD80(no_APD80>0);
need_to_remove_index = [no_APD80,no_APD80+1];
need_to_remove_index  = need_to_remove_index(need_to_remove_index<=length(upstroke_time_point)) ;


APD80(need_to_remove_index) = []; 
AP_area(need_to_remove_index) = []; 
peak_value(need_to_remove_index)= []; 

 %alternans analysis 

%taking the ratio

APD_ratio = APD80(2:end)./APD80(1:end-1); 
peak_ratio = peak_value(2:end)./peak_value(1:end-1);
AP_area_ratio = AP_area(2:end)./AP_area(1:end-1);


% only do alternanting test if we have more than 4 of action potential left
% allow 0.001 error here 
if length(APD_ratio)>=4
    
    %check for alternating pattern for each of these three paramters 
    if sum(find(APD_ratio(1:2:end)>1+0.001))==0&&sum(find(APD_ratio(2:2:end)<1-0.001))==0
        APD_alternans_index = 1; 
    elseif sum(find(APD_ratio(1:2:end)<1-0.001))==0&&sum(find(APD_ratio(2:2:end)>1+0.001))==0
        APD_alternans_index = 1; 
    else 
        APD_alternans_index = 0 ; 
    end 
    
    if sum(find(peak_ratio(1:2:end)>1))==0&&sum(find(peak_ratio(2:2:end)<1))==0
        peak_alternans_index = 1; 
    elseif sum(peak_ratio(1:2:end)<1)==0||sum(peak_ratio(2:2:end)>1)==0
        peak_alternans_index = 1; 
    else 
        peak_alternans_index = 0 ; 
    end 
    
    if sum(AP_area_ratio(1:2:end)>1+0.001)==0&&sum(AP_area_ratio(2:2:end)<1-0.001)==0
        area_alternans_index = 1; 
    elseif sum(AP_area_ratio(1:2:end)<1+0.001)==0&&sum(AP_area_ratio(2:2:end)>1-0.001)==0
        area_alternans_index = 1; 
    else 
        area_alternans_index = 0 ; 
    end 
    
    % only if all 3 are alternating, AP is alternating 
    if APD_alternans_index==1||area_alternans_index==1 %peak_alternans_index==1||
        alternans_index =1;
        alternans_score = sum(abs(APD_ratio-1))./length(APD_ratio); %find the level of alternans (if that's too small then may not be reliable alternans)
    else 
        alternans_index =0; 
        alternans_score =0; 
    end 
else 
    %disp('not enough AP for alternans analysis') 
    alternans_index=0; 
    alternans_score =0;
end 
        
        
        





     
    
    

