
function [APD_time_point,Threshold] = Get_CaTD_time_point (APD_number, upstroke_time_point,rapid_depolar_end,peakvalue, baseline, frameperiod,poly_coef) 
%% This function find the point for APD?? and APD?? (say APD50, APD20 or APD90)
% This function takes a single AP trace and the user need to input which
% APD he watns 
% For the fitted trace we need to know the intersection points of the
% fitted trace with the threshold line 
%% 
Threshold = peakvalue-APD_number/100*(peakvalue-baseline);

% too comput the intersection points of this fitted polynomial with the
% threshold 
poly_coef(end) = poly_coef(end)-Threshold;% compute the constant term 
if isnan(poly_coef)~=1
    root_vector = roots(poly_coef);  % find all the roots 
    root_vector = sort(root_vector);
    for i = 1 :length(root_vector) 
        if isreal(root_vector(i))==0
            root_vector(i)=0; 
        end 
    end 
    
    APD_time_index = find(root_vector >(rapid_depolar_end-upstroke_time_point)*frameperiod,1,'first');
    APD_time_point =root_vector(APD_time_index); 
        
else 
    APD_time_point = 0; 
end 

% This is designed for the cases where the APD80 or APD50 are no found because of funny cubic shape which bends away from the threshold line in the bad quality case 

% if isempty(APD_time_point)==1
%     APD_time_point = rapid_depolar_end-upstroke_time_point+find(row_trace(rapid_depolar_end-upstroke_time_point:end)>Threshold,1,'last'); % find the neareast( towardas +inifinity) where the filtered row trace meet the treshold 
% end 

