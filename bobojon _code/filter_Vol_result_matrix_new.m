
function filtered_signal_matrix = filter_Vol_result_matrix_new(signal_vol, signal_swing_matrix, upstroke_index) 
%% This function filter out the outliers of the median matrix 

%% code
filtered_signal_matrix = signal_vol; 
%length(find(isnan(filtered_median_matrix)~=1))
% figure; 
% contourf(filtered_median_matrix,10);colorbar 
% disp('before any filtering')
%filter with the signal swing 
filtered_signal_matrix(find(signal_swing_matrix)<=0.003)=NaN;
% length(find(isnan(filtered_median_matrix)~=1))
% figure; 
% contourf(filtered_median_matrix,10);colorbar 
% disp('signal_swing filtering')
% pause
% close all; 
%filter with AP number matrix (optional)
% length(find(isnan(filtered_median_matrix)~=1))
% figure; 
% contourf(filtered_median_matrix,10); colorbar
% disp('standard_AP filtering')
% pause
% close all
median_median = median(filtered_signal_matrix(find(isnan(filtered_signal_matrix)~=1)));

% length(find(isnan(filtered_median_matrix)~=1))
% figure; 
% contourf(filtered_median_matrix,10); colorbar
% disp('std over time filtering')
% pause
% close all

% filter with itself with quartile method 

upper_quartile =  quantile(filtered_signal_matrix(find(isnan(filtered_signal_matrix)~=1)),0.75);
lower_quartile = quantile(filtered_signal_matrix(find(isnan(filtered_signal_matrix)~=1)),0.25);
if upstroke_index ~=1
    upper_treshold = 2.5*(upper_quartile-lower_quartile)+median_median;
    lower_treshold = -2.5*(upper_quartile-lower_quartile)+median_median;
else 
    upper_treshold = 2*(upper_quartile-lower_quartile)+median_median;
    lower_treshold = -1.5*(upper_quartile-lower_quartile)+median_median;
end 

 filtered_signal_matrix(find(filtered_signal_matrix>upper_treshold))=NaN;
 filtered_signal_matrix(find(filtered_signal_matrix<lower_treshold))=NaN;
% length(find(isnan(filtered_median_matrix)~=1))
% figure; 
% contourf(filtered_median_matrix,10); colorbar
% disp('quartile filtering')
% pause
% close all

% if upstroke_index ==1
%    filtered_signal_matrix = filtered_signal_matrix-min(min(filtered_signal_matrix)); 
% end 

