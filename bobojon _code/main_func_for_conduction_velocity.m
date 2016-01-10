function conduction_struct=main_func_for_conduction_velocity(upstroke_matrix,mean_upstroke_matrix,save_figure_index,fname,file_name_start_index,...
    file_name_end_index,known_pacing, pacing_frequency,pos, frameperiod, total_row_number, total_col_number, result_path)
%% This is the main function for finding conduction velocity 

%% code 
%  filter upstroke time for each beat 

% filtered_upstroke_matrix = filter_upstroke_onset (upstroke_cell, number_upstroke_matrix, known_pacing, pacing_frequency, save_figure_index,fname, file_name_start_index,...
%  file_name_end_index,fibrillation_index,signal_swing_matrix,frameperiod,result_path); 
upstroke_odd= upstroke_matrix(:,:,1:2:end); 
upstroke_even= upstroke_matrix(:,:,2:2:end);
% analyze the variance of the upstroke time point 
std_matrix_odd = std(upstroke_odd,0,3); 
mean_std_odd = mean(std_matrix_odd(find(isnan(std_matrix_odd)~=1)));

std_matrix_even = std(upstroke_even,0,3); 
mean_std_even = mean(std_matrix_even(find(isnan(std_matrix_even)~=1)));

% use the simple method to estimate conduction velocity 
[conduction_velocity_matrix_odd, speed_matrix_odd, rough_angel_matrix_odd] = evaluate_conduction_velocity_simple(mean_upstroke_matrix(:,:,1), 1, 0);
rough_mean_speed_odd = mean(speed_matrix_odd([find(speed_matrix_odd>0);find(speed_matrix_odd<5000)])); 
[conduction_velocity_odd,fitting_error_matrix,conduction_speed_matrix_odd,angel_matrix_odd] = evaluate_conduction_velocity_Bayly(mean_upstroke_matrix(:,:,1), rough_mean_speed_odd,1, mean_std_odd, 3,0);

[conduction_velocity_matrix_even, speed_matrix_even, rough_angel_matrix_even] = evaluate_conduction_velocity_simple(mean_upstroke_matrix(:,:,2), 1, 0);
rough_mean_speed_even = mean(speed_matrix_odd([find(speed_matrix_even>0);find(speed_matrix_even<5000)])); 
[conduction_velocity_even,fitting_error_matrix,conduction_speed_matrix_even,angel_matrix_even] = evaluate_conduction_velocity_Bayly(mean_upstroke_matrix(:,:,2), rough_mean_speed_even,1, mean_std_even, 3,0);

plot_conduction_vector_map (mean_upstroke_matrix(:,:,1),conduction_velocity_odd,save_figure_index,fname,file_name_start_index,file_name_end_index, total_row_number, ...
    total_col_number, pos,result_path, 1)

plot_conduction_vector_map (mean_upstroke_matrix(:,:,2),conduction_velocity_even,save_figure_index,fname,file_name_start_index,file_name_end_index, total_row_number, ...
    total_col_number, pos,result_path, 0)
%plot_cleaned_conduction_velocity_map (filtered_upstroke_matrix(:,:,1),conduction_velocity,3,126)
% use Bayly's method to estimate conduction velocity 
%speed_mean = mean(speed_matrix(find(speed_matrix >0))); 

conduction_struct = struct('conduction_veolocity_odd',conduction_velocity_odd, 'conduction_velocity_even',conduction_velocity_even,'conduction_speed_odd',...
conduction_speed_matrix_odd,'conduction_speed_even',conduction_speed_matrix_even,'angel_odd',angel_matrix_odd,'angel_even',angel_matrix_even); 
