function [bined_image] = Taking_moving_average (image,num_pixels_to_bin)
%% This function takes a square image and taking the moving averaged of
% area as the signal value of the new pixel 

%% 
original_size_row = size(image,1);
original_size_col = size(image,2); 
original_time_length = size(image,3); 
bined_image = zeros(original_size_row-num_pixels_to_bin+1,original_size_col-num_pixels_to_bin+1, original_time_length); 
for i = 1: size(bined_image,1)
    for j = 1: size(bined_image,2) 
        bined_image(i,j,:) = mean(mean(image(i:i+num_pixels_to_bin-1,j:j+num_pixels_to_bin-1,:))); 
    end
end

        
