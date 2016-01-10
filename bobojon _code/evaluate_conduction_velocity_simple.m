function [conduction_velocity_matrix, speed_matrix, angel_matrix] = evaluate_conduction_velocity_simple(upstroke_matrix, spatial_callibration, plot_index) 
%% This function evaluates conduction velocity by simply adding time velocity along x axis and velocity along y axis 
% This function takes in: upstroke time matix and spatial callibration
% value (i.e.the real distance between two pixel) 

%% code 
conduction_velocity_matrix = zeros(size(upstroke_matrix,1), size(upstroke_matrix,2), 2);
speed_matrix = NaN(size(upstroke_matrix,1), size(upstroke_matrix,2));
angel_matrix = NaN(size(upstroke_matrix,1), size(upstroke_matrix,2));
 
if plot_index ==1
    figure; 
    hold on; 

    contourf( upstroke_matrix,10); colorbar
end 
for r = 2: size(upstroke_matrix,1)-1
    for c = 2:size(upstroke_matrix,2)-1
        %only find action potential if we have information about upstroke time at its four neighbours  
        
        if isnan(upstroke_matrix(r-1,c))~=1&&isnan(upstroke_matrix(r+1,c))~=1&&isnan(upstroke_matrix(r,c-1))~=1&&isnan(upstroke_matrix(r,c+1))~=1
           
           time_delay_x = upstroke_matrix(r,c+1)-upstroke_matrix(r,c-1);
           time_delay_y = upstroke_matrix(r+1,c)-upstroke_matrix(r-1,c);
           velocity_x = spatial_callibration/time_delay_x; 
           velocity_y = spatial_callibration/time_delay_y; 
           conduction_speed = sqrt(velocity_x^2+velocity_y^2); 
           angel = atan2(velocity_x, velocity_y);
        else 
           velocity_x = zeros; 
           velocity_y = zeros;
           conduction_speed = NaN; 
           angel= NaN;
        end 
        conduction_velocity_matrix(r,c,:) = [velocity_x,velocity_y]; 
        speed_matrix(r,c) = conduction_speed; 
        angel_matrix(r,c) = angel; 
        if plot_index==1&& isnan(velocity_x)~=1&&isnan(velocity_y)~=1
            [arrowx, arrowy] = plot_vector(c,r,velocity_x,velocity_y,1);
            plot(arrowx,arrowy); 
        end 
        
    end 
end 


% if plot_index ==1
%     quiver(
%     hold off;
% end

     
     
       

