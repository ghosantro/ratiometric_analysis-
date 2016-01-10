function [conduction_velocity,fitting_error_matrix, conduction_speed_matrix,angel_matrix] = evaluate_conduction_velocity_Bayly(upstroke_matrix, rough_mean_speed, spatial_callibration,  mean_std, num_bined,plot_index)
%% This function uses Bayly's method to find conduction velocity
% This function takes in: upstroke time matrix 
%                         average conudction speed estimated from the simple method
%                         rough conudction speed measured from the simple method 
%                         spatial callibraration value (i.e. distance between pixles)
                      
%% Code
if plot_index ==1
   figure; 
   hold on; 
   contourf(upstroke_matrix,10);
end 

row_number = size(upstroke_matrix,1); 
col_number = size(upstroke_matrix,2); 

conduction_velocity = zeros(row_number, col_number,2); 

% define the windows for x and y constrain for the fitting window 

t_window = max(2*mean_std,1); 
x_window_t_constrained = t_window*rough_mean_speed; % not sure if needs to split into x speed and y speed 
y_window_t_constrained = t_window*rough_mean_speed; 

if x_window_t_constrained >5000
   x_window_t_constrained=0; 
end 

if y_window_t_constrained >5000
   y_window_t_constrained=0; 
end 

x_window_min = num_bined*2; 
y_window_min = num_bined*2; 

x_window = round(max(x_window_t_constrained,x_window_min))
y_window = round(max(y_window_t_constrained,y_window_min)) 
fitting_error_matrix = zeros(size(upstroke_matrix,1), size(upstroke_matrix,2)); 
conduction_speed_matrix = zeros(size(upstroke_matrix,1),size(upstroke_matrix,2));
% for each point look for conduction velocity  
for r = 2:row_number-1
    for c = 2:col_number-1 
       
        % find points within x, y constrains 
        x_limit = [max(1,r-x_window),min(row_number, r+x_window)];
        y_limit = [max(1,c-y_window),min(col_number, c+y_window)] ;
        
        area_upstroke = upstroke_matrix(x_limit(1):x_limit(2),y_limit(1):y_limit(2)); 
        
        % check for time constrains 
        
        
        area_upstroke(find(abs(area_upstroke -upstroke_matrix(r,c))>t_window))=NaN; 

     
          
        % crate x, y, t vectors for polynomial fitting 
        time_vector = area_upstroke(find(isnan(area_upstroke)==0));
        % check wether we have enough points after all the constrains to do
        % fitting 
        if length(find(time_vector))>6             
           
            x_index_matrix = zeros(size(area_upstroke,1), size(area_upstroke,2));
            for i = 1:size(area_upstroke,1)
                x_index_matrix(i,:) = (x_limit(1)+i-1).*ones(1,size(area_upstroke,2)); 
            end 

            y_index_matrix = zeros(size(area_upstroke,1), size(area_upstroke,2));
            for i = 1:size(area_upstroke,2)
                y_index_matrix(:,i) = (y_limit(1)+i-1).*ones(size(area_upstroke,1),1); 
            end 
        
            x_vector = x_index_matrix (find(isnan(area_upstroke)==0));
            y_vector = y_index_matrix(find(isnan(area_upstroke)==0)); 

            % polynomial fitting: Aa = t i.e. use T= ax^2+by^2+cxy+dx+ey+f
            % create matrix A
            A = zeros(length(x_vector),6); 
            A(:,1) = x_vector.^2; 
            A(:,2) = y_vector.^2; 
            A(:,3) = x_vector.*y_vector; 
            A(:,4) = x_vector; 
            A(:,5) = y_vector; 
            A(:,6) = ones(length(x_vector),1); 

            [U,S,V] = svd(A,0); 
            a = V*inv(S)*U'*time_vector;
            fitting_error = sqrt(sum((A*a-time_vector).^2));
            fitting_error_matrix(r,c) = fitting_error;
            if fitting_error<30
               
            % estimate condunction velocity from the fitted polynomial surface
            % at point (r,c)

               Tx = 2*a(1)*r+a(3)*c+a(4); 
               Ty = 2*a(2)*c+a(3)*r+a(5); 
               normalise_factor = Tx^2+Ty^2; 
               conduction_velocity(r,c,:) = [Ty/normalise_factor,Tx/normalise_factor].*spatial_callibration;
            end
        end 
         if plot_index==1&& isnan(conduction_velocity(r,c,1))~=1&&isnan(conduction_velocity(r,c,2))~=1
            [arrowx, arrowy] = plot_vector(c,r,conduction_velocity(r,c,1),conduction_velocity(r,c,2),0.5);
            plot(arrowx,arrowy); 
         end 
        
         conduction_speed_matrix (r,c) = sqrt(conduction_velocity(r,c,1)^2+ conduction_velocity(r,c,2)^2);      
                
     end 
end 

% filter out outlier from the conduction speed
mean_conduction_speed = mean(conduction_speed_matrix(find(conduction_speed_matrix>0)))
std_conduction_speed= std(conduction_speed_matrix(find(conduction_speed_matrix>0)))
 
conduction_velocity_frame1 = conduction_velocity(:,:,1); 
conduction_velocity_frame2 = conduction_velocity(:,:,2); 
conduction_velocity_frame1(find(conduction_speed_matrix<mean_conduction_speed-2*std_conduction_speed))=0;       
conduction_velocity_frame1(find(conduction_speed_matrix>mean_conduction_speed+1.5*std_conduction_speed))=0;
conduction_velocity_frame2(find(conduction_speed_matrix>mean_conduction_speed+1.5*std_conduction_speed))=0;
conduction_velocity_frame2(find(conduction_speed_matrix<mean_conduction_speed-2*std_conduction_speed))=0;
conduction_velocity(:,:,1) = conduction_velocity_frame1; 
conduction_velocity(:,:,2) = conduction_velocity_frame2;

angel_matrix = atan2(conduction_velocity(:,:,2),conduction_velocity(:,:,1)); 
conduction_speed_matrix(find(conduction_speed_matrix>mean_conduction_speed+2*std_conduction_speed))=0; 
conduction_speed_matrix(find(conduction_speed_matrix<mean_conduction_speed-2*std_conduction_speed))=0;         