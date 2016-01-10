function  [median_matrix, mean_matrix, std_matrix] = calculate_stats (signal_matrix,alternans_matrix) 
median_matrix= NaN(size(signal_matrix,1),size(signal_matrix,2),2);
mean_matrix= NaN(size(signal_matrix,1),size(signal_matrix,2),2);
std_matrix = NaN(size(signal_matrix,1),size(signal_matrix,2),2);

for i = 1:size(signal_matrix,1)
    for j = 1:size(signal_matrix,2)
        if alternans_matrix(i,j)~=1
            signal_vector = signal_matrix(i,j,:);
            median_matrix(i,j,1) = median(signal_vector(isnan(signal_vector)~=1));
            median_matrix(i,j,2) = median(signal_vector(isnan(signal_vector)~=1));
            mean_matrix(i,j,1) = mean(signal_vector(isnan(signal_vector)~=1));
            mean_matrix(i,j,2) = mean(signal_vector(isnan(signal_vector)~=1));
            std_matrix(i,j,1) = mean(std(isnan(signal_vector)~=1));
            std_matrix(i,j,2) = mean(std(isnan(signal_vector)~=1));
        else
            signal_vector = signal_matrix(i,j,:);
            %median(signal_vector(isnan(signal_vector(1:2:end))~=1))
            median_calculation_even = signal_vector(1:2:end);
            median_calculation_odd = signal_vector(2:2:end);
            median_calculation_even(isnan(median_calculation_even)==1)=[];
            median_calculation_odd(isnan(median_calculation_odd)==1)=[];
             median(median_calculation_odd);
              median(median_calculation_even);
            median_matrix(i,j,1) = median(median_calculation_odd);
            mean_matrix(i,j,1) = mean(median_calculation_odd);
            std_matrix(i,j,1) = mean(median_calculation_odd);
           % median(signal_vector(isnan(signal_vector(2:2:end))~=1))
            median_matrix(i,j,2) = median(median_calculation_even);
            mean_matrix(i,j,2) = mean(median_calculation_even);
            std_matrix(i,j,2) = mean(median_calculation_even);
            
        end
    end 
end 


    