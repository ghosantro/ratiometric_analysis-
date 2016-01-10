function [APD50_matrix, APD80_matrix, area_under_AP_matrix, area_ratio_matrix, AP_skewness_matrix,upstroke_matrix_all,missing_out_first_beat_index] = align_signal (APD_infor_cell, Area_AP_cell,...
    AP_skewness_cell, number_upstroke_matrix, upstroke_cell, frameperiod, fibrillation_index,known_pacing, pacing_frequency, save_file_index,...
    fname,file_name_start_index,file_name_end_index,result_path)
% because the start end end frame are pretty arbitary, therefore in one
% slice some pixels may miss out the first or the last action potential
% therefore all the action potential signal needs to be aligned in a slice,
% this is particular important for alternating case 

%% code
% check the upstroke cell

if fibrillation_index ~=1 % if there's no fibrillation  
    standard_AP_number = mode(number_upstroke_matrix(find(number_upstroke_matrix>0)));
    upstroke_matrix_all = NaN(size(upstroke_cell,1), size(upstroke_cell,2),standard_AP_number);
     
    one_less_index = []; 
    one_more_index = []; 
    Mean_pacing_interval_matrix = NaN(size(upstroke_cell,1), size(upstroke_cell,2));
    for r = 1:size(upstroke_cell,1)
        for c = 1:size(upstroke_cell,2) 
            
            if number_upstroke_matrix(r,c)==standard_AP_number;
               upstroke_matrix_all(r,c,:) = upstroke_cell{r,c}(1:number_upstroke_matrix(r,c));
               mean_pacing_interval = mean(upstroke_cell{r,c}(2:number_upstroke_matrix(r,c))-upstroke_cell{r,c}(1:number_upstroke_matrix(r,c)-1)); 
               Mean_pacing_interval_matrix(r,c) = mean_pacing_interval; 
            elseif number_upstroke_matrix(r,c) == standard_AP_number-1 
                one_less_index= [one_less_index;r,c]; 
            elseif number_upstroke_matrix(r,c) == standard_AP_number+1
                one_more_index = [one_more_index;r,c]; 
            end 
                
        end 
    end 
    
%     figure; contourf(Mean_pacing_interval_matrix,20); colorbar; 
%     one_less_index
%     one_more_index
    first_frame = upstroke_matrix_all(:,:,1);
%     figure; contourf(first_frame,20); caxis([0,200]);colorbar
%     figure; hist(first_frame(first_frame>0),200); 
    last_frame = upstroke_matrix_all(:,:,end); 
%     figure; hist(last_frame(last_frame>0),200); 
%     figure; contourf(last_frame,20); colorbar
    median_upstroke_frame_one = median(first_frame(find(isnan(first_frame)~=1))); 
    median_upstroke_frame_last = median(last_frame(find(isnan(last_frame)~=1))) ;
    
    over_all_mean_pacing_interval = mean(Mean_pacing_interval_matrix(isnan(Mean_pacing_interval_matrix)~=1))*frameperiod;
    pacing_interval = 1000/pacing_frequency;
    if known_pacing ==1 
        % check whether this is kept up with pacing. If not, report it to a
        % file and record down the real pacing interval 
      if over_all_mean_pacing_interval > pacing_interval*1.1 && save_file_index~=1
        disp('did not manage to keep up with the pacing or pacing frequency recorded is wrong')
      elseif over_all_mean_pacing_interval < pacing_interval*0.9 && save_file_index~=1
        disp('pacing frequency recorded is wrong')            
      elseif over_all_mean_pacing_interval >  pacing_interval*1.1 && save_file_index==1
        fid = fopen(strcat(result_path,'\matlab_result_from_autanalysis','Not_keeping_up_pacing.txt'),'a'); 
        fprintf(fid,'%s \n',strcat(fname(file_name_start_index:file_name_end_index),'did not keep up with pacing or the pacing requency recorded is wrong')); 
        fclose(fid); 
     elseif over_all_mean_pacing_interval <  pacing_interval*0.9 && save_file_index==1
       fid = fopen(strcat(result_path,'\matlab_result_from_autanalysis','Not_keeping_up_pacing.txt'),'a'); 
       fprintf(fid,'%s \n',strcat(fname(file_name_start_index:file_name_end_index),'the pacing requency recorded is wrong')); 
       fclose(fid); 
      end  
    end    
        
        
    % if there's need to trunct the first one as in some pixles we miss the
    % first upstroke point 
   
    filtering_threshold = pacing_interval*0.8;
    
    
    missing_out_first_beat_index = zeros(size(upstroke_cell,1), size(upstroke_cell,2));
    having_extra_first_beat_index = zeros(size(upstroke_cell,1), size(upstroke_cell,2));
    miss_out_last_beat_index = zeros(size(upstroke_cell,1), size(upstroke_cell,2));
    having_extra_last_beat_index = zeros(size(upstroke_cell,1), size(upstroke_cell,2));
    
    if isempty(one_less_index)~=1
        for i =1: size(one_less_index,1)
             row_num = one_less_index(i,1);
             col_num= one_less_index(i,2);
             if isempty(upstroke_cell{row_num,col_num})~=1
                if abs(upstroke_cell{row_num,col_num}(1)-median_upstroke_frame_one)> filtering_threshold 
                   upstroke_vec = upstroke_cell{row_num,col_num};
                   upstroke_matrix_all(row_num,col_num,:) =[NaN,upstroke_vec(1:standard_AP_number-1)];
                    missing_out_first_beat_index(row_num,col_num)=1;
                    
                elseif abs(upstroke_cell{row_num,col_num}(end)-median_upstroke_frame_last)> filtering_threshold 
                   upstroke_vec = upstroke_cell{row_num,col_num};
                   upstroke_matrix_all(row_num,col_num,:) =[upstroke_vec(1:standard_AP_number-1),NaN];
                   miss_out_last_beat_index (row_num,col_num) = 1; 
                end 
             end
            %pause
        end 
    end 
    
    
    
 
    if isempty(one_more_index)~=1    
        for i =1: size(one_more_index,1)
            row_num = one_more_index(i,1);
            col_num= one_more_index(i,2); 
            if abs(upstroke_cell{row_num,col_num}(standard_AP_number+1)-median_upstroke_frame_last)> filtering_threshold 

               upstroke_matrix_all(row_num,col_num,:) =upstroke_cell{row_num,col_num}(1:standard_AP_number); 
               having_extra_last_beat_index(r,c) = 1;
            elseif abs(upstroke_cell{row_num,col_num}(1)-median_upstroke_frame_one)> filtering_threshold  
                upstroke_matrix_all(row_num,col_num,:) = upstroke_cell{row_num, col_num}(2:standard_AP_number+1);
                 having_extra_first_beat_index (row_num,col_num) =1; 

            end
        end 
    end 
    
   
  
   % align APD and other AP parameters 
   APD50_matrix = NaN(size(upstroke_cell,1), size(upstroke_cell,2),standard_AP_number);
   APD80_matrix = NaN(size(upstroke_cell,1), size(upstroke_cell,2),standard_AP_number);
   area_ratio_matrix = NaN(size(upstroke_cell,1), size(upstroke_cell,2),standard_AP_number);
   area_under_AP_matrix = NaN(size(upstroke_cell,1), size(upstroke_cell,2),standard_AP_number);
   AP_skewness_matrix = NaN(size(upstroke_cell,1), size(upstroke_cell,2),standard_AP_number);
   
   for r = 1:size(upstroke_cell,1)
        for c = 1:size(upstroke_cell,2) 
                                
            if length(find(APD_infor_cell{r,c}.APD50>0))==standard_AP_number
               % count = count+1;
               % disp('a')
               APD50_vect = APD_infor_cell{r,c}.APD50;
               APD50_matrix(r,c,:) = APD50_vect(APD50_vect>0); 
              % APD_infor_cell{r,c}.APD50
            elseif having_extra_first_beat_index(r,c)==1&&length(find(APD_infor_cell{r,c}.APD50>0))==standard_AP_number+1;
               %disp('b')
               APD50_vect = APD_infor_cell{r,c}.APD50(2:end);
               APD50_matrix(r,c,:) = APD50_vect(APD50_vect>0); 
               
              % APD_infor_cell{r,c}.APD50
            elseif  missing_out_first_beat_index(r,c) ==1&&length(find(APD_infor_cell{r,c}.APD50>0))==standard_AP_number-1;
               %disp('c')
               APD50_vect = APD_infor_cell{r,c}.APD50(1:end);
               APD50_matrix(r,c,:) =[NaN,APD50_vect(APD50_vect>0)]; 
               
               %APD_infor_cell{r,c}.APD50
            elseif miss_out_last_beat_index (r,c) == 1&&length(find(APD_infor_cell{r,c}.APD50>0))==standard_AP_number-1;
               % disp('d')
               APD50_vect = APD_infor_cell{r,c}.APD50(1:end);
               APD50_matrix(r,c,:) =[ APD50_vect(APD50_vect>0),NaN]; 
                
              % APD_infor_cell{r,c}.APD50
            elseif having_extra_last_beat_index(r,c) == 1&&length(find(APD_infor_cell{r,c}.APD50)>0)==standard_AP_number+1;
              %  disp('e')
               APD50_vect = APD_infor_cell{r,c}.APD50(1:end-1); 
               APD50_matrix(r,c,:) = APD50_vect(APD50_vect>0); 
               
              % APD_infor_cell{r,c}.APD50
            end
            
            if length(APD_infor_cell{r,c}.APD80)==standard_AP_number
                APD80_vect = APD_infor_cell{r,c}.APD80;
               APD80_matrix(r,c,:) = APD80_vect;                
            elseif having_extra_first_beat_index(r,c)==1&&length(APD_infor_cell{r,c}.APD80)==standard_AP_number+1;
                APD80_vect = APD_infor_cell{r,c}.APD80(2:end);
                APD80_matrix(r,c,:) = APD80_vect;
            elseif  missing_out_first_beat_index(r,c) ==1&&length(APD_infor_cell{r,c}.APD80)==standard_AP_number-1;
                APD80_vect = APD_infor_cell{r,c}.APD80;
                APD80_matrix(r,c,:) =[NaN, APD80_vect];                     
           elseif miss_out_last_beat_index (r,c) == 1&&length(APD_infor_cell{r,c}.APD80)==standard_AP_number-1;
                APD80_vect = APD_infor_cell{r,c}.APD80;
                APD80_matrix(r,c,:) = [APD80_vect,NaN];
            elseif length(APD_infor_cell{r,c}.APD80)==standard_AP_number-1&&length(APD_infor_cell{r,c}.APD50)==standard_AP_number
                APD80_vect = APD_infor_cell{r,c}.APD80;
                APD80_matrix(r,c,:) = [APD80_vect,NaN];
                  
           elseif having_extra_last_beat_index(r,c) == 1&&length(APD_infor_cell{r,c}.APD80)==standard_AP_number+1;
               APD80_vect = APD_infor_cell{r,c}.APD80(1:end-1);
               APD80_matrix(r,c,:) = APD80_vect;
               
           end
            
            if length(Area_AP_cell{r,c}.area_ratio_vect)==standard_AP_number
                
               area_ratio_matrix(r,c,:) = Area_AP_cell{r,c}.area_ratio_vect;                
            elseif having_extra_first_beat_index(r,c)==1&&length(Area_AP_cell{r,c}.area_ratio_vect)==standard_AP_number+1;
               area_ratio_matrix(r,c,:) = Area_AP_cell{r,c}.area_ratio_vect(2:end); 
            elseif  missing_out_first_beat_index(r,c) ==1&&length(Area_AP_cell{r,c}.area_ratio_vect)==standard_AP_number-1;
               area_ratio_matrix(r,c,:) = [NaN,Area_AP_cell{r,c}.area_ratio_vect(1:end)];               
           elseif miss_out_last_beat_index (r,c) == 1&&length(Area_AP_cell{r,c}.area_ratio_vect)==standard_AP_number-1;
               area_ratio_matrix(r,c,:) = [Area_AP_cell{r,c}.area_ratio_vect,NaN];
            elseif length(Area_AP_cell{r,c}.area_ratio_vect)==standard_AP_number-1&&length(find(APD_infor_cell{r,c}.APD50>0))==standard_AP_number
               area_ratio_matrix(r,c,:) = [Area_AP_cell{r,c}.area_ratio_vect,NaN]; 
           elseif having_extra_last_beat_index(r,c) == 1&&length(Area_AP_cell{r,c}.area_ratio_vect)==standard_AP_number+1;
               area_ratio_matrix(r,c,:) =Area_AP_cell{r,c}.area_ratio_vect(1:end-1); 
            end
            
            if length( Area_AP_cell{r,c}.area_under_AP_vect)==standard_AP_number
               area_under_AP_matrix(r,c,:) =  Area_AP_cell{r,c}.area_under_AP_vect;                
            elseif having_extra_first_beat_index(r,c)==1&&length( Area_AP_cell{r,c}.area_under_AP_vect)==standard_AP_number+1;
               area_under_AP_matrix(r,c,:) =  Area_AP_cell{r,c}.area_under_AP_vect(2:end); 
            elseif  missing_out_first_beat_index(r,c) ==1&&length( Area_AP_cell{r,c}.area_under_AP_vect)==standard_AP_number-1;
               area_under_AP_matrix(r,c,:) = [NaN, Area_AP_cell{r,c}.area_under_AP_vect(1:end)];               
           elseif miss_out_last_beat_index (r,c) == 1&&length( Area_AP_cell{r,c}.area_under_AP_vect)==standard_AP_number-1;
               area_under_AP_matrix(r,c,:) = [ Area_AP_cell{r,c}.area_under_AP_vect,NaN]; 
            elseif length( Area_AP_cell{r,c}.area_under_AP_vect)==standard_AP_number-1&&length(find(APD_infor_cell{r,c}.APD50>0))==standard_AP_number
               area_under_AP_matrix(r,c,:) = [ Area_AP_cell{r,c}.area_under_AP_vect,NaN]; 
           elseif having_extra_last_beat_index(r,c) == 1&&length( Area_AP_cell{r,c}.area_under_AP_vect)==standard_AP_number+1;
               area_under_AP_matrix(r,c,:)= Area_AP_cell{r,c}.area_under_AP_vect(1:end-1); 
               
            end
            
            if length(AP_skewness_cell{r,c})==standard_AP_number
               AP_skewness_matrix(r,c,:) = AP_skewness_cell{r,c};                
            elseif having_extra_first_beat_index(r,c)==1&&length(AP_skewness_cell{r,c})==standard_AP_number+1;
               AP_skewness_matrix(r,c,:) = AP_skewness_cell{r,c}(2:end); 
            elseif  missing_out_first_beat_index(r,c) ==1&&length(AP_skewness_cell{r,c})==standard_AP_number-1;
               AP_skewness_matrix(r,c,:) = [NaN,AP_skewness_cell{r,c}(1:end)];               
           elseif miss_out_last_beat_index (r,c) == 1&&length(AP_skewness_cell{r,c})==standard_AP_number-1;
               AP_skewness_matrix(r,c,:) = [AP_skewness_cell{r,c},NaN]; 
           elseif having_extra_last_beat_index(r,c) == 1&&length(AP_skewness_cell{r,c})==standard_AP_number+1;
               AP_skewness_matrix(r,c,:) = AP_skewness_cell{r,c}(1:end-1); 
           end
               
        end 
    end 
 
   
 
 
  length(find(number_upstroke_matrix==7))
else 
    error(' need to implement an algorithm for fibrillating tissue')
    
end 
