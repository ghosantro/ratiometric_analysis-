%% analysis for ratiometric


%% initiating:
Folder_path='J:\Abbirami Sathappan\cell culture\201510062\set7';
fname='2015-10-06_17h18m01s_PM1394Cam00.dat';
fname_full=strcat(Folder_path,'\',fname) ;
Save_Restult_Path='J:\Abbirami Sathappan\cell culture\201510062_analysis\set7';

pacing_frequency=0; %unknown/ 0/ 0/0.5
known_pacing=0; %unknown/ 0/ 0.5
pixels = 128;  
light_number = 2; %%%%%%%%%%%%need' to change here 
num_pixels_to_bin = 3; 
plot_index = 0; % we want to see contour maps 
fibrillation_index = 0; % no fibrillation in the data 
save_figure_index = 1; % want to save figure to file 

dual_vol_CaT_index=1;

xml_file_fname_end_index=strfind(fname,'_PM')-1; % appended_matrix
xml_file_fname=strcat(Folder_path,'\',fname(1:xml_file_fname_end_index),'.xml'); 
xml_struct=xml2struct( xml_file_fname ); 
frameperiod=str2num( xml_struct.XML.Camera.Recorded_Data.Mean_Time_Per_Frame.Text)

%% Draw Contour 
[imag1,imag2,column_number,row_number] = Dual_imaging_opread(fname_full,1, 3, 4);
imag1 = double(imag1); 
imag2 = double(imag2);

%rotating the picture to match cascade view of the camera.
figure
subplot(2,1,1)
imagesc(imag1(:,:,1))
subplot(2,1,2)
imagesc(imag2(:,:,1))

for k=1:size(imag1,3)
    imag1(:,:,k)=imag1(:,:,k)';
    
end

for k=1:size(imag2,3)
    imag2(:,:,k)= imag2(:,:,k)';
end

figure
subplot(2,1,1)
imagesc(imag1(:,:,1))
subplot(2,1,2)
imagesc(imag2(:,:,1))

%new variable for the contour 
imag_contour = imag1(:,:,1);

% spatial filtering the image before drawing the mask
bined_imag_contour = Taking_moving_average(imag_contour,num_pixels_to_bin)
close all;
figure;
imshow(bined_imag_contour./max(max(bined_imag_contour))) 
h=imfreehand;
pos=h.getPosition();
mask_matrix=createMask(h);

[imag1,imag2,column_number,row_number] = Dual_imaging_opread(fname_full,dual_vol_CaT_index, 700, 7400);
imag1 = double(imag1); 
imag2 = double(imag2);
imag1 = Taking_moving_average (imag1,num_pixels_to_bin);
imag2 = Taking_moving_average (imag2,num_pixels_to_bin);

%rotate the picture such that it matches up with the cascade view 
for k = 1:size(imag1,3)
     % to make movie
     %imshow(imag1(:,:,k)'./max(max(imag1(:,:,k)')));
     %  pause

     imag1(:,:,k)= imag1(:,:,k)'; 
end
% 
for k = 1:size(imag2,3)
     % to make movie
     %imshow(imag1(:,:,k)'./max(max(imag1(:,:,k)')));
     %  pause

     imag2(:,:,k)= imag2(:,:,k)'; 
end


% separating trace 1 and trace 2
tissue_index = find(mask_matrix>0);
    trace1= zeros(1,size(imag1,3));
    for u = 1:size(imag1,3)
       average_frame = imag1(:,:,u); 
       trace1(u) = mean(average_frame(tissue_index));        
    end 
    
    tissue_index = find(mask_matrix>0);
    trace2= zeros(1,size(imag2,3));
    for u = 1:size(imag2,3)
       average_frame = imag2(:,:,u); 
       trace2(u) = mean(average_frame(tissue_index));        
    end 

%% Finding derivatives and extreme pattern
derivative1= (trace1(2:end)-trace1(1:end-1))./frameperiod;

% Find a very extreme pattern

extreme_drop_derivative = find(derivative1<-45);
if isempty(extreme_drop_derivative)~=1
     drop_start_index = min(extreme_drop_derivative); 
     cut_index1 = drop_start_index-1; 
else
     cut_index1=length(trace1);
end 

derivative2= (trace2(2:end)-trace2(1:end-1))./frameperiod;

% Find a very extreme pattern 

extreme_drop_derivative = find(derivative2<-45);
  if isempty(extreme_drop_derivative)~=1
      drop_start_index = min(extreme_drop_derivative); 
      cut_index2 = drop_start_index-1; 
  else 
     cut_index2 =length(trace1);
  end 
 
  if abs(cut_index2-cut_index1>3)
    error('something is wrong with detecting light intensity swtich'); 
  else 
        imag1= imag1(:,:,1:cut_index1); 
        trace1=trace1(1:cut_index1); 
        imag2= imag2(:,:,1:cut_index2); 
        trace2=trace2(1:cut_index2); 
  end 

  trace3=trace1./trace2;
  trace3=-trace3;
figure; 
T1=0:frameperiod:frameperiod*(length(trace1)-1);
subplot(3,1,1); plot(T1,trace1);
xlabel('time in ms')
ylabel('trace')
title('415nm')
subplot(3,1,2); plot(0:frameperiod:frameperiod*(length(trace2)-1),trace2);
xlabel('time in ms')
ylabel('trace')
title('365nm')
subplot(3,1,3);
plot(0:frameperiod:frameperiod*(length(trace1)-1),trace3);
xlabel('time in ms')
ylabel('trace')
title('ratio')
pause
    
trace1_skewness=skewness(trace1);
trace2_skewness=skewness(trace2);     

%% calcium bound is 360nm and calcium unbound is 415 (in our case behaves like APD)nm
  if trace1_skewness>trace2_skewness
    imag_vol=imag2;
    imag_CaT=imag1;
    disp('first is 365nm, ca bound state')
  else 
    imag_vol=imag1;
    imag_CaT=imag2;
    disp('First is 415nm, ca unbound state')
  end 
% end

file_name_start_index = 1;
file_name_end_index = findstr(fname,'.mat')-1;   
     
        % %% Bin the Image and select the area of interest
bined_image_415 = imag_vol;
bined_image_365 = imag_CaT;

% mask the area of interest in voltage matrix
if isempty(bined_image_415)~=1
    for v=1:size(bined_image_415,3)
        bined_image_415(:,:,v) = bined_image_415(:,:,v); %.*mask_matrix;
    end
end 

if isempty(bined_image_365)~=1
    for v=1:size(bined_image_365,3)
        bined_image_365(:,:,v)=bined_image_365(:,:,v); %.*mask_matrix;
    end 
end 


ratiometric_image= bined_image_415; %./bined_image_365;
%% split into ca bound and ca unbound state: ca bound: 365:similar to ca. ca unbound 425, similar to voltage. 

if isempty(ratiometric_image)==0
     Type_of_analysis=0;

    disp('doing 415 analysis'); 

    [APD_infor_cell, Area_AP_cell, AP_skewness_cell, signal_swing_matrix, alternans_index_matrix,alternans_score_matrix, number_AP_matrix,number_upstroke_matrix, upstroke_cell]=voltage_analysis_new...
    (frameperiod,known_pacing,ratiometric_image, plot_index);    
            
    %align all the signals          
            
    [APD50_matrix, APD80_matrix, area_under_AP_matrix, area_ratio_matrix, AP_skewness_matrix,upstroke_matrix_all,missing_out_first_beat_index] = align_signal (APD_infor_cell, Area_AP_cell,...
    AP_skewness_cell, number_AP_matrix, upstroke_cell, frameperiod, fibrillation_index,known_pacing, pacing_frequency, 0,...
    fname,file_name_start_index,file_name_end_index,strcat(Save_Restult_Path,'\Vol'));            


   %filter these signals
    filtered_APD50_matrix = APD50_matrix; 
    filtered_APD80_matrix=APD80_matrix;
    filtered_area_under_AP_matrix=area_under_AP_matrix;
    filtered_area_ratio_matrix=area_ratio_matrix;
    filtered_AP_skewness_matrix=AP_skewness_matrix;
    filtered_upstroke_matrix_all=upstroke_matrix_all;
    for k = 1:size(APD50_matrix,3)
      filtered_APD50_matrix(:,:,k) = filter_Vol_result_matrix_new(APD50_matrix(:,:,k), signal_swing_matrix, 0);
      filtered_APD80_matrix(:,:,k)=filter_Vol_result_matrix_new(APD80_matrix(:,:,k), signal_swing_matrix, 0); 
      filtered_area_under_AP_matrix(:,:,k) = filter_Vol_result_matrix_new(filtered_area_under_AP_matrix(:,:,k), signal_swing_matrix, 0); 
      filtered_area_ratio_matrix(:,:,k)=filter_Vol_result_matrix_new(area_ratio_matrix(:,:,k),signal_swing_matrix, 0); 
      filtered_AP_skewness_matrix(:,:,k)=filter_Vol_result_matrix_new(AP_skewness_matrix(:,:,k),signal_swing_matrix, 0);
      filtered_upstroke_matrix_all(:,:,k)=filter_Vol_result_matrix_new(upstroke_matrix_all(:,:,k),signal_swing_matrix, 1);
    end
   
    % Calculate stats of these signals 
    [median_APD50_matrix, mean_APD50_matrix, std_APD50_matrix] = calculate_stats (filtered_APD50_matrix,alternans_index_matrix); 
    [median_APD80_matrix, mean_APD80_matrix, std_APD80_matrix] = calculate_stats (filtered_APD80_matrix,alternans_index_matrix);
    [median_area_under_AP_matrix, mean_area_under_AP_matrix, std_area_under_AP_matrix] = calculate_stats (filtered_area_under_AP_matrix,alternans_index_matrix);
    [median_area_ratio_matrix, mean_area_ratio_matrix, std_area_ratio_matrix] = calculate_stats (filtered_area_ratio_matrix,alternans_index_matrix);
    [median_AP_skewness_matrix, mean_AP_skewness_matrix, std_AP_skewness_matrix] = calculate_stats (filtered_AP_skewness_matrix,alternans_index_matrix);
    [median_upstroke_matrix, mean_upstroke_matrix, std_upstroke_matrix] = calculate_stats (filtered_upstroke_matrix_all,alternans_index_matrix);

    median_upstroke_matrix(:,:,1) = median_upstroke_matrix(:,:,1)-min(min(median_upstroke_matrix(:,:,1))); 
    median_upstroke_matrix(:,:,2) = median_upstroke_matrix(:,:,2)-min(min(median_upstroke_matrix(:,:,2))); 
    mean_upstroke_matrix(:,:,1) = median_upstroke_matrix(:,:,1)-min(min(mean_upstroke_matrix(:,:,1))); 
    mean_upstroke_matrix(:,:,2) = median_upstroke_matrix(:,:,2)-min(min(mean_upstroke_matrix(:,:,2))); 
    
     % volage plot
    plot_Vol_infor_new2(median_APD50_matrix,median_APD80_matrix,median_area_under_AP_matrix,mean_area_ratio_matrix,...
    mean_AP_skewness_matrix,mean_upstroke_matrix,signal_swing_matrix,alternans_index_matrix,alternans_score_matrix,...
    number_AP_matrix,pos, size(median_APD50_matrix,1),size(median_APD50_matrix,2), save_figure_index,strcat(Save_Restult_Path,'\Vol'),fname,file_name_start_index,file_name_end_index);
    
    %conduction velocity    
    conduction_struct=main_func_for_conduction_velocity(filtered_upstroke_matrix_all,mean_upstroke_matrix,save_figure_index,fname,file_name_start_index,...
    file_name_end_index,known_pacing, pacing_frequency,pos, frameperiod, size(median_APD50_matrix,1),size(median_APD50_matrix,2),strcat(Save_Restult_Path,'\Vol'));

            
    File_name_starting_index=min(strfind(fname,'h'));  
    File_name_end_index= strfind(fname,'.dat')-1; 
    save_file_name = strcat(fname(File_name_starting_index: File_name_end_index),'Vol'); 
    save_file_name(findstr(save_file_name,'-'))=''; %chage invalid naming 
    save_file_name(findstr(save_file_name,'.'))='';
    
%     [CaTD_infor_cell_415, Area_CaT_cell_415, CaT_skewness_cell_415, CaT_signal_swing_matrix_415, CaT_alternans_index_matrix_415,CaT_alternans_score_matrix_415, CaT_peak_alternans_index_matrix_415,CaT_peak_alternans_score_matrix_415, number_CaT_matrix_415, number_CaT_upstroke_matrix_415, upstroke_cell_CaT_415]=CaT_analysis_new2(frameperiod,known_pacing,ratiometric_image,plot_index);
    
    
%      [CaT50_matrix_415, CaT80_matrix_415, area_under_CaT_matrix_415, area_CaT_ratio_matrix_415, CaT_skewness_matrix_415,upstroke_matrix_all_CaT_415,missing_out_first_beat_index_CaT_415] = align_signal (CaTD_infor_cell_415, Area_CaT_cell_415,...
%      CaT_skewness_cell_415, number_CaT_matrix_415, upstroke_cell_CaT_415, frameperiod, fibrillation_index,known_pacing, pacing_frequency, 0,fname,file_name_start_index,file_name_end_index,strcat(Save_Restult_Path,'\CaT')); 
%      
%  
%  [CaTD_infor_cell, Area_CaT_cell, CaT_skewness_cell, CaT_signal_swing_matrix, CaT_alternans_index_matrix,CaT_alternans_score_matrix, CaT_peak_alternans_index_matrix, ...
%     CaT_peak_alternans_score_matrix, number_CaT_matrix, number_CaT_upstroke_matrix, upstroke_cell_CaT]=CaT_analysis_new2 (frameperiod,known_pacing,ratiometric_image,...
%     pacing_frequency,  plot_index);
%     
% 
%     [CaT50_matrix, CaT80_matrix, area_under_CaT_matrix, area_CaT_ratio_matrix, CaT_skewness_matrix,upstroke_matrix_all_CaT,missing_out_first_beat_index_CaT] = align_signal (CaTD_infor_cell, Area_CaT_cell,...
%     CaT_skewness_cell, number_CaT_matrix, upstroke_cell_CaT, frameperiod, fibrillation_index,known_pacing, pacing_frequency, 0,...
%     fname,file_name_start_index,file_name_end_index,strcat(Save_Restult_Path,'\CaT')); 



%      filtered_CaT50_matrix_415 = CaT50_matrix_415; 
%      filtered_CaT80_matrix_415 = CaT80_matrix_415;
%      filtered_area_under_CaT_matrix_415 = area_under_CaT_matrix_415;
%      filtered_area_ratio_CaT_matrix_415 =area_CaT_ratio_matrix_415;
%      filtered_CaT_skewness_matrix_415 =CaT_skewness_matrix_415;
%      filtered_upstroke_matrix_all_CaT_415 =upstroke_matrix_all_CaT_415;
%      filtered_upstroke_matrix_all_CaT_original_415 = filtered_upstroke_matrix_all_CaT_415;
end 

% [APD_infor_cell, Area_AP_cell, AP_skewness_cell, signal_swing_matrix, alternans_index_matrix,alternans_score_matrix, number_AP_matrix,number_upstroke_matrix, upstroke_cell]=voltage_analysis_new(frameperiod,known_pacing,bined_image_415, plot_index);    
%     





