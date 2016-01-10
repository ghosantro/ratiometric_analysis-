%% Draw the contours and analyze in one go

% set all the input by user 

Folder_path='J:\Abbirami Sathappan\cell culture\201510062\set5';
fname='2015-10-06_17h04m37s_PM1394Cam00.dat';
fname_full=strcat(Folder_path,'\',fname) ;
Save_Restult_Path='J:\Abbirami Sathappan\cell culture\201510062_analysis\set5';

pacing_frequency=0; 
known_pacing=0; 
pixels = 128;   
% %  y_pixel= 125; 
% %  x_pixel =77; 
light_number = 2; %%%%%%%%%%%%need' to change here 
num_pixels_to_bin = 3; 
plot_index = 0; % we want to see contour maps 
fibrillation_index = 0; % no fibrillation in the data 
save_figure_index = 1; % want to save figure to file 

dual_vol_CaT_index=1;


xml_file_fname_end_index=strfind(fname,'_PM')-1; % appended_matrix
xml_file_fname=strcat(Folder_path,'\',fname(1:xml_file_fname_end_index),'.xml'); 
xml_struct=xml2struct( xml_file_fname ); 
frameperiod =str2num( xml_struct.XML.Camera.Recorded_Data.Mean_Time_Per_Frame.Text)

% % % 
% %Draw Contour 
[imag1,imag2,column_number,row_number] = Dual_imaging_opread(fname_full,1, 3, 4);
imag1 = double(imag1); 
imag2 = double(imag2);


%         %rotate the picture such that it matches up with the cascade view 
for k = 1:size(imag1,3)
     % to make movie
     %imshow(imag1(:,:,k)'./max(max(imag1(:,:,k)')));
     %  pause

     imag1(:,:,k)= imag1(:,:,k)'; 
end
for k = 1:size(imag2,3)
    % to make movie
     %imshow(imag1(:,:,k)'./max(max(imag1(:,:,k)')));
     %  pause


     imag2(:,:,k)= imag2(:,:,k)'; 
end
imag_contour = imag1(:,:,1); 
%        spatially filter the image 
bined_imag_contour =Taking_moving_average (imag_contour,num_pixels_to_bin); 
close all;
figure;
imshow(bined_imag_contour./max(max(bined_imag_contour))); 
h = imfreehand; %draw a line to define the tissue bondary 
pos = h.getPosition(); %get the coordinates of the boundary 
mask_matrix = createMask(h);
        

[imag1,imag2,column_number,row_number] = Dual_imaging_opread(fname_full,dual_vol_CaT_index, 10, 3000);
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

if dual_vol_CaT_index==0
    
    tissue_index = find(mask_matrix>0);
    trace1= zeros(1,size(imag1,3));
    for u = 1:size(imag1,3)
        average_frame = imag1(:,:,u); 
        trace1(u) = mean(average_frame(tissue_index));        
    end 
        
    derivative1 = (trace1(2:end)-trace1(1:end-1))./frameperiod; 
    %find a very extreme pattern 
    extreme_drop_derivative = find(derivative1<-45);
    if isempty(extreme_drop_derivative)~=1
       drop_start_index = min(extreme_drop_derivative); 
       cut_index1 = drop_start_index-1; 
    else 
       cut_index1 =length(trace1);
    end   
    imag1= imag1(:,:,1:cut_index1); 
    trace1=trace1(1:cut_index1); 
    
    figure; plot(0:frameperiod:frameperiod*(length(trace1)-1),trace1)
    Vol_CaT_index=input('is this CaT or Vol, 0=CaT, 1=Vol?'); 

    if Vol_CaT_index ==0 % i.e.CaT              
        imag_CaT=imag1(:,:,2:end); 
        imag_CaT = Taking_moving_average (imag_CaT,num_pixels_to_bin);

        imag_vol=[];
        disp('CaT only'); 
    elseif Vol_CaT_index ==1             
        imag_vol=imag1(:,:,2:end); 
        imag_vol = Taking_moving_average (imag_vol,num_pixels_to_bin);
        imag_CaT=[]; 
        disp('vol only'); 
    else 
        problem_data_index_temp=1; 
    end 
     
elseif dual_vol_CaT_index==1 
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
    
    derivative1 = (trace1(2:end)-trace1(1:end-1))./frameperiod; 
    %find a very extreme pattern 
    extreme_drop_derivative = find(derivative1<-45);
    if isempty(extreme_drop_derivative)~=1
       drop_start_index = min(extreme_drop_derivative); 
       cut_index1 = drop_start_index-1; 
    else 
       cut_index1 =length(trace1);
    end           
                
    derivative2 = (trace2(2:end)-trace2(1:end-1))./frameperiod; 
    %find a very extreme pattern 
    extreme_drop_derivative = find(derivative2<-45);
    if isempty(extreme_drop_derivative)~=1
        drop_start_index = min(extreme_drop_derivative); 
        cut_index2 = drop_start_index-1; 
    else 
        cut_index2 =length(trace1);
    end           

    if abs(cut_index2-cut_index1)>3
       error('soemthing is wrong with detecting light intensity swtich'); 
    else 
        imag1= imag1(:,:,1:cut_index1); 
        trace1=trace1(1:cut_index1); 
        imag2= imag2(:,:,1:cut_index2); 
        trace2=trace2(1:cut_index2); 
    end 
%       

    figure; subplot(1,2,1); plot(0:frameperiod:frameperiod*(length(trace1)-1),trace1);
    subplot(1,2,2); plot(0:frameperiod:frameperiod*(length(trace2)-1),trace2);
    pause
    
    trace1_skewness=skewness(trace1);
    trace2_skewness=skewness(trace2);     

    if trace1_skewness>trace2_skewness                
        imag_vol=imag2;
        imag_CaT=imag1; 
        disp('first is CaT'); 
    else 
        imag_vol=imag1; 
        imag_CaT=imag2; 
        disp('first is Vol'); 
    end      
end    
    
file_name_start_index = 1;
file_name_end_index = findstr(fname,'.mat')-1;   
        

% %% Bin the Image and select the area of interest
bined_image_vol = imag_vol;
bined_image_cat = imag_CaT;
% if isempty(bined_image_vol)~=0 
%     image_size_x = size(bined_image_vol,1); % need to double check
%     image_size_y = size(bined_image_vol,2); % need to doubel check
% elseif isempty(bined_image_cat)~=0
%     image_size_x = size(bined_image_cat,1); % need to double check
%     image_size_y = size(bined_image_cat,2); % need to doubel check
% end 
% 
% clearvars imag1 imag2 imag_overall imag1_new imag2_new
% % select area of interest (CaT and Vol will have the same ROI)
% close all; 
% 
% if isempty(bined_image_vol)~=1
%     for v = 1:size(bined_image_vol,3)
%         bined_image_vol(:,:,v) = bined_image_vol(:,:,v).*mask_matrix; 
%     end
% end 
% 
% % mask the are of interest in CaT matrix 
% if isempty(bined_image_cat)~=1
%     for v = 1:size(bined_image_cat,3)
%         bined_image_cat(:,:,v) = bined_image_cat(:,:,v).*mask_matrix; 
%     end   
% 
% end 

    
%% splilt into voltage and Calcium analysis 
        
 if isempty(bined_image_vol)==0
     Type_of_analysis=0;
    [APD_infor_cell, Area_AP_cell, AP_skewness_cell, signal_swing_matrix, alternans_index_matrix,alternans_score_matrix, number_AP_matrix,number_upstroke_matrix, upstroke_cell]=voltage_analysis_new...
    (frameperiod,known_pacing,bined_image_vol, plot_index);    
            
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
    
    
     assignin('base',save_file_name,struct('APD80_matrix',filtered_APD80_matrix,'APD50_matrix',filtered_APD50_matrix,...
     'median_APD50_matrix',median_APD50_matrix,'median_APD80_matrix',median_APD50_matrix,'area_ratio_matrix',filtered_area_ratio_matrix,...
     'mean_area_ratio_matrix',mean_area_ratio_matrix,'AP_skewness_matrix',filtered_AP_skewness_matrix,'mean_AP_skewness_matrix',mean_AP_skewness_matrix,...
     'onset_depolar',filtered_upstroke_matrix_all,'mean_upstroke_matrix',mean_upstroke_matrix,'alternans_index',alternans_index_matrix,...
     'conduction_structure',conduction_struct,'number_AP_matrix',number_AP_matrix));
 
      
  
 
%      Including type of analysis if the pacing frequency is known by the
%      user.
    save(strcat(Save_Restult_Path, '\Vol\matlab_figure_from_autanalysis\',fname(1:end-4),Type_of_analysis, '.mat'), save_file_name);   


     assignin('base',strcat(save_file_name,'_APInfCell'),APD_infor_cell);
     save(strcat(Save_Restult_Path,'\Vol\matlab_figure_from_autanalysis\', fname(1:end-4), '_APInfCell', Type_of_analysis,'.mat'),strcat(save_file_name,'_APInfCell'));

     assignin('base',strcat(save_file_name,'_APSkewCell'),AP_skewness_cell);
     save(strcat(Save_Restult_Path,'\Vol\matlab_figure_from_autanalysis\', fname(1:end-4), '_APSkewCell',Type_of_analysis, '.mat'),strcat(save_file_name,'_APSkewCell'));


     assignin('base',strcat(save_file_name,'_APAreaCell'),Area_AP_cell);
     save(strcat(Save_Restult_Path,'\Vol\matlab_figure_from_autanalysis\', fname(1:end-4), '_APAreaCell', Type_of_analysis,'.mat'),strcat(save_file_name,'_APAreaCell'));
               
               
%      assignin('base',strcat(save_file_name, '_APPaceCell'), pacing_interval_cell_AP);
%      save(strcat(Save_Restult_Path,'\Vol\matlab_figure_from_autanalysis\', fname(1:end-4), '_APPaceCell', Type_of_analysis,'.mat'),strcat(save_file_name,'_APPaceCell'))
% 
%                

    %     if exist('imag2')==1
    %         Calcium_analysis (imag2,num_pixels_to_bin,frameperiod) 
    %     else 
    %         disp('there is no calcium analysis')
        %end 

 end       

        

     
 if isempty(bined_image_cat)==0



    
    disp('doing CaT analysis'); 
    
    [CaTD_infor_cell, Area_CaT_cell, CaT_skewness_cell, CaT_signal_swing_matrix, CaT_alternans_index_matrix,CaT_alternans_score_matrix, CaT_peak_alternans_index_matrix, ...
    CaT_peak_alternans_score_matrix, number_CaT_matrix, number_CaT_upstroke_matrix, upstroke_cell_CaT]=CaT_analysis_new2 (frameperiod,known_pacing,bined_image_cat,...
    pacing_frequency,  plot_index);
    

    [CaT50_matrix, CaT80_matrix, area_under_CaT_matrix, area_CaT_ratio_matrix, CaT_skewness_matrix,upstroke_matrix_all_CaT,missing_out_first_beat_index_CaT] = align_signal (CaTD_infor_cell, Area_CaT_cell,...
    CaT_skewness_cell, number_CaT_matrix, upstroke_cell_CaT, frameperiod, fibrillation_index,known_pacing, pacing_frequency, 0,...
    fname,file_name_start_index,file_name_end_index,strcat(Save_Restult_Path,'\CaT')); 

    filtered_CaT50_matrix = CaT50_matrix; 
    filtered_CaT80_matrix=CaT80_matrix;
    filtered_area_under_CaT_matrix=area_under_CaT_matrix;
    filtered_area_ratio_CaT_matrix=area_CaT_ratio_matrix;
    filtered_CaT_skewness_matrix=CaT_skewness_matrix;
    filtered_upstroke_matrix_all_CaT=upstroke_matrix_all_CaT;
    filtered_upstroke_matrix_all_CaT_original= filtered_upstroke_matrix_all_CaT;
    for k = 1:size(CaT50_matrix,3)
      filtered_CaT50_matrix(:,:,k) = filter_Vol_result_matrix_new(CaT50_matrix(:,:,k), CaT_signal_swing_matrix, 0);
      filtered_CaT80_matrix(:,:,k)=filter_Vol_result_matrix_new(CaT80_matrix(:,:,k), CaT_signal_swing_matrix, 0); 
      filtered_area_under_CaT_matrix(:,:,k) = filter_Vol_result_matrix_new(filtered_area_under_CaT_matrix(:,:,k), CaT_signal_swing_matrix, 0); 
      filtered_area_ratio_CaT_matrix(:,:,k)=filter_Vol_result_matrix_new(area_CaT_ratio_matrix(:,:,k),CaT_signal_swing_matrix, 0); 
      filtered_CaT_skewness_matrix(:,:,k)=filter_Vol_result_matrix_new(CaT_skewness_matrix(:,:,k),CaT_signal_swing_matrix, 0);
      filtered_upstroke_matrix_all_CaT(:,:,k)=filter_Vol_result_matrix_new(upstroke_matrix_all_CaT(:,:,k),CaT_signal_swing_matrix, 1);
      filtered_upstroke_matrix_all_CaT_original(:,:,k)=filtered_upstroke_matrix_all_CaT(:,:,k);
      frame_matrix_temp=filtered_upstroke_matrix_all_CaT(:,:,k);
      filtered_upstroke_matrix_all_CaT(:,:,k)=filtered_upstroke_matrix_all_CaT(:,:,k)-prctile(frame_matrix_temp(find(frame_matrix_temp>0)),1);
    end

    % filtere the sorted matrix    


    % Calculate stats of these signals 
    [median_CAT50_matrix, mean_CaT50_matrix, std_CaT50_matrix] = calculate_stats (filtered_CaT50_matrix,CaT_alternans_index_matrix); 
    [median_CaT80_matrix, mean_CaT80_matrix, std_CaT80_matrix] = calculate_stats (filtered_CaT80_matrix,CaT_alternans_index_matrix);
    [median_area_under_CaT_matrix, mean_area_under_CaT_matrix, std_area_under_CaT_matrix] = calculate_stats (filtered_area_under_CaT_matrix,CaT_alternans_index_matrix);
    [median_area_ratio_CaT_matrix, mean_area_ratio_CaT_matrix, std_area_ratio_CaT_matrix] = calculate_stats (filtered_area_ratio_CaT_matrix,CaT_alternans_index_matrix);
    [median_CaT_skewness_matrix, mean_CaT_skewness_matrix, std_CaT_skewness_matrix] = calculate_stats (filtered_CaT_skewness_matrix,CaT_alternans_index_matrix);
    [median_upstroke_matrix_CaT, mean_upstroke_matrix_CaT, std_upstroke_matrix_CaT] = calculate_stats (filtered_upstroke_matrix_all_CaT,CaT_alternans_index_matrix);

    median_upstroke_matrix_CaT(:,:,1) = median_upstroke_matrix_CaT(:,:,1)-min(min(median_upstroke_matrix_CaT(:,:,1))); 
    median_upstroke_matrix_CaT(:,:,2) = median_upstroke_matrix_CaT(:,:,2)-min(min(median_upstroke_matrix_CaT(:,:,2))); 
    mean_upstroke_matrix_CaT(:,:,1) = median_upstroke_matrix_CaT(:,:,1)-min(min(mean_upstroke_matrix_CaT(:,:,1))); 
    mean_upstroke_matrix_CaT(:,:,2) = median_upstroke_matrix_CaT(:,:,2)-min(min(mean_upstroke_matrix_CaT(:,:,2)));   

    file_name_start_index= min(strfind(fname,'h')); 
    file_name_end_index=strfind(fname,'.dat')-1;
    
    plot_CaT_infor_new2(median_CAT50_matrix,median_CaT80_matrix,median_area_under_CaT_matrix,median_area_ratio_CaT_matrix,...
    median_CaT_skewness_matrix,mean_upstroke_matrix_CaT,CaT_signal_swing_matrix,CaT_alternans_index_matrix,CaT_alternans_score_matrix,...
    CaT_peak_alternans_index_matrix,CaT_peak_alternans_score_matrix,number_CaT_matrix,pos, size(mean_CaT50_matrix,1),size(mean_CaT50_matrix,2),...
    save_figure_index,strcat(Save_Restult_Path,'\CaT'),fname,file_name_start_index,file_name_end_index);
% 

    File_name_starting_index=min(strfind(fname,'h'));  
    File_name_end_index= strfind(fname,'.dat')-1; 
    save_file_name = strcat(fname(File_name_starting_index: File_name_end_index),'CaT'); 
    save_file_name(findstr(save_file_name,'-'))=''; %chage invalid naming 
    save_file_name(findstr(save_file_name,'.'))='';

    assignin('base',save_file_name,struct('CaT80_matrix',filtered_CaT80_matrix,'CaT50_matrix',filtered_CaT50_matrix,'mean_CaT80_matrix',mean_CaT80_matrix,'mean_CaT50_matrix',mean_CaT50_matrix...
    ,'mean_CaT_area_ratio',mean_area_ratio_CaT_matrix,'mean_CaT_skewness',mean_CaT_skewness_matrix,'upstroke_all_matrix_CaT',filtered_upstroke_matrix_all_CaT,...
    'CaT_alternans_index',CaT_alternans_index_matrix,'CaT_alternans_score',CaT_alternans_score_matrix, 'CaT_peak_alternans_index',CaT_peak_alternans_index_matrix,...
    'CaT_peak_alternans_score',CaT_peak_alternans_score_matrix,'CaT_upstroke_original',...
    filtered_upstroke_matrix_all_CaT_original, 'original_CaT80',CaT80_matrix,'origianl_CaT50',CaT50_matrix,'number_CaT',number_CaT_matrix));

    save(strcat(Save_Restult_Path, fname(1:end-4),'CaT', Type_of_analysis, '.mat'), save_file_name);

    assignin('base',strcat(save_file_name,'_CaTInfCell'),CaTD_infor_cell);
    save(strcat(Save_Restult_Path,'\CaT\matlab_figure_from_autanalysis\', fname(1:end-4), '_CaTInfCell', Type_of_analysis,'.mat'),strcat(save_file_name,'_CaTInfCell'));

    assignin('base',strcat(save_file_name,'_CaTSkewCell'),CaT_skewness_cell);
    save(strcat(Save_Restult_Path,'\CaT\matlab_figure_from_autanalysis\', fname(1:end-4), '_CaTSkewCell',Type_of_analysis, '.mat'),strcat(save_file_name,'_CaTSkewCell'));


    assignin('base',strcat(save_file_name,'_CaTAreaCell'),Area_CaT_cell);
    save(strcat(Save_Restult_Path,'\CaT\matlab_figure_from_autanalysis\', fname(1:end-4), '_CaTAreaCell', Type_of_analysis,'.mat'),strcat(save_file_name,'_CaTAreaCell'))

%     assignin('base',strcat(save_file_name,'_CaTPaceCell'),pacing_interval_cell_CaT);
%     save(strcat(Save_Restult_Path,'\CaT\matlab_figure_from_autanalysis\', fname(1:end-4), '_CaTPaceCell', Type_of_analysis,'.mat'),strcat(save_file_name,'_CaTPaceCell'))
%        
 end 



  
 
     
     
     
 
 
 
 
