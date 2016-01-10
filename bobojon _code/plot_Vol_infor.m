function   [filtered_APD80_median_matrix_even, filtered_APD50_median_matrix_even,filtered_APD80_median_matrix_odd, filtered_APD50_median_matrix_odd,...
filtered_area_ratio_median_matrix_even,filtered_area_ratio_median_matrix_odd,filtered_AP_skeness_median_matrix_even,filtered_AP_skeness_median_matrix_odd]...
= plot_Vol_infor(APD_infor_cell, Area_AP_cell,...
AP_skewness_cell, signal_swing_matrix, alternans_index_matrix,alternans_score_matrix, number_AP_matrix,upstroke_cell,total_row_pix, total_col_pix,...
fibrillation_index,save_figure_index,fname,file_name_start_index,file_name_end_index,known_pacing,pacing_frequency,frameperiod, pos, result_path)

%% This function gets the output from the voltage analaysis function and plot the contour maps 
% The contour map will be created are: 
% alternans map 
% signal swing map 
% APD80 median map (if alternans occurs, there will be two maps created for even beats and odd beats) 
% APD50 median map (if alternans occurs, there will be two maps created for even beats and odd beats)
% area under curve map (if alternans occurs, there will be two maps created for even beats and odd beats)
% area ratio map (if alternans occurs, there will be two maps created for even beats and odd beats)
% mean AP skewness map
% number of AP map

%% code
% plot alternans map

% get the contour of tissue out 

x_cord = pos(:,1); 
y_cord = pos(:,2); 

% get rid of anything out of bounds 
x_cord(x_cord>total_row_pix) =total_row_pix; 
x_cord(x_cord<0) = 0;
y_cord(y_cord>total_col_pix) = total_col_pix; 
y_cord(y_cord<0) = 0;

if save_figure_index ==0 
    figure; 
    hold on; 
    contourf(alternans_index_matrix,10); colorbar; 
    plot(x_cord,y_cord,'k'); 
    title('alternans index map');
    hold off; 

    figure; 
    hold on; 
    contourf(signal_swing_matrix,10); colorbar; 
    plot(x_cord,y_cord,'k'); 
    title('signal swing map');
    hold off; 

    figure; 
    hold on; 
    contourf(alternans_score_matrix,10); colorbar; 
    plot(x_cord,y_cord,'k');
    title('alternans_score_map'); 
    hold off; 

    figure; 
    hold on; 
    contourf(number_AP_matrix,10); colorbar; 
    plot(x_cord,y_cord,'k');
    title('number of action potential map'); 
    hold off; 
else 
    figure; 
    hold on; 
    contourf(signal_swing_matrix,10); colorbar;
    plot(x_cord,y_cord,'k');
    set(gca,'FontSize',20,'FontName','Times');
    title('deltaF/F0','fontname','Times','fontsize',20); 
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'signal_swing',  '.fig'));
    close all
    
    figure;
    hold on; 
    contourf(number_AP_matrix,10); colorbar; 
    plot(x_cord,y_cord,'k');
    set(gca,'FontSize',20,'FontName','Times');
    title('number of action potential map','fontname','Times','fontsize',20);
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'num_AP',  '.fig'));
    close all
      
end    
 

   
matrix_APD80_median_plot_even = NaN(total_row_pix, total_col_pix); 
matrix_APD50_median_plot_even = NaN(total_row_pix, total_col_pix); 
matrix_APD80_median_plot_odd = NaN(total_row_pix, total_col_pix); 
matrix_APD50_median_plot_odd = NaN(total_row_pix, total_col_pix); 
matrix_APD80_std_plot_even = NaN(total_row_pix, total_col_pix); 
matrix_APD50_std_plot_even = NaN(total_row_pix, total_col_pix); 
matrix_APD80_std_plot_odd = NaN(total_row_pix, total_col_pix); 
matrix_APD50_std_plot_odd = NaN(total_row_pix, total_col_pix); 
matrix_AP_area_plot_even = NaN(total_row_pix, total_col_pix); 
matrix_AP_area_std_plot_even = NaN(total_row_pix, total_col_pix); 
matrix_area_ratio_plot_even = NaN(total_row_pix,total_col_pix); 
matrix_area_ratio_std_plot_even = NaN(total_row_pix,total_col_pix);
matrix_AP_area_plot_odd = NaN(total_row_pix, total_col_pix); 
matrix_AP_area_std_plot_odd = NaN(total_row_pix, total_col_pix); 
matrix_area_ratio_plot_odd = NaN(total_row_pix,total_col_pix); 
matrix_area_ratio_std_plot_odd = NaN(total_row_pix,total_col_pix); 
matrix_AP_skeness_median_even = NaN(total_row_pix,total_col_pix); 
matrix_AP_skeness_median_odd = NaN(total_row_pix, total_col_pix); 
matrix_AP_skeness_std_even = NaN(total_row_pix,total_col_pix); 
matrix_AP_skeness_std_odd = NaN(total_row_pix, total_col_pix); 
for r = 1: total_row_pix
  for c = 1:total_col_pix 

    if alternans_index_matrix(r,c)==0&&isempty(find(APD_infor_cell{r,c}.APD80_median>0))==0   
       matrix_APD80_median_plot_even(r,c) = APD_infor_cell{r,c}.APD80_median; 
       matrix_APD80_median_plot_odd(r,c) = APD_infor_cell{r,c}.APD80_median; 
    elseif alternans_index_matrix (r,c)==1&&isempty(find(APD_infor_cell{r,c}.APD80_median_even>0))==0

        matrix_APD80_median_plot_even(r,c) = APD_infor_cell{r,c}.APD80_median_even; 
        if isempty(find(APD_infor_cell{r,c}.APD80_median_odd>0))==0 
           matrix_APD80_median_plot_odd(r,c) = APD_infor_cell{r,c}.APD80_median_odd; 
        end
    end 

    if alternans_index_matrix(r,c)==0&&isempty(find(APD_infor_cell{r,c}.APD80_std>0))==0   
       matrix_APD80_std_plot_even(r,c) = APD_infor_cell{r,c}.APD80_std; 
       matrix_APD80_std_plot_odd(r,c) = APD_infor_cell{r,c}.APD80_std; 
    elseif alternans_index_matrix (r,c)==1&&isempty(find(APD_infor_cell{r,c}.APD80_std_even>0))==0

        matrix_APD80_std_plot_even(r,c) = APD_infor_cell{r,c}.APD80_std_even; 
        if isempty(find(APD_infor_cell{r,c}.APD80_std_odd>0))==0 
           matrix_APD80_std_plot_odd(r,c) = APD_infor_cell{r,c}.APD80_std_odd; 
        end
    end 


    if alternans_index_matrix(r,c)==0&& isempty(find(APD_infor_cell{r,c}.APD50_median>0))==0  
       matrix_APD50_median_plot_even(r,c) = APD_infor_cell{r,c}.APD50_median; 
       matrix_APD50_median_plot_odd(r,c) = APD_infor_cell{r,c}.APD50_median;
    elseif alternans_index_matrix(r,c)==1&&isempty(find(APD_infor_cell{r,c}.APD50_median_even>0))==0 
        matrix_APD50_median_plot_even(r,c) = APD_infor_cell{r,c}.APD50_median_even;
        if isempty(find(APD_infor_cell{r,c}.APD50_median_odd>0))==0 
           matrix_APD50_median_plot_odd(r,c) = APD_infor_cell{r,c}.APD50_median_odd;
        end
    end 

    if alternans_index_matrix(r,c)==0&&isempty(find(APD_infor_cell{r,c}.APD50_std>0))==0   
       matrix_APD50_std_plot_even(r,c) = APD_infor_cell{r,c}.APD50_std; 
       matrix_APD50_std_plot_odd(r,c) = APD_infor_cell{r,c}.APD50_std; 
    elseif alternans_index_matrix (r,c)==1&&isempty(find(APD_infor_cell{r,c}.APD50_std_even>0))==0

        matrix_APD50_std_plot_even(r,c) = APD_infor_cell{r,c}.APD50_std_even; 
        if isempty(find(APD_infor_cell{r,c}.APD50_std_odd>0))==0 
           matrix_APD50_std_plot_odd(r,c) = APD_infor_cell{r,c}.APD50_std_odd; 
        end
    end 

    


    if alternans_index_matrix(r,c)==0&&isempty(find(Area_AP_cell{r,c}.area_under_AP_mean>0))==0   
       matrix_AP_area_plot_even (r,c) = Area_AP_cell{r,c}.area_under_AP_mean; 
       matrix_AP_area_plot_odd(r,c) = Area_AP_cell{r,c}.area_under_AP_mean; 
    elseif alternans_index_matrix (r,c)==1&&isempty(find(Area_AP_cell{r,c}.area_under_AP_mean_even>0))==0

       matrix_AP_area_plot_even(r,c) = Area_AP_cell{r,c}.area_under_AP_mean_even; 
       if isempty(find(Area_AP_cell{r,c}.area_under_AP_mean_odd>0))==0 
          matrix_AP_area_plot_odd(r,c) = Area_AP_cell{r,c}.area_under_AP_mean_odd; 
       end
    end 

    if alternans_index_matrix(r,c)==0&&isempty(find(Area_AP_cell{r,c}.area_under_AP_std>0))==0   
       matrix_AP_area_std_plot_even (r,c) = Area_AP_cell{r,c}.area_under_AP_std; 
       matrix_AP_area_std_plot_odd(r,c) = Area_AP_cell{r,c}.area_under_AP_std; 
    elseif alternans_index_matrix (r,c)==1&&isempty(find(Area_AP_cell{r,c}.area_under_AP_std_even>0))==0

       matrix_AP_area_std_plot_even(r,c) = Area_AP_cell{r,c}.area_under_AP_std_even; 
       if isempty(find(Area_AP_cell{r,c}.area_under_AP_std_odd>0))==0 
          matrix_AP_area_std_plot_odd(r,c) = Area_AP_cell{r,c}.area_under_AP_std_odd; 
       end
    end 

    if alternans_index_matrix(r,c)==0&&isempty(find(Area_AP_cell{r,c}.area_ratio_mean>0))==0   
       matrix_area_ratio_plot_even (r,c) = Area_AP_cell{r,c}.area_ratio_mean; 
       matrix_area_ratio_plot_odd(r,c) = Area_AP_cell{r,c}.area_ratio_mean; 
    elseif alternans_index_matrix (r,c)==1&&isempty(find(Area_AP_cell{r,c}.area_ratio_mean_even>0))==0

       matrix_area_ratio_plot_even(r,c) = Area_AP_cell{r,c}.area_ratio_mean_even; 
       if isempty(find(Area_AP_cell{r,c}.area_ratio_mean_odd>0))==0 
          matrix_area_ratio_plot_odd(r,c) = Area_AP_cell{r,c}.area_ratio_mean_odd; 
       end
    end 

    if alternans_index_matrix(r,c)==0&&isempty(find(Area_AP_cell{r,c}.area_ratio_std>0))==0   
       matrix_area_ratio_std_plot_even (r,c) = Area_AP_cell{r,c}.area_ratio_std; 
       matrix_area_ratio_std_plot_odd(r,c) = Area_AP_cell{r,c}.area_ratio_std; 
    elseif alternans_index_matrix (r,c)==1&&isempty(find(Area_AP_cell{r,c}.area_ratio_std_even>0))==0

       matrix_area_ratio_std_plot_even(r,c) = Area_AP_cell{r,c}.area_ratio_std_even; 
       if isempty(find(Area_AP_cell{r,c}.area_ratio_std_odd>0))==0 
          matrix_area_ratio_std_plot_odd(r,c) = Area_AP_cell{r,c}.area_ratio_std_odd; 
       end
    end

     if alternans_index_matrix(r,c)==0&&isempty(AP_skewness_cell{r,c}.AP_skewness_median)==0   
       matrix_AP_skeness_median_even(r,c) = AP_skewness_cell{r,c}.AP_skewness_median; 
       matrix_AP_skeness_median_odd(r,c) = AP_skewness_cell{r,c}.AP_skewness_median; 
    elseif alternans_index_matrix (r,c)==1&&isempty(AP_skewness_cell{r,c}.AP_skewness_median_even)==0

      matrix_AP_skeness_median_even(r,c) = AP_skewness_cell{r,c}.AP_skewness_median_even; 
       if isempty(AP_skewness_cell{r,c}.AP_skewness_median_odd)==0 
          matrix_AP_skeness_median_odd(r,c) = AP_skewness_cell{r,c}.AP_skewness_median_odd; 
       end   
     end    

     if alternans_index_matrix(r,c)==0&&isempty(AP_skewness_cell{r,c}.AP_skewness_std)==0   
       matrix_AP_skeness_std_even(r,c) = AP_skewness_cell{r,c}.AP_skewness_std; 
       matrix_AP_skeness_std_odd(r,c) = AP_skewness_cell{r,c}.AP_skewness_std; 
    elseif alternans_index_matrix (r,c)==1&&isempty(AP_skewness_cell{r,c}.AP_skewness_std_even)==0

      matrix_AP_skeness_std_even(r,c) = AP_skewness_cell{r,c}.AP_skewness_std_even; 
       if isempty(AP_skewness_cell{r,c}.AP_skewness_std_odd)==0 
          matrix_AP_skeness_std_odd(r,c) = AP_skewness_cell{r,c}.AP_skewness_std_odd; 
       end   
     end    



  end 
end 


%remove outliers 

filtered_APD80_median_matrix_odd = filter_Vol_result_matrix( matrix_APD80_median_plot_odd, matrix_APD80_std_plot_odd, signal_swing_matrix, number_AP_matrix, fibrillation_index,0,0);

filtered_APD50_median_matrix_odd = filter_Vol_result_matrix( matrix_APD50_median_plot_odd, matrix_APD50_std_plot_odd, signal_swing_matrix, number_AP_matrix, fibrillation_index,0,0);
filtered_AP_area_median_matrix_odd = filter_Vol_result_matrix( matrix_AP_area_plot_odd, matrix_AP_area_std_plot_odd, signal_swing_matrix, number_AP_matrix, fibrillation_index,0,0);
filtered_area_ratio_median_matrix_odd = filter_Vol_result_matrix( matrix_area_ratio_plot_odd, matrix_area_ratio_std_plot_odd, signal_swing_matrix, number_AP_matrix, fibrillation_index,0,0);

filtered_APD80_median_matrix_even = filter_Vol_result_matrix( matrix_APD80_median_plot_even, matrix_APD80_std_plot_even, signal_swing_matrix, number_AP_matrix, fibrillation_index,0,0);
filtered_APD50_median_matrix_even = filter_Vol_result_matrix( matrix_APD50_median_plot_even, matrix_APD50_std_plot_even, signal_swing_matrix, number_AP_matrix, fibrillation_index,0,0);
filtered_AP_area_median_matrix_even = filter_Vol_result_matrix( matrix_AP_area_plot_even, matrix_AP_area_std_plot_even, signal_swing_matrix, number_AP_matrix, fibrillation_index,0,0);
filtered_area_ratio_median_matrix_even = filter_Vol_result_matrix( matrix_area_ratio_plot_even, matrix_area_ratio_std_plot_even, signal_swing_matrix, number_AP_matrix, fibrillation_index,0,0);
filtered_AP_skeness_median_matrix_even =filter_Vol_result_matrix(matrix_AP_skeness_median_even, matrix_AP_skeness_std_even, signal_swing_matrix, number_AP_matrix, fibrillation_index,0,1);
filtered_AP_skeness_median_matrix_odd = filter_Vol_result_matrix(matrix_AP_skeness_median_odd, matrix_AP_skeness_std_odd, signal_swing_matrix, number_AP_matrix, fibrillation_index,0,1);


filtered_upstroke_matrix = filter_upstroke_onset (upstroke_cell, number_AP_matrix, known_pacing, pacing_frequency, save_figure_index,fname, file_name_start_index,file_name_end_index,fibrillation_index,signal_swing_matrix,frameperiod,result_path); 

%filtered_upstroke_matrix = filter_upstroke_onset (upstroke_cell, number_AP_matrix, 0,0,0,0,0,0,0,signal_swing_matrix,frameperiod); 


if save_figure_index==0
    figure;
    hold on; 
    contourf(filtered_APD80_median_matrix_even,10); colorbar;
    plot(x_cord,y_cord,'k'); 
    title('APD80 median contour plot for alternating data:even beats '); 
    hold off; 
    figure; 
    hold on; 
    contourf(filtered_APD50_median_matrix_odd,10); colorbar;
    plot(x_cord,y_cord,'k'); 
    title('APD50 median contour plot for alternating data: odd beats');
    hold off; 
    figure; 
    hold on; 
    contourf(filtered_APD80_median_matrix_odd,10); colorbar;
    plot(x_cord,y_cord,'k'); 
    title('APD80 median contour plot for alternating data: odd beats'); 
    hold off; 
    figure; 
    hold on; 
    contourf(filtered_APD50_median_matrix_even,10); colorbar;
    plot(x_cord,y_cord,'k'); 
    title('APD50 median contour plot for alternating data: even beats');
    hold off; 
    figure;
    hold on; 
    contourf(filtered_area_ratio_median_matrix_odd,10); colorbar;
    plot(x_cord,y_cord,'k'); 
    title('area ratio map: odd beats');
    hold off; 
    figure; 
    hold on; 
    contourf(filtered_area_ratio_median_matrix_even,10); colorbar;
    plot(x_cord,y_cord,'k'); 
    title('area ratio map: even beats');
    hold off; 
    
    figure;
    hold on; 
    contourf(filtered_AP_area_median_matrix_even,10); colorbar;
    plot(x_cord,y_cord,'k'); 
    title('area under AP curve map:even beats '); 
    hold off; 
    figure; 
    hold on; 
    contourf(filtered_AP_area_median_matrix_odd,10); colorbar;
    plot(x_cord,y_cord,'k'); 
    title('area under AP curve map: odd beats');
    hold off; 

    figure; 
    hold on; 
    contourf(filtered_upstroke_matrix(:,:,1),10); colorbar;
    plot(x_cord,y_cord,'k'); 
    title('upstroke map')
    hold off; 

    figure; 
    hold on; 
    contourf(filtered_AP_skeness_median_matrix_even,10); colorbar;
    plot(x_cord,y_cord,'k'); 
    title('AP skeness plot: even beats'); 
    hold off; 

    figure; 
    hold on; 
    contourf(filtered_AP_skeness_median_matrix_odd,10); colorbar;
    plot(x_cord,y_cord,'k'); 
    title('AP skeness plot: odd beats'); 
    hold off; 
else 
    figure; 
    hold on;
    contourf(filtered_APD80_median_matrix_even,10); colorbar;
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('APD80 even beats ','fontname','Times','fontsize',20);
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'APD80_even',  '.fig'));
    close all;
    
    figure; 
    hold on; 
    contourf(filtered_APD50_median_matrix_odd,10); colorbar;
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('APD50 odd beats','fontname','Times','fontsize',20);
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'APD50_odd',  '.fig'));
    close all; 
    
    figure;
    hold on; 
    contourf(filtered_APD80_median_matrix_odd,10); colorbar;
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('APD80 odd beats','fontname','Times','fontsize',20);
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'APD80_odd',  '.fig'));
    close all; 
    
    
    figure; 
    hold on; 
    contourf(filtered_APD50_median_matrix_even,10); colorbar;
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('APD50 even beats','fontname','Times','fontsize',20);
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'APD50_even',  '.fig'));
    close all; 
    
    figure;
    hold on; 
    contourf(filtered_area_ratio_median_matrix_odd,10); colorbar;
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('area ratio:odd beats','fontname','Times','fontsize',20); 
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'area_ratio_odd',  '.fig'));
    close all
    
    figure;
    hold on; 
    contourf(filtered_area_ratio_median_matrix_even,10); colorbar;
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('area ratio: even beats','fontname','Times','fontsize',20);
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'area_ratio_even',  '.fig'));
    close all
    
    
    

    figure; 
    hold on; 
    contourf(filtered_upstroke_matrix(:,:,1),10); colorbar;
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('onset depolarisation','fontname','Times','fontsize',20); 
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'onset_upstroke',  '.fig'));
    close all; 
    
    figure; 
    hold on; 
    contourf(filtered_AP_skeness_median_matrix_even,10); colorbar;
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('AP skeness: even beats','fontname','Times','fontsize',20); 
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'AP_skewness_even',  '.fig'));
    close all; 
    
    figure; 
    hold on; 
    contourf(filtered_AP_skeness_median_matrix_odd,10); colorbar;
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('AP skeness: odd beats','fontname','Times','fontsize',20);
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'AP_skewness_odd',  '.fig'));
    close all;
end  
     

