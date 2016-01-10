function  plot_Vol_infor_new2(median_APD50_matrix,median_APD80_matrix,median_area_under_AP_matrix,mean_area_ratio_matrix,...
    mean_AP_skewness_matrix,mean_upstroke_matrix,signal_swing_matrix,alternans_index_matrix,alternans_score_matrix,number_AP_matrix,...
    pos, total_row_pix,total_col_pix,save_figure_index,result_path,fname,file_name_start_index,file_name_end_index)

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
    set(gca,'FontSize',20,'FontName','Times');
    title('alternans index map','fontname','Times','fontsize',20);
    hold off; 
       
    
    figure; 
    hold on; 
    contourf(signal_swing_matrix,10); colorbar; 
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('signal swing map','fontname','Times','fontsize',20);
    hold off; 

    figure; 
    hold on; 
    contourf(alternans_score_matrix,10); colorbar; 
    plot(x_cord,y_cord,'k');
    set(gca,'FontSize',20,'FontName','Times');
    title('alternans score map','fontname','Times','fontsize',20); 
    hold off; 

    figure; 
    hold on; 
    contourf(number_AP_matrix,10); colorbar; 
    plot(x_cord,y_cord,'k');
    set(gca,'FontSize',20,'FontName','Times');
    title('number of action potential map','fontname','Times','fontsize',20); 
    hold off; 

    figure;
    hold on; 
    contourf(median_APD80_matrix(:,:,2),10); colorbar;
    plot(x_cord,y_cord,'k');
    set(gca,'FontSize',20,'FontName','Times');
    title('APD80 median contour plot for alternating data:even beats ','fontname','Times','fontsize',20); 
    hold off; 
    
    
    figure; 
    hold on; 
    contourf(median_APD50_matrix(:,:,1),10); colorbar;
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('APD50 median contour plot for alternating data: odd beats','fontname','Times','fontsize',20); 
    hold off; 
    
    figure; 
    hold on; 
    contourf(median_APD80_matrix(:,:,1),10); colorbar;
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('APD80 median contour plot for alternating data: odd beats','fontname','Times','fontsize',20); 
    hold off;
    
    figure; 
    hold on; 
    contourf(median_APD50_matrix(:,:,2),10); colorbar;
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('APD50 median contour plot for alternating data: even beats','fontname','Times','fontsize',20);
    hold off; 
    
    figure;
    hold on; 
    contourf(median_area_under_AP_matrix(:,:,1),10); colorbar;
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('area under AP curve map:odd beats ','fontname','Times','fontsize',20); 
    hold off; 
    
    figure; 
    hold on; 
    contourf(median_area_under_AP_matrix(:,:,2),10); colorbar;
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('area under AP curve map: even beats','fontname','Times','fontsize',20);
    hold off; 
    
    figure;
    hold on; 
    contourf(mean_area_ratio_matrix(:,:,2),10); colorbar;
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('area ratio map: even beats','fontname','Times','fontsize',20); 
    hold off;
    
    figure; 
    hold on; 
    contourf(mean_area_ratio_matrix(:,:,1),10); colorbar;
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('area ratio map: odd beats','fontname','Times','fontsize',20);
    hold off; 

    figure; 
    hold on; 
    contourf(mean_upstroke_matrix(:,:,1),10); colorbar;
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('upstroke map: odd beats','fontname','Times','fontsize',20)
    hold off; 

    figure; 
    hold on; 
    contourf(mean_upstroke_matrix(:,:,2),10); colorbar;
    plot(x_cord,y_cord,'k');
    set(gca,'FontSize',20,'FontName','Times');
    title('upstroke map: even beats','fontname','Times','fontsize',20)
    hold off; 

    figure; 
    hold on; 
    contourf(mean_AP_skewness_matrix(:,:,2),10); colorbar;
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('AP skeness plot: even beats','fontname','Times','fontsize',20); 
    hold off; 

    figure; 
    hold on; 
    contourf(mean_AP_skewness_matrix(:,:,1),10); colorbar;
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('AP skeness plot: odd beats','fontname','Times','fontsize',20); 
    hold off; 
else 
    figure; 
    hold on; 
    [ch,ch]=contourf(signal_swing_matrix,10); colorbar;
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k');
    set(gca,'FontSize',20,'FontName','Times');
    title('deltaF/F0','fontname','Times','fontsize',20); 
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'signal_swing',  '.fig'));
    close all
    
    figure;
    hold on; 
    [ch,ch]=contourf(number_AP_matrix,10); colorbar; 
    plot(x_cord,y_cord,'k');
    set(gca,'FontSize',20,'FontName','Times');
    title('number of action potential map','fontname','Times','fontsize',20);
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'num_AP',  '.fig'));
    close all
    
    figure; 
    hold on; 
    [ch,ch]=contourf(alternans_index_matrix,1); colorbar; 
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('alternans index map','fontname','Times','fontsize',20);
    hold off;
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'alternans_index',  '.fig'));
    close all
      
 
     figure; 
    hold on; 
    [ch,ch]=contourf(alternans_score_matrix,20); colorbar; 
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('alternans index map','fontname','Times','fontsize',20);
    hold off;
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'alternans_score_index',  '.fig'));
    close all
    
    figure;
    hold on; 
    [ch,ch]=contourf(median_APD80_matrix(:,:,2),20); colorbar;
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k');
    set(gca,'FontSize',20,'FontName','Times');
    title('APD80 median contour plot for alternating data:even beats ','fontname','Times','fontsize',20); 
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'APD80_even',  '.fig'));
    close all
    
    
    figure; 
    hold on; 
    [ch,ch]=contourf(median_APD50_matrix(:,:,1),20); colorbar;
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('APD50 median contour plot for alternating data: odd beats','fontname','Times','fontsize',20); 
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'APD50_odd',  '.fig'));
    close all
    
    figure; 
    hold on; 
    [ch,ch]=contourf(median_APD80_matrix(:,:,1),20); colorbar;
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('APD80 median contour plot for alternating data: odd beats','fontname','Times','fontsize',20); 
    hold off;
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'APD80_odd',  '.fig'));
    close all
    
    figure; 
    hold on; 
    [ch,ch]=contourf(median_APD50_matrix(:,:,2),20); colorbar;
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('APD50 median contour plot for alternating data: even beats','fontname','Times','fontsize',20);
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'APD50_even',  '.fig'));
    close all
    
   
    
    figure;
    hold on; 
    [ch,ch]=contourf(mean_area_ratio_matrix(:,:,2),20); colorbar;
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('area ratio map: even beats','fontname','Times','fontsize',20); 
    hold off;
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'area_ratio_even',  '.fig'));
    close all
    
    figure; 
    hold on; 
    [ch,ch]=contourf(mean_area_ratio_matrix(:,:,1),20); colorbar;
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('area ratio map: odd beats','fontname','Times','fontsize',20);
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'area_ratio_odd',  '.fig'));
    close all

    figure; 
    hold on; 
    [ch,ch]=contourf(mean_upstroke_matrix(:,:,1),30); colorbar;
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('upstroke map: odd beats','fontname','Times','fontsize',20)
    hold off;
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'upstroke_map_odd',  '.fig'));
    close all

    figure; 
    hold on; 
    [ch,ch]=contourf(mean_upstroke_matrix(:,:,2),30); colorbar;
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k');
    set(gca,'FontSize',20,'FontName','Times');
    title('upstroke map: even beats','fontname','Times','fontsize',20)
    hold off;
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'upstroke_map_even',  '.fig'));
    close all

    figure; 
    hold on; 
    [ch,ch]=contourf(mean_AP_skewness_matrix(:,:,2),20); colorbar;
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('AP skeness plot: even beats','fontname','Times','fontsize',20); 
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'skewness_even',  '.fig'));
    close all

    figure; 
    hold on; 
    [ch,ch]=contourf(mean_AP_skewness_matrix(:,:,1),20); colorbar;
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('AP skeness plot: odd beats','fontname','Times','fontsize',20); 
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'skewness_odd',  '.fig'));
    close all
    
    
end  
     

