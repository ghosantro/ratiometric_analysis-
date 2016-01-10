function  plot_CaT_infor_new2(mean_CaT50_matrix,mean_CaT80_matrix,median_area_under_CaT_matrix,mean_area_ratio_CaT_matrix,...
    mean_CaT_skewness_matrix,mean_CaT_upstroke_matrix,CaT_signal_swing_matrix,CaT_alternans_index_matrix,CaT_alternans_score_matrix,...
    CaT_peak_alternans_index_matrix,CaT_peak_alternans_score_matrix,number_CaT_matrix,pos, total_row_pix,total_col_pix,save_figure_index,result_path,fname,file_name_start_index,file_name_end_index)

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
    [ch,ch]=contourf(CaT_alternans_index_matrix,10); colorbar; 
    set(ch,'edgecolor',none);
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT alternans index map','fontname','Times','fontsize',20);
    hold off; 
           
    
    figure; 
    hold on; 
    [ch,ch]=contourf(CaT_signal_swing_matrix,20); colorbar; 
    set(ch,'edgecolor',none);
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT signal swing map','fontname','Times','fontsize',20);
    hold off; 

    figure; 
    hold on; 
    [ch,ch]=contourf(CaT_alternans_score_matrix,10); colorbar; 
    set(ch,'edgecolor',none);
    plot(x_cord,y_cord,'k');
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT alternans score map','fontname','Times','fontsize',20); 
    hold off; 

    figure; 
    hold on; 
    [ch,ch]=contourf(number_CaT_matrix,20); colorbar; 
    set(ch,'edgecolor',none);
    plot(x_cord,y_cord,'k');
    set(gca,'FontSize',20,'FontName','Times');
    title('number of CaT map','fontname','Times','fontsize',20); 
    hold off; 

     figure; 
    hold on; 
    [ch,ch]=contourf(mean_CaT50_matrix(:,:,2),10); colorbar;
    set(ch,'edgecolor',none);
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT50 mean contour plot for alternating data: even beats','fontname','Times','fontsize',20); 
    hold off;  
    
    
    figure; 
    hold on; 
    [ch,ch]=contourf(mean_CaT50_matrix(:,:,1),10); colorbar;
    set(ch,'edgecolor',none);
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT50 mean contour plot for alternating data: odd beats','fontname','Times','fontsize',20); 
    hold off; 
    
    figure; 
    hold on; 
    [ch,ch]=contourf(mean_CaT80_matrix(:,:,1),10); colorbar;
    set(ch,'edgecolor',none);
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT80 median contour plot for alternating data: odd beats','fontname','Times','fontsize',20); 
    hold off;
    
    figure; 
    hold on; 
    [ch,ch]=contourf(mean_CaT50_matrix(:,:,2),10); colorbar;
    set(ch,'edgecolor',none);
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT50 median contour plot for alternating data: even beats','fontname','Times','fontsize',20);
    hold off; 
    
    figure;
    hold on; 
   [ch,ch]= contourf(median_area_under_CaT_matrix(:,:,1),10); colorbar;
    set(ch,'edgecolor',none);
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('area under CaT curve map:odd beats ','fontname','Times','fontsize',20); 
    hold off; 
    
    figure; 
    hold on; 
    [ch,ch]=contourf(median_area_under_CaT_matrix(:,:,2),10); colorbar;
    set(ch,'edgecolor',none);
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('area under CaT curve map: even beats','fontname','Times','fontsize',20);
    hold off; 
    
    figure;
    hold on; 
    [ch,ch]=contourf(mean_area_ratio_CaT_matrix(:,:,2),10); colorbar;
    set(ch,'edgecolor',none);
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT area ratio map: even beats','fontname','Times','fontsize',20); 
    hold off;
    
    figure; 
    hold on; 
    [ch,ch]=contourf(mean_area_ratio_CaT_matrix(:,:,1),10); colorbar;
    set(ch,'edgecolor',none);
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT area ratio map: odd beats','fontname','Times','fontsize',20);
    hold off; 

    figure; 
    hold on; 
    [ch,ch]=contourf(mean_CaT_upstroke_matrix(:,:,1),30); colorbar;
    set(ch,'edgecolor',none);
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT upstroke map: odd beats','fontname','Times','fontsize',20)
    hold off; 

    figure; 
    hold on; 
    [ch,ch]=contourf(mean_CaT_upstroke_matrix(:,:,2),30); colorbar;
    set(ch,'edgecolor',none);
    plot(x_cord,y_cord,'k');
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT upstroke map: even beats','fontname','Times','fontsize',20)
    hold off; 

    figure; 
    hold on; 
    [ch,ch]=contourf(mean_CaT_skewness_matrix(:,:,2),10); colorbar;
    set(ch,'edgecolor',none);
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('Cat skeness plot: even beats','fontname','Times','fontsize',20); 
    hold off; 

    figure; 
    hold on; 
    [ch,ch]=contourf(mean_CaT_skewness_matrix(:,:,1),10); colorbar;
    set(ch,'edgecolor',none);
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT skeness plot: odd beats','fontname','Times','fontsize',20); 
    hold off; 
else 
    figure; 
    hold on; 
    [ch,ch]=contourf(CaT_signal_swing_matrix,20); colorbar; 
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT signal swing map','fontname','Times','fontsize',20);
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'CaT_signal_swing',  '.fig'));
    close all
    
    figure; 
    hold on; 
    [ch,ch]=contourf(number_CaT_matrix,20); colorbar; 
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k');
    set(gca,'FontSize',20,'FontName','Times');
    title('number of CaT map','fontname','Times','fontsize',20); 
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'num_CaT',  '.fig'));
    close all
    
    figure; 
    hold on;     
    [ch,ch]=contourf(CaT_alternans_index_matrix,10); colorbar; 
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT alternans index map','fontname','Times','fontsize',20);
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'CaT_alternans_index',  '.fig'));
    close all
      
    figure; 
    hold on; 
    [ch,ch]=contourf(CaT_alternans_score_matrix,10); colorbar; 
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k');
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT alternans score map','fontname','Times','fontsize',20); 
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'CaT_alternans_score',  '.fig'));
    close all
    
    figure; 
    hold on;     
    [ch,ch]=contourf(CaT_peak_alternans_index_matrix,10); colorbar; 
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT peak alternans index map','fontname','Times','fontsize',20);
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'CaT_peak_alternans_index',  '.fig'));
    close all
      
    figure; 
    hold on; 
    [ch,ch]=contourf(CaT_peak_alternans_score_matrix,10); colorbar; 
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k');
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT peak alternans score map','fontname','Times','fontsize',20); 
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'CaT_peak_alternans_score',  '.fig'));
    close all

    
    figure;
    hold on; 
    [ch,ch]=contourf(mean_CaT80_matrix(:,:,2),20); colorbar;
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k');
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT80 mean contour plot for alternating data:even beats ','fontname','Times','fontsize',20); 
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'CaT80_even',  '.fig'));
    close all
    
    
    figure; 
    hold on; 
    [ch,ch]=contourf(mean_CaT80_matrix(:,:,1),20); colorbar;
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT80 mean contour plot for alternating data: odd beats','fontname','Times','fontsize',20); 
    hold off;
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'CaT80_odd',  '.fig'));
    close all
    
     figure; 
    hold on; 
    [ch,ch]=contourf(mean_CaT50_matrix(:,:,1),20); colorbar;
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT50 mean contour plot for alternating data: odd beats','fontname','Times','fontsize',20); 
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'CaT50_odd',  '.fig'));
    close all
    
    figure; 
    hold on; 
    [ch,ch]=contourf(mean_CaT50_matrix(:,:,2),20); colorbar;
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT50 mean contour plot for alternating data: even beats','fontname','Times','fontsize',20); 
    hold off;  
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'CaT50_even',  '.fig'));
    close all
    
    figure;
    hold on; 
    [ch,ch]=contourf(mean_area_ratio_CaT_matrix(:,:,2),20); colorbar;
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT area ratio map: even beats','fontname','Times','fontsize',20); 
    hold off;
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'CaT_area_ratio_even',  '.fig'));
    close all
    
        figure; 
    hold on; 
    [ch,ch]=contourf(mean_area_ratio_CaT_matrix(:,:,1),20); colorbar;
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT area ratio map: odd beats','fontname','Times','fontsize',20);
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'CaT_area_ratio_odd',  '.fig'));
    close all
    
        

    figure; 
    hold on; 
    [ch,ch]=contourf(mean_CaT_upstroke_matrix(:,:,1),30); colorbar;
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT upstroke map: odd beats','fontname','Times','fontsize',20)
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'CaT_upstroke_map_odd',  '.fig'));
    close all

    figure; 
    hold on; 
    [ch,ch]=contourf(mean_CaT_upstroke_matrix(:,:,2),30); colorbar;
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k');
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT upstroke map: even beats','fontname','Times','fontsize',20)
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'CaT_upstroke_map_even',  '.fig'));
    close all

    
    
    
    figure; 
    hold on; 
    [ch,ch]=contourf(mean_CaT_skewness_matrix(:,:,2),20); colorbar;
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('Cat skeness plot: even beats','fontname','Times','fontsize',20); 
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'CaT_skewness_even',  '.fig'));
    close all
    
    figure; 
    hold on; 
    [ch,ch]=contourf(mean_CaT_skewness_matrix(:,:,1),20); colorbar;
    set(ch,'edgecolor','none');
    plot(x_cord,y_cord,'k'); 
    set(gca,'FontSize',20,'FontName','Times');
    title('CaT skeness plot: odd beats','fontname','Times','fontsize',20); 
    hold off; 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'CaT_skewness_odd',  '.fig'));
    close all
    
    
end  
     

