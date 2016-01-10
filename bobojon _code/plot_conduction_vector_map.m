function plot_conduction_vector_map (upstroke_time_matrix,conduction_velocity_matrix,save_figure_index,fname,file_name_start_index,file_name_end_index,...
    total_row_number, total_col_number, pos,result_path,even_odd_index)

%plot the contour map for upstroke time matrix first 
if even_odd_index ==1
   plot_title= 'odd';
elseif even_odd_index ==0
   plot_title='even'
end 

figure; 
hold on; 

[ch,ch]=contourf(upstroke_time_matrix,30);colorbar;caxis([0,35]); 
set(ch,'edgecolor','none');

%plot vector field 

x_index_matirx = zeros(size(upstroke_time_matrix,1), size(upstroke_time_matrix,2));
for i = 1:size(upstroke_time_matrix,1)
    x_index_matrix(i,:) = i.*ones(1,size(upstroke_time_matrix,2)); 
end 
y_index_matirx = zeros(size(upstroke_time_matrix,1), size(upstroke_time_matrix,2));
for i = 1:size(upstroke_time_matrix,2)
    y_index_matrix(:,i) = i.*ones(size(upstroke_time_matrix,1),1); 
end 

u = conduction_velocity_matrix (:,:,1);
v = conduction_velocity_matrix(:,:,2);
[arrowx, arrowy] = plot_vector(y_index_matrix,x_index_matrix,u,v,1); 
grey = [0.4,0.4,0.4];
plot(arrowx,arrowy,'Color',grey,'Linewidth',1.5);
%title(strcat('conduction_vector',plot_title),'fontname','Times','fontsize',20); 
% get the contour of tissue out 

% x_cord = pos(:,1); 
% y_cord = pos(:,2); 
% 
% % get rid of anything out of bounds 
% x_cord(x_cord>total_row_number) = total_row_number; 
% x_cord(x_cord<0) = 0;
% y_cord(y_cord>total_col_number) = total_col_number; 
% y_cord(y_cord<0) = 0;
% 
% plot(x_cord,y_cord,'k'); 
hold off
if save_figure_index~=0 
    hgsave(gcf,strcat(result_path,'\matlab_figure_from_autanalysis\', fname(file_name_start_index:file_name_end_index),'conduction_vector',plot_title,'.fig'));
    close all
end    
