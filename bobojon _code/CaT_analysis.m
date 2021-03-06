function [CaT_infor_cell, Area_CaT_cell, CaT_skewness_cell, signal_swing_matrix, alternans_index_matrix,alternans_score_matrix, peak_alternans_index_matrix, peak_alternans_score_matrix,number_CaT_matrix,number_upstroke_matrix, upstroke_cell]=CaT_analysis (frameperiod,known_pacing,bined_image, plot_index)

total_row_pix = size(bined_image,1); %extracting row
total_col_pix = size(bined_image,2); %extracting column 


%% now we start with our analysis tool by going through each bined pixel and look for information 
%  some information are stored as cells of structure  
CaT_infor_cell = cell(total_row_pix, total_col_pix); 
Area_CaT_cell = cell(total_row_pix,total_col_pix); 
CaT_skewness_cell =  cell(total_row_pix,total_col_pix);

%upstroke information is stored in the a cell of vector 
upstroke_cell = cell(total_row_pix,total_col_pix); 

%some are store as matrix 
signal_swing_matrix = zeros(total_row_pix,total_col_pix); 
alternans_index_matrix = zeros(total_row_pix, total_col_pix); 
alternans_score_matrix = zeros(total_row_pix,total_col_pix); 
number_CaT_matrix = zeros(total_row_pix, total_col_pix); 
number_upstroke_matrix = zeros(total_row_pix, total_col_pix); 
number_to_bin=20

for r = 1: 1:total_row_pix 
    for c = 1:1:total_col_pix 
        %disp(strcat('r is ',num2str(r),' and  c is ',num2str(c)));
        trace1 = bined_image(r,c,:); 
        if sum(trace1) ~=0  % If this is within the tissue 
           trace1 = reshape(trace1,1,length(trace1));  
           trace1 = -trace1; 
          
           %disp(strcat('this pixel has row index',num2str(r),'and col index',num2str(c)))
            %% Get rough upstroke
            %Get upstroke first 
             %step1: filter the trace roughly to identify upstroke time point 
            [butter_coef_b, butter_coef_a] = Butter_filter_design(35,50, frameperiod); %% original value:Butter_filter_design(60,75, 1.96); changed for ratiometric calcium on 29.11.2015  
            row_filtered_sig1 = ButterFilter(trace1,butter_coef_b, butter_coef_a);
            filtered_sig1 = smooth(row_filtered_sig1,15);
            
            
            %% step2: idenfication of upstroke time point (moderate threshold applied)
         
            [upstroke_time_points,rapid_depolar_end] = Get_rough_upstroke_CaT(filtered_sig1',frameperiod,0,trace1,0.0045,0.20);  %Get_rough_upstroke_CaT(filtered_sig1',frameperiod,0,trace1,0.0045,0.2);changed for ratiometric calcium on 29.11.2015 
            wrong_upstroke_index = find(rapid_depolar_end-upstroke_time_points<1); 
            upstroke_time_points(wrong_upstroke_index)=[]; 
            rapid_depolar_end(wrong_upstroke_index)=[];
            if isempty(upstroke_time_points)~=1  % if we manage to find any action potential 
                %% Get baseline 
                % Step 3: get baseline for between each two upstrokes 
                baseline = zeros(length(upstroke_time_points),1); 
                for i = 1:length(upstroke_time_points)-1
                    
                   
                    %plotting histogram of the signal added: 30.11.2015
%                     number_to_bin= (max(filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1)))-...
%                     min(filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1))))./10;
%                      if number_to_bin<3
%                         number_to_bin= (max(filtered_sig1(upstroke_time_points(end):end))-...
%                         min(filtered_sig1(upstroke_time_points(end):end))).*10;
%                     end 
% 
%                     if number_to_bin<6
%                         number_to_bin=8; 
%                     end 
                     
                     % histogram
                    [accum_freq, bin_centre] = hist(filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1)),number_to_bin); 
                    total_freq = sum(accum_freq); 
                    sum_freq = 0 ; 
                    bin_num = 1;
                    while sum_freq <total_freq*0.3
                        sum_freq = sum_freq+accum_freq(bin_num); 
                        bin_num = bin_num+1; 
                    end 
                    
                    [~, max_freq_num] = max(accum_freq(1:bin_num-1));
                    
                    %range of bin to fit the normal distribution length
                    %(bin centre) max_freq_num 
                    
                    if length(bin_centre)>max_freq_num+4 && length(bin_centre)>max_freq_num*2
                        cut_off= max(bin_centre(max_freq_num*2),bin_centre(max_freq_num+4));
                    elseif length(bin_centre)>max_freq_num+4
                        cut_off= bin_centre(max_freq_num+4); 
                    elseif length(bin_centre)>max_freq_num*4
                        cut_off= bin_centre(max_freq_num*4); 
                    else 
                        cut_off = bin_centre(end); 
                    end 
                    
                    %fit normal distribution over the bins identifies
                    normal_fitting_vector = filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1));
                    normal_fitting_vector(normal_fitting_vector>cut_off)=[];
                    if length(normal_fitting_vector)>3
                        pd= fitdist(normal_fitting_vector,'normal'); 
                        baseline(i) = mean(pd);
                    else 
                        baseline(i) = filtered_sig1(upstroke_time_points(i)); 
                    end
                    if abs(baseline(i)-filtered_sig1(upstroke_time_points(i)))>0.4*abs(filtered_sig1(rapid_depolar_end(i))-filtered_sig1(upstroke_time_points(i)))
                       baseline(i) = (baseline(i)+filtered_sig1(upstroke_time_points(i)))/2;
                    end 
                end 
                
                
                
                %For the last CaT in the trace 
                
                % added on the 30.11.2015 for ratiometric calcium 
                
%                 number_to_bin= (max(filtered_sig1(upstroke_time_points(end):end))-...
%                 min(filtered_sig1(upstroke_time_points(end):end)))./10;
%             
%                  if number_to_bin<3
%                     number_to_bin= (max(filtered_sig1(upstroke_time_points(end):end))-...
%                     min(filtered_sig1(upstroke_time_points(end):end))).*10;
%                  end 
%         
%                  if number_to_bin<6
%                     number_to_bin=8; 
%                  end 
                 
                 
                [accum_freq, bin_centre] = hist(filtered_sig1(upstroke_time_points(end):end),number_to_bin); 
                total_freq = sum(accum_freq); 
                sum_freq = 0 ; 
                bin_num = 1;
                
                while sum_freq <total_freq*0.3
                      sum_freq = sum_freq+accum_freq(bin_num); 
                      bin_num = bin_num+1; 
                end 
                
                [~, max_freq_num] = max(accum_freq(1:bin_num-1)); 
                 if length(bin_centre)>max_freq_num+4 && length(bin_centre)>max_freq_num*2
                    cut_off= max(bin_centre(max_freq_num*2),bin_centre(max_freq_num+4));
                elseif length(bin_centre)>max_freq_num+4
                    cut_off= bin_centre(max_freq_num+4); 
                elseif length(bin_centre)>max_freq_num*4
                    cut_off= bin_centre(max_freq_num*4); 
                else 
                    cut_off = bin_centre(end); 
                 end 
                 
                normal_fitting_vector = filtered_sig1(upstroke_time_points(end):end); 
                normal_fitting_vector(normal_fitting_vector>cut_off)=[]; 
                
                if length(normal_fitting_vector)>3
                    pd= fitdist(normal_fitting_vector,'normal');
                    baseline(end) = mean(pd);
                else 
                    baseline(end) = filtered_sig1(upstroke_time_points(end));
                end 
                                               
                if abs(baseline(end)-filtered_sig1(upstroke_time_points(end)))>0.4*abs(filtered_sig1(rapid_depolar_end(end))-filtered_sig1(upstroke_time_points(end)))
                       baseline(end) = trace1(upstroke_time_points(end));
                end 
                %check with a plot 
%                 figure ; 
%                 time = 1:frameperiod:length(filtered_sig1).*frameperiod;
%                 hold on; 
%                 plot(time, trace1,'b'); 
%                 plot(time, filtered_sig1,'r')
% 
%                 for j = 1:length(upstroke_time_points)-1
%                     time_interval_base = frameperiod*upstroke_time_points(j):frameperiod:frameperiod*upstroke_time_points(j+1); 
%                     plot(time_interval_base, ones(1,length(time_interval_base)).*baseline(j),'k'); 
% 
%                     % check with the upstroke plot 
%                    % plot(upstroke_time_points(j)*frameperiod,trace1(upstroke_time_points(j)),'b*'); 
% 
%                 end 
% 
%                 time_interval_base = frameperiod*upstroke_time_points(end):frameperiod:frameperiod*length(trace1); 
%                 plot(time_interval_base, ones(1,length(time_interval_base)).*baseline(end),'k'); 
%                 %plot(upstroke_time_points(end)*frameperiod,trace1(upstroke_time_points(end)),'b*'); 
%                 xlabel('time (ms)'); 
%                 ylabel('vol sig'); 
%                 legend('row trace', 'filtered trace','baseline');

%                pause
                %% Estimate signal swing 

                % find out the 95 percentile value as the maximum signal value  
                signal_max = zeros(1,length(upstroke_time_points)); 
                for l = 1:length(upstroke_time_points)-1
                    AP_chunk = trace1(upstroke_time_points(l): upstroke_time_points(l+1));
                    quartile_value = quantile(AP_chunk, 0.95);
                    signal_max(l) = quartile_value;
                end 

                AP_chunk = trace1(upstroke_time_points(end): end);
                quartile_value = quantile(AP_chunk, 0.95);
                signal_max(end) = quartile_value; 

                % check the max signal with a plot 
                % figure; 
                % hold on; 
                % plot(time, trace1,'b'); 
                % plot(time, filtered_sig1,'r')
                % for m = 1:length(upstroke_time_points)-1
                %     time_interval_base = 1.96*upstroke_time_points(m):1.96:1.96*upstroke_time_points(m+1); 
                %     plot(time_interval_base, ones(1,length(time_interval_base)).*signal_max(m),'k'); 
                % end 
                % time_interval_base = 1.96*upstroke_time_points(end):1.96:1.96*length(trace1); 
                % plot(time_interval_base, ones(1,length(time_interval_base)).*signal_max(end),'k'); 
                % xlabel('time (ms)'); 
                % ylabel('vol sig'); 
                % legend('row trace','filtered trace','maximum sig'); 
                % hold off; 

                
                % step 5: calculate signal swing 
                signal_swing = abs((signal_max'-baseline)./baseline);
                average_swing = mean(signal_swing); 
                
           %% Based on the initial estimation of signal swing, classify these signals and recalculate the upstroke, basline ect

                % only process signal with signal swing bigger than 0.3%
                     
                
                if average_swing >0.003                   
                
                    %% case1: for weakish signal: Correct upstroke point again if the signal is too week 
                    if average_swing<0.006 
                        [upstroke_time_points,rapid_depolar_end] = Get_rough_upstroke_CaT(filtered_sig1',frameperiod,0,trace1,0.0045,0.25); % did not hve 0.0045 before 
                        wrong_upstroke_index = find(rapid_depolar_end-upstroke_time_points<1); 
                        upstroke_time_points(wrong_upstroke_index)=[]; 
                        if isempty(upstroke_time_points)~=1  % if we manage to find any action potential 
                            %% Get baseline 
                            %get baseline for between each two upstrokes 
                            baseline = zeros(length(upstroke_time_points),1); 
                            for i = 1:length(upstroke_time_points)-1
                                
                                
                                %number_to bin section added: on date 29th
                                %Novmber 2015
%                                 number_to_bin= (max(filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1)))-min(filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1))))./10;
%                                 if number_to_bin<3
%                                 number_to_bin= (max(filtered_sig1(upstroke_time_points(end):end))-min(filtered_sig1(upstroke_time_points(end):end))).*10;
%                                 end 
% 
%                                 if number_to_bin<6
%                                 number_to_bin=8; % replacing 8 by 10
%                                 end 
                                
                                % 20 replaced by number_to_bin_variable
                                [accum_freq, bin_centre] = hist(filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1)),number_to_bin); 
                                
                                % Identify the bin with peak value within the bottom 1/3 of the signal 
                   
                                total_freq = sum(accum_freq); 
                                sum_freq = 0 ; 
                                bin_num = 1;
                                
                                
                                while sum_freq <total_freq*0.3
                                    sum_freq = sum_freq+accum_freq(bin_num); 
                                    bin_num = bin_num+1; 
                                end 
                                
                                [~, max_freq_num] = max(accum_freq(1:bin_num-1));
                                
                                 %Identify range of bins to fit normal distribution 
                                 % length(bin_centre)
                                 % max_freq_num
                                
                           
                                if length(bin_centre)>max_freq_num+4 && length(bin_centre)>max_freq_num*2
                                    cut_off= max(bin_centre(max_freq_num*2),bin_centre(max_freq_num+4));
                                elseif length(bin_centre)>max_freq_num+4
                                    cut_off= bin_centre(max_freq_num+4); 
                                elseif length(bin_centre)>max_freq_num*4
                                    cut_off= bin_centre(max_freq_num*4); 
                                else 
                                    cut_off = bin_centre(end); 
                                end 
                                
                                
                                 % Fit a normal distribution over the bins identified (these bins should represent roughly the baseline of the signal)             
                   
                                normal_fitting_vector = filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1));
                                normal_fitting_vector(normal_fitting_vector>cut_off)=[];
                                if length(normal_fitting_vector)>3
                                    pd= fitdist(normal_fitting_vector,'normal'); 
                                    baseline(i) = mean(pd);
                                else 
                                    baseline(i) = filtered_sig1(upstroke_time_points(i)); 
                                end
                                
                                % if there is no significant 'baseline' data point in the signal, the the baseline can be identified to be too high, correct it using upstroke time point  
                   
                                if abs(baseline(i)-filtered_sig1(upstroke_time_points(i)))>0.4*abs(filtered_sig1(rapid_depolar_end(i))-filtered_sig1(upstroke_time_points(i)))
                                   baseline(i) = (baseline(i)+filtered_sig1(upstroke_time_points(i)))/2;
                                end 
                            end 
                
                             
                            %For the last CaT in the trace
                            
                            % number to bin section is added here as well
%                             number_to_bin= (max(filtered_sig1(upstroke_time_points(end):end))-...
%                              min(filtered_sig1(upstroke_time_points(end):end)))./10;
%             
%                             if number_to_bin<3
%                                number_to_bin= (max(filtered_sig1(upstroke_time_points(end):end))-...
%                                min(filtered_sig1(upstroke_time_points(end):end))).*10;
%                             end 
%         
%                             if number_to_bin<6
%                                number_to_bin=8; %replacing 10 by  8
%                             end 

                            [accum_freq, bin_centre] = hist(filtered_sig1(upstroke_time_points(end):end),number_to_bin); 
                            total_freq = sum(accum_freq); 
                            sum_freq = 0 ; 
                            bin_num = 1;
                            
                            
                            while sum_freq <total_freq*0.3
                                  sum_freq = sum_freq+accum_freq(bin_num); 
                                  bin_num = bin_num+1; 
                            end 
                            
                            
                            [~, max_freq_num] = max(accum_freq(1:bin_num-1)); 
                             if length(bin_centre)>max_freq_num+4 && length(bin_centre)>max_freq_num*2
                                cut_off= max(bin_centre(max_freq_num*2),bin_centre(max_freq_num+4));
                            elseif length(bin_centre)>max_freq_num+4
                                cut_off= bin_centre(max_freq_num+4); 
                            elseif length(bin_centre)>max_freq_num*4
                                cut_off= bin_centre(max_freq_num*4); 
                            else 
                                cut_off = bin_centre(end); 
                             end 
                            
                             
                             
                            normal_fitting_vector = filtered_sig1(upstroke_time_points(end):end); 
                            normal_fitting_vector(normal_fitting_vector>cut_off)=[]; 
                
                            if length(normal_fitting_vector)>3
                                pd= fitdist(normal_fitting_vector,'normal');
                                baseline(end) = mean(pd);
                            else 
                                baseline(end) = filtered_sig1(upstroke_time_points(end));
                            end 
                                               
                            if abs(baseline(end)-filtered_sig1(upstroke_time_points(end)))>0.4*abs(filtered_sig1(rapid_depolar_end(end))-filtered_sig1(upstroke_time_points(end)))
                                   baseline(end) = trace1(upstroke_time_points(end));
                            end 
%                             check with a plot 
%                             figure ; 
%                             time = 1:frameperiod:length(filtered_sig1).*frameperiod;
%                             hold on; 
%                             plot(time, trace1,'b'); 
%                             plot(time, filtered_sig1,'r')
%             
%                             for j = 1:length(upstroke_time_points)-1
%                                 time_interval_base = frameperiod*upstroke_time_points(j):frameperiod:frameperiod*upstroke_time_points(j+1); 
%                                 plot(time_interval_base, ones(1,length(time_interval_base)).*baseline(j),'k'); 
%             
%                                 % check with the upstroke plot 
%                                % plot(upstroke_time_points(j)*frameperiod,trace1(upstroke_time_points(j)),'b*'); 
%             
%                             end 
%             
%                             time_interval_base = frameperiod*upstroke_time_points(end):frameperiod:frameperiod*length(trace1); 
%                             plot(time_interval_base, ones(1,length(time_interval_base)).*baseline(end),'k'); 
%                             %plot(upstroke_time_points(end)*frameperiod,trace1(upstroke_time_points(end)),'b*'); 
%                             xlabel('time (ms)'); 
%                             ylabel('vol sig'); 
%                             legend('row trace', 'filtered trace','baseline');
% 
%             %                pause
                            %% Step 4: Estimate signal swing 

                            % find out the 90 percentile value as the maximum signal value  
                            signal_max = zeros(1,length(upstroke_time_points)); 
                            for l = 1:length(upstroke_time_points)-1
                                AP_chunk = trace1(upstroke_time_points(l): upstroke_time_points(l+1));
                                quartile_value = quantile(AP_chunk, 0.90);
                                signal_max(l) = quartile_value;
                            end 

                            AP_chunk = trace1(upstroke_time_points(end): end);
                            quartile_value = quantile(AP_chunk, 0.90);
                            signal_max(end) = quartile_value; 

                            % check the max signal with a plot 
%                             figure; 
%                             hold on; 
%                             plot(time, trace1,'b'); 
%                             plot(time, filtered_sig1,'r')
%                             for m = 1:length(upstroke_time_points)-1
%                                 time_interval_base = 1.96*upstroke_time_points(m):1.96:1.96*upstroke_time_points(m+1); 
%                                 plot(time_interval_base, ones(1,length(time_interval_base)).*signal_max(m),'k'); 
%                             end 
%                             time_interval_base = 1.96*upstroke_time_points(end):1.96:1.96*length(trace1); 
%                             plot(time_interval_base, ones(1,length(time_interval_base)).*signal_max(end),'k'); 
%                             xlabel('time (ms)'); 
%                             ylabel('vol sig'); 
%                             legend('row trace','filtered trace','maximum sig'); 
%                             hold off; 
%                                  

                            % Step 5: calculation of signal swing
                            signal_swing = abs((signal_max'-baseline)./baseline);
                            average_swing = mean(signal_swing); 
                        end 
                    end 
                    %% Correction on upstroke time points using baseline
                    corrected_upstroke_point = zeros(1,length(upstroke_time_points)); 
                    for n = 1: length(upstroke_time_points)
                        upstroke_chunk  = trace1(upstroke_time_points(n):rapid_depolar_end(n)); 
                        intersection_point = find(upstroke_chunk<baseline(n),1,'last');
                        if isempty(intersection_point)==1 && n==1 % if there's no intersection with the baseline in the upstroke phase and that's the first action potential (or only action potential) 
                            corrected_upstroke_point(n) = upstroke_time_points(n); % no correction of the onset of upstroke point 
                        elseif isempty(intersection_point)==1&& n>1 % if there's no intersection point but this is not the first action potential in the trace 
                            upstroke_chunk = trace1(rapid_depolar_end(n-1):upstroke_time_points(n)); % look for intersection in between the end of last depolar phase and onset of this depolarisation 
                            intersection_point = find(upstroke_chunk<baseline(n),1,'last');
                            if isempty(intersection_point)~=1 
                               intersection_time_point = intersection_point+rapid_depolar_end(n-1);
                                corrected_upstroke_point(n) = max(intersection_time_point,upstroke_time_points(n)); 
                            else 
                                corrected_upstroke_point(n) = upstroke_time_points(n);
                            end

                        else 

                            intersection_time_point = intersection_point+upstroke_time_points(n);
                            corrected_upstroke_point(n) = max(intersection_time_point,upstroke_time_points(n)); 
                        end 
                        %plot(ceil(corrected_upstroke_point(n))*frameperiod, trace1(ceil(corrected_upstroke_point(n))),'ro'); % check with a plot need to be used test baseline plot with upstroke plot uncommented 
                    end 
                    %plot(upstroke_time_points.*frameperiod,trace1(upstroke_time_points),'k*'); 
                    %plot(ceil(corrected_upstroke_point).*frameperiod, trace1(ceil(corrected_upstroke_point)),'ko');
                    corrected_upstroke_point = ceil(corrected_upstroke_point);


                    % need to test the corrected upstroke point as if the upstroke point is
                    % accurately identified but it's rounded up wrong then use floor instead of
                    % ceil 
                    difference = rapid_depolar_end-corrected_upstroke_point; 
                    index_wrong = find(difference<=0);
                    corrected_upstroke_point(index_wrong) = corrected_upstroke_point(index_wrong)-1; 


                    %% Correction on the end of upstrok time point 
                    rapid_depolar_end_corrected = zeros(1,length(rapid_depolar_end)); 
                    for s = 1: length(corrected_upstroke_point)
                        corrected_upstroke_end = corrected_upstroke_point(s)+find(row_filtered_sig1(corrected_upstroke_point(s):rapid_depolar_end(s))>signal_max(s),1,'first');    % use 90 percentile to find end of upstroke 
                        if isempty(corrected_upstroke_end)==1
                            rapid_depolar_end_corrected(s) = rapid_depolar_end(s);
                  %          disp('no corrected depolar end point')
                        else 
                            rapid_depolar_end_corrected(s) = corrected_upstroke_end; 
                        end 
                    end 
                  %  plot(rapid_depolar_end.*frameperiod,row_filtered_sig1(rapid_depolar_end),'k*'); 
                  % plot(rapid_depolar_end_corrected.*frameperiod, row_filtered_sig1(rapid_depolar_end_corrected),'ko');

                    % check with a plot 
                    
                     %% check error
          
                    Problem_index=find(rapid_depolar_end_corrected-corrected_upstroke_point<=0);
                    rapid_depolar_end(Problem_index)=[];
                    rapid_depolar_end_corrected(Problem_index)=[];
                    corrected_upstroke_point(Problem_index)=[];
                    upstroke_time_points(Problem_index)=[];
                    baseline(Problem_index)=[];
                    
                    if length(upstroke_time_points)>0
                        %% Curve Fitting (cubic) 
                        % upstroke fitting
                        further_corrected_upstroke = zeros(1,length(corrected_upstroke_point));
                        trace_upstroke_fitted = zeros(1,length(trace1)); 
                        linear_coefficient = zeros(2,length(corrected_upstroke_point)); 
                        for q = 1:length(upstroke_time_points)-1
                            [fitted_Vm,linear_coef,further_corrected_upstroke(q)] = fitt_AP_upstroke(row_filtered_sig1(corrected_upstroke_point(q):corrected_upstroke_point(q+1)),trace1(corrected_upstroke_point(q):corrected_upstroke_point(q+1)),corrected_upstroke_point(q),rapid_depolar_end_corrected(q), baseline(q),frameperiod, 0,signal_swing(q)); 
                            trace_upstroke_fitted(corrected_upstroke_point(q):corrected_upstroke_point(q+1)) = fitted_Vm; 
                            linear_coefficient(:,q) = linear_coef;  

                        end 
                        [fitted_Vm,fit_coef,further_corrected_upstroke(end)] =fitt_AP_upstroke(row_filtered_sig1(corrected_upstroke_point(end):end),trace1(corrected_upstroke_point(end):end),corrected_upstroke_point(end),rapid_depolar_end_corrected(end), baseline(end),frameperiod, 0,signal_swing(end));
                        trace_upstroke_fitted(corrected_upstroke_point(end):end) = fitted_Vm; 
                        linear_coefficient(:,end) = fit_coef; 

                        trace_upstroke_fitted(1:corrected_upstroke_point(1)) =row_filtered_sig1(1:corrected_upstroke_point(1)); 

                        wrong_corrected_index = find(further_corrected_upstroke- rapid_depolar_end_corrected>0); 
                        further_corrected_upstroke(wrong_corrected_index) = corrected_upstroke_point(wrong_corrected_index); 

                        %repolarisation fitting 
                        trace_curve_fitted = trace_upstroke_fitted;

                        poly_coefficient = zeros(4,length(corrected_upstroke_point)); 
                        baseline_intersection_time_point = zeros(1,length(further_corrected_upstroke)); % this uses further corrected upstroke point as an reference point
                        for q = 1:length(upstroke_time_points)-1
                           [fitted_Vm,fit_coef,baseline_intersection] = fitt_CaT_curve(row_filtered_sig1(further_corrected_upstroke(q):further_corrected_upstroke(q+1)),trace1(further_corrected_upstroke(q):further_corrected_upstroke(q+1)),further_corrected_upstroke(q),rapid_depolar_end_corrected(q), further_corrected_upstroke(q+1),baseline(q),frameperiod, 0,signal_swing(q));
                           if isempty(baseline_intersection) ~=1
                              trace_curve_fitted(further_corrected_upstroke(q):further_corrected_upstroke(q+1)) = fitted_Vm; 
                              poly_coefficient(:,q) = fit_coef; 
                              baseline_intersection_time_point(q)= baseline_intersection;
                           else 
                               poly_coefficient(:,q) = NaN; 
                               %disp('no fitting can be done'); 
                           end 

                        end 


                        [fitted_Vm,fit_coef,baseline_intersection] = fitt_CaT_curve(row_filtered_sig1(further_corrected_upstroke(end):end),trace1(further_corrected_upstroke(end):end),further_corrected_upstroke(end),rapid_depolar_end_corrected(end), length(row_filtered_sig1),baseline(end),frameperiod, 0,signal_swing(end));
                        if isempty(baseline_intersection)~=1
                            baseline_intersection_time_point(end) = baseline_intersection; 
                            trace_curve_fitted(further_corrected_upstroke(end):end) = fitted_Vm; 
                            poly_coefficient(:,end) = fit_coef; 
                        else 
                            poly_coefficient(:,end)=NaN; 
                            %disp('no fitting can be done'); 
                        end

                        %% check the length of the baseline to see whether we have a decent baseline 
                        % check this only if we have more upstroke than 1        
                        % for the first beat to the second last beat check wehther we get enough baseline for each beat 
                        if length(further_corrected_upstroke)>1

                            upstroke_interval = further_corrected_upstroke(2:end)-further_corrected_upstroke(1:end-1);
                            baseline_length_interval = upstroke_interval-baseline_intersection_time_point(1:end-1);

                            non_significant_baseline_index = find(baseline_length_interval<=0);
                            % the upstroke time point has the same value as the end of repolarisation which may caused by illy defined upstroke point or end of repolarisation point 

                            % only do motion correction if we have got a trace with more
                            % than two action potentials 
                            if length(further_corrected_upstroke)>2
                                motion_correction_start_end_points = [baseline_intersection_time_point(1:end-2)+further_corrected_upstroke(1:end-2);further_corrected_upstroke(3:end)]';
                                motion_correction_start_end_points = [0,0; motion_correction_start_end_points];
                                motion_correction_before_after_baseline = [baseline(1:end-1),baseline(2:end)];
                                motion_correction_before_after_baseline = [0,0; motion_correction_before_after_baseline]; 


                                binary_vect_motion_correction = zeros(1,length(further_corrected_upstroke)); 
                                binary_vect_motion_correction(1) = 1; 
                                binary_vect_motion_correction(end) = 1; 
                                binary_vect_motion_correction([non_significant_baseline_index,non_significant_baseline_index+1])=1;              
                                non_significant_baseline_index_full=[non_significant_baseline_index,non_significant_baseline_index+1];
                        %% motion correction:
                        % to find the drift line 
        %                         figure; 
        %                         hold on; 
        %                         plot(1:frameperiod:frameperiod*length(trace_curve_fitted),trace_curve_fitted,'b');
        %                         for i = 1:length(further_corrected_upstroke)
        %                             if binary_vect_motion_correction(i)==0 
        %                                plot(frameperiod.*[motion_correction_start_end_points(i,1),motion_correction_start_end_points(i,2)],[motion_correction_before_after_baseline(i,1),motion_correction_before_after_baseline(i,2)],'r');
        %                             end 
        %                         end 

                               % find the correction line and hence do the motion correction 
                                motion_corrected_trace = trace_curve_fitted;

                                for i = 1:length(further_corrected_upstroke) 
                                    if binary_vect_motion_correction(i) ==0 
                                      baseline_drift_gradient = (motion_correction_before_after_baseline(i,2)-motion_correction_before_after_baseline(i,1))./(motion_correction_start_end_points(i,2)-motion_correction_start_end_points(i,1))./frameperiod;
                                      correction_line = motion_correction_before_after_baseline(i,1)+baseline_drift_gradient.*(frameperiod:frameperiod:frameperiod*(motion_correction_start_end_points(i,2)-motion_correction_start_end_points(i,1)));

                                      motion_corrected_trace(further_corrected_upstroke(i):baseline_intersection_time_point(i)+further_corrected_upstroke(i)) = trace_curve_fitted(further_corrected_upstroke(i):baseline_intersection_time_point(i)+further_corrected_upstroke(i))...
                                      ./correction_line (further_corrected_upstroke(i)-motion_correction_start_end_points(i,1):further_corrected_upstroke(i)+baseline_intersection_time_point(i)-motion_correction_start_end_points(i,1));
                                      motion_corrected_trace(baseline_intersection_time_point(i)+1+further_corrected_upstroke(i):further_corrected_upstroke(i+1))=1;

                                    end
                                end

                                % this is to remove the wrongly fitted AP at the end of the trace

                                special_remove_index = find(non_significant_baseline_index_full+1-length(further_corrected_upstroke)>=0);
                                non_special_remove_index =  find(non_significant_baseline_index_full+1-length(further_corrected_upstroke)<0);
                                % normalize the other APs in the trace 
                                motion_corrected_trace (1:further_corrected_upstroke(2)) = trace_curve_fitted(1:further_corrected_upstroke(2))./baseline(1); 
                                motion_corrected_trace (further_corrected_upstroke(end)+1:end) = trace_curve_fitted(further_corrected_upstroke(end)+1:end)./baseline(end);
                                % remove the part which was found no significant baseline 
                                if isempty(non_special_remove_index)~=1  
                                    for u = 1:length(non_special_remove_index)

                                           motion_corrected_trace(further_corrected_upstroke(non_significant_baseline_index_full(non_special_remove_index(u))):further_corrected_upstroke(non_significant_baseline_index_full(non_special_remove_index(u))+1)-1)=1; 

                                    end 
                                end
                                if isempty(special_remove_index)~=1

                                   motion_corrected_trace(further_corrected_upstroke(non_significant_baseline_index_full(special_remove_index(1))):end)=1;
                                end 
                            else 
                                motion_corrected_trace = trace_curve_fitted;
                                motion_corrected_trace(1:further_corrected_upstroke(2)-1) = trace_curve_fitted(1:further_corrected_upstroke(2)-1)./baseline(1); 
                                motion_corrected_trace(further_corrected_upstroke(2):end) = trace_curve_fitted(further_corrected_upstroke(2):end)./baseline(1); 
                            end 
                        else
                            motion_corrected_trace = trace_curve_fitted./baseline(1); 
                            non_significant_baseline_index=[];
                        end

                         %% F/F0 for row trace and filtered trace 
                        % This can be done for row trace, filtered trace or fitted trace
                        trace1_norm = zeros(1,length(trace1)); 
                        row_filtered_norm = zeros(1,length(row_filtered_sig1));
                        %trace_fitted_norm = zeros(1,length(trace_curve_fitted)); 
                        if length(upstroke_time_points)>1
                            trace1_norm(1:upstroke_time_points(2))= trace1(1:upstroke_time_points(2))./baseline(1); 
                            row_filtered_norm(1:upstroke_time_points(2))= row_filtered_sig1(1:upstroke_time_points(2))./baseline(1); 
                            %trace_fitted_norm(1:upstroke_time_points(2))= -trace_curve_fitted(1:upstroke_time_points(2))./baseline(1); 
                            for p=2: length(upstroke_time_points)-1
                                trace1_norm(upstroke_time_points(p):upstroke_time_points(p+1))= trace1(upstroke_time_points(p):upstroke_time_points(p+1))./baseline(p); 
                                row_filtered_norm(upstroke_time_points(p):upstroke_time_points(p+1))= row_filtered_sig1(upstroke_time_points(p):upstroke_time_points(p+1))./baseline(p); 
                                %trace_fitted_norm(upstroke_time_points(p):upstroke_time_points(p+1))= -trace_curve_fitted(upstroke_time_points(p):upstroke_time_points(p+1))./baseline(p); 
                            end 

                            trace1_norm(upstroke_time_points(end):end) = trace1(upstroke_time_points(end):end)./baseline(end);
                            row_filtered_norm(upstroke_time_points(end):end)= row_filtered_sig1(upstroke_time_points(end):end)./baseline(end); 
                            %trace_fitted_norm(upstroke_time_points(end):end)= -trace_curve_fitted(upstroke_time_points(end):end)./baseline(end); 
                        else 


                            trace1_norm(1:end) = trace1(1:end)./baseline; 
                            row_filtered_norm(1:end) = row_filtered_sig1(1:end)./baseline; 
                        end
                        %normalized_trace = row_filtered_norm;
                        %row_filtered_norm(1:end) = row_filtered_sig1; 
                        %trace1_norm = trace1;
        %                 figure; 
        %                 hold on; 
        %                 plot(trace1_norm,'b') 
        %                 plot(row_filtered_norm,'r'); 
        %                 %plot(trace_fitted_norm,'k'); 
        %                 hold off
        %                 pause
                        % after this rescaling, baselibe becomes -1 
                        %baseline = ones(1,length(corrected_upstroke_point)).*-1; 
                         %% Delete the beats without significant portion of baseline %  these beats either have illy defined end of repolarisation point or wrongly identified upstroke time point  

                        further_corrected_upstroke([non_significant_baseline_index,non_significant_baseline_index+1]) = []; 
                        baseline_intersection_time_point([non_significant_baseline_index,non_significant_baseline_index+1])=[];    
                        baseline([non_significant_baseline_index,non_significant_baseline_index+1])= [];    
                        signal_swing([non_significant_baseline_index,non_significant_baseline_index+1])=[];
                        rapid_depolar_end_corrected([non_significant_baseline_index,non_significant_baseline_index+1])=[]; 
                        if ~isempty(further_corrected_upstroke)             
                            %% Refit repolarisation phase after motion correction and normalisation   
                            %repolarisation refitting to get the coefficient  
                            trace_curve_refitted = motion_corrected_trace;
                            baseline_norm = ones(1,length(corrected_upstroke_point)).*1;
                            poly_coefficient = zeros(4,length(corrected_upstroke_point)); 
                            baseline_intersection_time_point = zeros(1,length(further_corrected_upstroke)); % this uses further corrected upstroke point as an reference point
                            for q = 1:length(further_corrected_upstroke)-1
                               [fitted_Vm,fit_coef,baseline_intersection] = fitt_CaT_curve(row_filtered_norm(further_corrected_upstroke(q):further_corrected_upstroke(q+1)),motion_corrected_trace(further_corrected_upstroke(q):further_corrected_upstroke(q+1)),further_corrected_upstroke(q),rapid_depolar_end_corrected(q), further_corrected_upstroke(q+1),baseline_norm(q),frameperiod, 0,signal_swing(q));
                               if isempty(baseline_intersection) ~=1
                                  trace_curve_refitted(further_corrected_upstroke(q):further_corrected_upstroke(q+1)) = fitted_Vm; 
                                  poly_coefficient(:,q) = fit_coef; 
                                  baseline_intersection_time_point(q)= baseline_intersection;
                               else 
                                   poly_coefficient(:,q) = NaN; 
                                  % disp('no fitting can be done'); 
                               end 

                            end 


                            [fitted_Vm,fit_coef,baseline_intersection] = fitt_CaT_curve(row_filtered_norm(further_corrected_upstroke(end):end),motion_corrected_trace(further_corrected_upstroke(end):end),further_corrected_upstroke(end),rapid_depolar_end_corrected(end), length(motion_corrected_trace),baseline_norm(end),frameperiod, 0,signal_swing(end));
                            if isempty(baseline_intersection)~=1
                                baseline_intersection_time_point(end) = baseline_intersection; 
                                trace_curve_refitted(further_corrected_upstroke(end):end) = fitted_Vm; 
                                poly_coefficient(:,end) = fit_coef; 
                            else 
                                poly_coefficient(:,end)=NaN; 
                                %disp('no fitting can be done'); 
                            end                            

                            %% Caculate APD 
                            APD50 = zeros(1,length(further_corrected_upstroke)); 
                            APD80 = zeros(1,length(further_corrected_upstroke)); 
                            AP_peak_value = zeros(1,length(further_corrected_upstroke));
                            for r2 = 1:length(further_corrected_upstroke)-1
                                [APD,AP_peak_value(r2)] = Find_all_APD(trace_curve_refitted(further_corrected_upstroke(r2):further_corrected_upstroke(r2+1)),further_corrected_upstroke(r2), rapid_depolar_end_corrected(r2),baseline_norm(r2),frameperiod, poly_coefficient(:,r2),0); %the baseline here is -1
                                if isempty(APD.APD80)~=1
                                    APD80(r2)= APD.APD80; 
                                end
                                 if isempty(APD.APD50)~=1
                                    APD50(r2)= APD.APD50; 
                                 end      

                            end 
                            [APD,AP_peak_value(end)] = Find_all_APD(trace_curve_refitted(further_corrected_upstroke(end):end),further_corrected_upstroke(end),rapid_depolar_end_corrected(end),baseline_norm(end),frameperiod, poly_coefficient(:,end),0); %the baseline here is -1
                            if isempty(APD.APD80)~=1
                                APD80(end)= APD.APD80; 
                            end
                            if isempty(APD.APD50)~=1
                                APD50(end)= APD.APD50; 
                            end      

                            %% calculate number of AP found 
                            number_AP = length(find(APD50>0)); 
                            number_upstroke = length(further_corrected_upstroke);
                            % correct for wrongly identified upstroke  



                            %% calculate area under the curve and the ratio 
                            area_under_AP_vect = zeros(1,length(further_corrected_upstroke)); 
                            area_ratio_vect =  zeros(1,length(further_corrected_upstroke)); 
                            for s = 1:length(further_corrected_upstroke)-1
                                [area_under_AP,area_ratio]= find_area_APD(trace_curve_refitted(further_corrected_upstroke(s):further_corrected_upstroke(s+1)),baseline_norm(s),frameperiod,baseline_intersection_time_point(s),APD80(s),AP_peak_value(s), 0 );
                                area_under_AP_vect(s) = area_under_AP; 
                                area_ratio_vect(s) = area_ratio;

                            end 
                            [area_under_AP,area_ratio]= find_area_APD(trace_curve_refitted(further_corrected_upstroke(end):end),baseline_norm(end),frameperiod,baseline_intersection_time_point(end),APD80(end), AP_peak_value(end),0 );
                            area_under_AP_vect(end) = area_under_AP; 
                            area_ratio_vect(end) = area_ratio;



                            %% calculate Skewness of the AP trace (for fitted trace or row trace)
                            AP_skewness = zeros(1,length(further_corrected_upstroke)); 

                            for s = 1:length(further_corrected_upstroke)

                                AP_skewness(s) = skewness(trace_curve_refitted(further_corrected_upstroke(s):further_corrected_upstroke(s)+baseline_intersection_time_point(s)-1)); 

                            end 

                            % and correct for APD50 and APD80 area ratio and
                            % skewness, remove the wrongly ideentified beats 
                            further_corrected_upstroke (find(APD50==0))=[];
                            APD80(find(APD50==0))=[];
                            area_under_AP_vect(find(APD50==0))=[]; 
                            area_ratio_vect(find(APD50==0))=[]; 
                            AP_skewness(find(APD50==0))=[];
                            APD50(find(APD50==0))=[]; 


                            %% Test for alternans 
                            if isempty(APD80)~=1
                                [alternans_index,alternans_score,peak_alternans_index,peak_alternans_score] = test_alternans_CaT(AP_peak_value, further_corrected_upstroke, APD80, area_under_AP_vect, known_pacing);
                            else 
                                alternans_index= 0; 
                                alternans_score =0;
                            end 
                        else
                            APD80 = []; 
                            APD50 = []; 
                            area_under_AP_vect = []; 
                            area_ratio_vect = []; 
                            AP_skewness = []; 
                            alternans_index = 0;  
                            peak_alternans_score=0;
                            peak_alternans_index=0;
                            average_swing = 0; 
                            alternans_score =0; 
                            number_AP =0; 
                            number_upstroke = 0; 
                            further_corrected_upstroke = []; 
                        end 

                             
                    else 
                        APD80 = []; 
                        APD50 = []; 
                        area_under_AP_vect = []; 
                        area_ratio_vect = []; 
                        AP_skewness = []; 
                        alternans_index = 0;  
                        peak_alternans_score=0;
                        peak_alternans_index=0;
                        average_swing = 0; 
                        alternans_score =0; 
                        number_AP =0; 
                        number_upstroke = 0; 
                        further_corrected_upstroke = []; 
                    end                        
                        
                else 
                    APD80 = []; 
                    APD50 = []; 
                    area_under_AP_vect = []; 
                    area_ratio_vect = []; 
                    AP_skewness = []; 
                    alternans_index = 0;  
                    peak_alternans_score=0;
                    peak_alternans_index=0;
                    average_swing = 0; 
                    alternans_score =0; 
                    number_AP =0; 
                    number_upstroke = 0; 
                    further_corrected_upstroke = []; 
                end 
                    
           else
                APD80 = []; 
                APD50 = []; 
                area_under_AP_vect = []; 
                area_ratio_vect = []; 
                AP_skewness = []; 
                alternans_index = 0; 
                peak_alternans_score=0;
                peak_alternans_index=0;
                average_swing = 0; 
                alternans_score =0; 
                number_AP =0; 
                number_upstroke = 0; 
                further_corrected_upstroke = []; 
                
            end 
        else 
                               
                APD80 = []; 
                APD50 = []; 
                area_under_AP_vect = []; 
                area_ratio_vect = []; 
                AP_skewness = []; 
                alternans_index = 0;
                peak_alternans_score=0;
                peak_alternans_index=0;
                average_swing = 0;
                alternans_score=0; 
                number_AP=0;
                number_upstroke = 0; 
                further_corrected_upstroke = []; 
        end
        % start output here as a cell of structure 
        %create structures to store  in cells
        APD_structure= struct('APD50',APD50, 'APD80',APD80); 
        AP_area_structure = struct('area_under_AP_vect',area_under_AP_vect,'area_ratio_vect',area_ratio_vect);            
        CaT_infor_cell{r,c} = APD_structure; 
        Area_CaT_cell{r,c} = AP_area_structure; 
        CaT_skewness_cell{r,c} = AP_skewness;
        upstroke_cell{r,c}= further_corrected_upstroke;
        %some are store as matrix 
        signal_swing_matrix(r,c) = average_swing; 
        alternans_index_matrix(r,c) = alternans_index; 
        alternans_score_matrix(r,c) = alternans_score; 
        peak_alternans_index_matrix(r,c) = peak_alternans_index; 
        peak_alternans_score_matrix(r,c) = peak_alternans_score; 
        number_CaT_matrix(r,c) = number_AP; 
        number_upstroke_matrix(r,c)= number_upstroke;
     
        
    end 
end 

%% plot 

            

