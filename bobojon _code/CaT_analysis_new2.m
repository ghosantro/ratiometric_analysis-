function [CaT_infor_cell, Area_CaT_cell, CaT_skewness_cell, signal_swing_matrix, alternans_index_matrix,alternans_score_matrix, peak_alternans_index_matrix, peak_alternans_score_matrix,number_CaT_matrix,number_upstroke_matrix, upstroke_cell]=CaT_analysis_new2 (frameperiod,known_pacing,bined_image,pacing_frequency,  plot_index)

total_row_pix = size(bined_image,1); 
total_col_pix = size(bined_image,2); 


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
peak_alternans_index_matrix = zeros(total_row_pix,total_col_pix); 
peak_alternans_score_matrix = zeros(total_row_pix,total_col_pix); 
alternans_score_matrix = zeros(total_row_pix,total_col_pix); 
number_CaT_matrix = zeros(total_row_pix, total_col_pix); 
number_upstroke_matrix = zeros(total_row_pix, total_col_pix); 

%% need to check 
pacing_interval=1000/pacing_frequency; 
%         
for r = 1: 1:total_row_pix 
    for c = 1:1:total_col_pix 
        
        disp(strcat('r is ',num2str(r),' and  c is ',num2str(c)));
        trace1 = bined_image(r,c,:);   
         
        extreme_alternans_index=0;

        if sum(trace1) ~=0  % If this is within the tissue 
           trace1 = reshape(trace1,1,length(trace1));  
           trace1=-trace1;

            %% Get rough upstroke time point, baseline and signal maximal for singal swing estimation  
            %Get upstroke first 

            %stpe1: filter the trace roughly to identify upstroke time point 
            [butter_coef_b, butter_coef_a] = Butter_filter_design(35,50, frameperiod);
            row_filtered_sig1 = ButterFilter(trace1,butter_coef_b, butter_coef_a);
            filtered_sig1 = smooth(row_filtered_sig1,15);
            %check step 1 with an plot 
            if plot_index==1
                figure; 
                hold on; 
                plot(1:frameperiod:frameperiod*length(row_filtered_sig1),trace1,'r'); 
                plot(1:frameperiod:frameperiod*length(row_filtered_sig1),row_filtered_sig1,'b'); 
                xlabel('time(ms)'); 
                ylabel('Vm'); 
                legend('row trace','filtered trace');
                title(' Check whether filtering is correct');
                hold off;
                pause
            end 

            %step2: idenfication of upstroke time point (moderate threshold applied)
            [upstroke_time_points,rapid_depolar_end] = Get_rough_upstroke_CaT(filtered_sig1',frameperiod,plot_index,trace1,0.0030,0.40);
            wrong_upstroke_index = find(rapid_depolar_end-upstroke_time_points<1); 
            upstroke_time_points(wrong_upstroke_index)=[]; 
            rapid_depolar_end(wrong_upstroke_index)=[];
            if isempty(upstroke_time_points)~=1  % if we manage to find any action potential 

                %step3: Find baseline for each identified CaT (between each two upstrokes) 
                baseline = zeros(length(upstroke_time_points),1); 

                for i = 1:length(upstroke_time_points)-1

                    % plot histogram of the signal find how the signal distributes 
                    number_to_bin= (max(filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1)))-...
                    min(filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1))))./10;
                     if number_to_bin<3
                        number_to_bin= (max(filtered_sig1(upstroke_time_points(end):end))-...
                        min(filtered_sig1(upstroke_time_points(end):end))).*10;
                    end 

                    if number_to_bin<6
                        number_to_bin=8; 
                    end 
                    %hist(filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1)), number_to_bin) ;
                    [accum_freq, bin_centre] = hist(filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1)), number_to_bin);

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
%                     length(bin_centre)
%                     max_freq_num
                    if length(bin_centre)>max_freq_num+4 && length(bin_centre)>max_freq_num*2
                        cut_off= min(bin_centre(max_freq_num*2),bin_centre(max_freq_num+4));
                    elseif length(bin_centre)>max_freq_num+4
                        cut_off= bin_centre(max_freq_num+4); 
                    elseif length(bin_centre)>max_freq_num*4
                        cut_off= bin_centre(max_freq_num*4) ;
                    else 
                        cut_off = bin_centre(end);
                    end 

                    % Fit a normal distribution over the bins identified (these bins should represent roughly the baseline of the signal)             
                    normal_fitting_vector = filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1));
                    normal_fitting_vector(normal_fitting_vector>cut_off)=[];
                    if length(normal_fitting_vector)>3
                        pd= fitdist(normal_fitting_vector,'normal'); 
                        baseline(i) = mean(pd);  
                        % the mean value of this normal distribution will be used as the baseline value 

                        % only uncomment this if you need to check this par of the routine          
                        %figure; hold on; hist(filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1)), number_to_bin) 
                        %histfit(normal_fitting_vector,8,'normal')
                        %plot([baseline(i),baseline(i)],[0,50],'r-');
                        %hold off; 
                        %close all
                    else 
                        baseline(i) = filtered_sig1(upstroke_time_points(i)); 
                    end

                    % if there is no significant 'baseline' data point in the signal, the the baseline can be identified to be too high, correct it using upstroke time point  
                    if abs(baseline(i)-filtered_sig1(upstroke_time_points(i)))>0.25*abs(filtered_sig1(rapid_depolar_end(i))-filtered_sig1(upstroke_time_points(i)))
                       baseline(i) = (baseline(i)+filtered_sig1(upstroke_time_points(i)))/2;
                    end 
                end

                % For the last CaT in the trace 
                number_to_bin= (max(filtered_sig1(upstroke_time_points(end):end))-...
                min(filtered_sig1(upstroke_time_points(end):end)))./10;
            
                 if number_to_bin<3
                    number_to_bin= (max(filtered_sig1(upstroke_time_points(end):end))-...
                    min(filtered_sig1(upstroke_time_points(end):end))).*10;
                 end 
        
                 if number_to_bin<6
                    number_to_bin=8; 
                 end 
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

                % check step 3 with a plot
                if plot_index==1
                    figure ; 
                    time = 1:frameperiod:length(filtered_sig1)*frameperiod;
                    hold on; 
                    plot(time, trace1,'b'); 
                    for j = 1:length(upstroke_time_points)-1
                        time_interval_base = frameperiod*upstroke_time_points(j):frameperiod:frameperiod*upstroke_time_points(j+1); 
                        plot(time_interval_base, ones(1,length(time_interval_base)).*baseline(j),'k');                 
                    end 

                    time_interval_base = frameperiod*upstroke_time_points(end):frameperiod:frameperiod*length(trace1); 
                    plot(time_interval_base, ones(1,length(time_interval_base)).*baseline(end),'k'); 
                    xlabel('time (ms)'); 
                    ylabel('vol sig'); 
                    legend('row trace', 'baseline');
                    title(' Check whether baselines are correctly identified');
                    pause
                end


                % Step4: identify signal maximal: find out the 95 percentile value as the maximum signal value   
                signal_max = zeros(1,length(upstroke_time_points)); 
                for l = 1:length(upstroke_time_points)-1
                    CaT_chunk = trace1(upstroke_time_points(l): upstroke_time_points(l+1));
                    quartile_value = quantile(CaT_chunk, 0.95);
                    signal_max(l) = quartile_value;
                end 

                CaT_chunk = trace1(upstroke_time_points(end): end);
                quartile_value = quantile(CaT_chunk, 0.95);
                signal_max(end) = quartile_value; 

                % check step 4 with a plot 
                if plot_index==1
                    figure; 
                    hold on; 
                    plot(time, trace1,'b'); 
                    plot(time, filtered_sig1,'r')
                    for m = 1:length(upstroke_time_points)-1
                        time_interval_base = frameperiod*upstroke_time_points(m):frameperiod:frameperiod*upstroke_time_points(m+1); 
                        plot(time_interval_base, ones(1,length(time_interval_base)).*signal_max(m),'k'); 
                         plot(time_interval_base, ones(1,length(time_interval_base)).*baseline(m),'k'); 
                    end 
                    time_interval_base = frameperiod*upstroke_time_points(end):frameperiod:frameperiod*length(trace1); 
                    plot(time_interval_base, ones(1,length(time_interval_base)).*signal_max(end),'k'); 
                    xlabel('time (ms)'); 
                    ylabel('vol sig'); 
                    legend('row trace','filtered trace','maximum sig'); 
                    title('Check baseline and signal maixmal');
                    hold off; 
                end 

                 % Step5: Calculate signal swing 
                signal_swing = abs((signal_max'-baseline)./baseline);
                average_swing = mean(signal_swing);        


                %% Based on the initial estimation of signal swing, classify these signals and recalculate the upstroke, basline ect

                % only process signal with signal swing bigger than 0.3%
                if average_swing>0.003

                    %% Case1: For weakish signal, recalculate upstroke point and baseline with more strict threshold   
                    if average_swing<0.006
                        %Re-identification of upstroke time point based on more strict threshold applied (especially the upstroke amplitude)
                        [upstroke_time_points,rapid_depolar_end] = Get_rough_upstroke_CaT(filtered_sig1',frameperiod,plot_index,trace1,0.0045,0.25);
                         wrong_upstroke_index = find(rapid_depolar_end-upstroke_time_points<1); 
                         upstroke_time_points(wrong_upstroke_index)=[]; 
                         rapid_depolar_end(wrong_upstroke_index)=[];
                         if isempty(upstroke_time_points)~=1  % if we manage to find any CaT

                            % Re-identification fo baseline (same threshold as the general case)  
                            baseline = zeros(length(upstroke_time_points),1); 
                            for i = 1:length(upstroke_time_points)-1
                                 number_to_bin= (max(filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1)))-...
                                min(filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1))))./10;
                                 if number_to_bin<3
                                    number_to_bin= (max(filtered_sig1(upstroke_time_points(end):end))-...
                                     min(filtered_sig1(upstroke_time_points(end):end))).*10;
                                 end 
        
                                if number_to_bin<6
                                    number_to_bin=10; 
                                end 
                                [accum_freq, bin_centre] = hist(filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1)),number_to_bin); 
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
                                normal_fitting_vector = filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1));
                                normal_fitting_vector(normal_fitting_vector>cut_off)=[];
                                if length(normal_fitting_vector)>3
                                    pd= fitdist(normal_fitting_vector,'normal'); 
                                    baseline(i) = mean(pd);
                                else 
                                    baseline(i) = filtered_sig1(upstroke_time_points(i)); 
                                end
                                if abs(baseline(i)-filtered_sig1(upstroke_time_points(i)))>0.25*abs(filtered_sig1(rapid_depolar_end(i))-filtered_sig1(upstroke_time_points(i)))
                                   baseline(i) = (baseline(i)+filtered_sig1(upstroke_time_points(i)))/2;
                                end 
                            end 
                             number_to_bin= (max(filtered_sig1(upstroke_time_points(end):end))-...
                             min(filtered_sig1(upstroke_time_points(end):end)))./10;
                         
                             if number_to_bin<3
                                number_to_bin= (max(filtered_sig1(upstroke_time_points(end):end))-...
                                min(filtered_sig1(upstroke_time_points(end):end))).*10;
                            end 
        
                            if number_to_bin<6 %was 3, changed to 6 on 30.11.2015 
                                number_to_bin=10; 
                            end 
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
                            % check re-identified baseline  
                            if plot_index==1
                                figure ; 
                                time = 1:frameperiod:length(filtered_sig1)*frameperiod;
                                hold on; 
                                plot(time, trace1,'b');                         
                                for j = 1:length(upstroke_time_points)-1
                                    time_interval_base = frameperiod*upstroke_time_points(j):frameperiod:frameperiod*upstroke_time_points(j+1); 
                                    plot(time_interval_base, ones(1,length(time_interval_base)).*baseline(j),'k'); 
                                end 

                                time_interval_base = frameperiod*upstroke_time_points(end):frameperiod:frameperiod*length(trace1); 
                                plot(time_interval_base, ones(1,length(time_interval_base)).*baseline(end),'k'); 
                                xlabel('time (ms)'); 
                                ylabel('vol sig'); 
                                legend('row trace', 'baseline');
                                title('Check reidentification of baseline in case of weak signal')
                                pause
                            end                    

                            % re-identification of CaT the 90 percentile value will be used as the 'maximum'signal value  
                            signal_max = zeros(1,length(upstroke_time_points)); 
                            for l = 1:length(upstroke_time_points)-1
                                CaT_chunk = trace1(upstroke_time_points(l): upstroke_time_points(l+1));
                                quartile_value = quantile(CaT_chunk, 0.9);
                                signal_max(l) = quartile_value;
                            end 
                            CaT_chunk = trace1(upstroke_time_points(end): end);
                            quartile_value = quantile(CaT_chunk, 0.9);
                            signal_max(end) = quartile_value; 

                            % check reidentification of maximal signal with a plot 
                            if plot_index==1
                                figure; 
                                hold on; 
                                plot(time, trace1,'b'); 
                                plot(time, filtered_sig1,'r')
                                for m = 1:length(upstroke_time_points)-1
                                    time_interval_base = frameperiod*upstroke_time_points(m):frameperiod:frameperiod*upstroke_time_points(m+1); 
                                    plot(time_interval_base, ones(1,length(time_interval_base)).*signal_max(m),'k'); 
                                     plot(time_interval_base, ones(1,length(time_interval_base)).*baseline(m),'k'); 
                                end 
                                time_interval_base = frameperiod*upstroke_time_points(end):frameperiod:frameperiod*length(trace1); 
                                plot(time_interval_base, ones(1,length(time_interval_base)).*signal_max(end),'k'); 
                                xlabel('time (ms)'); 
                                ylabel('vol sig'); 
                                legend('row trace','filtered trace','maximum sig'); 
                                title('Check reidentification of baseline and maximal signal')
                                hold off; 
                            end

                            %recalculate signal swing 
                            signal_swing = abs((signal_max'-baseline)./baseline);
                            average_swing = mean(signal_swing); 
                         end 
                    end

                    %% Case2: For good signal, recalculate upstroke point and baseline with more lose threshold and identify extreme CaT alternans for tailored processing    
                    if average_swing>0.2 % ???? May need to modify the threshold here 

                        % Re-idntification of upstroke time point with particular lose threshold for upstroke amplitude 
                        [upstroke_time_points,rapid_depolar_end] = Get_rough_upstroke_CaT(filtered_sig1',frameperiod,plot_index,trace1,0.0025,0.02);
                        wrong_upstroke_index = find(rapid_depolar_end-upstroke_time_points<1); 
                        upstroke_time_points(wrong_upstroke_index)=[]; 
                        rapid_depolar_end(wrong_upstroke_index)=[];

                        if isempty(upstroke_time_points)~=1  % if we manage to find any CAT
                            baseline = zeros(length(upstroke_time_points),1); 
                            % seperate the normal case and the case wtih extreme CaT alternans      
                            % if a pttern of alternating low/high upstroke is observed, this is  likely to be a case of extreme CaT alternans 
                            if length(upstroke_time_points>6) % only make sense to look for the extreme alternans case if there are more than 6 CaT identified 
                                upstroke_CaT_value=filtered_sig1(upstroke_time_points);
                                eval_temp=upstroke_CaT_value(2:2:end-1)-upstroke_CaT_value (3:2:end);
                                if length(find(eval_temp>0))==length(eval_temp)||length(find(eval_temp<0))==length(eval_temp)
                                    %alterns_extend = abs(1-upstroke_CaT_value(2:2:end)./upstroke_CaT_value(1:2:end-1));
                                    extreme_alternans_index=1; 

                                else 
                                    extreme_alternans_index=0; 
                                end                    

                                % Need to check that the small CaT is 'riding' on the big CaT release by looking for incomplete decay for the big CaT 
                                pass_test_index = []; 

                                if length(find(eval_temp>0))==length(eval_temp)% this is when the big realease is on the odd beat 

                                     for po=3:2:(length(upstroke_time_points)-1)      % then start looking from beat one for big CaT decay 
                                        CaT_decay_chunk = trace1(rapid_depolar_end(po):upstroke_time_points(po+1)); 
                                        CaT_decay_chunk = smooth(CaT_decay_chunk,20); 
                                        FD_CaT_decay = (CaT_decay_chunk(2:2:end)-CaT_decay_chunk(1:2:end-1))./frameperiod;
                                        FD_CaT_decay_90prct= prctile(FD_CaT_decay,90); 
                                        FD_CaT_decay(find(FD_CaT_decay>FD_CaT_decay_90prct))=[];

                                        if length(find(FD_CaT_decay>0))==0
                                            pass_test_index = [pass_test_index,1];   %if the the trace was just decay before the next upstroke (without reach resting level) then it is  the case of incomplete big CaT 
                                        else
                                            pass_test_index=[pass_test_index,0];
                                        end 
                                    end 

                                    if length(find(pass_test_index==0))~=0;
                                       extreme_alternans_index=0;
                                    end 
                                elseif length(find(eval_temp<0))==length(eval_temp) % this is when the big realease is on the even beat

                                    for po=4:2:(length(upstroke_time_points)-1)      % then start looking from beat two for big CaT decay 
                                        CaT_decay_chunk = trace1(rapid_depolar_end(po):upstroke_time_points(po+1)); 
                                        CaT_decay_chunk = smooth(CaT_decay_chunk,20); 
                                        FD_CaT_decay = (CaT_decay_chunk(2:2:end)-CaT_decay_chunk(1:2:end-1))./frameperiod;
                                        FD_CaT_decay_90prct= prctile(FD_CaT_decay,90); 
                                        FD_CaT_decay(find(FD_CaT_decay>FD_CaT_decay_90prct))=[];

                                        if length(find(FD_CaT_decay>0))==0
                                            pass_test_index = [pass_test_index,1] ;   %if the the trace was just decay before the next upstroke (without reach resting level) then it is  the case of incomplete big CaT 
                                        else
                                            pass_test_index=[pass_test_index,0];
                                        end 

                                    end                         
                                    if length(find(pass_test_index==0))~=0;
                                       extreme_alternans_index=0;
                                    end 

                                end 

                             else
                                extreme_alternans_index=0;
                            end 
                            % With special case identified perform a completely different routine compared to the other cases  
                            if extreme_alternans_index==1
                                %the first beat might not be a complete beat and
                                %therefore the value should be replaced with the
                                %the third one 
                                upstroke_CaT_value(1)=upstroke_CaT_value(3); 
                                if length(find(eval_temp>0))==length(eval_temp)  % if an alternans in upstroke CaT value is identified                   
                                   % baseline should be the upstroke value 
                                   upstroke_CaT_value (2:2:end)=upstroke_CaT_value(1:2:end-1);
                                elseif length(find(eval_temp<0))==length(eval_temp)
                                   upstroke_CaT_value (1:2:end-1)=upstroke_CaT_value(2:2:end);
                                end 


                                if length(upstroke_time_points)==length(upstroke_CaT_value)
                                    baseline = upstroke_CaT_value;
                                elseif length(upstroke_time_points) ==length(upstroke_CaT_value)+1
                                    baseline(1:end-1)= upstroke_CaT_value;
                                    baseline(end)=upstroke_CaT_value(end);
                                end 

                                if plot_index==1
                                    figure ; 
                                    time = 1:frameperiod:length(filtered_sig1)*frameperiod;
                                    hold on; 
                                    plot(time, trace1,'b');                         
                                    for j = 1:length(upstroke_time_points)-1
                                        time_interval_base = frameperiod*upstroke_time_points(j):frameperiod:frameperiod*upstroke_time_points(j+1); 
                                        plot(time_interval_base, ones(1,length(time_interval_base)).*baseline(j),'k'); 
                                    end 
                                    time_interval_base = frameperiod*upstroke_time_points(end):frameperiod:frameperiod*length(trace1); 
                                    plot(time_interval_base, ones(1,length(time_interval_base)).*baseline(end),'k'); 
                                    xlabel('time (ms)'); 
                                    ylabel('vol sig'); 
                                    legend('row trace', 'filtered trace','baseline');
                                    title('Check re-identification of baseline in extreme alternans case (good signal)');
                                    pause
                                end 

                                % Get maximual signal and re-calculate signal swing
                                signal_max = zeros(1,length(upstroke_time_points)); 
                                for l = 1:length(upstroke_time_points)-1
                                    CaT_chunk = trace1(upstroke_time_points(l): upstroke_time_points(l+1));
                                    quartile_value = quantile(CaT_chunk, 0.975);
                                    signal_max(l) = quartile_value;
                                end 

                                CaT_chunk = trace1(upstroke_time_points(end): end);
                                quartile_value = quantile(CaT_chunk, 0.975);
                                signal_max(end) = quartile_value; 

                                % check the max signal with a plot 
                                if plot_index==1 
                                    figure; 
                                    hold on; 
                                    plot(time, trace1,'b'); 
                                    plot(time, filtered_sig1,'r')
                                    for m = 1:length(upstroke_time_points)-1
                                        time_interval_base = frameperiod*upstroke_time_points(m):frameperiod:frameperiod*upstroke_time_points(m+1); 
                                        plot(time_interval_base, ones(1,length(time_interval_base)).*signal_max(m),'k'); 
                                         plot(time_interval_base, ones(1,length(time_interval_base)).*baseline(m),'k'); 
                                    end 
                                    time_interval_base = frameperiod*upstroke_time_points(end):frameperiod:frameperiod*length(trace1); 
                                    plot(time_interval_base, ones(1,length(time_interval_base)).*signal_max(end),'k'); 
                                    xlabel('time (ms)'); 
                                    ylabel('vol sig'); 
                                    legend('row trace','filtered trace','maximum sig'); 
                                    title('Re-identification of baseline and signal maximal in case of extreme CaT alternans (good signal)');
                                    hold off; 
                                end 

                                % recalculate signal swing 
                                signal_swing = abs((signal_max'-baseline)./baseline);
                                average_swing = mean(signal_swing); 

                                % Correction on upstroke time points using baseline

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
                                        if intersection_time_point <upstroke_time_points(n) 
                                           corrected_upstroke_point(n)= upstroke_time_points(n);
                                        else 
                                           corrected_upstroke_point(n) = max(intersection_time_point,upstroke_time_points(n)); 
                                        end 
                                    end 
                                   % plot(ceil(corrected_upstroke_point(n))*frameperiod, trace1(ceil(corrected_upstroke_point(n))),'ro'); % check with a plot need to be used test baseline plot with upstroke plot uncommented 
                                end 
                                
                
                                corrected_upstroke_point = ceil(corrected_upstroke_point);
                                if plot_index==1
                                    figure; 
                                    hold on; 
                                    plot(1:frameperiod:frameperiod*length(trace1),trace1,'k'); 
                                    plot(upstroke_time_points.*frameperiod,trace1(upstroke_time_points),'k*'); 
                                    plot(ceil(corrected_upstroke_point).*frameperiod, trace1(ceil(corrected_upstroke_point)),'ko');
                                    legend('raw trace','upstroke time', 'corrected upstroke time');
                                    title('Check correction on upstroke time point in extreme case of CaT alternans (good signal)'); 
                                    hold off; 
                                end 

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

                                if plot_index==1
                                    figure; 
                                    hold on; 
                                    plot(1:frameperiod:frameperiod*length(trace1), trace1,'k');
                                    plot(rapid_depolar_end.*frameperiod,trace1(rapid_depolar_end),'b*'); 
                                    plot(rapid_depolar_end_corrected.*frameperiod, trace1(rapid_depolar_end_corrected),'bo');
                                    hold off; 
                                    legend('raw trace','rapid upstroke end','corrected rapid upstroke end'); 
                                    % check with a plot 
                                    pause  
                                end 

                                % checke whether the correction induces error, if
                                % so delete the problematic part 
                                Problem_index=find(rapid_depolar_end_corrected-corrected_upstroke_point<=0);
                                rapid_depolar_end(Problem_index)=[];
                                rapid_depolar_end_corrected(Problem_index)=[];
                                corrected_upstroke_point(Problem_index)=[];
                                upstroke_time_points(Problem_index)=[];
                                baseline(Problem_index)=[];


                                % F/F0 for row trace,filtered trace 

                                trace1_norm = zeros(1,length(trace1)); 
                                row_filtered_norm = zeros(1,length(row_filtered_sig1));
                                %trace_fitted_norm = zeros(1,length(trace_curve_fitted)); 
                                if length(upstroke_time_points)>1
                                    trace1_norm(1:upstroke_time_points(2))= trace1(1:upstroke_time_points(2))./baseline(1); 
                                    row_filtered_norm(1:upstroke_time_points(2))= row_filtered_sig1(1:upstroke_time_points(2))./baseline(1); 
                                    %trace_fitted_norm(1:upstroke_time_points(2))= trace_curve_fitted(1:upstroke_time_points(2))./baseline(1); 
                                    for p=2: length(upstroke_time_points)-1
                                        trace1_norm(upstroke_time_points(p):upstroke_time_points(p+1))= trace1(upstroke_time_points(p):upstroke_time_points(p+1))./baseline(p); 
                                        row_filtered_norm(upstroke_time_points(p):upstroke_time_points(p+1))= row_filtered_sig1(upstroke_time_points(p):upstroke_time_points(p+1))./baseline(p); 
                                        %trace_fitted_norm(upstroke_time_points(p):upstroke_time_points(p+1))= trace_curve_fitted(upstroke_time_points(p):upstroke_time_points(p+1))./baseline(p); 
                                    end 

                                    trace1_norm(upstroke_time_points(end):end) = trace1(upstroke_time_points(end):end)./baseline(end);
                                    row_filtered_norm(upstroke_time_points(end):end)= row_filtered_sig1(upstroke_time_points(end):end)./baseline(end); 
                                    %trace_fitted_norm(upstroke_time_points(end):end)= trace_curve_fitted(upstroke_time_points(end):end)./baseline(end); 
                                else 


                                    trace1_norm(1:end) = trace1(1:end)./baseline; 
                                    %row_filtered_norm(1:end) = row_filtered_sig1(1:end)./baseline; 
                                end

                                % check F/F0 with a plot 
                                if plot_index==1               
                                    figure; 
                                    hold on; 
                                    plot(trace1_norm,'b') 
                                    plot(row_filtered_norm,'r'); 
                                    %plot(trace_fitted_norm,'k'); 
                                    legend('raw trace','filtered trace');
                                    title('check fitted trace after F/F0 in case of extreme CaT alternans (good signal)');
                                    hold off
                                    pause
                                end 

                                %refit after normalisation 
                                if length(upstroke_time_points)>0
                                   % fit the trace with straight line and polynomial curve 
                                   % upstroke fitting
                                   further_corrected_upstroke = zeros(1,length(corrected_upstroke_point));
                                   trace_upstroke_fitted = zeros(1,length(trace1)); 
                                   linear_coefficient = zeros(2,length(corrected_upstroke_point)); 
                                   for q = 1:length(upstroke_time_points)-1
                                       [fitted_Vm,linear_coef,further_corrected_upstroke(q)] = fitt_AP_upstroke(row_filtered_norm(corrected_upstroke_point(q):corrected_upstroke_point(q+1)),trace1_norm(corrected_upstroke_point(q):corrected_upstroke_point(q+1)),corrected_upstroke_point(q),rapid_depolar_end_corrected(q), 1,frameperiod, plot_index,signal_swing(q)); 
                                       trace_upstroke_fitted(corrected_upstroke_point(q):corrected_upstroke_point(q+1)) = fitted_Vm; 
                                       linear_coefficient(:,q) = linear_coef;  
                                   end 
                                    [fitted_Vm,fit_coef,further_corrected_upstroke(end)] =fitt_AP_upstroke(row_filtered_norm(corrected_upstroke_point(end):end),trace1_norm(corrected_upstroke_point(end):end),corrected_upstroke_point(end),rapid_depolar_end_corrected(end), 1,frameperiod, plot_index,signal_swing(end));
                                    trace_upstroke_fitted(corrected_upstroke_point(end):end) = fitted_Vm; 
                                    linear_coefficient(:,end) = fit_coef; 
                                    trace_upstroke_fitted(1:corrected_upstroke_point(1)) =row_filtered_sig1(1:corrected_upstroke_point(1)); 
                                    wrong_corrected_index = find(further_corrected_upstroke- rapid_depolar_end_corrected>0); 
                                    further_corrected_upstroke(wrong_corrected_index) = corrected_upstroke_point(wrong_corrected_index); 
                                    if plot_index==1
                                        pause
                                        close all
                                    end

                                    %repolarisation fitting
                                    trace_curve_fitted = trace_upstroke_fitted;

                                    poly_coefficient = zeros(8,length(corrected_upstroke_point)); 
                                    baseline_intersection_time_point = zeros(1,length(further_corrected_upstroke)); % this uses further corrected upstroke point as an reference point
                                    for q = 1:length(upstroke_time_points)-1
                                       [fitted_Vm,fit_coef,baseline_intersection] = fitt_CaT_curve(row_filtered_norm(further_corrected_upstroke(q):further_corrected_upstroke(q+1)),trace1_norm(further_corrected_upstroke(q):further_corrected_upstroke(q+1)),further_corrected_upstroke(q),rapid_depolar_end_corrected(q), further_corrected_upstroke(q+1),1,frameperiod, plot_index,signal_swing(q));
                                       trace_curve_fitted(further_corrected_upstroke(q):further_corrected_upstroke(q+1)) = fitted_Vm; 
                                       poly_coefficient(:,q) = fit_coef;                                                             

                                    end 

                                    [fitted_Vm,fit_coef,baseline_intersection] = fitt_CaT_curve(row_filtered_norm(further_corrected_upstroke(end):end),trace1_norm(further_corrected_upstroke(end):end),further_corrected_upstroke(end),rapid_depolar_end_corrected(end), length(row_filtered_sig1),1,frameperiod, plot_index,signal_swing(end));
                                    trace_curve_fitted(further_corrected_upstroke(end):end) = fitted_Vm; 
                                    poly_coefficient(:,end) = fit_coef; 

                                    if plot_index==1
                                        pause
                                        close all;
                                        figure;
                                        hold on; 
                                        plot(1:frameperiod:frameperiod*length(trace_curve_fitted),trace1_norm(1:end),'b') 
                                        plot(further_corrected_upstroke*frameperiod,trace1_norm(further_corrected_upstroke),'c*') 
                                        %plot(1:frameperiod:frameperiod*(length(trace_curve_fitted)-further_corrected_upstroke(1)+1),trace_curve_fitted(further_corrected_upstroke(1):end),'r');
                                        xlabel('time (ms)'); 
                                        ylabel('F/F0'); 
                                        legend('raw trace','further_corrected_upstroke')
                                        title('Check further corrected upstroke via upstroke fitting in case of extreme CaT alternans (good signal)');
                                        hold off;
                                        pause
                                        close all
                                        pause

                                        figure; 
                                        hold on; 
                                        plot(1:frameperiod:frameperiod*(length(trace_curve_fitted)-further_corrected_upstroke(1)+1),trace1_norm(further_corrected_upstroke(1):end),'b') 
                                        plot(1:frameperiod:frameperiod*(length(trace_curve_fitted)-further_corrected_upstroke(1)+1),trace_curve_fitted(further_corrected_upstroke(1):end),'r');
                                        xlabel('time (ms)'); 
                                        ylabel('F/F0'); 
                                        legend('raw trace','fitted trace');
                                        title('Check fitted trace in case of extreme CaT alternans (good signal)');
                                        hold off;
                                        pause
                                        close all
                                        pause

                                    end 

                                end




                                baseline_norm=ones(1,length(baseline));
                                APD50 = zeros(1,length(further_corrected_upstroke)); 
                                APD80 = zeros(1,length(further_corrected_upstroke)); 
                                AP_peak_value = zeros(1,length(further_corrected_upstroke));
                                for r2 = 1:length(further_corrected_upstroke)-1

                                    [APD,AP_peak_value(r2)] = Find_all_CaTD(trace_curve_fitted(further_corrected_upstroke(r2):further_corrected_upstroke(r2+1)),trace1_norm(further_corrected_upstroke(r2):further_corrected_upstroke(r2+1)),further_corrected_upstroke(r2), rapid_depolar_end_corrected(r2),baseline_norm(r2),frameperiod, poly_coefficient(:,r2),pacing_interval,plot_index); %the baseline here is -1
                                    if isempty(APD.APD80)~=1
                                        APD80(r2)= APD.APD80; 
                                    end
                                    if isempty(APD.APD50)~=1
                                       APD50(r2)= APD.APD50; 
                                    end      

                                end 
                                [APD,AP_peak_value(end)] = Find_all_CaTD(trace_curve_fitted(further_corrected_upstroke(end):end),trace1_norm(further_corrected_upstroke(end):end),further_corrected_upstroke(end),rapid_depolar_end_corrected(end),baseline_norm(end),frameperiod, poly_coefficient(:,end),pacing_interval,plot_index); %the baseline here is -1
                                if isempty(APD.APD80)~=1
                                    if APD.APD80<pacing_interval
                                       APD80(end)= APD.APD80;
                                    else 
                                       APD80(end)=0;
                                    end 
                                end
                                if isempty(APD.APD50)~=1
                                    if APD.APD50<pacing_interval
                                        APD50(end)= APD.APD50; 
                                    else 
                                        APD50(end)=0;
                                    end
                                end      

                                 if length(find(eval_temp>0))==length(eval_temp)%big realease on the odd beat 
                                    APD80_odd=APD80(1:2:end);
                                    index_long_APD80=find(APD80_odd>pacing_interval);
                                    APD80_odd(index_long_APD80)=pacing_interval; 
                                    APD80(1:2:end)=APD80_odd; 
                                    if sum(abs(APD80(1:2:end)-APD80_odd))~=0
                                       APD80_saturation_index=1; 
                                    end 
                                    APD50_odd=APD50(1:2:end);
                                    index_long_APD50=find(APD50_odd>pacing_interval);
                                    APD50_odd(find(APD50_odd>pacing_interval))=pacing_interval; 
                                    APD50(1:2:end)=APD50_odd; 
                                    if sum(abs(APD50(1:2:end)-APD50_odd))~=0
                                       APD50_saturation_index=1; 
                                    end
                                 elseif length(find(eval_temp<0))==length(eval_temp)%big realease on the even beat 
                                    APD80_even=APD80(2:2:end);
                                    index_long_APD80=find(APD80_even>pacing_interval);
                                    APD80_even(index_long_APD80)=pacing_interval; 
                                    APD80(2:2:end)=APD80_even; 
                                    if sum(abs(APD80(2:2:end)-APD80_even))~=0
                                       APD80_saturation_index=1; 
                                    end 

                                    APD50_even=APD50(2:2:end);
                                    index_long_APD50=find(APD50_even>pacing_interval);
                                    APD50_even(find(APD50_even>pacing_interval))=pacing_interval; 
                                    APD50(2:2:end)=APD50_even;
                                    if sum(abs(APD50(2:2:end)-APD50_even))~=0
                                       APD50_saturation_index=1; 
                                    end 

                                 end

                                area_under_AP_vect = []; 
                                area_ratio_vect = []; 
                                AP_skewness = []; 

                                APD80_median = median(APD80); 
                                APD50_median = median(APD50); 
                                APD50_std = std(APD50); 
                                APD80_std =std(APD80); 
                                AP_skewness_median = [];
                                AP_skewness_std = []; 
                                area_under_AP_mean = [];
                                area_under_AP_std = []; 
                                area_ratio_mean = []; 
                                area_ratio_std = [];                     

                                number_AP =length(APD80); 
                                number_upstroke = length(corrected_upstroke_point); 
                                %since the last CaT can be incomplete and result
                                %in a problematic CaT duration, therefore don't use
                                %it for alternans analysis 

                                [alternans_index,alternans_score,peak_alternans_index,peak_alternans_score] = test_alternans_CaT_special(AP_peak_value(1:end-1), further_corrected_upstroke(1:end-1), APD80(1:end-1), APD50(1:end-1), known_pacing)


                            else
                                for i = 1:length(upstroke_time_points)-1
                                    
                                    number_to_bin= (max(filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1)))-...
                                    min(filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1))))./10;
                            
                                    if number_to_bin<3
                                       number_to_bin= (max(filtered_sig1(upstroke_time_points(end):end))-...
                                       min(filtered_sig1(upstroke_time_points(end):end))).*10;
                                    end 

                                    if number_to_bin<3
                                        number_to_bin=10; 
                                    end 
                         
                                    [accum_freq, bin_centre] = hist(filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1)),number_to_bin); 
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
                                    normal_fitting_vector = filtered_sig1(upstroke_time_points(i):upstroke_time_points(i+1));
                                    normal_fitting_vector(normal_fitting_vector>cut_off)=[];
                                    if length(normal_fitting_vector)>3
                                        pd= fitdist(normal_fitting_vector,'normal'); 
                                        baseline(i) = mean(pd);
                                    else 
                                        baseline(i) = filtered_sig1(upstroke_time_points(i)); 
                                    end
                                    if abs(baseline(i)-filtered_sig1(upstroke_time_points(i)))>0.1*abs(filtered_sig1(rapid_depolar_end(i))-filtered_sig1(upstroke_time_points(i)))
                                       baseline(i) = (baseline(i)+filtered_sig1(upstroke_time_points(i)))/2;
                                    end 
                                end 
                                
                                number_to_bin= (max(filtered_sig1(upstroke_time_points(end):end))-...
                                min(filtered_sig1(upstroke_time_points(end):end)))./10;
                    
                                if number_to_bin<3
                                   number_to_bin= (max(filtered_sig1(upstroke_time_points(end):end))-...
                                    min(filtered_sig1(upstroke_time_points(end):end))).*10;
                                end 

                                if number_to_bin<3
                                    number_to_bin=10; 
                                end 
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

                                if abs(baseline(end)-filtered_sig1(upstroke_time_points(end)))>0.1*abs(filtered_sig1(rapid_depolar_end(end))-filtered_sig1(upstroke_time_points(end)))
                                       baseline(end) = trace1(upstroke_time_points(end));
                                end 
                                %check with a plot 
                                if plot_index==1
                                    figure ; 
                                    time = 1:frameperiod:length(filtered_sig1)*frameperiod;
                                    hold on; 
                                    plot(time, trace1,'b'); 
                                    %plot(time, filtered_sig1,'r')

                                    for j = 1:length(upstroke_time_points)-1
                                        time_interval_base = frameperiod*upstroke_time_points(j):frameperiod:frameperiod*upstroke_time_points(j+1); 
                                        plot(time_interval_base, ones(1,length(time_interval_base)).*baseline(j),'k'); 

                                        % check with the upstroke plot 
                                        % plot(upstroke_time_points(j)*frameperiod,trace1(upstroke_time_points(j)),'b*'); 
                                    end 

                                    time_interval_base = frameperiod*upstroke_time_points(end):frameperiod:frameperiod*length(trace1); 
                                    plot(time_interval_base, ones(1,length(time_interval_base)).*baseline(end),'k'); 
                                    %plot(upstroke_time_points(end)*frameperiod,trace1(upstroke_time_points(end)),'b*'); 
                                    xlabel('time (ms)'); 
                                    ylabel('vol sig'); 
                                    legend('row trace', 'filtered trace','baseline');
                                    pause
                                end 
                                %% Estimate signal swing 

                                % find out the 95 percentile value as the maximum signal value  
                                signal_max = zeros(1,length(upstroke_time_points)); 
                                for l = 1:length(upstroke_time_points)-1
                                    AP_chunk = trace1(upstroke_time_points(l): upstroke_time_points(l+1));
                                    quartile_value = quantile(AP_chunk, 0.975);
                                    signal_max(l) = quartile_value;
                                end 

                                AP_chunk = trace1(upstroke_time_points(end): end);
                                quartile_value = quantile(AP_chunk, 0.975);
                                signal_max(end) = quartile_value; 

                                % check the max signal with a plot 
                                if plot_index==1
                                    figure; 
                                    hold on; 
                                    plot(time, trace1,'b'); 
                                    %plot(time, filtered_sig1,'r')
                                    for m = 1:length(upstroke_time_points)-1
                                        time_interval_base = frameperiod*upstroke_time_points(m):frameperiod:frameperiod*upstroke_time_points(m+1); 
                                        plot(time_interval_base, ones(1,length(time_interval_base)).*signal_max(m),'k'); 
                                         plot(time_interval_base, ones(1,length(time_interval_base)).*baseline(m),'k'); 
                                    end 
                                    time_interval_base = frameperiod*upstroke_time_points(end):frameperiod:frameperiod*length(trace1); 
                                    plot(time_interval_base, ones(1,length(time_interval_base)).*signal_max(end),'k'); 
                                    xlabel('time (ms)'); 
                                    ylabel('vol sig'); 
                                    legend('row trace','filtered trace','maximum sig'); 
                                    hold off; 
                                    pause
                                end 

                                signal_swing = abs((signal_max'-baseline)./baseline);
                                average_swing = mean(signal_swing); 
                             end 
                        end 
                    end


                    %% Correction on upstroke time points using baseline

                    if extreme_alternans_index==0

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
                                if intersection_time_point <upstroke_time_points(n) 
                                   corrected_upstroke_point(n)= upstroke_time_points(n);
                                else 
                                   corrected_upstroke_point(n) = max(intersection_time_point,upstroke_time_points(n)); 
                                end 
                            end 
                           % plot(ceil(corrected_upstroke_point(n))*frameperiod, trace1(ceil(corrected_upstroke_point(n))),'ro'); % check with a plot need to be used test baseline plot with upstroke plot uncommented 
                        end 
                        corrected_upstroke_point = ceil(corrected_upstroke_point);

                        if plot_index==1
                            figure; 
                            hold on; 
                            plot(1:frameperiod:frameperiod*length(trace1),trace1,'k'); 
                            plot(1:frameperiod:frameperiod*length(trace1), filtered_sig1,'r')
                            plot(upstroke_time_points.*frameperiod,trace1(upstroke_time_points),'k*'); 
                            plot(ceil(corrected_upstroke_point).*frameperiod, trace1(ceil(corrected_upstroke_point)),'ko');
                            legend('raw trace','filtered_trace','upstroke time', 'corrected upstroke time');
                            title('Check upstroke time point correction')
                            hold off; 
                            pause
                        end


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
                        if plot_index==1
                            figure; 
                            hold on; 
                            plot(1:frameperiod:frameperiod*length(trace1), trace1,'k');
                            plot(1:frameperiod:frameperiod*length(trace1),filtered_sig1,'r');
                            plot(rapid_depolar_end.*frameperiod,trace1(rapid_depolar_end),'b*'); 
                            plot(rapid_depolar_end_corrected.*frameperiod, trace1(rapid_depolar_end_corrected),'bo');
                            legend('raw trace','filtered_trace','end of upstroke', 'corrected end of upstroke');
                            title('Check end of upstroke correction')
                            hold off;                   
                            pause  
                            close all
                        end 

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
                                [fitted_Vm,linear_coef,further_corrected_upstroke(q)] = fitt_AP_upstroke(row_filtered_sig1(corrected_upstroke_point(q):corrected_upstroke_point(q+1)),trace1(corrected_upstroke_point(q):corrected_upstroke_point(q+1)),corrected_upstroke_point(q),rapid_depolar_end_corrected(q), baseline(q),frameperiod, plot_index,signal_swing(q)); 
                                trace_upstroke_fitted(corrected_upstroke_point(q):corrected_upstroke_point(q+1)) = fitted_Vm; 
                                linear_coefficient(:,q) = linear_coef;  

                            end 
                            [fitted_Vm,fit_coef,further_corrected_upstroke(end)] =fitt_AP_upstroke(row_filtered_sig1(corrected_upstroke_point(end):end),trace1(corrected_upstroke_point(end):end),corrected_upstroke_point(end),rapid_depolar_end_corrected(end), baseline(end),frameperiod, plot_index,signal_swing(end));
                            trace_upstroke_fitted(corrected_upstroke_point(end):end) = fitted_Vm; 
                            linear_coefficient(:,end) = fit_coef; 

                            trace_upstroke_fitted(1:corrected_upstroke_point(1)) =row_filtered_sig1(1:corrected_upstroke_point(1)); 

                            wrong_corrected_index = find(further_corrected_upstroke- rapid_depolar_end_corrected>0); 
                            further_corrected_upstroke(wrong_corrected_index) = corrected_upstroke_point(wrong_corrected_index); 

                            %repolarisation fitting 
                            trace_curve_fitted = trace_upstroke_fitted;

                            poly_coefficient = zeros(8,length(corrected_upstroke_point)); 
                            baseline_intersection_time_point = zeros(1,length(further_corrected_upstroke)); % this uses further corrected upstroke point as an reference point
                            for q = 1:length(upstroke_time_points)-1
                               [fitted_Vm,fit_coef,baseline_intersection] = fitt_CaT_curve(row_filtered_sig1(further_corrected_upstroke(q):further_corrected_upstroke(q+1)),trace1(further_corrected_upstroke(q):further_corrected_upstroke(q+1)),further_corrected_upstroke(q),rapid_depolar_end_corrected(q), further_corrected_upstroke(q+1),baseline(q),frameperiod, plot_index,signal_swing(q));
                               if isempty(baseline_intersection) ~=1
                                  trace_curve_fitted(further_corrected_upstroke(q):further_corrected_upstroke(q+1)) = fitted_Vm; 
                                  poly_coefficient(:,q) = fit_coef; 
                                  baseline_intersection_time_point(q)= baseline_intersection;
                               else 
                                   poly_coefficient(:,q) = NaN; 
                                   %disp('no fitting can be done'); 
                               end 

                            end 


                            [fitted_Vm,fit_coef,baseline_intersection] = fitt_CaT_curve(row_filtered_sig1(further_corrected_upstroke(end):end),trace1(further_corrected_upstroke(end):end),further_corrected_upstroke(end),rapid_depolar_end_corrected(end), length(row_filtered_sig1),baseline(end),frameperiod, plot_index,signal_swing(end));
                            if isempty(baseline_intersection)~=1
                                baseline_intersection_time_point(end) = baseline_intersection; 
                                trace_curve_fitted(further_corrected_upstroke(end):end) = fitted_Vm; 
                                poly_coefficient(:,end) = fit_coef; 
                            else 
                                poly_coefficient(:,end)=NaN; 
                                %disp('no fitting can be done'); 
                            end

                            if plot_index==1
                                figure; 
                                hold on; 

                                plot(1:frameperiod:frameperiod*length(trace_curve_fitted),trace1(1:end),'b') 
                                plot(further_corrected_upstroke*frameperiod,trace1(further_corrected_upstroke),'c*') 
                                %plot(1:frameperiod:frameperiod*(length(trace_curve_fitted)-further_corrected_upstroke(1)+1),trace_curve_fitted(further_corrected_upstroke(1):end),'r');
                                xlabel('time (ms)'); 
                                ylabel('F/F0'); 
                                legend('raw trace','filtered trace');
                                hold off;
                                pause
                                close all
                                pause

                                figure; 
                                hold on; 

                                plot(1:frameperiod:frameperiod*(length(trace_curve_fitted)-further_corrected_upstroke(1)+1),trace1(further_corrected_upstroke(1):end),'b') 
                                plot(1:frameperiod:frameperiod*(length(trace_curve_fitted)-further_corrected_upstroke(1)+1),trace_curve_fitted(further_corrected_upstroke(1):end),'r');
                                xlabel('time (ms)'); 
                                ylabel('F/F0'); 
                                legend('raw trace','fitted trace');
                                hold off;
                                pause
                                close all
                                pause
                            end


                            %% check the length of the baseline to see whether we have a decent baseline 
                            % check this only if we have more upstroke than 1        
                            % for the first beat to the second last beat check wehther we get enough baseline for each beat 
                            if length(further_corrected_upstroke)>1

                                upstroke_interval = further_corrected_upstroke(2:end)-further_corrected_upstroke(1:end-1);
                                baseline_length_interval = upstroke_interval-baseline_intersection_time_point(1:end-1);

                                non_significant_baseline_index = find(baseline_length_interval<=0);
                                if average_swing>0.015


                                   baseline_intersection_time_point(non_significant_baseline_index)=further_corrected_upstroke(non_significant_baseline_index+1)-further_corrected_upstroke(non_significant_baseline_index)-1;

                                   non_significant_baseline_index=[];
                                end 
                                
                                motion_corrected_trace = trace_curve_fitted; 
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
                                    if plot_index==1
                                        figure; 
                                        hold on; 
                                        plot(1:frameperiod:frameperiod*length(trace_curve_fitted),trace_curve_fitted,'b');
                                        for i = 1:length(further_corrected_upstroke)
                                            if binary_vect_motion_correction(i)==0 
                                               plot(frameperiod.*[motion_correction_start_end_points(i,1),motion_correction_start_end_points(i,2)],[motion_correction_before_after_baseline(i,1),motion_correction_before_after_baseline(i,2)],'r');
                                            end 
                                        end
                                        title('motion correction');
                                        pause 
                                    end 

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

                                        special_remove_index = find(non_significant_baseline_index_full+1-length(further_corrected_upstroke)>=0) ;
                                        non_special_remove_index =  find(non_significant_baseline_index_full+1-length(further_corrected_upstroke)<0) ;
        %                                normalize the other APs in the trace 
                                        motion_corrected_trace (1:further_corrected_upstroke(2)) = trace_curve_fitted(1:further_corrected_upstroke(2))./baseline(1); 
                                        motion_corrected_trace (further_corrected_upstroke(end)+1:end) = trace_curve_fitted(further_corrected_upstroke(end)+1:end)./baseline(end);
           %                             remove the part which was found no significant baseline 
                                        if isempty(non_special_remove_index)~=1  
                                            for u = 1:length(non_special_remove_index)
                                               motion_corrected_trace(further_corrected_upstroke(non_significant_baseline_index_full(non_special_remove_index(u))):further_corrected_upstroke(non_significant_baseline_index_full(non_special_remove_index(u))+1)-1)=1; 

                                            end 
                                        end
                                        
                                        
                                                                     
                                
                        
                            
                                    
                                        if isempty(special_remove_index)~=1
                                           motion_corrected_trace(further_corrected_upstroke(non_significant_baseline_index_full(special_remove_index(1))):end)=1;
                                          
                                        end 
                                else 
                                    motion_corrected_trace(1:further_corrected_upstroke(2)-1) = trace_curve_fitted(1:further_corrected_upstroke(2)-1)./baseline(1); 
                                    motion_corrected_trace(further_corrected_upstroke(2):end)=1; 
                                end 
                            else
                                motion_corrected_trace = trace_curve_fitted./baseline(1); 
                                non_significant_baseline_index = [];
                            end
                            if plot_index==1
                               figure; 
                               plot(motion_corrected_trace)
                               title('motioned_corrected_trace');
                               pause
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
                            if plot_index==1 
                                figure; 
                                 hold on; 
                                 plot(trace1_norm,'b') 
                                 plot(row_filtered_norm,'r'); 
                    %                 %plot(trace_fitted_norm,'k'); 
                                 hold off
                                 pause
                            end 
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
                                poly_coefficient = zeros(8,length(corrected_upstroke_point)); 
                                baseline_intersection_time_point = zeros(1,length(further_corrected_upstroke)); % this uses further corrected upstroke point as an reference point
                                                           
                                
                           
                                for q = 1:length(further_corrected_upstroke)-1
                                    
                                   [fitted_Vm,fit_coef,baseline_intersection] = fitt_CaT_curve(row_filtered_norm(further_corrected_upstroke(q):further_corrected_upstroke(q+1)),motion_corrected_trace(further_corrected_upstroke(q):further_corrected_upstroke(q+1)),further_corrected_upstroke(q),rapid_depolar_end_corrected(q),further_corrected_upstroke(q+1), baseline_norm(q),frameperiod, 0,signal_swing(q));
                                   if isempty(baseline_intersection) ~=1
                                      trace_curve_refitted(further_corrected_upstroke(q):further_corrected_upstroke(q+1)) = fitted_Vm; 
                                      poly_coefficient(:,q) = fit_coef; 
                                      baseline_intersection_time_point(q)= baseline_intersection;
                                   else 
                                       poly_coefficient(:,q) = NaN; 
                                      % disp('no fitting can be done'); 
                                   end 

                                end 


                                [fitted_Vm,fit_coef,baseline_intersection] = fitt_CaT_curve(row_filtered_norm(further_corrected_upstroke(end):end),motion_corrected_trace(further_corrected_upstroke(end):end),further_corrected_upstroke(end),rapid_depolar_end_corrected(end),length(motion_corrected_trace), baseline_norm(end),frameperiod, plot_index,signal_swing(end));
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
                                    [APD,AP_peak_value(r2)] = Find_all_CaTD(trace_curve_refitted(further_corrected_upstroke(r2):further_corrected_upstroke(r2+1)),trace1_norm(further_corrected_upstroke(r2):further_corrected_upstroke(r2+1)),further_corrected_upstroke(r2), rapid_depolar_end_corrected(r2),baseline_norm(r2),frameperiod, poly_coefficient(:,r2),pacing_interval,plot_index); %the baseline here is -1
                                    if isempty(APD.APD80)~=1
                                        APD80(r2)= APD.APD80; 
                                    end
                                    if isempty(APD.APD50)~=1
                                       APD50(r2)= APD.APD50; 
                                    end      

                                end 
                                [APD,AP_peak_value(end)] = Find_all_CaTD(trace_curve_refitted(further_corrected_upstroke(end):end),trace1_norm(further_corrected_upstroke(end):end),further_corrected_upstroke(end),rapid_depolar_end_corrected(end),baseline_norm(end),frameperiod, poly_coefficient(:,end),pacing_interval,plot_index); %the baseline here is -1
                                if isempty(APD.APD80)~=1
                                    if APD.APD80<pacing_interval
                                        APD80(end)= APD.APD80; 
                                    else 
                                        APD80(end)=0; 
                                    end 

                                end
                                if isempty(APD.APD50)~=1
                                    if APD.APD50<pacing_interval 
                                       APD50(end)= APD.APD50; 
                                    else 
                                        APD50(end)=0; 
                                    end 
                                end      

                                %% calculate number of AP found 
                                number_AP = length(find(APD50>0)); 
                                number_upstroke = length(further_corrected_upstroke); 

                                %% calculate area under the curve and the ratio 
                                area_under_AP_vect = zeros(1,length(further_corrected_upstroke)); 
                                area_ratio_vect =  zeros(1,length(further_corrected_upstroke)); 
                                for s = 1:length(further_corrected_upstroke)-1
                                    [area_under_AP,area_ratio]= find_area_APD(trace_curve_refitted(further_corrected_upstroke(s):further_corrected_upstroke(s+1)),baseline_norm(s),frameperiod,baseline_intersection_time_point(s),APD80(s),AP_peak_value(s),plot_index );
                                    area_under_AP_vect(s) = area_under_AP; 
                                    area_ratio_vect(s) = area_ratio;

                                end 
                                [area_under_AP,area_ratio]= find_area_APD(trace_curve_refitted(further_corrected_upstroke(end):end),baseline_norm(end),frameperiod,baseline_intersection_time_point(end),APD80(end), AP_peak_value(end),plot_index );
                                area_under_AP_vect(end) = area_under_AP; 
                                area_ratio_vect(end) = area_ratio;
                                if plot_index==1
                                    pause
                                    close all
                                end 

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
                                
                                if length(APD80)>1
                                    [alternans_index,alternans_score,peak_alternans_index,peak_alternans_score] = test_alternans_CaT(AP_peak_value(1:end-1), further_corrected_upstroke(1:end-1), APD80(1:end-1), area_under_AP_vect(1:end-1), known_pacing);
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

            


