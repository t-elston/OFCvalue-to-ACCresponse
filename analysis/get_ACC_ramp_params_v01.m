function [xx] = get_ACC_ramp_params_v01(s_contra_norm_full,s_ipsi_norm_full, choice_ts, sel_monkey_ix, thresh)
xx=[];

[~, w_start] = min(abs(choice_ts - -1000));
[~, w_end] = min(abs(choice_ts - 1000));
win_times = choice_ts(w_start:w_end);


s_ipsi_norm = s_ipsi_norm_full(:,w_start:w_end);
s_contra_norm = s_contra_norm_full(:,w_start:w_end);

[n_units, n_times] = size(s_ipsi_norm);

contra_over_thresh_ix = double(s_contra_norm >=thresh);
ipsi_over_thresh_ix = double(s_ipsi_norm >=thresh);

for u = 1:n_units
    
    contra_thresh_crossing = min(find(contra_over_thresh_ix(u,:)));
    ipsi_thresh_crossing = min(find(ipsi_over_thresh_ix(u,:)));
    
    % get max values and times
    [contra_max_val,contra_max_ix] = max(s_contra_norm(u,:));
    [ipsi_max_val,ipsi_max_ix] = max(s_ipsi_norm(u,:));


    try
    contra_thresh_time(u,1) = win_times(contra_thresh_crossing);
    ipsi_thresh_time(u,1) = win_times(ipsi_thresh_crossing);
    
    contra_max_time(u,1) = win_times(contra_max_ix(1));
    ipsi_max_time(u,1) = win_times(ipsi_max_ix(1));
    
    contra_max_z(u,1) = contra_max_val;
    ipsi_max_z(u,1) = ipsi_max_val;
    
    monkey_ix(u) = sel_monkey_ix(u);
    
    catch  
        contra_thresh_time(u,1) = NaN;
        ipsi_thresh_time(u,1) = NaN;
        contra_max_time(u,1) = NaN;
        ipsi_max_time(u,1) = NaN;
        contra_max_z(u,1) = NaN;
        ipsi_max_z(u,1) = NaN;
    end
        
end % of cylcing over units

% now sort the cells
[~,contra_sort_ix] = sort(contra_thresh_time);
[~,ipsi_sort_ix] = sort(ipsi_thresh_time);

contraM_contraS = s_contra_norm(contra_sort_ix,:);




plot(contra_thresh_time(monkey_ix==0), ipsi_thresh_time(monkey_ix==0),'.')





end % of function