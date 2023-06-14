function [down_pps] = boxcar_downsample_PPs(pp, decoder_ts, unit_ts)
down_pps=[];

decoder_ts = decoder_ts+1;

[n_trials, n_original_times] = size(pp);
n_unit_times = numel(unit_ts);

for ut = 1:n_unit_times
    
    u_time = unit_ts(ut);
    [~,u_win_start] = min(abs((u_time-50) - decoder_ts));
    [~,u_win_end] = min(abs((u_time+50) - decoder_ts));

    
    down_pps(:,ut) = nanmean(pp(:,u_win_start:u_win_end),2);
    
    
end % of cylcing over unit times



end % of function