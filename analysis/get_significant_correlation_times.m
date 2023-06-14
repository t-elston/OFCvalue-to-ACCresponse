function [sig_times] = get_significant_correlation_times(indata, len_thresh, t_mids)


sig_times = zeros(size(t_mids));

time_step = mean(diff(t_mids));
% find moments that were significantly greater than 0
candidate_sig_times = ttest(indata,0,"Tail","right",'Alpha',0.01); 
candidate_sig_times(isnan(candidate_sig_times))=0;

% find stretches of significance
[start, len, k1] = ZeroOnesCount(candidate_sig_times);

selectivity_lens = len*time_step;

valid_moments = find(selectivity_lens > len_thresh);

for m = 1:numel(valid_moments)
    
    m_start = start(valid_moments(m));
    m_len = len(valid_moments(m));
    
    sig_times(m_start:m_start+m_len) = 1;
    
end

sig_times = find(sig_times);

end % of function