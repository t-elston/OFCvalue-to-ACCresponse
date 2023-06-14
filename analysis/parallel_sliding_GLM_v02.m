function [u_sig, u_beta] = parallel_sliding_GLM_v02(reg_tbl, CI_choice_dir, z_units, choice_ts)

[n_trials, n_times] = size(z_units);

% get indices of when the first 500ms prior to choice were
[~, w_start] = min(abs(choice_ts - -500));
[~, w_end] = min(abs(choice_ts - 0));
% [~, w_start] = min(abs(choice_ts - 0));
% [~, w_end] = min(abs(choice_ts - 500));

win_fr = nanmean(z_units(:,w_start:w_end),2);


reg_tbl.CIdir = CI_choice_dir;
reg_tbl.win_fr = win_fr;

sig_times = NaN(1, n_times, 3);
t_betas = NaN(1, n_times, 3);


u_sig = NaN(1,3);
u_beta = NaN(1,3);




for t = 1 : n_times
    
    reg_tbl.firing_rate = z_units(:,t);
    
    % now fit the GLM
    u_t_mdl = fitlm(reg_tbl,'firing_rate ~ max_val + trial_type + CIdir + t_num');
    
    sig_times(1,t,1) = u_t_mdl.Coefficients.pValue(5) < .01; % choice direction
    
    t_betas(1,t,1) = u_t_mdl.Coefficients.tStat(5); % choice direction
    
    sig_times(1,t,2) = u_t_mdl.Coefficients.pValue(3) < .01; % maxval
    
    t_betas(1,t,2) = u_t_mdl.Coefficients.tStat(3); % maxval
    
    sig_times(1,t,3) = u_t_mdl.Coefficients.pValue(2) < .01; % trial type
    
    t_betas(1,t,3) = u_t_mdl.Coefficients.tStat(1); % trial type

end % of looping over times


% now ask whether this unit was significant with significance defined as at least 100ms of selectivity in the 500ms preceding
% choice

[~, dir_len, ~] = ZeroOnesCount(sig_times(1,w_start:w_end,1)); % choice direction
[~, val_len, ~] = ZeroOnesCount(sig_times(1,w_start:w_end,2)); % max val
[~, trialtype_len, ~] = ZeroOnesCount(sig_times(1,w_start:w_end,3)); % trial type


u_sig(1,1) = any(dir_len >=4);
u_sig(1,2) = any(val_len >=4);
u_sig(1,3) = any(trialtype_len >=4);

u_beta(1,1) = nanmean(t_betas(1,w_start:w_end,1));
u_beta(1,2) = nanmean(t_betas(1,w_start:w_end,2));
u_beta(1,3) = nanmean(t_betas(1,w_start:w_end,3));


end % of function