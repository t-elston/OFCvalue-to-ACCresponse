function [sig_times, t_betas, ipsi_frs, contra_frs, u_sig, u_beta, CPD] = parallel_sliding_GLM_v01(reg_tbl, CI_choice_dir, z_units, choice_ts)

[n_trials, n_times] = size(z_units);

% get indices of when the first 500ms prior to choice were
[~, w_start] = min(abs(choice_ts - -500));
[~, w_end] = min(abs(choice_ts - 0));

win_fr = nanmean(z_units(:,w_start:w_end),2);


reg_tbl.CIdir = CI_choice_dir;
reg_tbl.win_fr = win_fr;

sig_times = NaN(1, n_times, 2);
t_betas = NaN(1, n_times, 2);





u_win_mdl = fitlm(reg_tbl,'win_fr ~ max_val + trial_type + CIdir + t_num');
% u_sig(1,1) = u_win_mdl.Coefficients.pValue(5) < .01; % choice direction
u_beta(1,1) = u_win_mdl.Coefficients.Estimate(5); % choice direction

u_tbl = anova(u_win_mdl);
u_beta(1,2) = GetPartialEtaSquaredFromFAndDFs_v01(u_tbl.F(2),u_tbl.DF(2),u_tbl.DF(5)); % value


u_beta(1,3) = u_win_mdl.Coefficients.Estimate(3); % max val


for t = 1 : n_times
    
    reg_tbl.firing_rate = z_units(:,t);
    
    % now fit the GLM
    u_t_mdl = fitlm(reg_tbl,'firing_rate ~ max_val + trial_type + CIdir + t_num');
    
    sig_times(1,t,1) = u_t_mdl.Coefficients.pValue(5) < .01; % choice direction
    
    t_betas(1,t,1) = u_t_mdl.Coefficients.tStat(5); % choice direction
    
    sig_times(1,t,2) = u_t_mdl.Coefficients.pValue(3) < .01; % maxval
    
    t_betas(1,t,2) = u_t_mdl.Coefficients.tStat(3); % maxval
    
    % calculate the CPD
    tbl = anova(u_t_mdl);
    CPD(1,t,1) = GetPartialEtaSquaredFromFAndDFs_v01(tbl.F(4),tbl.DF(4),tbl.DF(5)); % choice dir
    CPD(1,t,2) = GetPartialEtaSquaredFromFAndDFs_v01(tbl.F(2),tbl.DF(2),tbl.DF(5)); % max val

    
end % of looping over times


% now ask whether this unit was significant with significance defined as at least 100ms of selectivity in the 500ms preceding
% choice

[~, dir_len, ~] = ZeroOnesCount(sig_times(1,w_start:w_end,1));
[~, val_len, ~] = ZeroOnesCount(sig_times(1,w_start:w_end,2));


% u_sig(1,1) = sum(sig_times(1,w_start:w_end,1)) >= 4;
% u_sig(1,2) = sum(sig_times(1,w_start:w_end,2)) >= 4;

u_sig(1,1) = any(dir_len >=4);
u_sig(1,2) = any(val_len >=4);

ipsi_frs(1,:) = nanmean(z_units(reg_tbl.CIdir==1,:));
contra_frs(1,:) = nanmean(z_units(reg_tbl.CIdir==-1,:));


end % of function