function [match_p, nonmatch_p, LR_betas] = regress_FRs_against_state_val(ch_FRs, unch_FRs, ch_details, unch_details, xt)


% format of ch_details and unch_details is:
% column 1 = trial state occured in
% column 2 = the state within the trial (first, second, third, etc)
% column 3 = the value of the state
% column 4 = the screen side of the state (-1 = left, 1 = right)
% column 5 = the time of state onset
% column 6 = the value of the alternative option from that of the state
% column 7 = whether it was a free or forced choice


all_FRs = [ch_FRs ; unch_FRs];
all_details = [ch_details ; unch_details];

% only use free choices with different values
states2use = all_details(:,7) == 2 & (all_details(:,3) ~= all_details(:,6));

% subselect the relevant trials
FRs = all_FRs(states2use,:);
match_vals    = all_details(states2use,3);
nonmatch_vals = all_details(states2use,6);

% construct the regressors
match_regressors = [ones(size( match_vals)), match_vals];
nonmatch_regressors = [ones(size(nonmatch_vals)), nonmatch_vals];

% loop over timesteps
n_times = length(xt);
for t = 1:n_times
    
    [~,~,~,~, match_stats] = regress(FRs(:,t), match_regressors);
    [~,~,~,~, nonmatch_stats] = regress(FRs(:,t), nonmatch_regressors);
    
    match_p(1,t) = match_stats(3);
    nonmatch_p(1,t) = nonmatch_stats(3);
    
end % of looping over times

% now compute the betas for when the states were to the left and right
ch_win_FRs = mean(ch_FRs(:, 21:26),2);

left_ix = ch_details(:,4) == -1 & ch_details(:,7) == 2;
right_ix = ch_details(:,4) == 1 & ch_details(:,7) == 2;

left_vals = ch_details(left_ix,3);
right_vals = ch_details(right_ix,3);

left_b = regress(ch_win_FRs(left_ix), [ones(size(left_vals)), left_vals]);
right_b = regress(ch_win_FRs(right_ix), [ones(size(right_vals)), right_vals]);

LR_betas = [left_b(2), right_b(2)];

end % of function