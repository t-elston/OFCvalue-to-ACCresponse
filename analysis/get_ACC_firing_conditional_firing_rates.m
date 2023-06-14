function [hv_pref_mean, hv_pref_sem, lv_pref_mean, lv_pref_sem,...
    hv_nonpref_mean, hv_nonpref_sem, lv_nonpref_mean, lv_nonpref_sem] =...
    get_ACC_firing_conditional_firing_rates(ch_details, unch_details, dir_betas, ch_FRs, unch_FRs, xt)

unch_t=[];
ch_t=[];

% format of ch_details and unch_details is:
% column 1 = trial state occured in
% column 2 = the state within the trial (first, second, third, etc)
% column 3 = the value of the state
% column 4 = the screen side of the state (-1 = left, 1 = right)
% column 5 = the time of state onset

n_units = numel(dir_betas);

% get some indices of when chosen and unchosen states were on the left/right side of the screen
left_Chstate_ix = ch_details(:,4) ==-1;
right_Chstate_ix = ch_details(:,4) ==1;
left_unChstate_ix = unch_details(:,4) ==-1;
right_unChstate_ix = unch_details(:,4) ==1;


for u = 1:n_units
    
    % extract this unit's firing rates
    u_ch_FRs = squeeze(ch_FRs(:,:,u));
    u_unch_FRs = squeeze(unch_FRs(:,:,u));


    % figure out which individual states were congruent with this neuron's preferred direction
    if dir_betas(u) < 0 % this is a left-preferring neuron
        
        % pref side is ultimately chosen
        % OFC represents that side
        hv_pref_FR = u_ch_FRs(left_Chstate_ix,:);
        
        % OFC represents the other side
        lv_pref_FR = u_unch_FRs(right_unChstate_ix,:);


        %nonpref side is ultimately chosen
        % OFC represents the nonpref side
        hv_nonpref_FR = u_ch_FRs(right_Chstate_ix,:);


        % OFC represents the preferred side
        lv_nonpref_FR = u_unch_FRs(left_unChstate_ix,:);
        
        
    else % it's a right-preferring neuron
        
        % pref side is ultimately chosen
        % OFC represents that side
        hv_pref_FR = u_ch_FRs(right_Chstate_ix,:);

      
        % OFC represents the other side
        lv_pref_FR = u_unch_FRs(left_unChstate_ix,:);

        
        %nonpref side is ultimately chosen
        % OFC represents the nonpref side
        hv_nonpref_FR = u_ch_FRs(left_Chstate_ix,:);

        % OFC represents the preferred side
        lv_nonpref_FR = u_unch_FRs(right_unChstate_ix,:);


    end % of determining which side this neuron prefers
    
    
[hv_pref_mean(u,:), hv_pref_sem(u,:)] = GetMeanCI(hv_pref_FR,'sem');
[lv_pref_mean(u,:), lv_pref_sem(u,:)] = GetMeanCI(lv_pref_FR,'sem');
[hv_nonpref_mean(u,:), hv_nonpref_sem(u,:)] = GetMeanCI(hv_nonpref_FR,'sem');
[lv_nonpref_mean(u,:), lv_nonpref_sem(u,:)] = GetMeanCI(lv_nonpref_FR,'sem');
    
end % of looping over units





end % of function
    
    
    