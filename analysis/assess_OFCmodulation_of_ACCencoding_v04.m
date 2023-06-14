function [unch_p, ch_p, cong_p, cong_t] = assess_OFCmodulation_of_ACCencoding_v04(win_ch, win_unch, ch_details, unch_details, dir_betas, ch_FRs, unch_FRs, xt)

unch_p=[];
ch_p=[];
cong_p=[];
cong_t=[];

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
    u_ch_win = squeeze(win_ch(:,:,u));
    u_unch_win = squeeze(win_unch(:,:,u));
    
    u_ch_FRs = squeeze(ch_FRs(:,:,u));
    u_unch_FRs = squeeze(unch_FRs(:,:,u));


    % figure out which individual states were congruent with this neuron's preferred direction
    if dir_betas(u) < 0 % this is a left-preferring neuron
        
        % pref side is ultimately chosen
        % OFC represents that side
        ch_pref_win = u_ch_win(left_Chstate_ix,:);
        ch_pref_FR = u_ch_FRs(left_Chstate_ix,:);
        
        % OFC represents the other side
        unch_pref_win = u_unch_win(right_unChstate_ix,:);
        unch_pref_FR = u_unch_FRs(right_unChstate_ix,:);


        %nonpref side is ultimately chosen
        % OFC represents the nonpref side
        ch_nonpref_win = u_ch_win(right_Chstate_ix,:);
        ch_nonpref_FR = u_ch_FRs(right_Chstate_ix,:);


        % OFC represents the preferred side
        unch_nonpref_win = u_unch_win(left_unChstate_ix,:);
        unch_nonpref_FR = u_unch_FRs(left_unChstate_ix,:);
        
        
    else % it's a right-preferring neuron
        
        % pref side is ultimately chosen
        % OFC represents that side
        ch_pref_win = u_ch_win(right_Chstate_ix,:);
        ch_pref_FR = u_ch_FRs(right_Chstate_ix,:);

      
        % OFC represents the other side
        unch_pref_win = u_unch_win(left_unChstate_ix,:);
        unch_pref_FR = u_unch_FRs(left_unChstate_ix,:);

        
        %nonpref side is ultimately chosen
        % OFC represents the nonpref side
        ch_nonpref_win = u_ch_win(left_Chstate_ix,:);
        ch_nonpref_FR = u_ch_FRs(left_Chstate_ix,:);


        % OFC represents the preferred side
        unch_nonpref_win = u_unch_win(right_unChstate_ix,:);
        unch_nonpref_FR = u_unch_FRs(right_unChstate_ix,:);


    end % of determining which side this neuron prefers
    
    % these compare the effect of congruency when states do/don't match
    [~,pref_p,~,pref_stats] = ttest2(ch_pref_win, unch_pref_win);
    [~,nonpref_p,~,nonpref_stats] = ttest2(ch_nonpref_win, unch_nonpref_win);
    
    % these compare the effect of congruency, generally
    [~, congruency_p, ~, congruency_stats] = ttest2([ch_pref_win ; ch_nonpref_win], [unch_pref_win ; unch_nonpref_win]);
    
    
    cong_p(u,:) = congruency_p;
    cong_t(u,:) = congruency_stats.tstat;
    
    ch_p(u,:) = pref_p;
    unch_p(u,:) = nonpref_p;

    

    
    %*******
    % verify that differences make sense with the full firing rate traces
    %*******
%     CT2 = cbrewer('qual','Paired',12);
%     
%     if pref_p(2) < .01
%         figure; 
%         hold on
%         plot(xt, mean(ch_pref_FR),'color',CT2(6,:), 'LineWidth',2);
%         plot(xt, mean(unch_pref_FR),'color',CT2(5,:), 'LineWidth',2);
%         plot(xt, mean(ch_nonpref_FR),'color',CT2(2,:), 'LineWidth',2);
%         plot(xt, mean(unch_nonpref_FR),'color',CT2(1,:), 'LineWidth',2);
%     end


    
end % of looping over units



end % of function