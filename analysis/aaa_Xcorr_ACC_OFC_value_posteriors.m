% aaa_Xcorr_ACC_OFC_value_posteriors.m

load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\ACC_valdec_output.mat');

ACCval = valdata_pics;
% ACCval = valdata_choice;
t_mids = valdata_pics.t_mids;

load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\OFC_valdec_output.mat');
OFCval = valdata_pics;
% OFCval = valdata_choice;


% get the session IDs
OFC_s_ids = unique(OFCval.session);
ACC_s_ids = unique(ACCval.session);


% get indices of a window in which to detect states
[~,win_start] = min(abs(t_mids - 0));
[~,win_end] = min(abs(t_mids - 300));

all_ch_lags =[];
all_unch_lags =[];
animal_id=[];

all_s_ch_lags=[];
all_s_unch_lags=[];

% loop over the sessions with simultaneous recordings
simultaneous_recs = intersect(OFC_s_ids, ACC_s_ids);
for s = 1:numel(simultaneous_recs)
    
    s_trialtype = bhvdata.trialtype(contains(bhvdata.session,simultaneous_recs{s}));
    
    % extract this session's indices of when states occur
    s_OFC_ch_state = OFCval.ch_state(contains(OFCval.session,simultaneous_recs{s}),:);
    s_ACC_ch_state = ACCval.ch_state(contains(ACCval.session,simultaneous_recs{s}),:);
    s_OFC_unch_state = OFCval.unch_state(contains(OFCval.session,simultaneous_recs{s}),:);
    s_ACC_unch_state = ACCval.unch_state(contains(ACCval.session,simultaneous_recs{s}),:);
    
    % identify trials with valid states
    valid_OFC_ch_trials = nansum(s_OFC_ch_state(:,win_start:win_end),2)> 0;
    valid_ACC_ch_trials = nansum(s_ACC_ch_state(:,win_start:win_end),2)> 0;
    valid_OFC_unch_trials = nansum(s_OFC_unch_state(:,win_start:win_end),2)> 0;
    valid_ACC_unch_trials = nansum(s_ACC_unch_state(:,win_start:win_end),2)> 0;
    
    % which trials had both valid ACC and OFC states?
    valid_ch_trials = valid_OFC_ch_trials & valid_ACC_ch_trials & s_trialtype ==2;
    valid_unch_trials = valid_OFC_unch_trials & valid_ACC_unch_trials & s_trialtype ==2;

  
    % extract this session's value posteriors
    s_OFC_chval_pps = OFCval.ch_ppd(contains(OFCval.session,simultaneous_recs{s}),:);
    s_ACC_chval_pps = ACCval.ch_ppd(contains(ACCval.session,simultaneous_recs{s}),:);
    s_OFC_unchval_pps = OFCval.unch_ppd(contains(OFCval.session,simultaneous_recs{s}),:);
    s_ACC_unchval_pps = ACCval.unch_ppd(contains(ACCval.session,simultaneous_recs{s}),:);
    
    % now loop over trials and do the cross-correlation
    n_trials = numel(valid_ch_trials);
    
    ch_max_lag = [];
    unch_max_lag = [];
    
    shuffled_ch_max_lag=[];
    shuffled_unch_max_lag=[];
    
    for t = 1:n_trials
        
        % check if this trial had a valid set of ACC and OFC ch states
        if valid_ch_trials(t)
            
            t_OFC_ch_pps = s_OFC_chval_pps(t,win_start:win_end);
            t_ACC_ch_pps = s_ACC_chval_pps(t,win_start:win_end);
            
            [c,lags] = xcorr(t_OFC_ch_pps,t_ACC_ch_pps);
            
            % get max lag for real data
            [max_corr, max_ix] = max(c);
            ch_max_lag(t,1) = lags(max_ix) * mean(diff(t_mids));
            
            % now do a shuffle
            [shuff_c, shuff_lags] = xcorr(shuffle(t_OFC_ch_pps), shuffle(t_ACC_ch_pps));
            
            % get max lag for shuffled data
            [~, shuffled_max_ix] = max(shuff_c);
            shuffled_ch_max_lag(t,1) = lags(shuffled_max_ix) * mean(diff(t_mids));
            

            
        else     
            ch_max_lag(t,1) = NaN;
        end
        
        % check if this trial had a valid set of ACC and OFC unch states
        if valid_unch_trials(t)
            
            t_OFC_unch_pps = s_OFC_unchval_pps(t,win_start:win_end);
            t_ACC_unch_pps = s_ACC_unchval_pps(t,win_start:win_end);
            
            [c,lags] = xcorr(t_OFC_unch_pps,t_ACC_unch_pps);
            
            % get max lag for real data
            [max_corr, max_ix] = max(c);
            unch_max_lag(t,1) = lags(max_ix) * mean(diff(t_mids));
            
            % now do a shuffle
            [shuff_c, shuff_lags] = xcorr(shuffle(t_OFC_unch_pps), shuffle(t_ACC_unch_pps));
            
            % get max lag for shuffled data
            [~, shuffled_max_ix] = max(shuff_c);
            shuffled_unch_max_lag(t,1) = lags(shuffled_max_ix) * mean(diff(t_mids));
            

            
        else     
            unch_max_lag(t,1) = NaN;
        end
        
    end % of looping over trials
    
    % accumulate the data
    all_ch_lags = [all_ch_lags; ch_max_lag];
    all_unch_lags = [all_unch_lags; unch_max_lag];
    
    all_s_ch_lags = [all_s_ch_lags; shuffled_ch_max_lag];
    all_s_unch_lags = [all_s_unch_lags; shuffled_unch_max_lag];
    
    % make an index of which animal this was from
    animal_id = [animal_id; ones(n_trials, 1) * contains(simultaneous_recs{s},'George')];


end % of looping over sessions


    figure; 
    subplot(1,2,1);
    hold on
    histogram(all_ch_lags(animal_id==0),'Normalization','Probability','BinWidth',20); 
    histogram(all_unch_lags(animal_id==0),'Normalization','Probability','BinWidth',20); 
    h = gca;
    h.YAxis.Visible = 'off';
    plot([0,0],ylim,'k','LineWidth',2);
    xlabel('Trial Max Lag (ms)');
    title('Animal C');
    legend('Chosen','Unchosen');
    axis square
    
    
    
    [h,p,~, stats] = ttest(all_ch_lags(animal_id==0));
%         [h,p,~, stats] = ttest(all_unch_lags(animal_id==0));


    subplot(1,2,2);
    hold on
    histogram(all_ch_lags(animal_id==1),'Normalization','Probability','BinWidth',20); 
    histogram(all_unch_lags(animal_id==1),'Normalization','Probability','BinWidth',20);
    h = gca;
    h.YAxis.Visible = 'off';
    plot([0,0],ylim,'k','LineWidth',2);
    xlabel('Trial Max Lag (ms)');
    title('Animal G');
    axis square

    [h,p,~, stats] = ttest(all_ch_lags(animal_id==1));

