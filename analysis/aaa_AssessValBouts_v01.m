
% aaa_AssessValBouts_v01
%---------------------------
% Does the difference in number and duration of representational bouts covary with the value difference and/or mean value (of the
% choice options)?
%---------------------------

% load the states and find how long each one was and see whether duration was a function of value difference
% find out how many states were associated with the chosen option on each trial 
%---------------------------

load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\OFC_valdec_output.mat');
% load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\ACC_valdec_output.mat');


t_mids = valdata_choice.t_mids;

f_names = unique(valdata_pics.session);

n_trials = numel(valdata_pics.ch_state(:,1));

chosen_side = bhvdata.lever+1;
chosen_side(chosen_side==0) = 1;
unchosen_side = chosen_side;
unchosen_side(chosen_side==1) = 2;
unchosen_side(chosen_side==2) = 1;

[~,pics_on_time_ix] = min(abs(t_mids - 0));
[~,max_time_ix] = min(abs(t_mids - 500));

ch_dur=[];
ch_dur_val=[];
ch_dur_val_diff=[];
ch_dur_animal = [];

unch_dur=[];
unch_dur_val=[];
unch_dur_val_diff=[];
unch_dur_animal = [];

na_dur=[];
na_dur_animal = [];


n_ch=[];
n_unch=[];
n_na = [];
n_animal=[];
n_ch_val=[];
n_unch_val=[];
n_valdiff=[];
n_ch_total_time=[];
n_unch_total_time=[];



% go through each trial and find the start indices of each state and their associated values
for t = 1:n_trials
      
    ch_durations=[];
    unch_durations=[];

    % what were the values for this trial?
    t_ch_val(t,1) = bhvdata.valbin_expval(t,chosen_side(t));
    t_unch_val(t,1) = bhvdata.valbin_expval(t,unchosen_side(t));
    
    if contains(bhvdata.subject(t),'George')
        t_animal(t,1) = 1;
    else
        t_animal(t,1) = 0;
    end
    
    % get this trial's chosen state
    ch_state = valdata_pics.ch_state(t,:);
    
    % get this trial's unchosen state
    unch_state = valdata_pics.unch_state(t,:);
    
    % get this trial's unavailable state
    na_state = valdata_pics.na_state(t,:,1);

    % find when the chosen state was 'up'
    [~, ch_state_starts, ch_state_widths] = findpeaks(ch_state);
    
    % find when the unchosen state was 'up'
    [~, unch_state_starts, unch_state_widths] = findpeaks(unch_state);
    
    % find when the not-avaiable state was 'up'
    [~, na_state_starts, na_state_widths] = findpeaks(na_state);
    
    % now find the valid state starts and durations
    ch_state_starts = ch_state_starts((ch_state_starts>pics_on_time_ix) & (ch_state_starts<max_time_ix));
    ch_state_widths = ch_state_widths((ch_state_starts>pics_on_time_ix) & (ch_state_starts<max_time_ix));

    unch_state_starts = unch_state_starts((unch_state_starts>pics_on_time_ix) & (unch_state_starts<max_time_ix));
    unch_state_widths = unch_state_widths((unch_state_starts>pics_on_time_ix) & (unch_state_starts<max_time_ix));
    
    na_state_starts = na_state_starts((na_state_starts>pics_on_time_ix) & (na_state_starts<max_time_ix));
    na_state_widths = na_state_widths((na_state_starts>pics_on_time_ix) & (na_state_starts<max_time_ix));
    
    ch_durations = mean(diff(t_mids))*ch_state_widths;
    unch_durations = mean(diff(t_mids))*unch_state_widths;
    na_durations = mean(diff(t_mids))*na_state_widths;

    
    % now wrangle the individual state data for aggregate analysis later   
    if bhvdata.trialtype(t)==2 & ~isnan(t_ch_val(t,1)) & numel(ch_durations)>0
        ch_dur = [ch_dur ; ch_durations'];
        ch_dur_val = [ch_dur_val ; ones(size(ch_durations'))*t_ch_val(t)];
        ch_dur_val_diff = [ch_dur_val_diff ; ones(size(ch_durations'))*(t_ch_val(t) - t_unch_val(t))];
        ch_dur_animal = [ch_dur_animal ; ones(size(ch_durations'))*t_animal(t)];
        
        unch_dur = [unch_dur; unch_durations'];
        unch_dur_val = [unch_dur_val ; ones(size(unch_durations'))*t_unch_val(t)];
        unch_dur_val_diff = [unch_dur_val_diff ; ones(size(unch_durations'))*(t_ch_val(t) - t_unch_val(t))];
        unch_dur_animal = [unch_dur_animal; ones(size(unch_durations'))*t_animal(t)];
        
        na_dur = [na_dur; na_durations'];
        na_dur_animal = [na_dur_animal; ones(size(na_durations'))*t_animal(t)];
        
        %-------------------
        
        % number of states on a trial
        n_ch = [n_ch; numel(ch_state_widths)]; 
        n_ch_total_time = [n_ch_total_time; sum(ch_durations)];
        n_unch_total_time = [n_unch_total_time; sum(unch_durations)];
        n_unch = [n_unch; numel(unch_state_widths)];
        n_animal = [n_animal ; t_animal(t)];
        n_ch_val = [n_ch_val; t_ch_val(t)];
        n_unch_val = [n_unch_val; t_unch_val(t)];
        n_valdiff = [n_valdiff ; t_ch_val(t) - t_unch_val(t)];
        
        n_na = [n_na ; numel(na_state_widths)];
        
        
    end


end % of cycling over trials
%---------------------------------------

%---
% REGRESSION MODELS to assess state dynamics
n_ch_tbl = table; 
n_ch_tbl.n_ch_states = n_ch; 
n_ch_tbl.n_unch_states = n_unch; 
n_ch_tbl.max_val = n_ch_val;
n_ch_tbl.val_diff = n_valdiff; 


n_ch_mdl = fitglm(n_ch_tbl(n_animal==0,:), 'n_ch_states ~ max_val + val_diff');
n_unch_mdl = fitglm(n_ch_tbl(n_animal==1,:), 'n_unch_states ~ max_val + val_diff');

ch_duration_tbl = table; 
ch_duration_tbl.ch_dur = ch_dur; 
ch_duration_tbl.max_val = ch_dur_val; 
ch_duration_tbl.val_diff = ch_dur_val_diff; 

ch_dur_mdl = fitglm(ch_duration_tbl(ch_dur_animal==1,:), 'ch_dur ~ max_val + val_diff');

unch_duration_tbl = table; 
unch_duration_tbl.unch_dur = unch_dur; 
unch_duration_tbl.max_val = unch_dur_val; 
unch_duration_tbl.val_diff = unch_dur_val_diff; 

unch_dur_mdl = fitglm(unch_duration_tbl(unch_dur_animal==1,:), 'unch_dur ~ max_val + val_diff');
%---














