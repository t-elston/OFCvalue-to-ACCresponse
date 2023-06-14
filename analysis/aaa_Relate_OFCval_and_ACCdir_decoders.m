% aaa_Relate_OFCval_and_ACCdir_decoders

load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\OFC_valdec_output.mat');
load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\ACC_dirdec_output.mat');

% only keep the trials that were the same across both recording types
simultaneous_sessions = intersect(valdata_choice.session, dirdata_choice.session);



t_mids = valdata_choice.t_mids+1;
n_trials = numel(valdata_choice.ch_state(:,1));

[~,choice_on_time_ix] = min(abs(t_mids - -400));
[~,max_time_ix] = min(abs(t_mids - -50));



ACC_ch=[];
ACC_unch=[];
ACC_na=[];

ch_animal_id = [];
unch_animal_id = [];
na_animal_id = [];



for s = 1:numel(simultaneous_sessions)
    
    val_s_ix = contains(bhvdata.session, simultaneous_sessions{s});
    dir_s_ix = contains(dirdata_choice.session, simultaneous_sessions{s});
        
    chosen_side = bhvdata.lever(val_s_ix)+1;
    chosen_side(chosen_side==0) = 1;
    unchosen_side = chosen_side;
    unchosen_side(chosen_side==1) = 2;
    unchosen_side(chosen_side==2) = 1;
    
    n_trials = numel(chosen_side); 
    
    trialtype = bhvdata.trialtype(val_s_ix);
    
    s_values = bhvdata.valbin_expval(val_s_ix,:);
    
    val_ch_state = valdata_choice.ch_state(val_s_ix,:);
    val_unch_state = valdata_choice.unch_state(val_s_ix,:);
    val_na_states = valdata_choice.na_state(val_s_ix,:,:);
    
    ACC_ch_pp = dirdata_choice.postprob_ch(dir_s_ix,:);
    
    animal_id = contains(simultaneous_sessions(s),'George');


    % go through each trial and find the start indices of each state and their associated values
    for t = 1:n_trials
        
        % what were the values for this trial?
        t_ch_val(t,1) = s_values(t,chosen_side(t));
        t_unch_val(t,1) = s_values(t,unchosen_side(t));
        
        ch_state = val_ch_state(t,:);
        
        % get this trial's unchosen state
        unch_state = val_unch_state(t,:);
        
        % get this trial's non-available state
        na_state = val_na_states(t,:,1);
        
        % find when the chosen state was 'up'
        [~, ch_state_starts] = findpeaks(ch_state);
        
        % find when the unchosen state was 'up'
        [~, unch_state_starts] = findpeaks(unch_state);
        
        % find when the non-available state was 'up'
        [~, na_state_starts] = findpeaks(na_state);
        
        % get the indices of the state starts that happened prior to the choice
        ch_state_ix = ch_state_starts((ch_state_starts>=choice_on_time_ix) & (ch_state_starts<=max_time_ix));
        unch_state_ix = unch_state_starts((unch_state_starts>=choice_on_time_ix) & (unch_state_starts<=max_time_ix));
        na_state_ix = na_state_starts((na_state_starts>=choice_on_time_ix) & (na_state_starts<=max_time_ix));
        
        % now get the ACC ch_pps for the 3 OFC alignments
        ch_pp=[];
        unch_pp=[];
        na_pp=[];
        
        % only include trials where the options had different values and it was free choices
        if trialtype(t) == 2 & t_ch_val(t,1) ~=4
            
            for ch_ix = 1:numel(ch_state_ix)
                
                ch_pp(ch_ix,:) = ACC_ch_pp(t, ch_state_ix(ch_ix)-80 : ch_state_ix(ch_ix)+ 80);
                
            end
            
            for unch_ix = 1:numel(unch_state_ix)
                
                unch_pp(unch_ix,:) = ACC_ch_pp(t, unch_state_ix(unch_ix)-80 : unch_state_ix(unch_ix)+ 80);
                
            end
            
            for na_ix = 1:numel(na_state_ix)
                
                na_pp(na_ix,:) = ACC_ch_pp(t, na_state_ix(na_ix)-80 : na_state_ix(na_ix)+ 80);
                
            end
            
            % now save the ACC posteriors
            ACC_ch = [ACC_ch ; ch_pp];
            ACC_unch = [ACC_unch ; unch_pp];
            ACC_na = [ACC_na ; na_pp];
            
            % now keep track of which animal this was
            if ~isempty(ch_pp)
                ch_animal_id = [ch_animal_id; ones(numel(ch_pp(:,1)),1) * animal_id];
            end
            if ~isempty(unch_pp)
                unch_animal_id = [unch_animal_id; ones(numel(unch_pp(:,1)),1) * animal_id];
            end
            if ~isempty(na_pp)
                na_animal_id = [na_animal_id; ones(numel(na_pp(:,1)),1) * animal_id];
            end

        end
        
        
    end % of cycling over trials
    
end % of cycling over sessions
%---------------------------------------

xt = [1:numel(ACC_ch(1,:))]*mean(diff(t_mids)) - (round(numel(ACC_ch(1,:))/2) * mean(diff(t_mids)));

[chap_ch_mean, chap_ch_sem] = GetMeanCI(ACC_ch(ch_animal_id==0,:), 'sem');
[chap_unch_mean, chap_unch_sem] = GetMeanCI(ACC_unch(unch_animal_id==0,:), 'sem');
[chap_na_mean, chap_na_sem] = GetMeanCI(ACC_na(na_animal_id==0,:), 'sem');

[george_ch_mean, george_ch_sem] = GetMeanCI(ACC_ch(ch_animal_id==1,:), 'sem');
[george_unch_mean, george_unch_sem] = GetMeanCI(ACC_unch(unch_animal_id==1,:), 'sem');
[george_na_mean, george_na_sem] = GetMeanCI(ACC_na(na_animal_id==1,:), 'sem');


CT = cbrewer('qual','Set1',9);

[~,xt_start] = min(abs(xt-0));
[~,xt_end] = min(abs(xt-100));

ch_window_mean = mean(ACC_ch(:,xt_start:xt_end),2);
unch_window_mean = mean(ACC_unch(:,xt_start:xt_end),2);
na_window_mean = mean(ACC_na(:,xt_start:xt_end),2);





figure; 
subplot(1,2,1); 
hold on
shadedErrorBar(xt,chap_ch_mean, chap_ch_sem,'LineProps',{'LineWidth',2,'color',CT(4,:)});
shadedErrorBar(xt,chap_unch_mean, chap_unch_sem,'LineProps',{'LineWidth',2,'color',CT(2,:)});
shadedErrorBar(xt,chap_na_mean, chap_na_sem,'LineProps',{'LineWidth',2,'color',CT(9,:)});
plot([0,0],ylim,'k','LineWidth',1);
xlabel('Time from OFC State Onset (ms)')
ylabel('DIRacc post. prob.')

axis square
title('Animal C'); 

subplot(1,2,2); 
hold on
shadedErrorBar(xt,george_ch_mean, george_ch_sem,'LineProps',{'LineWidth',2,'color',CT(4,:)});
shadedErrorBar(xt,george_unch_mean, george_unch_sem,'LineProps',{'LineWidth',2,'color',CT(2,:)});
shadedErrorBar(xt,george_na_mean, george_na_sem,'LineProps',{'LineWidth',2,'color',CT(9,:)});
plot([0,0],ylim,'k','LineWidth',1);
axis square
title('Animal G'); 






