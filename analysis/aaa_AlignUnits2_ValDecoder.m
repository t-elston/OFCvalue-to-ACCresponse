% aaa_AlignUnits2_ValDecoder

% load the value decoder output and see which files were used. We'll use those same files
load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\OFC_valdec_output.mat');
% load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\ACC_valdec_output.mat');

bhv_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\bhv\';
rec_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\recdata\';
alignment = 'pics';
brain_area = 'OFC';
offsets=[-500, 1000];

files_for_decoding = unique(bhvdata.session);

n_files = numel(files_for_decoding);

match_p=[];
nonmatch_p=[];
LR_betas=[];
u_animal=[];

% now loop over each file
for f_ix = 1:n_files
    
    current_file = files_for_decoding(f_ix);
    animal_id = contains(current_file,'George'); % 0 = Chap; 1 = George
    
    fprintf('\n')
    disp(['***' current_file{1}]);
    
    % pull out this file's behavior
    [trialinfo] = load_behavior(bhv_dir, current_file);
    
    choice_dirs = trialinfo.lever;
    
    % get the chosen and unchosen values for this session
    trial_values = trialinfo.valbin_expval;
    
    % pull out this file's neurons
    [units, u_hemi_ix, ts] = load_neurons(rec_dir, current_file, alignment, brain_area, offsets);
    
    [n_trials, n_times, n_units] = size(units);
    
    f_match_p=NaN(n_units, 41);
    f_nonmatch_p =NaN(n_units, 41);
    f_LR_betas=NaN(n_units, 2);
    
    
    % now loop over the individual units; hold out one and then train the decoder
    pw = PoolWaitbar(n_units, 'assessing units...');
    parfor u = 1:n_units
        increment(pw);
        
        %  disp(['u: ' num2str(u) ' / ' num2str(n_units)]);
        
        % get the rest of the units, excluding the held out one
        other_units = units;
        other_units(:,:,u)=[];
        
        % run the value decoder and get the chosen and unchosen states
        [ch_states, unch_states] = run_value_LDA_and_get_states(other_units, trialinfo, ts);
        
        % collect the details for the states
        [ch_val, unch_val, ch_nStates, unch_nStates, ch_times, unch_times] = get_trialState_details_v02(trial_values,choice_dirs,...
            animal_id,ch_states, unch_states,ts);
        
        % get long-form details about the states
        [ch_details, unch_details] = get_longform_state_details(ch_nStates, unch_nStates, ch_times, unch_times,...
            ch_val, unch_val, choice_dirs, trialinfo.trialtype);
        
        % now align this unit's FRs to the chosen and unchosen states
        [ch_FRs, unch_FRs, xt] = align_unit_to_states(units(:,:,u), ts, ch_details, unch_details);
        
        % now regress the held out neuron's firing rates agains the value of the states
        [f_match_p(u,:), f_nonmatch_p(u,:), f_LR_betas(u,:)] = regress_FRs_against_state_val(ch_FRs, unch_FRs, ch_details, unch_details, xt);
        
    end % of looping over units
    delete(pw);
    
    
    % now acculumulate the results
    
    match_p = [match_p ; f_match_p];
    nonmatch_p = [nonmatch_p ; f_nonmatch_p];
    LR_betas = [LR_betas ; f_LR_betas];
    u_animal = [u_animal ; ones(n_units,1)*animal_id];
    
    
end % of looping over files

xt = [1:41]*mean(diff(ts)) - 500;

% now plot
CT = cbrewer('qual','Set1',9);
CT2 = cbrewer('qual','Dark2',8);
cmap_c = 3;

n_g = sum(u_animal ==1);
n_c = sum(u_animal ==0);


% run chi squares on each time point for each animal
for xi = 1:numel(xt)
    

    [chap_p(xi),~] = chisquarecont([ sum(match_p(u_animal==0,xi)<.01), n_c ; 
                                     sum(nonmatch_p(u_animal==0,xi)<.01), n_c]);
                                 
    [george_p(xi),~] = chisquarecont([ sum(match_p(u_animal==1,xi)<.01), n_g ; 
                                     sum(nonmatch_p(u_animal==1,xi)<.01), n_g]);                              

end

c_X2_sig = NaN(xi,1);
c_X2_sig(chap_p<.01) = 1;

g_X2_sig = NaN(xi,1);
g_X2_sig(george_p<.01) = 1;


figure;
subplot(2,2,1);
n_units = numel(match_p(:,1));
hold on
bar(xt,mean(match_p(u_animal==0,:)<.01), 'EdgeColor','none', 'FaceColor',CT2(cmap_c,:),'BarWidth', 1);
bar(xt,mean(nonmatch_p(u_animal==0,:)<.01), 'EdgeColor','none', 'FaceColor',CT(9,:),'BarWidth', 1);
plot(xt, c_X2_sig*.2, 'LineWidth',2,'color','k')
ylim([0, .25]);
yticks([0,.1,.2, .3, .4]);
plot([0,0],ylim,'k','LineStyle',':');
xlim([-400, 400]);
xlabel('Time from State Onset (ms)');
ylabel('Prop. Units p<.01');
set(gca, 'FontSize',10, 'LineWidth',1);
axis square
title('Animal C');


subplot(2,2,2);
hold on
bar(xt,mean(match_p(u_animal==1,:)<.01), 'EdgeColor','none', 'FaceColor',CT2(cmap_c,:),'BarWidth', 1);
bar(xt,mean(nonmatch_p(u_animal==1,:)<.01), 'EdgeColor','none', 'FaceColor',CT(9,:),'BarWidth', 1);
plot(xt, g_X2_sig*.2, 'LineWidth',2,'color','k')
ylim([0, .25]);
yticks([0,.1,.2]);
plot([0,0],ylim,'k','LineStyle',':');
xlim([-400, 400]);
xlabel('Time from State Onset (ms)');
ylabel('Prop. Units p<.01');
set(gca, 'FontSize',10, 'LineWidth',1);
axis square
title('Animal G');



subplot(2,2,3); 
hold on
plot(LR_betas(u_animal==0,1),LR_betas(u_animal==0,2),'.','color',CT2(cmap_c,:));
xlim([-.5,.5]);
ylim([-.5,.5]);
plot(xlim,[0,0],'k','LineWidth',1);
plot([0,0],ylim,'k','LineWidth',1);


l = lsline; 
l.Color = [.5,.5,.5];
l.LineWidth=1;

% xlabel('\beta left | left state');
% ylabel('\beta right | right state');
set(gca, 'FontSize',10, 'LineWidth',1);
[chap_r, chap_p] = corr(LR_betas(u_animal==0,1),LR_betas(u_animal==0,2));
axis square


subplot(2,2,4); 
hold on
plot(LR_betas(u_animal==1,1),LR_betas(u_animal==1,2),'.','color',CT2(cmap_c,:));
xlim([-.5,.5]);
ylim([-.5,.5]);
plot(xlim,[0,0],'k','LineWidth',1);
plot([0,0],ylim,'k','LineWidth',1);


l = lsline; 
l.Color = [.5,.5,.5];
l.LineWidth=1;

% xlabel('\beta left | left state');
% ylabel('\beta right | right state');
set(gca, 'FontSize',10, 'LineWidth',1);
[george_r, george_p] = corr(LR_betas(u_animal==1,1),LR_betas(u_animal==1,2));
axis square





