%%% aaa_Assess_OFC_DIRdecoding_UnitEncoding

%-------------------------------
bhv_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\bhv\';
rec_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\recdata\';
alignment = 'choice';
brain_area = 'ACC';
offsets=[-600, 0];
n_boots = 1000;
%-------------------------------

% find names of bhv files and rec files
bhv_folder_names = get_subfolder_names(bhv_dir);

% load the OFC and ACC decoder output
load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\OFC_valdec_output.mat');
OFC_sessions = unique(valdata_choice.session);
t_mids = valdata_choice.t_mids;

load('C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\decoder_output\ACC_valdec_output.mat');
ACC_sessions = unique(valdata_choice.session);

% now list the indices of the peak representation strengths (see aaa_Assess_OFC_DIRencoding_UnitEncoding.m)
% these are the times, in milliseconds relative to choice, that value decoding peaks
ACC_val_peak_times = [-230, -165]; % first animal C, then animal G
OFC_val_peak_times = [-240, -145]; % first animal C, then animal G

% these are the times, in milliseconds relative to choice, that direction decoding peaks
ACC_dir_peak_times = [-150, -115]; % first animal C, then animal G
OFC_dir_peak_times = [-250, -200]; % first animal C, then animal G

files_for_decoding = unique([OFC_sessions ; ACC_sessions]);

n_files = numel(files_for_decoding);


% initialize decoder output
all_ch_dir = [];
all_dir_animal = [];
ctr=0;

n_neurons = 5:5:130;


% now loop over each file
for f_ix = 1:n_files
    
    current_file = files_for_decoding(f_ix);
    monkey_id = contains(current_file,'George'); % 0 = Chap; 1 = George
    
    fprintf('\n')
    disp([current_file{1}]);
    
    % pull out this file's behavior
    [trialinfo] = load_behavior(bhv_dir, current_file);
    
    % pull out this file's neurons
    [ACC_units, ACC_hemi_ix, choice_ts] = load_neurons(rec_dir, current_file, alignment, 'ACC', offsets);
    [OFC_units, OFC_hemi_ix, choice_ts] = load_neurons(rec_dir, current_file, alignment, 'OFC', offsets);
    
    % get the peak decoding times for this animal
    ACC_val_time = ACC_val_peak_times(monkey_id+1);
    ACC_dir_time = ACC_dir_peak_times(monkey_id+1);
    
    OFC_val_time = OFC_val_peak_times(monkey_id+1);
    OFC_dir_time = OFC_dir_peak_times(monkey_id+1);
    
    % get the indices from choice_ts for the beginning/end of fixed windows centered on the decoding peaks
    [~, ACC_val_start_ix] = min(abs(choice_ts - (ACC_val_time-100)));
    [~, ACC_val_end_ix] = min(abs(choice_ts - (ACC_val_time+100)));

    [~, ACC_dir_start_ix] = min(abs(choice_ts - (ACC_dir_time-100)));
    [~, ACC_dir_end_ix] = min(abs(choice_ts - (ACC_dir_time+100)));

    [~, OFC_val_start_ix] = min(abs(choice_ts - (OFC_val_time - 100)));
    [~, OFC_val_end_ix] = min(abs(choice_ts - (OFC_val_time + 100)));

    [~, OFC_dir_start_ix] = min(abs(choice_ts - (OFC_dir_time - 100)));
    [~, OFC_dir_end_ix] = min(abs(choice_ts - (OFC_dir_time + 100)));

    
    % check whether there are ACC / OFC neurons for this file
    has_ACC_units = ~isempty(ACC_units);
    has_OFC_units = ~isempty(OFC_units);
    
    %------------------
    % begin neuron-dropping / bootstrapping for ACC
    if has_ACC_units & any(contains(ACC_sessions, current_file))
        
        disp('- analyzing ACC neurons');
        
        % get mean firing rates for the fixed window
        ACC_val_win_FRs = squeeze(nanmean(ACC_units(:, ACC_val_start_ix: ACC_val_end_ix,:),2));
        ACC_dir_win_FRs = squeeze(nanmean(ACC_units(:, ACC_dir_start_ix: ACC_dir_end_ix,:),2));
       
        [~, n_ACC_units] = size(ACC_dir_win_FRs);
        
        % initialize output of parallel looping
        ACC_val_output = NaN(n_boots, 26);
        ACC_dir_output = NaN(n_boots, 26);
        
        parfor b = 1 : n_boots
            [ACC_val_output(b,:), ACC_dir_output(b,:), n_neurons] = neuron_drop_value_and_direction_v01(ACC_val_win_FRs, ACC_dir_win_FRs, trialinfo);
        end
        
        ACC_val_acc(f_ix,:) = nanmean(ACC_val_output);
        ACC_dir_acc(f_ix,:) = nanmean(ACC_dir_output);
        
        % compute the value and direction slopes for this session
        %         ACC_val_reg = regress(nanmean(ACC_val_output)', [ones(size(n_neurons')), n_neurons']);
        %         ACC_dir_reg = regress(nanmean(ACC_dir_output)', [ones(size(n_neurons')), n_neurons']);
        %
        %         ACC_slopes(f_ix,1) = ACC_val_reg(2);
        %         ACC_slopes(f_ix,2) = ACC_dir_reg(2);
        
        ACC_slopes(f_ix,1) = nanmean(diff(nanmean(ACC_val_output)));
        ACC_slopes(f_ix,2) = nanmean(diff(nanmean(ACC_dir_output)));

    else
        ACC_val_acc(f_ix,:) = NaN(1, 26);
        ACC_dir_acc(f_ix,:) = NaN(1, 26);
        ACC_slopes(f_ix,:) = NaN(1,2);

        
    end % of doing neuron dropping for ACC units
    %------------------
    
    %------------------
    % begin neuron-dropping / bootstrapping for OFC
    if has_OFC_units & any(contains(OFC_sessions, current_file))
        disp('- analyzing OFC neurons');
        
        % get mean firing rates for the fixed window
        OFC_val_win_FRs = squeeze(nanmean(OFC_units(:, OFC_val_start_ix: OFC_val_end_ix,:),2));
        OFC_dir_win_FRs = squeeze(nanmean(OFC_units(:, OFC_dir_start_ix: OFC_dir_end_ix,:),2));
        
        % initialize output of parallel looping
        OFC_val_output = NaN(n_boots, 26);
        OFC_dir_output = NaN(n_boots, 26);
        
        parfor b = 1 : n_boots
            [OFC_val_output(b,:), OFC_dir_output(b,:), n_neurons] = neuron_drop_value_and_direction_v01(OFC_val_win_FRs, OFC_dir_win_FRs, trialinfo);
        end
        
        OFC_val_acc(f_ix,:) = nanmean(OFC_val_output);
        OFC_dir_acc(f_ix,:) = nanmean(OFC_dir_output);
        
        % compute the value and direction slopes for this session
%         OFC_val_reg = regress(nanmean(OFC_val_output)', [ones(size(n_neurons')), n_neurons']);
%         OFC_dir_reg = regress(nanmean(OFC_dir_output)', [ones(size(n_neurons')), n_neurons']);
%         
%         OFC_slopes(f_ix,1) = OFC_val_reg(2);
%         OFC_slopes(f_ix,2) = OFC_dir_reg(2);
        
        OFC_slopes(f_ix,1) = nanmean(diff(nanmean(OFC_val_output)));
        OFC_slopes(f_ix,2) = nanmean(diff(nanmean(OFC_dir_output)));
        
                

    else
        OFC_val_acc(f_ix,:) = NaN(1, 26);
        OFC_dir_acc(f_ix,:) = NaN(1, 26);
        OFC_slopes(f_ix,:) = NaN(1,2);
          
    end % of doing neuron dropping for OFC units
    %------------------
    
    animal_id(f_ix,1) = monkey_id;
    
    
end % of looping over files


CT = cbrewer('qual','Dark2',8);

% get means and SEMs for each condition
[c_OFC_val_mean, c_OFC_val_sem] = GetMeanCI(OFC_val_acc(animal_id==0,:),'sem');
[c_ACC_val_mean, c_ACC_val_sem] = GetMeanCI(ACC_val_acc(animal_id==0,:),'sem');
[c_OFC_dir_mean, c_OFC_dir_sem] = GetMeanCI(OFC_dir_acc(animal_id==0,:),'sem');
[c_ACC_dir_mean, c_ACC_dir_sem] = GetMeanCI(ACC_dir_acc(animal_id==0,:),'sem');

[g_OFC_val_mean, g_OFC_val_sem] = GetMeanCI(OFC_val_acc(animal_id==1,:),'sem');
[g_ACC_val_mean, g_ACC_val_sem] = GetMeanCI(ACC_val_acc(animal_id==1,:),'sem');
[g_OFC_dir_mean, g_OFC_dir_sem] = GetMeanCI(OFC_dir_acc(animal_id==1,:),'sem');
[g_ACC_dir_mean, g_ACC_dir_sem] = GetMeanCI(ACC_dir_acc(animal_id==1,:),'sem');

[c_OFC_b_means, c_OFC_b_sems] = GetMeanCI(OFC_slopes(animal_id==0,:),'sem');
[c_ACC_b_means, c_ACC_b_sems] = GetMeanCI(ACC_slopes(animal_id==0,:),'sem');
[g_OFC_b_means, g_OFC_b_sems] = GetMeanCI(OFC_slopes(animal_id==1,:),'sem');
[g_ACC_b_means, g_ACC_b_sems] = GetMeanCI(ACC_slopes(animal_id==1,:),'sem');


Fig = figure;
set(Fig,'renderer','Painters');
set(Fig,'Position',[300 300 800 500]);
chap_val_ax       = axes('Position',[.05, .55, .25, .3]); 
chap_val_slope_ax = axes('Position',[.35, .55, .1, .3]); 

george_val_ax       = axes('Position',[.55, .55, .25, .3]); 
george_val_slope_ax = axes('Position',[.85, .55, .1, .3]); 

chap_dir_ax       = axes('Position',[.05, .1, .25, .3]); 
chap_dir_slope_ax = axes('Position',[.35, .1, .1, .3]); 

george_dir_ax       = axes('Position',[.55, .1, .25, .3]); 
george_dir_slope_ax = axes('Position',[.85, .1, .1, .3]); 

% let's compute the slopes





axes(chap_val_ax);
hold on
plot(n_neurons,OFC_val_acc(animal_id==0,:)','color',CT(1,:));
plot(n_neurons,ACC_val_acc(animal_id==0,:)','color',CT(3,:));
shadedErrorBar(n_neurons, c_OFC_val_mean, c_OFC_val_sem,'LineProps',{'color',CT(1,:), 'LineWidth',2});
shadedErrorBar(n_neurons, c_ACC_val_mean, c_ACC_val_sem,'LineProps',{'color',CT(3,:), 'LineWidth',2});
xlabel('Number of Neurons');
ylabel('Value Accuracy');
title('Animal C')

axes(chap_val_slope_ax);
hold on
plot(ones(size(OFC_slopes(animal_id==0,1))),OFC_slopes(animal_id==0,1),'.','color',CT(1,:),'MarkerSize',15);
plot(ones(size(ACC_slopes(animal_id==0,1)))*1.5,ACC_slopes(animal_id==0,1),'.','color',CT(3,:),'MarkerSize',15);
errorbar([1,1.5],[c_OFC_b_means(1), c_ACC_b_means(1)], [c_OFC_b_sems(1), c_ACC_b_sems(1)],...
    'Marker','.','MarkerSize',20, 'color',CT(8,:),'CapSize',0,'LineWidth',2);
xlim([.8,1.7]);
xticklabels({'OFC','ACC'});
title('info gain / unit');

axes(george_val_ax);
hold on
plot(n_neurons,OFC_val_acc(animal_id==1,:)','color',CT(1,:));
plot(n_neurons,ACC_val_acc(animal_id==1,:)','color',CT(3,:));
shadedErrorBar(n_neurons, g_OFC_val_mean, g_OFC_val_sem,'LineProps',{'color',CT(1,:), 'LineWidth',2});
shadedErrorBar(n_neurons, g_ACC_val_mean, g_ACC_val_sem,'LineProps',{'color',CT(3,:), 'LineWidth',2});
title('Animal G')

axes(george_val_slope_ax);
hold on
plot(ones(size(OFC_slopes(animal_id==1,1))),OFC_slopes(animal_id==1,1),'.','color',CT(1,:),'MarkerSize',15);
plot(ones(size(ACC_slopes(animal_id==1,1)))*1.5,ACC_slopes(animal_id==1,1),'.','color',CT(3,:),'MarkerSize',15);
errorbar([1,1.5],[g_OFC_b_means(1), g_ACC_b_means(1)], [g_OFC_b_sems(1), g_ACC_b_sems(1)],...
    'Marker','.','MarkerSize',20, 'color',CT(8,:),'CapSize',0,'LineWidth',2);
xlim([.8,1.7]);
xticklabels({'OFC','ACC'});
title('info gain / unit');

axes(chap_dir_ax);
hold on
plot(n_neurons,OFC_dir_acc(animal_id==0,:)','color',CT(1,:));
plot(n_neurons,ACC_dir_acc(animal_id==0,:)','color',CT(3,:));
shadedErrorBar(n_neurons, c_OFC_dir_mean, c_OFC_dir_sem,'LineProps',{'color',CT(1,:), 'LineWidth',2});
shadedErrorBar(n_neurons, c_ACC_dir_mean, c_ACC_dir_sem,'LineProps',{'color',CT(3,:), 'LineWidth',2});
yticks([.5,.6,.7]);
ylabel('Direction Accuracy')

axes(chap_dir_slope_ax);
hold on
plot(ones(size(OFC_slopes(animal_id==0,2))),OFC_slopes(animal_id==0,2),'.','color',CT(1,:),'MarkerSize',15);
plot(ones(size(ACC_slopes(animal_id==0,2)))*1.5,ACC_slopes(animal_id==0,2),'.','color',CT(3,:),'MarkerSize',15);
errorbar([1,1.5],[c_OFC_b_means(2), c_ACC_b_means(2)], [c_OFC_b_sems(2), c_ACC_b_sems(2)],...
    'Marker','.','MarkerSize',20, 'color',CT(8,:),'CapSize',0,'LineWidth',2);
xlim([.8,1.7]);


axes(george_dir_ax);
hold on
plot(n_neurons,OFC_dir_acc(animal_id==1,:)','color',CT(1,:));
plot(n_neurons,ACC_dir_acc(animal_id==1,:)','color',CT(3,:));
shadedErrorBar(n_neurons, g_OFC_dir_mean, g_OFC_dir_sem,'LineProps',{'color',CT(1,:), 'LineWidth',2});
shadedErrorBar(n_neurons, g_ACC_dir_mean, g_ACC_dir_sem,'LineProps',{'color',CT(3,:), 'LineWidth',2});

axes(george_dir_slope_ax);
hold on
plot(ones(size(OFC_slopes(animal_id==1,2))),OFC_slopes(animal_id==1,2),'.','color',CT(1,:),'MarkerSize',15);
plot(ones(size(ACC_slopes(animal_id==1,2)))*1.5,ACC_slopes(animal_id==1,2),'.','color',CT(3,:),'MarkerSize',15);
errorbar([1,1.5],[g_OFC_b_means(2), g_ACC_b_means(2)], [g_OFC_b_sems(2), g_ACC_b_sems(2)],...
    'Marker','.','MarkerSize',20, 'color',CT(8,:),'CapSize',0,'LineWidth',2);
xlim([.8,1.7]);
