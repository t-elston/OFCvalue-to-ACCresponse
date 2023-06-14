%%% aaa_Localize_val_and_dir_in_ACC_v01.m

%--------------------------------
bhv_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\bhv\';
rec_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\recdata\';
coords_dir = 'C:\Users\Thomas Elston\Documents\MATLAB\Projects\OFC-ACC\rec_info\';
alignment = 'choice';
brain_area = 'ACC';
offsets = [-1000,100];
%--------------------------------


% find names of bhv files and rec files
bhv_folder_names = get_subfolder_names(bhv_dir);

n_files = numel(bhv_folder_names);


pvals=[];
betas=[];
MLAPDV_coords=[];
animal_id=[];


% now loop over each file
for f_ix = 1:n_files
    
    current_file = bhv_folder_names(f_ix);
    monkey_id = contains(current_file,'George'); % 0 = Chap; 1 = George
    
    fprintf('\n') 
    disp(['***' current_file{1}]);
    
    % pull out this file's behavior
    [trialinfo] = load_behavior(bhv_dir, current_file);
    
    % pull out this file's neurons
    [units, u_hemi_ix, choice_ts] = load_neurons(rec_dir, current_file, alignment, brain_area, offsets);
    
    [n_trials, n_times, n_units] = size(units);
    
    % pull out the coordinates associated with these ACC units
    f_MLAPDV_coords = get_ACC_unit_coords(coords_dir, current_file);
    
    
    % get the selectivity associated with the units 
%     [f_pvals, f_betas] = fixed_window_rel2choice_GLMs_v02(units, choice_ts, trialinfo);

    % now loop over each unit and check out direction encoding
    reg_tbl = table;
    reg_tbl.trial_type = trialinfo.trialtype;
    reg_tbl.trial_type(reg_tbl.trial_type==2) = -1;
    reg_tbl.max_val = nanmax(trialinfo.subjval_expval ,[],2);
    reg_tbl.t_num = trialinfo.TrialNumber;
    
    
    % initialize output of GLMs
    sig_times = NaN(n_units, n_times,2);
    t_betas = NaN(n_units, n_times,2);
    t_CPD = NaN(n_units, n_times,2);
    
    ipsi_frs = NaN(n_units, n_times);
    contra_frs = NaN(n_units, n_times);
    
    u_sig=NaN(n_units,2);
    u_beta=NaN(n_units,3);
    
    
    % now loop over the units
    pw = PoolWaitbar(n_units, 'Assessing file units...');
    parfor u = 1:n_units
        increment(pw);
        
        [~, ~, ~, ~, u_sig(u,:), u_beta(u,:), ~] =...
            parallel_sliding_GLM_v01(reg_tbl, trialinfo.lever, units(:,:,u), choice_ts);
        
    end % of looping over units
    delete(pw);
    
    pvals = [pvals; u_sig];
    betas = [betas ; u_beta];
    MLAPDV_coords = [MLAPDV_coords ; f_MLAPDV_coords];
    animal_id = [animal_id ; ones(n_units,1)*monkey_id];

end % of looping over files

g_ix = animal_id ==1; 
c_ix = animal_id ==0;
MLAPDV_coords(g_ix,3) = MLAPDV_coords(g_ix,3)+5.5;


c_sulcus_estimate = 27; % should be 27
g_sulcus_estimate = 27.5; % should be 27.5

dir_units = pvals(:,1)==1;

g_total = sum(dir_units & g_ix);
c_total = sum(dir_units & c_ix);



all_c_coords = MLAPDV_coords(c_ix,3) - c_sulcus_estimate;
all_g_coords = MLAPDV_coords(g_ix,3) - g_sulcus_estimate;

c_dir_coords = MLAPDV_coords(c_ix & dir_units,3) - c_sulcus_estimate;
g_dir_coords = MLAPDV_coords(g_ix & dir_units ,3) - g_sulcus_estimate;

CT = cbrewer('qual','Set1',9);

figure; 
subplot(1,3,1);
hold on
histogram(all_c_coords, 'BinWidth',.2,'EdgeColor','none','Normalization','Probability', 'FaceColor', [0 0.4470 0.7410]);
plot([0,0],ylim,'k');
xlim([-4,4]);
xlabel('DV Distance from ACC Fundus (mm)')
ylabel('p(ACC units)');
title('All Animal C Units');
axis square

subplot(1,3,2);
hold on
histogram(all_g_coords, 'BinWidth',.2,'EdgeColor','none','Normalization','Probability', 'FaceColor', [0.8500 0.3250 0.0980]);
plot([0,0],ylim,'k');
xlim([-4,4]);
xlabel('DV Distance from ACC Fundus (mm)')
ylabel('p(ACC units)');
title('All Animal G Units');
axis square

subplot(1,3,3);
hold on
histogram(c_dir_coords, 'BinWidth',.2,'EdgeColor','none','Normalization','Probability');
histogram(g_dir_coords, 'BinWidth',.2,'EdgeColor','none','Normalization','Probability');

plot([0,0],ylim,'k');
xlim([-4,4]);
xlabel('DV Distance from ACC Fundus (mm)')
ylabel('p(Dir-Encoding Units)');
axis square
title('Direction-Encoding Units')
legend({'Animal C','Animal G'});




CT = cbrewer('qual','Set1',9);
CT2 = cbrewer('qual','Dark2',8);



figure; 
subplot(1,2,1);
hold on
plot(MLAPDV_coords(~dir_units & c_ix,1), MLAPDV_coords(~dir_units & c_ix,3),'o','MarkerSize',5, 'color',CT2(8,:),'LineWidth',1);
plot(MLAPDV_coords(dir_units & c_ix,1), MLAPDV_coords(dir_units & c_ix,3),'.','MarkerSize',15, 'color',CT2(1,:));
xlim([-4,4]);
ylim([24, 30]);
ylim([23,30]);
title('Animal C');
legend('Direction p>.01','Direction p<.01'); 
xlabel('ML Coordinate (mm)'); 
ylabel('DV Coordinate (mm)');
axis square


subplot(1,2,2); 
hold on 
plot(MLAPDV_coords(~dir_units & g_ix,1), MLAPDV_coords(~dir_units & g_ix,3),'o','MarkerSize',5, 'color',CT2(8,:),'LineWidth',1);
plot(MLAPDV_coords(dir_units & g_ix,1), MLAPDV_coords(dir_units & g_ix,3),'.','MarkerSize',15, 'color',CT2(1,:));
xlim([-4,4]);
ylim([24, 30]);
title('Animal G');
xlabel('ML Coordinate (mm)'); 
ylabel('DV Coordinate (mm)');
axis square









