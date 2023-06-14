function [all_times, all_betas, all_hemi_ix, all_monkey_ix, all_ipsi, all_contra, choice_ts, all_pval, all_beta,all_meanFRs, all_CPD] = assessing_ACC_ramping_v01(bhv_dir, rec_dir, alignment, brain_area)

outdata = [];

offsets=[-1500, 1000];

% find names of bhv files and rec files
bhv_folder_names = get_subfolder_names(bhv_dir);

n_files = numel(bhv_folder_names);

all_times=[];
all_betas=[];
all_hemi_ix=[];
all_monkey_ix=[];
all_ipsi=[];
all_contra=[];
all_pval=[];
all_beta=[];
all_meanFRs=[];
all_CPD = [];

% now loop over each file
for f_ix = 1:n_files
    
    current_file = bhv_folder_names(f_ix);
    monkey_id = contains(current_file,'George'); % 0 = Cap; 1 = George
    
    fprintf('\n') 
    disp(['***' current_file{1}]);
    
    % pull out this file's behavior
    [trialinfo] = load_behavior(bhv_dir, current_file);
    
    % pull out this file's neurons
    [units, u_hemi_ix, choice_ts] = load_neurons(rec_dir, current_file, alignment, brain_area, offsets);
    
    if ~isempty(units)
        
        % get mean and std of the fixation period for future z-scoring
        [u_means, u_stds] = get_epoch_mean_std(rec_dir, current_file, brain_area, 'pics', [-600,-100]);
        [choice_FR, ~] = get_epoch_mean_std(rec_dir, current_file, brain_area, 'choice', [-500,0]);
        
        
        % put the fixation and choice firing rates together
        u_mean_FR = [u_means,choice_FR];
        
        % go ahead and z-score the units
        %     [z_units] = zscore_units_from_mean_and_std(units, u_means, u_stds);
        z_units = units;
        %
        % now loop over each unit and check out direction encoding
        reg_tbl = table;
        reg_tbl.trial_type = trialinfo.trialtype;
        reg_tbl.trial_type(reg_tbl.trial_type==2) = -1;
        reg_tbl.max_val = nanmax(trialinfo.subjval_expval ,[],2);
        reg_tbl.t_num = trialinfo.TrialNumber;
        
        % now loop over times and units
        [n_trials, n_times, n_units] = size(z_units);
        
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
            
            [sig_times(u,:,:), t_betas(u,:,:), ipsi_frs(u,:,:), contra_frs(u,:), u_sig(u,:), u_beta(u,:), t_CPD(u,:,:)] =...
                parallel_sliding_GLM_v01(reg_tbl, trialinfo.lever*u_hemi_ix(u), z_units(:,:,u), choice_ts);
            
        end % of looping over units
        delete(pw);
        
        
        all_times = [all_times; sig_times];
        all_betas = [all_betas; t_betas];
        all_hemi_ix = [all_hemi_ix; u_hemi_ix];
        all_monkey_ix = [all_monkey_ix ; ones(n_units,1)*monkey_id];
        all_ipsi = [all_ipsi; ipsi_frs];
        all_contra = [all_contra; contra_frs];
        all_pval = [all_pval;u_sig];
        all_beta = [all_beta;u_beta];
        all_meanFRs = [all_meanFRs; u_mean_FR];
        
        all_CPD = [all_CPD; t_CPD];
    end
    
end % of looping over files

end % of function