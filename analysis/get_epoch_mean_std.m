function [u_means, u_stds] = get_epoch_mean_std(rec_dir, current_file, brain_area, align_event, offsets)

units = []; % unit data: n_trials x n_times x n_units
u_hemi_ix = []; % which hemi are the units in? n_units x 1; -1 = left hemi; 1 = right hemi
ts = []; % timesteps - when, relative to the aligment

u_means = []; % mean of units during fixation
u_stds  = []; % std of units during fixation

folder_contents = dir([rec_dir, current_file{1}]);
% look for the non-directory contents
sub_directory_ix = [folder_contents(:).isdir]; % returns logical vector
file_names = {folder_contents(~sub_directory_ix).name}';

% now find the files with the proper alignment
aligned_ix = contains(file_names,align_event);

% now find files with the proper brain area
brain_area_ix = contains(file_names,brain_area);

% use these indices to get the names of the relevant file(s)
relevant_files = file_names(brain_area_ix & aligned_ix);

% initialize some empty arrays for the data from each hemisphere (if
% applicable)
left_hemi_data = [];
left_hemi_ix = [];
right_hemi_data = [];
right_hemi_ix = [];

% now load the neurons

% do we have left hemisphere data?
if any(contains(relevant_files,'L'))
    left_hemi_data = cell2mat(struct2cell(load([rec_dir, current_file{1},'/',relevant_files{contains(relevant_files,'L')}], 'SPKfr_units')));
    ts = cell2mat(struct2cell(load([rec_dir, current_file{1},'/',relevant_files{contains(relevant_files,'L')}], 't_mids')));
    left_hemi_ix = ones(numel(left_hemi_data(1,1,:)),1)*-1;
    
end

% do we have any right hemisphere data?
if any(contains(relevant_files,'R'))
    right_hemi_data = cell2mat(struct2cell(load([rec_dir, current_file{1},'/',relevant_files{contains(relevant_files,'R')}], 'SPKfr_units')));
    right_hemi_ix = ones(numel(right_hemi_data(1,1,:)),1);
end

units = cat(3, left_hemi_data, right_hemi_data);
u_hemi_ix = [left_hemi_ix ; right_hemi_ix];

ts = ts+1;

% get window indices for taking means
[~,w_start] = min(abs(ts - offsets(1)));
[~,w_end] = min(abs(ts - -offsets(2)));

% units(units==0) = NaN;

u_means = nanmean(squeeze(nanmean(units(:,w_start:w_end,:),2)))';
u_stds = nanstd(squeeze(nanmean(units(:,w_start:w_end,:),2)))';

end % of function