function [trialinfo] = load_behavior(bhv_dir, current_file)

trialinfo = [];
folder_contents = dir([bhv_dir, current_file{1}]);
% look for the non-directory contents
sub_directory_ix = [folder_contents(:).isdir]; % returns logical vector
file_names = {folder_contents(~sub_directory_ix).name}';

% load the behavior
behavioral_data = load([bhv_dir, current_file{1},'\', file_names{1}]);

trialinfo = behavioral_data.trialinfo;

end % of function 