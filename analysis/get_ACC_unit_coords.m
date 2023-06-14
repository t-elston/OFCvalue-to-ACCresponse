function [f_coords] = get_ACC_unit_coords(coords_dir, current_file)
f_coords=[];

folder_contents = dir([coords_dir, current_file{1}]);
% look for the non-directory contents
sub_directory_ix = [folder_contents(:).isdir]; % returns logical vector
file_names = {folder_contents(~sub_directory_ix).name}';

% find the .mat file
relevant_file = file_names(contains(file_names,'.mat'));

unit_data = load([coords_dir, current_file{1},'/',relevant_file{1}]);

f_coords = unit_data.coords_ML_AP_DV(contains(unit_data.region,'ACC'),:);

end % of function