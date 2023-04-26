function plot_folder = findPlotFolder()
% function plot_folder = findPlotFolder()
%
% Function that makes the path to the plot based on where the script is located
%
% Input:
%       - script_path: < 1 x _ > char array  
%
% Output:
%       - plot_folder: < 1 x _ > char array of what the appripriate folder name is for the plot to be saved
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 10-13-2022


%% Identify what (script or function) called findPlotFolder

[paths_info,I] = dbstack('-completenames', 1);

% If a script did not call findPlotFolder, then return empty string
if isempty(paths_info)
    plot_folder = '';
    return
end

% Only use the script/function that called findPlotFolder
path_info = paths_info(I);


%% Separate path of file that called findPlotFolder

% Split path apart
path_split = split(path_info.file,filesep);


% Find where script (path_info.name) is located
folder_ind = find(contains(path_split,path_info.name))-1;


% Use folder that contains script as folder for plotting
plot_folder = path_split{folder_ind};


end % End makePlotPath



