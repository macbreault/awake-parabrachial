% plot_photometry_percent_change.m
%
% Script that plots percent change using behavioral and activity table
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Modified: 10-10-2022

clear
clc


%% Set data

% ~~~~~~~~~~~ Set mouse information ~~~~~~~~~~~



% Other variables
save_it = 1;




% ~~~~~~~~~~~~~~~~~~~~~~ DO NOT CHANGE BELOW ~~~~~~~~~~~~~~~~~~~~~~

%%%%%% Set paths %%%%%%
paths = struct('dropbox',[]);

if contains(pwd,'/mac/')
    paths.dropbox = '/Users/mac/Dropbox (MIT)/Jesse';
else
    paths.dropbox = '/Users/Jessesmith 1/Dropbox/Jesse';
end

addpath(genpath(fullfile(paths.dropbox,'code'))); % Add libraries to directory

table_path  = fullfile(paths.dropbox,'tables');
if ~exist(table_path,'dir') && save_it; mkdir(table_path); end





%% Load data

% Get file names
behavior_file = fullfile(paths.dropbox,'tables','Photometry Behavior.xlsx');
trace_file = fullfile(paths.dropbox,'tables','Photometry Traces.xlsx');

% Get types of traces to plot
trace_types = cellstr(sheetnames(trace_file));
trace_types = trace_types(~strcmp(sheetnames(trace_file),'TEST'));
trace_types_field = cellfun(@(type) strrep(type,' ','_'), trace_types,'un',0);

% Load data
behavior_all = readtable(behavior_file,'ReadRowNames',1,'PreserveVariableNames',1);
traces_all = struct();
for i = 1:numel(trace_types)
    traces_all.(trace_types_field{i}) = readtable(trace_file,'Sheet',trace_types{i},'ReadRowNames',1,'PreserveVariableNames',1);
end



%% Prepare for precent change

mouse_ids = behavior_all.Properties.VariableNames;


% Make spreadsheet
days = unique(cell2mat([cellfun(@str2num, behavior_all.Properties.RowNames); cellfun(@(trace_type) cellfun(@str2num, traces_all.(trace_type).Properties.RowNames), strrep(trace_types,' ','_'),'un',0)]));
days = days(days > 0);

table_temp = array2table(NaN(numel(days)+1,numel(mouse_ids)));
table_temp.Properties.VariableNames = mouse_ids;
table_temp.Properties.RowNames = ['Baseline';arrayfun(@num2str, days,'un',0)];

% Make structure array
field_names = reshape(['Behavior'; trace_types_field],1,[]);
spreadsheet = struct();
for field = field_names
    spreadsheet.(field{:}) = table_temp;
end



%% Make behavior % change

% Loop through all mice...
for mouse = reshape(mouse_ids,1,[])

    plot_time = cellfun(@str2num, behavior_all(:,mouse).Properties.RowNames);
    plot_data = behavior_all{:,mouse};

    pre = plot_time <= 0;
    post = plot_time > 0;

    baseline = mean(plot_data(pre),'omitnan');
    percent_change_tmp = ((plot_data - baseline) / baseline) * 100;

    % Add to table
    spreadsheet.Behavior{'Baseline', mouse} = baseline;
    spreadsheet.Behavior{arrayfun(@num2str, plot_time(post), 'un', 0), mouse} = percent_change_tmp(post);

end



%% Make activity % change

for i = 1:numel(trace_types)

    trace_type = trace_types{i};
    trace_type_field = strrep(trace_type,' ','_');
    trace = traces_all.(trace_type_field);

    % Loop through all mice...
    for mouse = reshape(mouse_ids,1,[])

        plot_time = cellfun(@str2num, trace(:,mouse).Properties.RowNames);
        plot_data = trace{:,mouse};

        pre = plot_time <= 0;
        post = plot_time > 0;

        baseline = mean(plot_data(pre),'omitnan');
        percent_change_tmp = ((plot_data - baseline) / baseline) * 100;

        % Add to table
        spreadsheet.(trace_type_field){'Baseline', mouse} = baseline;
        spreadsheet.(trace_type_field){arrayfun(@num2str, plot_time(post), 'un', 0), mouse} = percent_change_tmp(post);

    end % Loop through mice

end % Loop through traces



%% Save table

if save_it
    
    % Delete existing table
    if exist(fullfile(table_path, 'Photometry Percent Change.xlsx'),'file')
        delete(fullfile(table_path, 'Photometry Percent Change.xlsx'))
    end

    % Save new table
    for field = reshape(field_names,1,[])
        writetable(spreadsheet.(field{:}), fullfile(table_path,'Photometry Percent Change.xlsx'),...
            'Sheet',strrep(field{:},'_',' '),'WriteRowNames',true,'WriteVariableNames',true,'WriteMode','overwritesheet',...
            'PreserveFormat',true)
    end

end

