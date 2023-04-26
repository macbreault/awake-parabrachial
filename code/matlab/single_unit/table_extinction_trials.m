% plot_laser_extinction.m
%
% Script that plots the number of trials it took for animals to reach extinction during laser condition
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 10-15-2022

clear
clc
set(0,'DefaultFigureWindowStyle','docked') %docks the figures onto the Matlab window


CI = @(x) 1.96 * std(x,[],1,'omitnan') / sqrt(size(x,2)); % Function



%% Set data

% ~~~~~~~~~~~ Set mouse information ~~~~~~~~~~~




% Other variables to change
save_it = 1; % Whether to save plots (1) or not (0)





remove_blink = 1; % Remove blink (1) or not (0)
rate_type = 'binned'; % Use 'binned' or 'sliding' window for firing rate
bin_width = 1; % seconds (MUST BE LARGER THAN 1/FRAME RATE=0.03333)
norm_it = 1; % Normalize pupil size by largest size (1) or not (0)

pupil_type = 'raw'; % 'raw' or 'filter'
neural_type = 'raw'; % 'raw' or 'filter'
T = 2; % Trigger number to plot in example

time_before = 5; % Seconds before trigger
time_after = 5; % Seconds after trigger
run_all = 1;  % Boolean as to whether to run all conditions (1) or not (0)





% ~~~~~~~~~~~~~~~~~~~~~~ DO NOT CHANGE BELOW ~~~~~~~~~~~~~~~~~~~~~~

% Figure configurations
fig_params = struct('Position',[0 0 8 5], 'Units','Inches', 'WindowStyle','normal', 'FontSize', 10);

%%%%%% Set paths %%%%%%

paths = struct('dropbox',[],'video',[],'neural',[]);

if contains(pwd,'/mac/')
    paths.dropbox = '/Users/mac/Dropbox (MIT)/Jesse';
else
    paths.dropbox = '/Users/Jessesmith 1/Dropbox/Jesse';
end

red_driver_path = '/Volumes/T7';
seagate_driver_path = '/Volumes/Seagate Backup Plus Drive';

% Check which driver is plugged in
if exist(red_driver_path,'dir')
    driver_path = red_driver_path;
elseif exist(seagate_driver_path,'dir')
    driver_path = seagate_driver_path;
else
    driver_path = [];
end

% Find path to neural data by first checking SEAGATE, then checking RED
paths.neural = fullfile(seagate_driver_path, 'Electrophysiology','Awake Recordings','Conditioning With pupil');
if ~exist(paths.neural,'dir') % Check seagate
    paths.neural = replace(paths.neural,seagate_driver_path,red_driver_path);
    if ~exist(paths.neural,'dir')% Check red drive
        paths.neural = fullfile(paths.dropbox,'data','neural'); % Else, assign dropbox
    end
end

% Find path to video/arduino data
paths.video = fullfile(driver_path,'Jesse','data');
if ~exist(paths.video,'dir')
    paths.video = fullfile(paths.dropbox,'data');
end

% Check that both paths exist
if ~exist(paths.dropbox,'dir'); error('Error in %s: Dropbox path not found', mfilename); end
if ~exist(paths.neural,'dir'); error('Error in %s: Neural hard drive not found', mfilename); end
if ~exist(paths.video,'dir'); error('Error in %s: Video hard drive not found', mfilename); end


addpath(genpath(fullfile(paths.dropbox,'code'))); % Add libraries to directory

table_path  = fullfile(paths.dropbox,'tables');
if ~exist(table_path,'dir') && save_it; mkdir(table_path); end




%% Collect data from all trials for plotting

condition_file = fullfile(paths.dropbox,'tables','Conditioning Recordings.xlsx');
mouse_ids = cellstr(sheetnames(condition_file));
exper_conds = {'test','retest'}; % What experiment were you recording
paths = rmfield(paths,'video'); % Delete pupil field so that loadData knows to skip it


if ~run_all
    mouse_ids = mouse_ids(ismember(mouse_ids,mouse_id));
end


% Initialize to collect data for each animal
trials = struct();
for i = 1:numel(exper_conds)
    trials.(exper_conds{i}) = cell(size(mouse_ids));
end


% Loop through all mice...
for mouse = 1:numel(mouse_ids)

    mouse_id = mouse_ids{mouse};

    disp(['Running mouse ',mouse_id,'...'])

    %%%%%% Import their table %%%%%%

    condition_table = readtable(condition_file, 'Sheet', mouse_id,'PreserveVariableNames',true); %'VariableNamingRule','preserve')

    % Filter out table to only use dates with exper_conds only
    dates = unique(condition_table.Date);
    keep_dates = arrayfun(@(d) all(ismember(exper_conds,lower(condition_table{condition_table.Date == dates(d),'Condition'}))), 1:numel(dates));
    date_inds = ismember(condition_table.Date, dates(keep_dates));

    condition_inds = any(cell2mat(cellfun(@(exper_cond) ismember(lower(condition_table.Condition), exper_cond),exper_conds,'un',0)),2);

    % Filter out table
    condition_table(~(date_inds & condition_inds),:) = [];

    % Filter out the conditions that do not match specifications
    if ~run_all

        % ~~~~~~~ Check which conditions should be run (if not empty) ~~~~~~~
        date_inds = true(size(condition_table,1),1);
        file_inds = true(size(condition_table,1),1);
        condition_inds = true(size(condition_table,1),1);

        % Date
        if ~isempty(date_str)
            datetime_str = datetime(date_str,'InputFormat','M-d-yy');
            date_inds = condition_table.Date == datetime_str;
        end

        % File
        if ~isempty(file_num)
            file_inds = strcmp(condition_table.File,file_num);
        end

        % Condition
        if ~isempty(exper_conds)
            condition_inds = ismember(lower(condition_table.Condition),exper_conds);
        end

        % Remove rows from condition_table that do not match conditions
        exper_inds = date_inds & file_inds & condition_inds;
        condition_table(~exper_inds,:) = [];

    end % Only keep conditions requested



    %% Loop through all dates

    dates = unique(condition_table.Date);

    for d = 1:numel(dates)

        % Extract conditions (if running all files in condition table)
        try
            date_str = strjoin(cellfun(@(num) num2str(str2double(num)),{datestr(dates(d),'mm'),datestr(dates(d),'dd'),datestr(dates(d),'yy')},'un',0),'-');

            % Check to make sure there is only 1 file
            if numel(unique(condition_table{condition_table.Date == dates(d),'File'})) == 1
                file_num = condition_table{condition_table.Date == dates(d),'File'}{1};
            else
                error('Error in %s: Found more than 1 File for this date. Skip!',mfilename)
            end

        catch
            warning('Invalid date in excel table. Skip!')
            continue
        end


        disp(' ')
        disp(['             Running ',date_str,', ',file_num, '...'])



        %% Load data

        % Set params
        params = struct('remove_blink',remove_blink,'rate_type',rate_type,'bin_width',bin_width,'norm_it',norm_it,'pupil_type',pupil_type,'neural_type',neural_type,'time_before',time_before,'time_after',time_after);


        % Try to load data from each condition
        try
            for exper_cond = exper_conds

                [~, ~, trigger_data] = loadData(mouse_id, date_str, file_num, exper_cond{:}, paths, params);

                trials.(exper_cond{:}){mouse} = [trials.(exper_cond{:}){mouse}; trigger_data.trigger_total];
            
            end
        catch ME
            if contains(ME.identifier,'MyComponent:')
                switch ME.identifier
                    case 'MyComponent:exit'
                        rethrow(ME)
                    otherwise
                        disp([ME.message, '. SKIP!'])
                        continue
                end
            else
                rethrow(ME)
            end
        end


    end % Loop through conditions

end % Loop through mice



%% Make table for saving

field_names = fieldnames(trials);
spreadsheet = struct();

% Run for each experiment type
for field = reshape(fieldnames(trials),1,[])

    table_temp = array2table(NaN(max(cellfun(@length,trials.(field{:}))), numel(mouse_ids)));
    table_temp.Properties.VariableNames = mouse_ids;
    
    for mouse = 1:numel(mouse_ids)
        table_temp{1:length(trials.(field{:}){mouse}),mouse_ids{mouse}} = trials.(field{:}){mouse};
    end

    spreadsheet.(field{:}) = table_temp;

end



%% Save table
if save_it
    
    % Delete existing table
    if exist(fullfile(table_path, 'Single-unit extinction trials.xlsx'),'file')
        delete(fullfile(table_path, 'Single-unit extinction trials.xlsx'))
    end

    % Save new table
    for field = reshape(field_names,1,[])
        writetable(spreadsheet.(field{:}), fullfile(table_path,'Single-unit extinction trials.xlsx'),...
            'Sheet',field{:},'WriteVariableNames',true,'WriteMode','overwritesheet',...
            'PreserveFormat',true)
    end

end




disp(' ')