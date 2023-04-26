% run_evoked_neural.m
%
% Script that plots firing rate comparing conditions
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 09-16-2022 - split from run_evoked.m


clear
clc
set(0,'DefaultFigureWindowStyle','docked') %docks the figures onto the Matlab window


%CI = @(x) 1.96 * std(x,[],1,'omitnan') / sqrt(size(x,2)); % Function
step_plot = @(x) reshape([x;x],1,2*numel(x));



%% Set data

% ~~~~~~~~~~~ Set mouse information ~~~~~~~~~~~
mouse_id = 'C7';                     % Mouse ID
date_str = '6-29-22';                     % Date of recording
file_num = 'File1';                     % File number ('File1')



% Other variables to change
run_all = 1;  % Boolean as to whether to run all conditions (1) or not (0)
save_it = 1; % Whether to save plots (1) or not (0)

remove_blink = 1; % Remove blink (1) or not (0)
rate_type = 'binned'; % Use 'binned' or 'sliding' window for firing rate
bin_width = 0.1; % seconds (MUST BE LARGER THAN 1/FRAME RATE=0.03333)
norm_it = 1; % Normalize pupil size by largest size (1) or not (0)

pupil_type = 'raw'; % 'raw' or 'filter'
neural_type = 'raw'; % 'raw' or 'filter'
T = 2; % Trigger number to plot in example

time_before = 5; % Seconds before trigger
time_after = 5; % Seconds after trigger




% ~~~~~~~~~~~~~~~~~~~~~~ DO NOT CHANGE BELOW ~~~~~~~~~~~~~~~~~~~~~~

% Figure configurations
fig_params = struct('Position',[0 0 8 3.5], 'Units','Inches', 'WindowStyle','normal', 'FontSize', 10);

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



%% Loop through all conditions

condition_file = fullfile(paths.dropbox,'tables','Conditioning Recordings.xlsx');
mouse_ids = sheetnames(condition_file);
exper_conds = {'test','tone alone','retest','extinction'}; % All conditions that are needed to complete this plot
paths = rmfield(paths,'video'); % Delete pupil field so that loadData knows to skip it


if ~run_all
    mouse_ids = mouse_ids(ismember(mouse_ids,mouse_id));
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

        % Initialize structure that hold onto data
        exper_data = struct();

        % Try to load data from each condition
        try
            for exper_cond = exper_conds
                exper_cond_field = replace(exper_cond{:},' ','_');
                exper_data.(exper_cond_field) = struct('pupil_data',struct(),'neural_data',struct(),'trigger_data',struct());
                [exper_data.(exper_cond_field).pupil_data, exper_data.(exper_cond_field).neural_data, exper_data.(exper_cond_field).trigger_data] = loadData(mouse_id, date_str, file_num, exper_cond_field, paths, params);
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


        %%%%%% Prepare for plotting %%%%%%

        % Create experiment title for all plots
        exper_title = {['Mouse ',num2str(mouse_id),', ',date_str,', ',file_num, ', {', strjoin(exper_conds,', '),'}'], ...
            ['rate type=',rate_type,', bin width=',num2str(bin_width),' (s), remove blink=',num2str(remove_blink)]};

        params_folder = strjoin({...
                    [rate_type,'_',num2str(bin_width)],...
                    ['time_(-',num2str(abs(time_before)),',',num2str(time_after),')'],...
                    },'|');

        plot_path  = fullfile(paths.dropbox,'plots',findPlotFolder,'evoked', params_folder, mouse_id, date_str, file_num);
        if ~exist(plot_path,'dir') && save_it; mkdir(plot_path); end



        %% Plot 4 panels

        % ~~~~~~~~~ Get data ~~~~~~~~~

        plot_time = step_plot(exper_data.test.trigger_data.time_bins);


        % ~~~~~~~~~ Plot ~~~~~~~~~
        % Inialize figure
        fig = figure(2);
        clf(fig)
        fig.WindowStyle = fig_params.WindowStyle;
        fig.Units = fig_params.Units;
        fig.Position = fig_params.Position;

        layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

        title(layout, exper_title)


        %%%% Plot firing rate (test vs. extinction) %%%%
        nexttile(layout)

        % Plot
        hold on
        exper_cond = 'tone_alone';
        plot_data = step_plot(mean(exper_data.(exper_cond).trigger_data.neural.data,1,'omitnan'));
        plot_CI = step_plot(calculateCR(exper_data.(exper_cond).trigger_data.neural.data')');
        thresh = max(calculateCI(reshape(exper_data.(exper_cond).trigger_data.neural.data(:,plot_time < 0),1,[])));
        H = shadedErrorBar(plot_time(2:end), plot_data(1:end-1), plot_CI(1:end-1),'lineProps',{'k'});
        H.mainLine.DisplayName = strrep(exper_cond,'_',' ');
        yline(thresh,'k')

        exper_cond = 'test';
        plot_data = step_plot(mean(exper_data.(exper_cond).trigger_data.neural.data,1,'omitnan'));
        plot_CI = step_plot(calculateCR(exper_data.(exper_cond).trigger_data.neural.data')');
        thresh = max(calculateCI(reshape(exper_data.(exper_cond).trigger_data.neural.data(:,plot_time < 0),1,[])));
        H = shadedErrorBar(plot_time(2:end), plot_data(1:end-1), plot_CI(1:end-1),'lineProps',{'r'});
        H.mainLine.DisplayName = strrep(exper_cond,'_',' ');
        keep_legend{2} = H.mainLine;
        yline(thresh,'r')
        hold off
        box on
        xlim([min(plot_time), max(plot_time)])
        %xline(0,'Tag','trigger','DisplayName','Trigger')
        xtickangle(0)
        xlabel('Time since trigger (s)')
        ylabel('Spike rate (Hz)')
        legend('Location','southoutside','Orientation','horizontal')


        % Plot firing rate (retest)
        nexttile(layout)

        hold on
        exper_cond = 'extinction';
        plot_data = step_plot(mean(exper_data.(exper_cond).trigger_data.neural.data,1,'omitnan'));
        plot_CI = step_plot(calculateCR(exper_data.(exper_cond).trigger_data.neural.data')');
        thresh = max(calculateCI(reshape(exper_data.(exper_cond).trigger_data.neural.data(:,plot_time < 0),1,[])));
        H = shadedErrorBar(plot_time(2:end), plot_data(1:end-1), plot_CI(1:end-1),'lineProps',{'k'});
        H.mainLine.DisplayName = strrep(exper_cond,'_',' ');
        yline(thresh,'k')

        exper_cond = 'retest';
        plot_data = step_plot(mean(exper_data.(exper_cond).trigger_data.neural.data,1,'omitnan'));
        plot_CI = step_plot(calculateCR(exper_data.(exper_cond).trigger_data.neural.data')');
        thresh = max(calculateCI(reshape(exper_data.(exper_cond).trigger_data.neural.data(:,plot_time < 0),1,[])));
        H = shadedErrorBar(plot_time(2:end), plot_data(1:end-1), plot_CI(1:end-1),'lineProps',{'r'});
        H.mainLine.DisplayName = strrep(exper_cond,'_',' ');
        keep_legend{2} = H.mainLine;
        yline(thresh,'r')
        hold off
        box on
        xlim([min(plot_time), max(plot_time)])
        %xline(0,'Tag','trigger','DisplayName','Trigger')
        xtickangle(0)
        xlabel('Time since trigger (s)')
        ylabel('Spike rate (Hz)')
        legend('Location','southoutside','Orientation','horizontal')



        % Save
        sameLimits(fig)
        uniformFormat(fig)
        if save_it
            exportgraphics(fig,fullfile(plot_path,'2_evoked_neural.pdf'),'Resolution',300,'ContentType','vector')
        end


        disp(['             Done running ',date_str,', ',file_num])

    end % End loop through conditions

    disp(' ')

end % End loop through mice

disp(' ')