% plot_conditions_laser_corr.m
%
% Script that plots STUFF using NEW conditioning experiment
%
% Uses info from video files and arduino for timing of tone
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 08-04-2022
% Modified: 09-15-2022 - Replaced anonymous function of CI with defined CI function, which uses a different method for calculating

clear
clc
set(0,'DefaultFigureWindowStyle','docked') %docks the figures onto the Matlab window


CI = @(x) 1.96 * std(x,[],1,'omitnan') / sqrt(size(x,2)); % Function
step_plot = @(x) reshape([x;x],1,2*numel(x));



%% Set data

% ~~~~~~~~~~~ Set mouse information ~~~~~~~~~~~
mouse_id = 'C3';              % Mouse ID
date_str = '9-6-2022';         % Date of recording
file_num = 'File1';           % File number ('File1')

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
time_after = 10; % Seconds after trigger




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



%% Loop through all conditions

condition_file = fullfile(paths.dropbox,'tables','Conditioning Recordings.xlsx');
mouse_ids = sheetnames(condition_file);
exper_cond = 'Laser'; % What experiment were you recording

if ~run_all
    mouse_ids = mouse_ids(ismember(mouse_ids,mouse_id));
end


% Loop through all mice...
for mouse = 1:numel(mouse_ids)

    mouse_id = mouse_ids{mouse};

    disp(['Running mouse ',mouse_id,'...'])

    %%%%%% Import their table %%%%%%

    condition_table = readtable(condition_file, 'Sheet', mouse_id,'PreserveVariableNames',true); %'VariableNamingRule','preserve')

    % Filter out table to LASER experiments only
    condition_table(~strcmpi(condition_table.Condition, exper_cond),:) = [];

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
        if ~isempty(exper_cond)
            condition_inds = strcmpi(condition_table.Condition, exper_cond);
        end

        % Remove rows from condition_table that do not match conditions
        exper_inds = date_inds & file_inds & condition_inds;
        condition_table(~exper_inds,:) = [];

    end % Only keep conditions requested



    %% Loop through all conditions

    for i = 1:size(condition_table,1)

        % Extract conditions (if running all files in condition table)
        try
            date_str = strjoin(cellfun(@(num) num2str(str2double(num)),{datestr(condition_table(i,:).Date,'mm'),datestr(condition_table(i,:).Date,'dd'),datestr(condition_table(i,:).Date,'yy')},'un',0),'-');
            file_num = condition_table(i,:).File{1};
            exper_cond = condition_table(i,:).Condition{1};
        catch
            warning('Invalid date in excel table. Skip!')
            continue
        end

        if ~strcmpi(exper_cond,'Laser')
            continue
        end


        disp(' ')
        disp(['             Running ',date_str,', ',file_num, ', ',exper_cond,'...'])



        %% Load data

        % Set params
        params = struct('remove_blink',remove_blink,'rate_type',rate_type,'bin_width',bin_width,'norm_it',norm_it,'pupil_type',pupil_type,'neural_type',neural_type,'time_before',time_before,'time_after',time_after);


        try
            [pupil_data, neural_data, trigger_data] = loadData(mouse_id, date_str, file_num, exper_cond, paths, params);
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
        exper_title = ['Mouse ',num2str(mouse_id),', ',date_str,', ',file_num, ', ', exper_cond, ...
            ', rate type = ',rate_type,', bin width = ',num2str(bin_width),' (s), remove blink=',num2str(remove_blink)];

        params_folder = strjoin({...
                    [rate_type,'_',num2str(bin_width)],...
                    ['time_(-',num2str(abs(time_before)),',',num2str(time_after),')'],...
                    },'|');

        plot_path  = fullfile(paths.dropbox,'plots',findPlotFolder,'laser', params_folder, mouse_id, date_str, file_num);
        if ~exist(plot_path,'dir') && save_it; mkdir(plot_path); end




        %% Plot raw data

        run('plot_raw_data.m')



        %% Plot split data

        run('plot_split_data.m')

        

        %% Plot average data

        run('plot_average_data.m')



        %% Plot correlation between pupil size at t=0 (or baseline) and MAGNITUDE firing rate after stimulation

        % Get data (for each trial)

        time_before_laser = 1; % (seconds) before t=0 to average pupil data over

        binned_time = trigger_data.time_bins;
        
        pupil_before = NaN(trigger_data.trigger_total,1);
        rate_after = NaN(trigger_data.trigger_total,1);

        for t = 1:trigger_data.trigger_total

            % Pupil at t=0
            pupil_ind = find((round(-1* time_before_laser,2) <= binned_time) & (binned_time < 0));
            shift = 0;
            while all(isnan(trigger_data.pupil.data(t, pupil_ind)))
                pupil_ind = find((round(-1* time_before_laser,2) <= binned_time) & (binned_time < 0))-shift;
                shift = shift + 1;
            end
            pupil_before(t) = mean(trigger_data.pupil.data(t, pupil_ind),'omitnan');


            % Peak rate after t=0
            neural_ind = 0 < binned_time;
            rate_after(t) = max(trigger_data.neural.data(t, neural_ind));

        end

        [r,p] = corr(pupil_before,rate_after, 'type', 'Pearson','rows','pairwise');


        % Inialize figure
        fig = figure(2);
        clf(fig)
        fig.WindowStyle = fig_params.WindowStyle;
        fig.Units = fig_params.Units;
        fig.Position = fig_params.Position;

        layout = tiledlayout(1,1,'TileSpacing','normal','Padding','normal');

        title(layout, exper_title)


        % Plot
        nexttile(layout)
        plot(pupil_before, rate_after, 'o', 'MarkerSize',15,'MarkerFaceColor','none','MarkerEdgeColor',[0.5 0.5 0.5],'DisplayName',['r=',num2str(r)])
        ls = lsline; % Make least square line
        ls.DisplayName = ['p=',num2str(p)];
        if p <= 0.05
            ls.Color = 'm'; % If significant, then lsline turns red
            ls.LineWidth = 2;
        else
            ls.Color = 'k';
        end

        % Format
        axis square
        box on

        % Label
        ylabel({'Peak spiking rate (Hz)'; ['between t=(0,',num2str(time_after),'] s']})
        if norm_it
            x_label = {'Average pupil area (norm pixel^2)'; ['between t=[-',num2str(time_before_laser),',0] s']};
        else
            x_label = {'Average pupil area (pixel^2)'; ['between t=[-',num2str(time_before_laser),',0] s']};
        end
        xlabel(x_label)
        legend('location','eastoutside')
        title([num2str(sum(~isnan(rate_after) & ~isnan(pupil_before))),' triggers total'])
        labelpoints(pupil_before, rate_after,trigger_data.trigger_num,'C',0,1)



        % Save
        if save_it
            savefig(fig,fullfile(plot_path,'2_magnitude.fig'))
            exportgraphics(fig,fullfile(plot_path,'2_magnitude.pdf'),'Resolution',300,'ContentType','vector')
        end




        %% Plot correlation between pupil size at t=0 (or baseline) and DURATION of firing rate response after stimulation

        % ~~~~~~ Get data (for each trial) ~~~~~~

        time_before_laser = 1; % (seconds) before t=0 to average pupil data over

        pupil_before = NaN(trigger_data.trigger_total,1);
        duration_after = NaN(trigger_data.trigger_total,1);

        for t = 1:trigger_data.trigger_total

            % Pupil at t=0
            pupil_ind = find((round(-1* time_before_laser,2) <= binned_time) & (binned_time < 0));
            shift = 0;
            while all(isnan(trigger_data.pupil.data(t, pupil_ind)))
                pupil_ind = find((round(-1* time_before_laser,2) <= binned_time) & (binned_time < 0))-shift;
                shift = shift + 1;
            end
            pupil_before(t) = mean(trigger_data.pupil.data(t, pupil_ind),'omitnan');



            % Get duration of response after t=0:

            % 1. Find average firing rate and upper bound of confidence interval before t=0
            neural_ind_before = find(binned_time < 0);
            neural_data_before = trigger_data.neural.data(t, neural_ind_before);
            threshold = max(calculateCI(neural_data_before));

            % 2. Any firing rate after t=0 that is above 95% confidence interval before t=0 is SIGNIFICANT?
            neural_ind_after = find(0 < binned_time);
            neural_data_after = trigger_data.neural.data(t, neural_ind_after);

            % 3. Onset must have 2 bins on ...
            try
                onset = threshold <= neural_data_after;
                labeledVector = bwlabel(onset);
                measurements = regionprops(labeledVector, neural_data_after, 'Area', 'PixelValues','PixelIdxList');
                measurements([measurements.Area] < 2) = []; % Remove indices that have less than 2 bins
                min_inds = cellfun(@min,{measurements.PixelIdxList}); % Get max index for each group
                onset_ind = min(measurements(min(min_inds) == min_inds).PixelIdxList);
    
    
                % 4. and offset must have 2 bins off
                offset = threshold > neural_data_after;
                labeledVector = bwlabel(offset);
                measurements = regionprops(labeledVector, neural_data_after, 'Area', 'PixelValues','PixelIdxList');
                measurements([measurements.Area] < 2) = []; % Remove indices that have less than 2 bins
                max_inds = cellfun(@max,{measurements.PixelIdxList}); % Get max index for each group
                offset_ind = min(measurements(max(max_inds) == max_inds).PixelIdxList);
    
                % 5. Record duration
                duration_after(t) = binned_time(neural_ind_after(offset_ind)) - binned_time(neural_ind_after(onset_ind));
            end

        end

        [r,p] = corr(pupil_before, duration_after, 'type', 'Pearson','rows','pairwise');


        % Inialize figure
        fig = figure(3);
        clf(fig)
        fig.WindowStyle = fig_params.WindowStyle;
        fig.Units = fig_params.Units;
        fig.Position = fig_params.Position;

        layout = tiledlayout(1,1,'TileSpacing','normal','Padding','normal');

        title(layout, exper_title)


        % Plot
        nexttile(layout)
        plot(pupil_before, duration_after, 'o', 'MarkerSize',15,'MarkerFaceColor','none','MarkerEdgeColor',[0.5 0.5 0.5],'DisplayName',['r=',num2str(r)])
        ls = lsline; % Make least square line
        ls.DisplayName = ['p=',num2str(p)];
        if p <= 0.05
            ls.Color = 'm'; % If significant, then lsline turns red
            ls.LineWidth = 2;
        else
            ls.Color = 'k';
        end

        % Format
        axis square
        box on

        % Label
        ylabel({'Duration of response (s)'; ['between t=(0,',num2str(time_after),'] s']})
        if norm_it
            x_label = {'Average pupil area (norm pixel^2)'; ['between t=[-',num2str(time_before_laser),',0] s']};
        else
            x_label = {'Average pupil area (pixel^2)'; ['between t=[-',num2str(time_before_laser),',0] s']};
        end
        xlabel(x_label)
        legend('location','eastoutside')
        title([num2str(sum(~isnan(duration_after) & ~isnan(pupil_before))),' triggers total'])
        labelpoints(pupil_before, duration_after,trigger_data.trigger_num,'C',0,1)



        % Save
        if save_it
            savefig(fig,fullfile(plot_path,'3_duration.fig'))
            exportgraphics(fig,fullfile(plot_path,'3_duration.pdf'),'Resolution',300,'ContentType','vector')
        end



        disp(['             Done running ',date_str,', ',file_num, ', ',exper_cond])



    end % End loop through conditions

    disp(' ')

end % End loop through mice

disp(' ')