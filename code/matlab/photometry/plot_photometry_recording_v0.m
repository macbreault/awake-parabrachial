% process_photometry_v0.m
%
% Script to loads and processes photometry data for Jesse
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 09-12-2022 - v0

clear
clc
set(0,'DefaultFigureWindowStyle','docked') %docks the figures onto the Matlab window

warning('TODO: 1. Handle remove triggers, 2. Combine overlapping significant response durations')

%% Set data

% ~~~~~~~~~~~ Set mouse information ~~~~~~~~~~~
mouse_id = 'M3';                     % Mouse ID
date_str = '09192022';               % Date of recording
cond_type = 'CFA';                   % Pain type ('pre','cfa','cci')
cond_day  = 15;                       % Day since pain type (0 for 'pre')
stim_type = 'heat';                  % Type of stimulation applied ('heat' or 'tail')



% Other variables to change
run_all = 1;  % Boolean as to whether to run all conditions (1) or not (0)
save_it = 1; % Whether to save plots (1) or not (0)


min_response_time_on  = 0.5; % seconds
min_response_time_off = 0.5; % seconds





% ~~~~~~~~~~~~~~~~~~~~~~ DO NOT CHANGE BELOW ~~~~~~~~~~~~~~~~~~~~~~
% Figure configurations
fig_params = struct('Position',[0 0 8 5], 'Units','Inches', 'WindowStyle','normal', 'FontSize', 10);




%%%%%% Set paths %%%%%%
paths = struct('dropbox',[],'neural',[]);

if contains(pwd,'/mac/')
    paths.dropbox = '/Users/mac/Dropbox (MIT)/Jesse';
else
    paths.dropbox = '/Users/Jessesmith 1/Dropbox/Jesse';
end

red_driver_path = '/Volumes/T7';
seagate_driver_path = '/Volumes/Seagate Backup Plus Drive';

% Find path to neural data by first checking SEAGATE, then checking RED, then assigning DROPBOX
paths.neural = fullfile(seagate_driver_path, 'Photometry');
if ~exist(paths.neural,'dir') % Check seagate
    paths.neural = replace(paths.neural,seagate_driver_path,red_driver_path);
    if ~exist(paths.neural,'dir')% Check red drive
        paths.neural = fullfile(paths.dropbox,'data','photometry'); % Else, assign dropbox
    end
end


% Check that both paths exist
if ~exist(paths.dropbox,'dir'); error('Error in %s: Dropbox path not found', mfilename); end
if ~exist(paths.neural,'dir'); error('Error in %s: Neural hard drive not found', mfilename); end

addpath(genpath(fullfile(paths.dropbox,'code'))); % Add libraries to directory




%% Loop through all conditions

condition_file = fullfile(paths.dropbox,'Photometry Recordings.xlsx');
mouse_ids = sheetnames(condition_file);

if ~run_all
    mouse_ids = mouse_ids(ismember(mouse_ids,mouse_id));
end

% Loop through all mice...
for mouse = 1:numel(mouse_ids)

    mouse_id = mouse_ids{mouse};

    disp(['Running mouse ',mouse_id,'...'])

    %%%%%% Import their table %%%%%%

    condition_table = readtable(condition_file, 'Sheet', mouse_id,'PreserveVariableNames',true); %'VariableNamingRule','preserve')

    % Filter out the conditions that do not match
    if ~run_all

        % ~~~~~~~ Check which conditions should be run (if not empty) ~~~~~~~
        date_inds = true(size(condition_table,1),1);
        cond_inds = true(size(condition_table,1),1);
        day_inds = true(size(condition_table,1),1);
        stim_ind = true(size(condition_table,1),1);

        % Date
        if ~isempty(date_str)
            datetime_str = datetime(date_str,'InputFormat','MMddyyyy');
            date_inds = condition_table.Date == datetime_str;
        end

        % Condition
        if ~isempty(cond_inds)
            cond_inds = strcmpi(condition_table.Condition,cond_type);
        end

        % Pain day
        if ~isempty(day_inds)
            if ~strcmpi(cond_type,'pre')
                day_inds = condition_table.('Days Post CFA') == cond_day;
            else
                day_inds = isnan(condition_table.('Days Post CFA')) | (condition_table.('Days Post CFA') == 0);
            end
        end

        % Stimulation type
        if ~isempty(stim_ind)
            stim_ind = contains(condition_table.('Stimulation type'), stim_type,'IgnoreCase',1);
        end


        % Remove rows from condition_table that do not match conditions
        exper_inds = date_inds & cond_inds & day_inds & stim_ind;
        condition_table(~exper_inds,:) = [];

    end % Only keep conditions requested



    %% Loop through all conditions

    for i = 1:size(condition_table,1)

        % Extract conditions (if running all files in condition table) to make it match folders/files
        try
            date_str = datestr(condition_table(i,:).Date,'mmddyyyy');
            cond_type = condition_table(i,:).Condition{1};

            if strcmpi(cond_type,'pre')
                cond_day = '0';
            else
                cond_day = num2str(condition_table(i,:).('Days Post CFA'));
            end

            stim_type = strrep(lower(condition_table(i,:).('Stimulation type'){1}),' ','-');
            stim_type(1) = upper(stim_type(1));
        catch
            warning('Invalid date in excel table. Skip!')
            continue
        end

        disp(' ')
        disp(['             Running ',strjoin({date_str,cond_type,cond_day,stim_type},', '),'...'])

        
        %% Load data

        % Load CSV file
        folder_name = strjoin({['#',mouse_id],cond_type,cond_day,stim_type},'-');
        data_path_file = fullfile(paths.neural, date_str, folder_name, 'Data','Tws_1_TrialTraceData.csv');
        if ~exist(data_path_file,'file')
            warning('Could not find file. Skip!')
            continue
        end
        raw_table = readtable(data_path_file); % Load table


        % Pre stimulation time
        pre_time = condition_table(i,:).('Pre Stimulus Time (sec)');


        % Post stimulation time
        post_time = condition_table(i,:).('Post Stimulus Time (sec)');


        % Triggers to ignore
        remove_triggers = condition_table(i,:).('Triggers removed?');
        

        %%%%%% Prepare for plotting %%%%%%

        % Create experiment title for all plots
        fig_title = ['Mouse ',strjoin({mouse_id,date_str,cond_type,cond_day,stim_type},', ')];

        plot_path  = fullfile(paths.dropbox,'plots','photometry',date_str, folder_name);
        if ~exist(plot_path,'dir') && save_it; mkdir(plot_path); end



        %% Extract data

        % Trial information
        total_trials = size(raw_table,1);
        good_trials = true(total_trials,1);
        good_trials(remove_triggers(~isnan(remove_triggers))) = 0;

        % Time information
        total_bins = size(raw_table,2);
        time_bins  = linspace(-1*abs(pre_time), abs(post_time), total_bins);

        % Photometry data
        data = table2array(raw_table); % trials x time bins



        %% Plot trials - SPLIT

        % Prepare figure
        fig = figure(1);
        clf(fig)
        fig.WindowStyle = fig_params.WindowStyle;
        fig.Units = fig_params.Units;
        fig.Position = fig_params.Position;
        colormap(fig,'prism')

        num = numSubplots(sum(good_trials));
        layout = tiledlayout(num(1),num(2),'TileSpacing','compact','Padding','none');
        title(layout, fig_title)

        % Loop through all trials
        for t = 1:total_trials

            if ~isnan(remove_triggers) && any(t == remove_triggers)
                continue
            end

            nexttile(layout)

            % Plot
            hold on
            plot(time_bins, data(t,:))
            xline(0,'DisplayName','Stimulation')
            hold off

            % Format
            axis tight

            % Label
            title(['Trial ',num2str(t)])
            
        end

        xlabel(layout,'Time since stimulation [s]')
        ylabel(layout,'Z-score')

        % ~~~~~~~~~~~ Uniform formatting ~~~~~~~~~~~
        AX = findobj(fig,'Type','Axes');
        y_lim = [min(min(cell2mat(get(AX,'YLim')))), max(max(cell2mat(get(AX,'YLim'))))];
        set(AX,'YLim',y_lim)

        % Save
        if save_it

            % Save figure
            exportgraphics(fig,fullfile(plot_path,'1_trials.pdf'),'Resolution',300,'ContentType','vector')
            savefig(fig,fullfile(plot_path,'1_trials.fig'))

        end



        %% Plot trials - AVERAGED

        % Get data for plotting
        plot_time = time_bins;
        plot_data = mean(data(good_trials,:),1,'omitnan');
        plot_CI =  calculateCR(data(good_trials,:)');

%         test_data = (1:3)'.*repmat(1:10,3,1);
%         plot_time = 1:10;
%         plot_data = mean(test_data,1,'omitnan');
%         plot_CI =  calculateCR(test_data');


        % ~~~ Calculate duration of significant response(s) ~~~

        duration = NaN;

        ind_before = find(time_bins < 0);
        ind_after  = find(time_bins > 0);

        data_before = plot_data(ind_before);
        data_after  = plot_data(ind_after);

        % 1. Use upper bound of CI from baseline as threshold
        thresh = max(calculateCI(data_before));

        % 2. Refine responses, so consecutive bins after threshold above/below threshold to be considered on/off
        min_response_bin_on = round(min_response_time_on / unique(round(diff(time_bins),5))); % Number of bins corresponding to min_response_time
        min_response_bin_off = round(min_response_time_off / unique(round(diff(time_bins),5))); % Number of bins corresponding to min_response_time

        % Onset
        onset = thresh <= data_after;
        labeledVector = bwlabel(onset);
        measurements = regionprops(labeledVector, data_after, 'Area', 'PixelValues','PixelIdxList');
        measurements([measurements.Area] < min_response_bin_on) = []; % Remove indices that are too small
        onset_inds = arrayfun(@(ind) min(measurements(ind).PixelIdxList), 1:numel(measurements));


        % Offset
        offset = thresh > data_after;
        labeledVector = bwlabel(offset);
        measurements = regionprops(labeledVector, data_after, 'Area', 'PixelValues','PixelIdxList');
        measurements([measurements.Area] < min_response_bin_off) = []; % Remove indices that are too small
        offset_inds = arrayfun(@(ind) min(measurements(ind).PixelIdxList), 1:numel(measurements));

        
        % Get duration for every onset...
        duration = zeros(size(onset_inds));
        plot_duration = cell(size(onset_inds));
        for j = 1:numel(onset_inds)
            
            onset_time = plot_time(ind_after(onset_inds(j)));
            
            % Find offset index that matches onset index
            offset_ind = min(offset_inds(onset_inds(j) < offset_inds));
            if isempty(offset_ind)
                offset_ind = numel(ind_after);
            end
            offset_time = plot_time(ind_after(offset_ind));

            duration(j) = offset_time - onset_time;
            plot_duration{j} = ind_after(onset_inds(j)):ind_after(offset_ind);

        end

        % Check that there is no overlap between durations
        overlap = NaN(numel(duration));
        for j1 = 1:numel(duration)
            for j2 = 1:numel(duration)
                overlap(j1,j2) = all(ismember(plot_duration{j1},plot_duration{j2}));
            end
        end
        actual_overlap = triu(overlap,1) | tril(overlap,-1);
        [row,col] = ind2sub(size(actual_overlap), find(actual_overlap));
        for j = 1:numel(row)
            % Keep the index with the largest duration
            if duration(row(j)) > duration(col(j))
                duration(col(j)) = NaN;
                plot_duration{col(j)} = {};
            elseif duration(row(j)) < duration(col(j))
                duration(row(j)) = NaN;
                plot_duration{row(j)} = {};
            end
        end

        % Remove empty
        duration(isnan(duration)) = [];
        plot_duration(cellfun(@isempty,plot_duration)) = [];


%         % Check to make sure durations do not overlap. If so, then combine
%         if numel(plot_duration) > 1
%             % Do something
%             old_plot_duration = plot_duration;
%             plot_duration = {};
%             plot_duration{1} = old_plot_duration{1};
%             new_plot_duration = {};
%             
%             for j = 2:numel(old_plot_duration)
%                 if any(intersect(plot_duration{j-1}, old_plot_duration{j}))
%                     plot_duration{j-1} = unique([plot_duration{j-1}, old_plot_duration{j}]);
%                 end
%             end
%         end



        % Prepare figure
        fig = figure(2);
        clf(fig)
        fig.WindowStyle = fig_params.WindowStyle;
        fig.Units = fig_params.Units;
        fig.Position = fig_params.Position;
        colormap(fig,'prism')

        layout = tiledlayout(1,1,'TileSpacing','normal','Padding','normal');
        title(layout, {fig_title; ['Minimum response time (onset, offset) = (',num2str(min_response_time_on),', ',num2str(min_response_time_off),') s']})

        % Plot
        nexttile(layout)
        H = shadedErrorBar(plot_time, plot_data, plot_CI);
        delete(H.edge)
        hold on
        arrayfun(@(j) plot(plot_time(plot_duration{j}),plot_data(plot_duration{j}),'.', 'DisplayName', ['Duration ',num2str(j),' = ',num2str(duration(j)),' s']), 1:numel(duration))
        xline(0,'DisplayName','Stimulation')
        yline(thresh,'r','DisplayName',['Threshold (',num2str(thresh),')'])
        hold off

        % Format
        axis tight

        % Label
        title(['Trials averaged (',num2str(sum(good_trials)),')'])
        xlabel('Time since stimulation [s]')
        ylabel('Z-score')
        H.mainLine.DisplayName = 'Average';
        legend('location','southoutside','Orientation','horizontal','NumColumns',3)

        
        % Save
        if save_it

            % Save plot
            exportgraphics(fig,fullfile(plot_path,'2_average.pdf'),'Resolution',300,'ContentType','vector')
            savefig(fig,fullfile(plot_path,'2_average.fig'))

            % Save duration into condition table
            warning('Save duration to spreadsheet')
            %condition_table.("Duration of Significatn Response") = 

        end



    end % End loop through all conditions
end % End loop through mice



