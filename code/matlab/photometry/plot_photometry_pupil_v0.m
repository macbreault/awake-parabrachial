% plot_photometry_pupil_v0.m
%
% Script to loads and processes photometry AND pupil data for Jesse
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 02-22-2022 - v0

clear
clc
close all
set(0,'DefaultFigureWindowStyle','docked') %docks the figures onto the Matlab window



%% Set data

% ~~~~~~~~~~~ Set mouse information ~~~~~~~~~~~
mouse_id = 'F4';                     % Mouse ID
date_str = '12152022';               % Date of recording
cond_type = 'cfa';                   % Pain type ('pre','cfa','cci')
cond_day  = 4;                      % Day since pain type (0 for 'pre')
stim_type = 'heat';                  % Type of stimulation applied ('heat' or 'tail')


% Other variables to change
run_all = 0;  % Boolean as to whether to run all conditions (1) or not (0)
save_it = 0; % Whether to save plots (1) or not (0)


% Data to run
types = {'zscore','delta'};

remove_blink = 1; % Remove blink (1) or not (0)
rate_type = 'binned'; % Use 'binned' or 'sliding' window for firing rate
norm_it = 1; % Normalize pupil size by largest size (1) or not (0)
pupil_type = 'raw'; % 'raw' or 'filter'

min_response_time_on  = 0.5; % seconds
min_response_time_off = 0.5; % seconds





% ~~~~~~~~~~~~~~~~~~~~~~ DO NOT CHANGE BELOW ~~~~~~~~~~~~~~~~~~~~~~
% Figure configurations
fig_params = struct('Position',[0 0 8 5], 'Units','Inches', 'WindowStyle','normal', 'FontSize', 10);




%%%%%% Set paths %%%%%%
paths = struct('dropbox',[],'video',[],'neural',[],'code',[],'plots',[]);

if ~nate_it

    if contains(pwd,'/mac/')
        paths.dropbox = '/Users/mac/Dropbox (MIT)/Jesse';
    else
        paths.dropbox = '/Users/Jessesmith 1/Dropbox/Jesse';
    end

    red_driver_path = '/Volumes/T7';
    seagate_driver_path = '/Volumes/Seagate Backup Plus Drive';
    nate_driver_path = 'ADD NATE HERE';


    % Find path to neural data
    paths.neural = fullfile(seagate_driver_path, 'Photometry');
    if ~exist(paths.neural,'dir')% Check red drive
        paths.neural = fullfile(paths.dropbox,'data','photometry'); % Else, assign dropbox
    end

    % Find path to video/arduino data
    paths.video = fullfile(paths.dropbox,'data');
    if ~exist(paths.video,'dir')
        paths.video = fullfile(red_driver_path,'Jesse','data');
    end

    paths.code = fullfile(paths.dropbox,'code');
    paths.plots = fullfile(paths.dropbox,'plots');

else

    % To Nate: Add your paths here!
    paths = struct('dropbox', 'ADD MAIN PATH HERE', ...
                   'video',   'ADD VIDEO PATH HERE', ...
                   'neural',  'ADD NEURAL PATH HERE', ...
                   'code',    'ADD CODE PATH HERE', ...
                   'plot',    'ADD CODE PATH HERE');

end


% Check that both paths exist
if ~exist(paths.dropbox,'dir'); error('Error in %s: Dropbox path not found', mfilename); end
if ~exist(paths.neural,'dir'); error('Error in %s: Neural hard drive not found', mfilename); end
if ~exist(paths.video,'dir'); error('Error in %s: Video hard drive not found', mfilename); end


% Add libraries to directory
addpath(genpath(paths.code));




%% Loop through all conditions

condition_file = fullfile(paths.dropbox,'tables','Photometry Recordings.xlsx');
mouse_ids = sheetnames(condition_file);

if ~run_all
    mouse_ids = mouse_ids(ismember(mouse_ids,mouse_id));
end

% Loop through all mice...
for mouse = 1:numel(mouse_ids)

    mouse_id = mouse_ids{mouse};

    disp(['Running mouse ',mouse_id,'...'])

    %%%%%% Import their table %%%%%%

    opts = detectImportOptions(condition_file, 'Sheet', mouse_id,'PreserveVariableNames',true);
    opts.VariableTypes(strcmp(opts.VariableNames,'Triggers removed?')) = {'char'};
    condition_table = readtable(condition_file, opts);

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

            if strcmp(stim_type,'pinch')
                warning('RENAME!')
            end

        catch
            warning('Invalid date in excel table. Skip!')
            continue
        end

        disp(' ')
        disp(['             Running ',strjoin({date_str,cond_type,cond_day,stim_type},', '),'...'])



        %% Load/process pupil data

        %makeME = @(errID, msgtext, strcell) MException(['MyComponent:',errID],msgtext,strcell{:});

        if any(strcmp(fieldnames(paths),'video'))
            video_time = num2str(condition_table{i,'Video File'});       % Time for video file
            arduino_time = num2str(condition_table{i,'Arduino File'});     % Time for arduino file

            % Convert date to string for video folder
            video_date_str = strjoin(cellfun(@(num) num2str(str2double(num)),{datestr(condition_table(i,:).Date,'mm'),datestr(condition_table(i,:).Date,'dd'),datestr(condition_table(i,:).Date,'yy')},'un',0),'-');

            if numel(video_time) ~= 9
                video_time = ['0',video_time];
            end

            if numel(arduino_time) ~= 12
                arduino_time = ['0',arduino_time];
            end

            arduino_file = fullfile(paths.video,'arduino', video_date_str, [arduino_time,'.csv']);
            video_folder = fullfile(paths.video,'video',   video_date_str, video_time);


            %%%%%% Check that ALL files exist!!! %%%%%%

            % video.mp4
            if ~exist(fullfile(video_folder, 'video.mp4'),'file')
                warning('Could not find video.mp4: %s... SKIP!', fullfile(video_folder, 'video.mp4'))
                continue
            end

            % video_log.mat
            if ~exist(fullfile(video_folder, 'video_log.mat'),'file')
                warning('Could not find video_log.mat: %s... SKIP!', fullfile(video_folder, 'video_log.mat'))
                continue
            end

            % video_proc.m
            if ~exist(fullfile(video_folder, 'video_proc.mat'),'file')
                warning('Could not find video_proc.mat: %s... SKIP!', fullfile(video_folder, 'video_proc.mat'))
                continue
            end

            % arduino_file
            if ~exist(arduino_file,'file')
                warning('Could not find arduino data: %s... SKIP!', arduino_file)
                continue
            end

        end % End checking video data




        %~~~~~~~~~~~ Initialize structure array where neural and pupil data will be saved ~~~~~~~~~~~
        pupil_data = [];


        % Get bad triggers that were manually removed from NEURAL DATA
        bad_triggers = condition_table{i,'Triggers removed?'};
        if iscell(bad_triggers)
            bad_triggers = str2double(split(bad_triggers{:},','));
        end
        if isnan(bad_triggers)
            bad_triggers = [];
        end


        %%%%%% Load pupil data %%%%%%
        if any(strcmp(fieldnames(paths),'video'))

            % ~~~~~~~~~~~ Prepare pupil data ~~~~~~~~~~~
            clearvars proc

            % Load data
            load(fullfile(video_folder, 'video_proc.mat'),'proc') % PUPIL % See https://github.com/MouseLand/facemap/blob/main/docs/svd_matlab_tutorial.md for outputs
            vidObj = VideoReader(fullfile(video_folder,'video.mp4')); % VIDEO

            if proc.nframes ~= vidObj.NumFrames % Check that pupil and video have the same number of frames
                warning('Number of frames does not match between pupil and video. Choosing smaller of the two...')
                numFrames = min([proc.nframes,vidObj.NumFrames]);
            else
                numFrames = vidObj.NumFrames;
            end

            % Get relative time of each frame (starting at time=0)

            FrameRate = vidObj.FrameRate;
            pupil_data.time = ((0:(double(numFrames)-1))./ FrameRate)';

            % Save pupil and blink data
            pupil_data.raw = double(proc.pupil(1).area);
            
            if numel(proc.pupil) ~= 1
                warning('Found more than 1 pupil file... SKIP!')
                continue
            end

            % Add blink
            if remove_blink

                if ~any(strcmp(fieldnames(proc),'blink'))
                    warning('Blink not in proc data')
                end

                pupil_data.blink.raw = double(proc.blink.area);
                pupil_data.blink.thresh = mean(pupil_data.blink.raw,'omitnan') - 3*std(pupil_data.blink.raw,'omitnan');
                pupil_data.blink.ind = (pupil_data.blink.raw < pupil_data.blink.thresh);

                pupil_data.raw(pupil_data.blink.ind) = NaN;

            end

            % Normalize pupil size
            if norm_it
                pupil_data.raw = normalize(pupil_data.raw,'range');
            end

            % Filter pupil data using moving average
            pupil_data.filter = smoothdata(pupil_data.raw,'movmean',100,'omitnan'); % https://www.mathworks.com/help/matlab/ref/smoothdata.html#bvhejau-method


            % Align triggers with pupil video
            try
                pupil_data.trigger.frames = alignTriggers(arduino_file, video_folder, video_date_str, vidObj);
            end

            switch lastwarn
                case 'no_arduino'
                    warning('Could not find arduino file')
                case 'no_video'
                    warning('Could not find video file')
                case 'bad_video_time'
                    warning('Video time broken in log')
                case 'empty_arduino'
                    warning('Arduino did not record')
                case 'no_triggers'
                    warning('No triggers found')
                case 'check_table'
                    error('Error in %s: video and arduino file names in table do not align',mfilename)
            end

            pupil_data.trigger.num = (1:numel(pupil_data.trigger.frames))';
            pupil_data.trigger.total = numel(pupil_data.trigger.frames);
            pupil_data.trigger.times = [pupil_data.time(cellfun(@(trigger) trigger(1), pupil_data.trigger.frames)), pupil_data.time(cellfun(@(trigger) trigger(end), pupil_data.trigger.frames))];


            % ~~~~~~~~~~~ Remove triggers that were removed manually ~~~~~~~~~~~
            pupil_data.trigger.num(bad_triggers) = [];
            pupil_data.trigger.frames(bad_triggers) = [];
            pupil_data.trigger.times(bad_triggers,:) = [];
            pupil_data.trigger.total = pupil_data.trigger.total - numel(bad_triggers);

        end




        %% Load data for plotting Z-Score and Delta versus pupil
        for type = reshape(types,1,[])

            % Pick file name based on type of data to process
            switch type{:}
                case 'zscore'
                    file_name = '_ZScorePETH';
                case 'delta'
                    file_name = '_DeltaFPETH';
            end


            % Load .fig of Z-Score data
            folder_name = strjoin({['#',mouse_id],cond_type,cond_day,stim_type},'-');
            data_path = fullfile(paths.neural, date_str, folder_name, 'Figures');
            listings = dir(data_path);
            file_ind = reshape(find(contains({listings.name}, file_name) & contains({listings.name},'.fig')),1,[]);
            if ~exist(data_path,'dir') || ~any(file_ind)
                warning('Could not find Figure folder or .fig file. Skip!')
                continue
            end



            %% Extract data for every .fig that meets the condition
            for k = file_ind

                % Load figure as invisible
                data_fig = openfig(fullfile(listings(k).folder, listings(k).name), 'invisible');

                % Photometry data (as image object from data_fig)
                img_obj = findobj(data_fig,'Type','Image');

                if numel(img_obj) ~= 1
                    figure(data_fig) % Make figure visible because of error
                    error('Error in %s: Found more than one image plotted in figure that contains data. Check!',mfilename)
                end

                data = img_obj.CData; % trials x time bins

                % Extract y-label from axes that has average plot on it
                lines = findobj(data_fig,'Type','Line');
                y_label = lines(1).Parent.YLabel.String;

                % Trial information
                total_trials = size(data,1);
                good_trials = true(total_trials,1);
                remove_triggers = condition_table(i,:).('Triggers removed?'); % Get triggers to ignore

                if iscell(remove_triggers) && (numel(remove_triggers) == 1)% Is a cell of size 1
                    if (numel(remove_triggers) == 1)
                        remove_triggers = cell2mat(remove_triggers);
                    else
                        error('Error in %s: Triggers to removed are not formatted for matlab. MSB must add new condition to this line', mfilename)
                    end
                end

                if isa(remove_triggers,'char')
                    remove_triggers = cellfun(@str2double, split(remove_triggers,','))';
                end

                if isa(remove_triggers,'double') % If a double
                    good_trials(remove_triggers(total_trials >= remove_triggers)) = 0;
                else
                    error('Error in %s: Triggers to removed are not formatted for matlab. MSB must add new condition to this line', mfilename)
                end

                % Time information
                pre_time = img_obj.Parent.XLim(1); % -5
                post_time = img_obj.Parent.XLim(2); % 10
                total_bins = size(data,2);
                time_bins  = linspace(pre_time, post_time, total_bins);
                bin_width = mode(diff(time_bins));


                



                %%%%%% Prepare for plotting %%%%%%

                % Close figure
                close(data_fig)

                % Create experiment title for all plots
                fig_title = ['Mouse ',strjoin({mouse_id,date_str,cond_type,cond_day,stim_type},', ')];

                extra = erase(cell2mat(regexp(listings(k).name,[file_name,'.*fig'],'match')),{file_name,'.fig'});
                params_folder = strjoin({...
                    ['time_(',num2str(pre_time),',',num2str(post_time),')',extra],...
                    },'|');

                plot_path  = fullfile(paths.plots,findPlotFolder,'pupil',date_str, folder_name, params_folder);
                if ~exist(plot_path,'dir') && save_it; mkdir(plot_path); end





                %% Process video data so it aligns with photometry data

                % Time lock using respective triggers
                if ~isempty(pupil_data)
                    trigger_num = pupil_data.trigger.num;
                    trigger_total = pupil_data.trigger.total;
                    pupil_2_corr = pupil_data.(pupil_type);
                else
                    pupil_2_corr = NaN;
                end

                time_bins_pupil = -1*abs(pre_time) : bin_width : post_time;
                total_bins_pupil = numel(time_bins_pupil); % Check that each trigger has this many bins in it


                pupil2corr = [];
                neural2corr = [];
                pupil2corr_trial = cell(trigger_total,1);
                neural2corr_trial = cell(trigger_total,1);
                pupil_start_ind = NaN(trigger_total,1);
                pupil_trig_ind = NaN(trigger_total,1);
                pupil_stop_ind = NaN(trigger_total,1);
                neural_start_ind = NaN(trigger_total,1);
                neural_trig_ind = NaN(trigger_total,1);
                neural_stop_ind = NaN(trigger_total,1);
                neural_ind = cell(trigger_total,1);
                pupil_ind = cell(trigger_total,1);

                pupil_start_time = NaN(trigger_total,1);
                pupil_trig_time = NaN(trigger_total,1);
                pupil_stop_time = NaN(trigger_total,1);
                neural_start_time = NaN(trigger_total,1);
                neural_trig_time = NaN(trigger_total,1);
                neural_stop_time = NaN(trigger_total,1);
                neural_time = cell(trigger_total,1);
                pupil_time = cell(trigger_total,1);

                pupil_trigger_times = cell(trigger_total,1);
                pupil_trigger_data = cell(trigger_total,1);
                neural_trigger_times = cell(trigger_total,1);
                neural_trigger_data = cell(trigger_total,1);

                trigger_loop = 1:trigger_total; % Triggers to actually loop through


                for t = trigger_loop

                    %%%%%%%% Pupil stuff %%%%%%%%
                    if ~isempty(pupil_data)
                        % Get exact times of events for pupil data
                        pupil_start_time(t) = pupil_data.trigger.times(t,1) - abs(pre_time);
                        pupil_trig_time(t)  = pupil_data.trigger.times(t,1);
                        pupil_stop_time(t)  = pupil_data.trigger.times(t,1) + post_time;

                        % Extract pupil data during trigger
                        pupil_trigger_ind = (pupil_start_time(t) <= pupil_data.time) & (pupil_data.time <= pupil_stop_time(t));
                        pupil_trigger_times{t} = pupil_data.time(pupil_trigger_ind);
                        pupil_trigger_data{t} = pupil_data.raw(pupil_trigger_ind);

                        % Get average pupil size per bin (using similar method as neural data)
                        pupil_binned = NaN(1,total_bins_pupil);

                        switch rate_type
                            case 'binned'

                                % Take median of pupil within bin
                                for bin = 1:length(time_bins_pupil)
                                    pupil_binned(bin) = median(pupil_trigger_data{t}((time_bins_pupil(bin) <= (pupil_trigger_times{t} - pupil_trig_time(t))) & ((pupil_trigger_times{t} - pupil_trig_time(t)) <= time_bins_pupil(bin)+bin_width)),'omitnan');
                                end

                            otherwise
                                error('Error in %s: Did not make method for calculating firing rate for time-locked trigger pupil data',mfilename)
                        end
                    else
                        pupil_binned = NaN(1,total_bins_pupil);
                    end

                    % Extract the data into cell arrays. Trials may not be the same size,
                    % but that will be fixed later
                    pupil2corr_trial{t} = pupil_binned;

                end % Loop through triggers

                pupil = struct(...
                    'type',         pupil_type,...
                    'rate_type',    rate_type,...
                    'raw',          struct(... % times, index, raw data, ect
                    'time',         struct('start',pupil_start_time,'trig',pupil_trig_time,'stop',pupil_stop_time,'time',{pupil_trigger_times}),...
                    'data',         {pupil_trigger_data}), ...
                    'data',         cell2mat(pupil2corr_trial)); % 10 x 16 double

                trigger_data = struct( ...
                    'pupil',         pupil, ...
                    'trigger_num',   trigger_num,...
                    'trigger_total', trigger_total,...
                    'time_before',   abs(pre_time),...
                    'bin_width',     bin_width,...
                    'time_after',    post_time,...
                    'time_bins',     time_bins_pupil,...
                    'total_bins',    total_bins_pupil);


                %% Plot trials - SPLIT

                % Prepare figure
                fig = figure(sub2ind([numel(types),2],1,find(strcmp(type,types))));
                clf(fig)
                fig.WindowStyle = fig_params.WindowStyle;
                fig.Units = fig_params.Units;
                fig.Position = fig_params.Position;

                P = numSubplots(sum(good_trials));
                layout = tiledlayout(P(1),P(2),'TileSpacing','compact','Padding','none');
                title(layout, fig_title)

                % Loop through all trials
                for t = 1:total_trials

                    if any(~isnan(remove_triggers) & (t == remove_triggers))
                        continue
                    end

                    nexttile(layout)
                    
                    trigger_trial_ind = trigger_data.trigger_num == t;
                    [y,x] = ind2sub(flip(P),find(trigger_trial_ind));

                    % Plot and format LEFT
                    yyaxis left
                    hold on
                    plot(time_bins, data(t,:),'Color','g','DisplayName',['Binned ', y_label],'LineWidth',1.5)
                    xline(0,'Tag','trigger')
                    hold off
                    axis tight
                    set(gca,'YColor','g')
                    if y == 1
                        ylabel(y_label)
                    end
                    title(['Trial ',num2str(t)])


                    % Plot and format RIGHT
                    try
                        pupil_time_bins = linspace(-1*trigger_data.time_before, trigger_data.time_after, numel(trigger_data.pupil.raw.data{trigger_trial_ind}));
                    catch
                        continue
                    end
                    yyaxis right
                    hold on
                    plot(pupil_time_bins, trigger_data.pupil.raw.data{trigger_trial_ind},'Color','b','DisplayName','Binned pupil','LineWidth', 1.5)
                    xline(0,'Tag','trigger')
                    hold off
                    axis tight
                    set(gca,'YColor','b')
                    if y == P(2)
                        if norm_it
                            ylabel({'Normalized','pupil area'})
                        else
                            ylabel({'Pupil area'})
                        end
                    end

                    % Format
                    axis tight

                    % Label
                    title(['Trial ',num2str(t)])



                end

                xlabel(layout,'Time since stimulation [s]')

                % ~~~~~~~~~~~ Uniform formatting ~~~~~~~~~~~
                %sameLimits(fig,'y')
                AX = findobj(fig,'Type','Axes');
                for a = 1:numel(AX)
                    yyaxis(AX(a),'left')
                end
                uniformFormat(fig)
                AX = findobj(fig,'Type','Axes');
                for a = 1:numel(AX)
                    yyaxis(AX(a),'right')
                end
                uniformFormat(fig)

                % Save
                if save_it

                    % Save figure
                    exportgraphics(fig,fullfile(plot_path,[type{:},'_1_trials.pdf']),'Resolution',300,'ContentType','vector')
                    savefig(fig,fullfile(plot_path,[type{:},'_1_trials.fig']))

                end



                %% Plot trials - AVERAGED

                % Get data for plotting
                plot_time = time_bins;
                plot_data = mean(data(good_trials,:),1,'omitnan');
                plot_CI =  calculateCR(data(good_trials,:)');


                % ~~~ Calculate duration of significant response(s) ~~~

                duration = NaN;

                ind_before = find(time_bins < 0);
                ind_after  = find(time_bins > 0);

                data_before = plot_data(ind_before);
                data_after  = plot_data(ind_after);

                % 1. Use upper bound of CI from baseline as threshold
                thresh_both = calculateCI(data_before);
                thresh = max(thresh_both);

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



                % ~~~ Calculate AREA under significant response(s) ~~~
                AUC = NaN(size(duration));
                plot_AUC = cell(size(duration));

                for j = 1:numel(AUC)

                    % Calculate AUC between these curves:
                    x  = plot_time(plot_duration{j});
                    y1 = plot_data(plot_duration{j}); % Response
                    y2 = repmat(thresh ,size(x));     % Threshold

                    %                     % Test cases
                    %                     thresh = 0.5;
                    %                     x = [0,1];
                    %                     y1 = [0,1];
                    %                     y2 = repmat(thresh,size(x));

                    % 1. Find all points where these curves intersect (not including ends)
                    [x0,y0] = intersections(x,y1,x,y2);

                    % 2. Add end points
                    x_intersect = unique([x(1); x0; x(end)]);

                    % 3. Find y1 value for each intersect
                    x_with_intersects =  [x, x0'];
                    y1_with_intersects = [y1,y0'];
                    [x_with_intersects_sorted, order] = sort(x_with_intersects);
                    y1_with_intersects_sorted = y1_with_intersects(order);

                    % 4. Integrate between each point of intersection
                    AUC_temp = NaN(numel(x_intersect)-1,1);
                    for jj = 1:(numel(x_intersect)-1)
                        mask = (x_intersect(jj) <= x_with_intersects_sorted) & (x_with_intersects_sorted <= x_intersect(jj+1));
                        AUC_temp(jj) = trapz(x_with_intersects_sorted(mask), abs(y1_with_intersects_sorted(mask)-thresh));
                    end

                    % 5. Sum all areas together
                    AUC(j) = sum(AUC_temp); % trapz(x,abs(y1-thresh)) % Can't use this one because of areas that go below threshold

                    % Make mask of AUC for plotting
                    plot_AUC{j}(1,:) = [x, fliplr(x)]; % x-value
                    plot_AUC{j}(2,:) = [y1, fliplr(y2)]; % y-value

                end

                
                % ~~~ Repeat for PUPIL ~~~
                pupil_plot_time = linspace(-1*trigger_data.time_before, trigger_data.time_after, mode(cellfun(@numel, trigger_data.pupil.raw.data)));
                pupil_data_mat = NaN(trigger_data.trigger_total, numel(pupil_plot_time));
                for t = 1:size(pupil_data_mat,1)
                    pupil_data_mat(t,1:numel(trigger_data.pupil.raw.data{t})) = trigger_data.pupil.raw.data{t};
                end
                pupil_plot_mean = mean(pupil_data_mat,1,'omitnan');
                pupil_plot_CI = calculateCR(pupil_data_mat');





                % ~~~ Figure ~~~

                % Prepare figure
                fig = figure(sub2ind([numel(types),2],2,find(strcmp(type,types))));
                clf(fig)
                fig.WindowStyle = fig_params.WindowStyle;
                fig.Units = fig_params.Units;
                fig.Position = fig_params.Position;

                sig_color = [0,0,0];

                layout = tiledlayout(1,1,'TileSpacing','normal','Padding','normal');
                title(layout, {fig_title; ['Minimum response time (onset, offset) = (',num2str(min_response_time_on),', ',num2str(min_response_time_off),') s']})

                %%%%%%%% Plot photometry %%%%%%%%
                nexttile(layout)
                yyaxis left
                H = shadedErrorBar(plot_time, plot_data, plot_CI);
                delete(H.edge)
                H.mainLine.Color = 'g';
                H.patch.FaceColor = 'g';
                hold on
                arrayfun(@(j) fill(plot_AUC{j}(1,:), plot_AUC{j}(2,:), sig_color,'EdgeColor','none','FaceAlpha',0.25,'DisplayName', ['AUC ',num2str(j),' = ',num2str(AUC(j))]), 1:numel(AUC))
                arrayfun(@(j) plot(plot_time(plot_duration{j}), plot_data(plot_duration{j}),'.-','Color',sig_color,'DisplayName', ['Duration ',num2str(j),' = ',num2str(duration(j)),' s']), 1:numel(duration))
                xline(0,'Tag','trigger')
                arrayfun(@(y) yline(y,'Tag','thresh'), thresh_both)
                hold off

                set(gca,'YColor','g')
                axis tight

                title(['Trials averaged (',num2str(sum(good_trials)),')'])
                xlabel('Time since stimulation [s]')
                ylabel(y_label)
                H.mainLine.DisplayName = 'Average';
                H.mainLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
                if ~isempty(AUC)
                    legend('location','southoutside','Orientation','horizontal','NumColumns',numel(AUC))
                end


                %%%%%%%% Plot PUPIL %%%%%%%%
                yyaxis right
                H = shadedErrorBar(pupil_plot_time, pupil_plot_mean, pupil_plot_CI);
                delete(H.edge)
                H.mainLine.Color = 'b';
                H.mainLine.DisplayName = 'Pupil';
                H.patch.FaceColor = 'b';
                set(gca,'YColor','b')
                axis tight

                ylabel(y_label)

                if norm_it
                    ylabel({'Normalized','pupil area'})
                else
                    ylabel({'Pupil area'})
                end



                % ~~~~~~~~~~~ Uniform formatting ~~~~~~~~~~~
                yyaxis left
                uniformFormat(fig)
                yyaxis right
                uniformFormat(fig)

                % Save
                if save_it

                    % Save plot
                    exportgraphics(fig,fullfile(plot_path,[type{:},'_2_average.pdf']),'Resolution',300,'ContentType','vector')
                    savefig(fig,fullfile(plot_path,[type{:},'_2_average.fig']))

                end

            end % Loop through files


        end % Loop through plot types


        disp(['             Done running ',strjoin({date_str,cond_type,cond_day,stim_type},', '),'.'])



    end % End loop through all conditions

    if size(condition_table,1) == 0
        disp(['No data found for mouse ',mouse_id,' and/or conditions. Skip!'])
    else
        disp(['Done running mouse ',mouse_id,'.'])
    end

    disp(' ')

end % End loop through mice



