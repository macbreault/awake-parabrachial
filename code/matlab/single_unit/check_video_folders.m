% check_video_folders.m
%
% Script that checks the timing of video files to make sure video_log lines up with file name
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 07-02-2022

error('OUTDATED')


%% Initalize folders
% ~~~~~~~~~~~~~~~~~~~~~~ DO NOT CHANGE BELOW ~~~~~~~~~~~~~~~~~~~~~~

%%%%%% Set paths %%%%%%

if contains(pwd,'/mac/')
    dropbox_path = '/Users/mac/Dropbox (MIT)/Jesse';
    drive_path = fullfile(dropbox_path, 'Neural data');
    %drive_path = '/Volumes/Seagate Backup Plus Drive/Electrophysiology/Awake Recordings/Conditioning With pupil';
    plot_path  = fullfile(dropbox_path,'matlab_code/plots', mouse_id, date_str, file_num, experiment_condition);
else
    dropbox_path = '/Users/Jessesmith 1/Dropbox/Jesse';
    drive_path = '/Volumes/Seagate Backup Plus Drive/Electrophysiology/Awake Recordings/Conditioning With pupil';
    plot_path = fullfile('/Users/Jessesmith 1/Dropbox/Jesse/matlab_code/plots', mouse_id, date_str, file_num, experiment_condition);
end

% Check that both paths exist
if ~exist(dropbox_path,'dir'); error('Error in %s: Dropbox path not found', mfilename); end
if ~exist(drive_path,'dir'); warning('Warning in %s: Hardrive path not found', mfilename); end


addpath(genpath(erase(mfilename('fullpath'),mfilename))); % Add libraries to directory

%%%%%% Import excel %%%%%%
excel_file = fullfile(dropbox_path,'Conditioning Recordings.xlsx');
[~,mouse_ids] = xlsfinfo(excel_file);

time_elapse = cell(numel(mouse_ids),1);
time_conds  = cell(numel(mouse_ids),1);


%% Loop through mice
for n = 1:numel(mouse_ids)

    mouse_id = mouse_ids{n};
    condition_table = readtable(excel_file, 'Sheet', mouse_id,'PreserveVariableNames',true); %'VariableNamingRule','preserve')

    % Initialize vector to keep time
    time_elapse{n} = NaN(size(condition_table,1),1);
    time_conds{n}  = cell(size(condition_table,1),1);


    %% Loop through data
    for i = 1:size(condition_table,1)

        try

            date_str = datestr(condition_table.Date(i));                % Date of recording
            file_num = condition_table.File{i};                         % File number
            experiment_condition = condition_table.Condition{i};        % What experiment were you recording


            %%%%%% Convert dates to different versions that are used in files %%%%%%
            dates = struct('video','','arduino','','neural','');

            % Video
            dates.video = split(date_str,'-');
            dates.video{2} = dates.video{2}(1:3);
            dates.video = strjoin(dates.video,'-');

            % Arduino
            dates.arduino = date_str;

            % Neural
            [y,m,d] = ymd(datetime(date_str));
            dates.neural = [num2str(m),'-',num2str(d),'-',num2str(y-2000)]; % 6-8-22


            %%%%%% Check that mouse_id, date_str, file_num, and experiment_condition exists %%%%%%

            date_rows = condition_table.Date == datetime(date_str);
            file_rows = strcmp(condition_table.File, file_num);
            condition_rows = strcmpi(condition_table.Condition, experiment_condition);

            % Find the row that matches all the information
            exper_row = date_rows & file_rows & condition_rows;

            % Check to make sure there is only 1 row that matches experiment
            if sum(exper_row) == 0
                error('Error in %s: No experiments were found that meet the conditions', mfilename)
            elseif sum(exper_row) > 1
                error('Error in %s: Too many experiments were found that meet the conditions', mfilename)
            end


            %%%%%% Get file names for video and arduino %%%%%%
            video_time = num2str(condition_table{exper_row,'Video File'});       % Time for video file
            arduino_time = num2str(condition_table{exper_row,'Arduino File'});     % Time for arduino file

            if isnan(str2double(video_time))
                continue
            elseif numel(video_time) ~= 9
                video_time = ['0',video_time];
            end

            if isnan(str2double(arduino_time))
                continue
            elseif numel(arduino_time) ~= 12
                arduino_time = ['0',arduino_time];
            end

            arduino_file = fullfile(dropbox_path,'arduino_logs',dates.arduino, [arduino_time,'.csv']);
            video_path = fullfile(dropbox_path,'video_logs',dates.video, video_time);



            %% Load video_log

            load(fullfile(video_path,'video_log.mat')); % Loads video_log variable
            start_ind = find(strcmp({video_log.Type},'Trigger'),1,'first'); % Used in eventlog


            % Calculate difference in time from  Trigger (from video_log) and file.
            startTime_log = video_log(start_ind).Data.AbsTime;
            startTime_file = datevec([dates.video, ' ',video_time],'dd-mmm-yyyy HHMMSSFFF');

            time_elapse{n}(i) = etime(startTime_log, startTime_file);
            time_conds{n}{i} = condition_table(i,:);
            time_conds{n}{i}.File = append(mouse_id,'_', time_conds{n}{i}.File);

            if (time_elapse{n}(i) < 0) && ~strcmpi(experiment_condition,'laser')
                disp('hi')
            end


        catch
            continue
        end


    end % Loop through files

end % Loop through mice



%% Plot results

time_all = vertcat(time_elapse{:});
conds_all = vertcat(time_conds{:});



% Good etimes must be:
bad_ind = false(size(time_all));

% 1) Positive
bad_ind(time_all < 0) = 1;

% 2) Not NaN
bad_ind(isnan(time_all)) = 1;

% 3) Not more than 300 seconds?
bad_ind(abs(time_all) > 300) = 1;


figure(100)
clf
histogram(time_all(~bad_ind))
xlabel('Elapsed time (s)')
ylabel('Count')

disp(['Average time between folder created and video_trigger = ',num2str(mean(time_all(~bad_ind))),' seconds'])



conds_bad_table = conds_all(time_all < 0);
conds_bad_cell = cellfun(@table2cell, conds_bad_table,'un',0);
conds_bad = vertcat(conds_bad_cell{:});

