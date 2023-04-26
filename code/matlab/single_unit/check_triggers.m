% check_triggers.m
%
% Script that checks whether there are the same number of triggers from video and neural
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 07-06-2022

save_it = 1; % Whether to save trigger problems to condition table (1) or not (0)


%% Initialize parameters

%%%%%% Set paths %%%%%%

if contains(pwd,'/mac/')
    dropbox_path = '/Users/mac/Dropbox (MIT)/Jesse';
else
    dropbox_path = '/Users/Jessesmith 1/Dropbox/Jesse';
end

red_driver_path = '/Volumes/T7';
seagate_driver_path = '/Volumes/Seagate Backup Plus Drive';

% Check which driver is plugged in
if exist(red_driver_path,'dir')
    driver_path = red_driver_path;
elseif exist(seagate_driver_path,'dir')
    driver_path = seagate_driver_path;
end

% Find path to neural data
neural_path = fullfile(driver_path, 'Electrophysiology','Awake Recordings','Conditioning With pupil');
if ~exist(neural_path,'dir')
    neural_path = fullfile(dropbox_path,'data','neural');
end

% Find path to video/arduino data
video_path = fullfile(driver_path,'Jesse','data');
if ~exist(video_path,'dir')
    video_path = fullfile(dropbox_path,'data');
end

% Check that both paths exist
if ~exist(dropbox_path,'dir'); error('Error in %s: Dropbox path not found', mfilename); end
if ~exist(neural_path,'dir'); error('Error in %s: Neural hard drive not found', mfilename); end
if ~exist(video_path,'dir'); error('Error in %s: Video hard drive not found', mfilename); end


addpath(genpath(fullfile(dropbox_path,'code','matlab'))); % Add libraries to directory


%% Loop through all conditions

condition_file = fullfile(dropbox_path,'Conditioning Recordings.xlsx');
mouse_ids = sheetnames(condition_file);

if ~save_it
    mouse_ids = mouse_ids(5);
end

% Initialize thing that is going to keep track of which data is correct
bad_exper =  cell(1,numel(mouse_ids));

% Loop through all mice...
for mouse = 6:numel(mouse_ids)%1:numel(mouse_ids)

    mouse_id = mouse_ids{mouse};

    disp(['Checking mouse ',mouse_id,'...'])

    %%%%%% Import their table %%%%%%

    condition_table = readtable(condition_file, 'Sheet', mouse_id,'PreserveVariableNames',true); %'VariableNamingRule','preserve')

    bad_exper{mouse} = NaN(size(condition_table,1),1);


    % Loop through all conditions
    for i = 1:size(condition_table,1)

        % Extract conditions (if running all files in condition table)
        try
            date_str = strjoin({num2str(month(datestr(condition_table(i,:).Date))), ...
                num2str(day(datestr(condition_table(i,:).Date))), ...
                datestr(condition_table(i,:).Date,'yy')},'-');
            file_num = condition_table(i,:).File{1};
            experiment_condition = condition_table(i,:).Condition{1};
        catch
            warning('Invalid date in excel table. Skip!')
            bad_exper{mouse}(i) = 1;
            continue
        end
        

        disp(' ')
        disp(['             Checking ',date_str,', ',file_num, ', ',experiment_condition,'...'])
        
        

        %% Check data


        %%%%%% Check that mouse_id, date_str, file_num, and experiment_condition exists %%%%%%

        date_inds = condition_table.Date == datetime(date_str,'InputFormat','MM-d-yy');
        file_inds = strcmp(condition_table.File, file_num);
        condition_inds = strcmpi(condition_table.Condition, experiment_condition);

        % Find the row that matches all the information
        exper_ind = date_inds & file_inds & condition_inds;

        % Check to make sure there is only 1 row that matches experiment
        if sum(exper_ind) == 0
            %error('Error in %s: No experiments were found that meet the conditions', mfilename)
            warning('Invalid date in excel table. Skip!')
            bad_exper{mouse}(i) = 1;
            continue
        elseif sum(exper_ind) > 1
            %error('Error in %s: Too many experiments were found that meet the conditions', mfilename)
            warning('Invalid date in excel table. Skip!')
            bad_exper{mouse}(i) = 1;
            continue
        end

        % Extract row that meets conditions
        condition_row = condition_table(exper_ind,:);


        %%%%%% Get file names for video and arduino %%%%%%
        video_time = num2str(condition_table{exper_ind,'Video File'});       % Time for video file
        arduino_time = num2str(condition_table{exper_ind,'Arduino File'});     % Time for arduino file

        if numel(video_time) ~= 9
            video_time = ['0',video_time];
        end

        if numel(arduino_time) ~= 12
            arduino_time = ['0',arduino_time];
        end

        arduino_file = fullfile(video_path,'arduino',date_str, [arduino_time,'.csv']);
        video_folder = fullfile(video_path,'video',date_str, video_time);



        %% Load data

        % Initialize structure array where neural and pupil data will be saved
        pupil = struct();
        neural = struct();


        % ~~~~~~~~~~~ Prepare pupil data ~~~~~~~~~~~

        % Check data
        if ~exist(arduino_file,'file') ||  ~exist(video_folder,'dir')
            warning('Could not find arduino or video file. Skip!')
            bad_exper{mouse}(i) = 2.1;
            continue
        end

        if ~exist(fullfile(video_folder, 'video.mp4'),'file') || ~exist(fullfile(video_folder, 'video_log.mat'),'file')
            warning('MP4 and/or log does not exist. Skip!')
            bad_exper{mouse}(i) = 2.2;
            continue
        end

%         if ~exist(fullfile(video_folder, 'video_proc.mat'),'file')
%             warning('video_proc.mat does not exist. Skip!')
%             bad_exper{mouse}(i) = 2.3;
%             continue
%         end


        % Load data
        vidObj = VideoReader(fullfile(video_folder,'video.mp4')); % VIDEO
        FrameRate = vidObj.FrameRate;
        pupil.time = ((0:(double(vidObj.NumFrames)-1))./ FrameRate)';

        % Align triggers with raw video
        pupil.trigger.frames = alignTriggers(arduino_file, video_folder, date_str, vidObj);

        if strcmp(lastwarn,'bad_video_time')
            warning('Error with times in video_log.mat during alignTriggers. Skip!')
            bad_exper{mouse}(i) = 3;
            continue
        end

        if strcmp(lastwarn,'empty_arduino')
            warning('Arduino did not record. Skip!')
            bad_exper{mouse}(i) = 4.1;
            continue
        end

        if strcmp(lastwarn,'no_triggers')
            warning('No triggers found. Skip!')
            bad_exper{mouse}(i) = 4.2;
            continue
        end

        if strcmp(lastwarn,'check_table')
            warning('Arduino and video times from table do not match. Skip!')
            bad_exper{mouse}(i) = 4.3;
            continue
        end

        if strcmp(lastwarn,'idk')
            warning('Something went wrong in alignTriggers? Skip!')
            bad_exper{mouse}(i) = 4.4;
            continue
        end

        pupil.trigger.num = numel(pupil.trigger.frames);
        pupil.trigger.times = [pupil.time(cellfun(@(trigger) trigger(1), pupil.trigger.frames)), pupil.time(cellfun(@(trigger) trigger(end), pupil.trigger.frames))];




        % ~~~~~~~~~~~ Load spike data ~~~~~~~~~~~

        if exist(fullfile(neural_path,num2str(mouse_id),date_str,'Sorted',[file_num,' ',experiment_condition],'MATLAB','spike.mat'),'file')
            filename_spike = fullfile(neural_path,num2str(mouse_id),date_str,'Sorted',[file_num,' ',experiment_condition],'MATLAB','spike.mat');
        elseif exist(fullfile(neural_path,num2str(mouse_id),date_str,'Sorted',[file_num,' ',lower(experiment_condition)],'MATLAB','spike.mat'),'file')
            filename_spike = fullfile(neural_path,num2str(mouse_id),date_str,'Sorted',[file_num,' ',lower(experiment_condition)],'MATLAB','spike.mat');
        else
            warning('Could not find spike file. Skip!')
            bad_exper{mouse}(i) = 5;
            continue
        end
        neural_struct = load(filename_spike);

        % Extract spike data (regardless of the name of the spike file)
        neural_name = fieldnames(neural_struct);

        if numel(neural_name) == 1
            spike_data = neural_struct.(neural_name{:}); % Time each spike is occuring
        else
            warning('Found more than 1 spike file. Skip!')
            bad_exper{mouse}(i) = 6;
            continue
        end


        % ~~~~~~~~~~~ Load neural trigger time since recording started (in seconds) ~~~~~~~~~~~

        trigger_neural_struct = load(replace(filename_spike,'spike.mat','event.mat'));

        % Extract neural video data (regardless of the name of the event file)
        trigger_neural_name = fieldnames(trigger_neural_struct);

        if numel(trigger_neural_name) == 1
            neural.trigger.time = trigger_neural_struct.(trigger_neural_name{:});
            neural.trigger.num = numel(neural.trigger.time);
        else
            warning('Found more than 1 neural file. Skip!')
            bad_exper{mouse}(i) = 6;
            continue
        end


        % ~~~~~~~~~~~ Remove triggers that were removed manually ~~~~~~~~~~~
        bad_triggers = condition_row{:,'Triggers removed?'};
        if iscell(bad_triggers)
            bad_triggers = str2double(split(bad_triggers{:},','));
        end

        if isnan(bad_triggers)
            bad_triggers = [];
        end

%         if (numel(bad_triggers) ~= 0) && any(~isnan(bad_triggers))
%             error('Idk what to do!')
%         end

        pupil.trigger.frames(bad_triggers) = [];
        pupil.trigger.times(bad_triggers,:) = [];
        pupil.trigger.num = pupil.trigger.num - numel(bad_triggers);



        % ~~~~~~~~~~~  Check that there the same number of triggers for pupil and neural data ~~~~~~~~~~~

        if pupil.trigger.num ~= neural.trigger.num
            warning('Found more neural triggers than video triggers. Skip!')
            if pupil.trigger.num < neural.trigger.num
                bad_exper{mouse}(i) = 8;
            elseif pupil.trigger.num > neural.trigger.num
                bad_exper{mouse}(i) = 7;
            end
            continue
        end


         % ~~~~~~~~~~~  Check if proc file exists ~~~~~~~~~~~
        if ~exist(fullfile(video_folder, 'video_proc.mat'),'file')
            warning('video_proc.mat does not exist. Skip!')
            bad_exper{mouse}(i) = 2.3;
            continue
        else
            proc = load(fullfile(video_folder, 'video_proc.mat'));
            if ~any(strcmp(fieldnames(proc.proc),'blink'))
                warning('Blink does not exist in video_proc.mat. Skip!')
                bad_exper{mouse}(i) = 2.4;
                continue
            end
        end


        bad_exper{mouse}(i) = 0;
        disp(['             Done checking ',date_str,', ',file_num, ', ',experiment_condition])

    end

    %% Save mouse to spreadsheet
    if save_it
        condition_table{:,'Neural Notes'} = bad_exper{mouse};
        writetable(condition_table,condition_file,'Sheet',mouse_id,'WriteVariableNames',true);
    end

end
