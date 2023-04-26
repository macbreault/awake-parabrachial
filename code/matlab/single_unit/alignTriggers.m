function trigger_frames = alignTriggers(arduino_file, video_path, date_str, vidObj)
% function trigger_frames = alignTriggers(arduino_file, video_path, date_str)
%
% Function that aligns the times of the triggers from the arduino to the video to find the frames during which the trigger was on
%
% Inputs:
%       - arduino_file: string of where arduino file is located to import
%       - video_path: string of path to where video file are located to import
%       - date_str: string of the date used for the experiment
%       - vidObj: video object imported from video_path
%
% Outputs:
%       - trigger_frames: < T x 1 > cell array of the frame numbers of the video taken during trigger T, which should match pupil data
%
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 06-15-2022


%% Check that paths exist

if ~exist(arduino_file,'file')
    warning('no_arduino')
    error('Error in %s: Arduino file called %s does not exist', mfilename, arduino_file)
end

if ~exist(video_path,'dir')
    warning('no_video')
    error('Error in %s: Video folder called %s does not exist', mfilename, video_path)
end



%% 1. Load arudino log
%
arduino_path_split = split(arduino_file,filesep);
arduino_path_datetime_str = erase(strjoin(arduino_path_split(end-1:end)),'.csv');
arduino_path_datetime = datetime(arduino_path_datetime_str,'InputFormat','M-d-yy HHmmssSSSSSSSSS');

arduino_log = readtable(arduino_file); % Loads arduino_log variable

arduino_time_raw = arduino_log.time;
arduino_value = logical(arduino_log.value);

% Arduino time format: %H%M%S%f
%   %H {00, 01, ..., 23}
%   %M {00, 01, ..., 59}
%   %S {00, 01, ..., 59}
%   %f {000000, ..., 999999}

arduino_time_short = floor(arduino_time_raw/ 1000); % Shorten so that there are only 3 milliseconds
                     % datetime('144818323884','InputFormat','HHmmssSSSSSS')
% Convert arduino time (string) into date vector (like in video_log)
arduino_date_time = arrayfun(@(time) [date_str, ' ', num2str(time)], arduino_time_short,'un',0);
%arduino_time = datevec(arduino_date_time,'dd-mmmm-yyyy HHMMSSFFF'); %************* Takes the most time (12 seconds), to get to 1.0e+03 * 2.0220    0.0050    0.0230    0.0110    0.0200    0.0086

%}

%% 2. Load video log

video_path_split = split(video_path,filesep);
video_path_datetime_str = strjoin(video_path_split(end-1:end));
video_path_datetime = datetime(video_path_datetime_str,'InputFormat','M-d-yy HHmmssSSS');

% Load events (video_log)
load(fullfile(video_path,'video_log.mat')); % Loads video_log variable

% % Load video
% vidObj = VideoReader(fullfile(video_path,'video.mp4'));

% Get video time (actual time for each frame)
if strcmp({video_log.Type},'Stop')
    error('STOP and add Stop as a possible better way to get frame time')
end
start_ind = find(strcmp({video_log.Type},'Trigger'),1,'first'); % Used in eventlog
startTime = video_log(start_ind).Data.AbsTime;

frameNum =  1:vidObj.NumFrames;
frameRate = vidObj.FrameRate;
frameTime = (frameNum ./ frameRate) - (1/frameRate); % 30 Frames per second -> SECONDS
warning('off','MATLAB:addtodate:NonIntegerValue')
video_time = cell2mat(arrayfun(@(t) datevec(addtodate(datenum(startTime),t*1000,'millisecond')), frameTime,'un',0)');

% Video time format: %H%M%S%f
%   %H {00, 01, ..., 23}
%   %M {00, 01, ..., 59}
%   %S {00, 01, ..., 59}
%   %f {000000, ..., 999999}



%% 3. Get frame number for each trigger event

% Find trigger start time of each trial (+1)
trigger_start = datevec(arduino_date_time([false;diff(arduino_value) > 0]),'mm-dd-yy HHMMSSFFF');  % Time when trigger turns on using arduino data 
trigger_end = datevec(arduino_date_time([diff(arduino_value) < 0;false]),'mm-dd-yy HHMMSSFFF');    % Time when trigger turns off using arduino data

% Check that there are equal number of trigger starts and ends
if size(trigger_start,1) ~= size(trigger_end,1)
    error('Error in %s: Trigger start and end are not the same size',mfilename)
end


% Remove triggers that are less than 2.5 seconds (5 second is the shortest stimuli)
bad_trigger = etime(trigger_end,trigger_start) <= 2.5;
trigger_start(bad_trigger,:) = [];
trigger_end(bad_trigger,:) = [];


% Get frames for each trigger
trigger_frames = cell(size(trigger_start,1),1);
for t = 1:numel(trigger_frames)

    [~,start_ind] = min(abs(etime(video_time,trigger_start(t,:)))); %min(abs(datenum(video_time - trigger_start(t,:))));
    [~,end_ind] = min(abs(etime(video_time,trigger_end(t,:)))); %min(abs(datenum(video_time - trigger_end(t,:))));

    trigger_frames{t} = frameNum(start_ind:end_ind); % Give me the frame number when the trigger is on

end


%% Run some conditions to check for things that can go wrong during alignTrigger

if isempty(arduino_value)
    warning('empty_arduino')
end

if (numel(trigger_frames) > 1) && isequal(trigger_frames{:})

    if abs(etime(datevec(video_path_datetime), datevec(arduino_path_datetime))) > 300 % More than 5 minutes between video and arduino starting
         warning('check_table')
    elseif etime(datevec(arduino_date_time{1},'mm-dd-yy HHMMSSFFF'),startTime) > 36000 % (more than 1 hour between video log and arduino)
        warning('bad_video_time')
    else
        warning('idk')
    end

end

if isempty(trigger_frames)
    warning('no_triggers')
end





% Check to make sure each trigger is unique
% if isequal(trigger_frames{:})
%     
%     error('Error in %s: Video time was not saved correctly', mfilename)
%     test = dir(fullfile(video_path,'video_log.mat'));
%     video_log(1).Data.AbsTime
%     etime(datevec(datenum(video_log(1).Data.AbsTime)),datevec(test.datenum))
% end


% % % Print example of trigger video (for testing)
% % trigger_num = 6;
% % vidframes = read(vidObj,[trigger_frames{trigger_num}(1), trigger_frames{trigger_num}(end)] );
% % montage(vidframes)





end % end alignTriggers



