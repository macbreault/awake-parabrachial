% load_data.m
%
% Script that loads arduino and video log data for comparison
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 05-10-2022


%% Set location

save_path = '/Users/Jessesmith 1/Dropbox/Jesse';
save_date = '12-May-2022';
save_time = '084803548';



%% Load arudino log

arduino_log = load(fullfile(save_path,'arduino_logs',save_date,[save_time,'.mat'])); % Loads arduino_log variable

arduino_time_raw = arduino_log.arduino_time(2:end);
arduino_value = arduino_log.arduino_value(2:end);

% Convert arduino time (string) into date vector (like in video_log)
arduino_date_time = cellfun(@(time) [save_date, ' ', time], arduino_time_raw,'un',0);
arduino_time = datevec(arduino_date_time,'dd-mmmm-yyyy HHMMSSFFF');



%% Load video log

% Load events (video_log)
load(fullfile(save_path,'video_logs',save_date,save_time,'video_log.mat')); % Loads video_log variable

% Load video
vidObj = VideoReader(fullfile(save_path,'video_logs',save_date,save_time,'video.mp4'));

% Get video time (actual time for each frame)
start_ind = find(strcmp({video_log.Type},'Trigger'),1,'first'); % Used in eventlog

frameNum =  1:vidObj.NumFrames;
frameRate = vidObj.FrameRate;
video_time = cell2mat(arrayfun(@(t) datevec(addtodate(datenum(video_log(start_ind).Data.AbsTime), t * 1000, 'millisecond')), (frameNum ./ frameRate), 'un', 0)');




%% Compare times


% Convert arduino time (string) into date vector (like in video_log)
arduino_date_time = cellfun(@(time) [save_date, ' ', time], arduino_time_raw,'un',0);
arduino_time_vec = datevec(arduino_date_time,'dd-mmmm-yyyy HHMMSSFFF');


% Find closest time in ARDUINO that matches when video STARTED



% Convert frames into seconds
disp(['Video start time: ', datestr(video_time(1,:),'HH:MM:SS:FFF')])
disp(['Arduino start time: ', datestr(arduino_time(1,:),'HH:MM:SS:FFF')])

disp(['Video stop time: ', datestr(video_time(end,:),'HH:MM:SS:FFF')])
disp(['Arduino stop time: ', datestr(arduino_time(end,:),'HH:MM:SS:FFF')])
