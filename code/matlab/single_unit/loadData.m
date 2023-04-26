function [pupil_data, neural_data, trigger_data] = loadData(mouse_id, date_str, file_num, exper_cond_field, paths, params)
% function [pupil_data, neural_data, trigger_data] = loadData_v1(mouse_id, date_str, file_num, strrep(exper_cond,'_',' '), paths, params)
%
% Script loads neural and pupil data for condition experiments
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 08-04-2022
% Modified: 08-10-2022 - v1 - Now uses firing rate found using smoothing window instead of binned count
% Modified: 08-18-2022 - v2 - Now finds time-locked firing rate seperate from neural_data.raw so time is more accurate at time 0 + keep
% spike data for trigger_data

% Make custom error message that will trigger termination of script ('exit') OR skip ('ignore') depending on errorID
makeME = @(errID, msgtext, strcell) MException(['MyComponent:',errID],msgtext,strcell{:});

disp(['Loading data for ',strjoin({mouse_id, date_str, file_num, exper_cond_field},', '),'...'])


%% Initialize

% Load condition table
condition_table = readtable(fullfile(paths.dropbox,'tables','Conditioning Recordings.xlsx'), 'Sheet', mouse_id,'PreserveVariableNames',true); %'VariableNamingRule','preserve')


%%%%%% Check that mouse_id, date_str, file_num, and experiment_condition exists %%%%%%

exper_cond = strrep(exper_cond_field,'_',' '); % Change formatting for uniform use

date_inds = condition_table.Date == datetime(date_str,'InputFormat','M-d-yy');
file_inds = strcmp(condition_table.File, file_num);
condition_inds = strcmpi(condition_table.Condition, exper_cond);

% Find the row that matches all the information
exper_ind = date_inds & file_inds & condition_inds;

% Check to make sure there is only 1 row that matches experiment
if sum(exper_ind) ~= 1
    throw(makeME('ignore','Error in %s: Problem finding experiment in Conditioning Recordings.xlsx', {mfilename}))
end

% Extract row that meets conditions
condition_row = condition_table(exper_ind,:);


%% Make and check path to video data

if any(strcmp(fieldnames(paths),'video'))
    video_time = num2str(condition_row{:,'Video File'});       % Time for video file
    arduino_time = num2str(condition_row{:,'Arduino File'});     % Time for arduino file

    if numel(video_time) ~= 9
        video_time = ['0',video_time];
    end

    if numel(arduino_time) ~= 12
        arduino_time = ['0',arduino_time];
    end

    arduino_file = fullfile(paths.video,'arduino',date_str, [arduino_time,'.csv']);
    video_folder = fullfile(paths.video,'video',date_str, video_time);


    %%%%%% Check that ALL files exist!!! %%%%%%

    % video.mp4
    if ~exist(fullfile(video_folder, 'video.mp4'),'file')
        throw(makeME('ignore','Error in %s: Could not find video.mp4: %s', {mfilename, fullfile(video_folder, 'video.mp4')}))
    end

    % video_log.mat
    if ~exist(fullfile(video_folder, 'video_log.mat'),'file')
        throw(makeME('ignore','Error in %s: Could not find video_log.mat: %s', {mfilename, fullfile(video_folder, 'video_log.mat')}))
    end

    % video_proc.m
    if ~exist(fullfile(video_folder, 'video_proc.mat'),'file')
        throw(makeME('ignore','Error in %s: Could not find video_proc.mat: %s', {mfilename, fullfile(video_folder, 'video_proc.mat')}))
    end

    % arduino_file
    if ~exist(arduino_file,'file')
        throw(makeME('ignore','Error in %s: Could not find arduino data: %s', {mfilename, arduino_file}))
    end

end % End checking video data



%% Make and check path to neural data

if any(strcmp(fieldnames(paths),'neural'))

    % Check that "Sorted" exists
    neural_folder = fullfile(paths.neural,num2str(mouse_id),date_str,'Sorted');

    if ~exist(neural_folder,'dir')
        throw(makeME('ignore','Error in %s: Could not find path to folders of neural data: %s', {mfilename, neural_folder}))
    end

    % List experimental folders under "Sorted" to find the one that matches experimental condition
    exper_folder_name = [file_num,' ',exper_cond];
    listings = dir(neural_folder);
    
    exper_folder_ind = strcmpi({listings.name},exper_folder_name);
    if ~any(exper_folder_ind)
        throw(makeME('ignore','Error in %s: Could not find experimental folders of neural data matching: %s',{mfilename, exper_folder_name}))
    elseif sum(exper_folder_ind) > 1
        throw(makeME('ignore','Error in %s: Found more than 1 experimental folder of neural data matching: %s',{mfilename, exper_folder_name}))
    end

    neural_file_path = fullfile(listings(exper_folder_ind).folder, listings(exper_folder_ind).name,'MATLAB');

    % Do final check to make sure MATLAB exists under experimental condition folder
    if ~exist(neural_file_path,'dir')
        throw(makeME('ignore','Error in %s: Could not find path to neural data: %s',{mfilename, neural_file_path}))
    end


    %%%%%% Check that ALL files exist!!! %%%%%%

    % Neural - spike.mat
    if ~exist(fullfile(neural_file_path,'spike.mat'),'file')
        throw(makeME('ignore','Error in %s: Could not find spike.mat: %s', {mfilename, fullfile(neural_file_path,'spike.mat')}))
    end

    % Neural - event.mat
    if ~exist(fullfile(neural_file_path,'event.mat'),'file')
        throw(makeME('ignore','Error in %s: Could not find event.mat: %s', {mfilename, fullfile(neural_file_path,'event.mat')}))
    end

end % End checking neural data



%% Initialize variables before loading data

% Initialize structure array where neural and pupil data will be saved
pupil_data = [];
neural_data = [];


% Get bad triggers that were manually removed from NEURAL DATA
bad_triggers = condition_row{:,'Triggers removed?'};
if iscell(bad_triggers)
    bad_triggers = str2double(split(bad_triggers{:},','));
end
if isnan(bad_triggers)
    bad_triggers = [];
end



%% Load pupil data

if any(strcmp(fieldnames(paths),'video'))

    % ~~~~~~~~~~~ Prepare pupil data ~~~~~~~~~~~

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

    % Check that bin width of neural data is big enough so pupil_data is not repeated
    if (1/FrameRate) > params.bin_width
        throw(makeME('ignore','Error in %s: Bin width is to small. Must be larger than 1 / frame rate!',{mfilename}))
    end

    % Save pupil and blink data
    pupil_data.raw = double(proc.pupil(1).area);

    % Add blink
    if params.remove_blink

        if ~any(strcmp(fieldnames(proc),'blink'))
            throw(makeME('ignore','Error in %s: Blink not in proc data',{mfilename}))
        end

        pupil_data.blink.raw = double(proc.blink.area);
        pupil_data.blink.thresh = mean(pupil_data.blink.raw,'omitnan') - 3*std(pupil_data.blink.raw,'omitnan');
        pupil_data.blink.ind = (pupil_data.blink.raw < pupil_data.blink.thresh);

        pupil_data.raw(pupil_data.blink.ind) = NaN;

    end

    % Normalize pupil size
    if params.norm_it
        pupil_data.raw = normalize(pupil_data.raw,'range');
    end

    % Filter pupil data using moving average
    pupil_data.filter = smoothdata(pupil_data.raw,'movmean',100,'omitnan'); % https://www.mathworks.com/help/matlab/ref/smoothdata.html#bvhejau-method


    % Align triggers with pupil video
    try
        pupil_data.trigger.frames = alignTriggers(arduino_file, video_folder, date_str, vidObj);
    end

    switch lastwarn
        case 'no_arduino'
            throw(makeME('ignore','Error in %s: Could not find arduino file',{mfilename}))
        case 'no_video'
            throw(makeME('ignore','Error in %s: Could not find video file',{mfilename}))
        case 'bad_video_time'
            throw(makeME('ignore','Error in %s: Video time broken in log',{mfilename}))
        case 'empty_arduino'
            throw(makeME('ignore','Error in %s: Arduino did not record',{mfilename}))
        case 'no_triggers'
            throw(makeME('ignore','Error in %s: No triggers found',{mfilename}))
        case 'check_table'
            throw(makeME('exit','Error in %s: video and arduino file names in table do not align',{mfilename}))
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



%% Load neural data

if any(strcmp(fieldnames(paths),'neural'))
    % ~~~~~~~~~~~ Load spike data ~~~~~~~~~~~

    neural_struct = load(fullfile(neural_file_path,'spike.mat'));

    % Extract spike data (regardless of the name of the spike file)
    neural_name = fieldnames(neural_struct);

    if numel(neural_name) == 1
        spike_data = neural_struct.(neural_name{:}); % Time each spike is occuring
    else
        throw(makeME('ignore','Error in %s: Found more than 1 spike file',{mfilename}))
    end

    % Round based on smpRate
    smpRate = 2.5e-05; % 40kHz -> 40000 Hz -> 1/40000 Rate
    spike_data = round(spike_data,6);
    neural_data.spike = spike_data;

    % % Re-create spike train from paper, and see if methods mimic their results.
    % smpRate = 0.0025;
    % spike_data = round([0+[1,2.5,3.5,8.5]/100,0.1+[1.5,3,4,5,6.5,7,7.5,8,9.5]/100,0.2+[0.5,1.5,2,9]/100,0.3+[1,3.5]/100,0.4+[6.5,8]/100,0.5+[0,1.5,4,4.5,5,6,7.5,9.5]/100,...
    %                            0.6+[0,1,2,6.5,9]/100,0.7+[6,6.5,9]/100,0.8+[3.5,5.5,6,8.5,9]/100,0.9+[0.5,2,3,7,9.5]/100,1+4.5/100,...
    %                            1.2+[6,9.5]/100,1.6+[5]/100,1.7+[2,8]/100,1.8+[1,4.5,6,7.5]/100,1.9+[3,4,4.5,7.5]/100,...
    %                            2+[1.5,3.5,6,8,9]/100,2.1+[0,4,7,8.5]/100,2.2+[1.5,3.5,5.5,6.5,8]/100,2.3+[1.5,5,7,7.5,9]/100,2.4+[0,0.5,3.5]/100,...
    %                            2.6+[1]/100,2.7+[8,9.5]/100,2.8+[2,3,4,6]/100,2.9+[0.5,3,5,7,7.5,8,9.5,10]/100],6);


    % Calculate firing rate (and assigns neural_data.time and neural_data.raw)
    switch params.rate_type
        case 'binned'

            % Use binned count
            time_end = ceil(spike_data(end)); % Time of last spike
            neural_data.time = (0:params.bin_width:time_end)';
            neural_data.raw = NaN(size(neural_data.time)); % Initalize vector that keeps track of how many spikes are in each 1 second bin
            for t = 1:length(neural_data.time)
                neural_data.raw(t) = sum(  (neural_data.time(t) <= spike_data) & (spike_data < (neural_data.time(t)+params.bin_width)) ) / params.bin_width; % Convert to per second (aka Hz)
            end


        case 'sliding'

            % Use sliding window
            neural_data.time = round((0:smpRate:spike_data(end))',6);

            % Match most of the spike_data
            spike_train = ismember(neural_data.time,spike_data);

            if numel(spike_data) ~= sum(spike_train)

                throw(makeME('exit','Error in %s: Did not match some spikes',{mfilename}))

                % Find matches for the remaining spike_data
                not_matched = spike_data(~ismember(spike_data,neural_data.time));

            end


            % Take sliding window
            k = rectwin(find(neural_data.time == params.bin_width)-1); % Put bin width here!
            firing_rate = conv(spike_train,k);
            firing_rate = firing_rate(ceil(length(k)/2):end-floor(length(k)/2));

            neural_data.raw = firing_rate ./ params.bin_width;

        otherwise
            throw(makeME('exit','Error in %s: Method for calculating firing rate was not recognized',{mfilename}))
    end



    % Check that neural data is not missing data (for more than 2 minutes)
    min_time = 120; % In seconds = 2 minutes
    if max(diff(neural_data.time(neural_data.raw ~= 0))) >= min_time % In seconds
        warning('Missing large chunk of time in neural data. Check!')
    end


    % Filter spike data
    neural_data.filter = smoothdata(neural_data.raw,'movmean',5);




    % ~~~~~~~~~~~ Load neural trigger time since recording started (in seconds) ~~~~~~~~~~~

    trigger_neural_struct = load(fullfile(neural_file_path,'event.mat'));

    % Extract neural video data (regardless of the name of the event file)
    trigger_neural_name = fieldnames(trigger_neural_struct);
    
    neural_data.trigger = struct('num',[],'time',[],'total',[]);
    if numel(trigger_neural_name) == 1
        neural_data.trigger.time = trigger_neural_struct.(trigger_neural_name{:});
        neural_data.trigger.total = numel(neural_data.trigger.time);
        neural_data.trigger.num = (1:neural_data.trigger.total)';
    else
        throw(makeME('ignore','Error in %s: Found more than 1 laser neural file',{mfilename}))
    end

    % ~~~~~~~~~~~ Remove triggers that were removed manually ~~~~~~~~~~~
    if ~isempty(bad_triggers)
        if ~isempty(pupil_data)
            neural_data.trigger.num = pupil_data.trigger.num;
            neural_data.trigger.total = length(pupil_data.trigger.num);
            if numel(neural_data.trigger.time) ~= neural_data.trigger.total 
                neural_data.trigger.time = neural_data.trigger.time(pupil_data.trigger.num);
            end
%         else
%             neural_data.trigger.num(bad_triggers) = [];
        end
    end

end


% ~~~~~~~~~~~  Check that there the same number of triggers for pupil and neural data ~~~~~~~~~~~
if ~isempty(neural_data) && ~isempty(pupil_data)

    if pupil_data.trigger.total ~= neural_data.trigger.total
        throw(makeME('exit','Error in %s: Found more neural triggers than video triggers',{mfilename}))
    end

end



%% Extract pupil and neural information for each trigger


% Time lock using respective triggers
if ~isempty(pupil_data)
    trigger_num = pupil_data.trigger.num;
    trigger_total = pupil_data.trigger.total;
    pupil_2_corr = pupil_data.(params.pupil_type);
else
    pupil_2_corr = NaN;
end

if ~isempty(neural_data)
    trigger_num = neural_data.trigger.num;
    trigger_total = neural_data.trigger.total;
    neural_2_corr = neural_data.(params.neural_type);
else
    neural_2_corr = NaN;
end
time_bins = -1*params.time_before : params.bin_width : params.time_after;
total_bins = numel(time_bins); % Check that each trigger has this many bins in it


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
neural_trigger_data = cell(trigger_total,1);
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
        pupil_start_time(t) = pupil_data.trigger.times(t,1) - params.time_before;
        pupil_trig_time(t)  = pupil_data.trigger.times(t,1);
        pupil_stop_time(t)  = pupil_data.trigger.times(t,1) + params.time_after;

        % Extract pupil data during trigger
        pupil_trigger_ind = (pupil_start_time(t) <= pupil_data.time) & (pupil_data.time <= pupil_stop_time(t));
        pupil_trigger_times{t} = pupil_data.time(pupil_trigger_ind);
        pupil_trigger_data{t} = pupil_data.raw(pupil_trigger_ind);

        % Get average pupil size per bin (using similar method as neural data)
        pupil_binned = NaN(1,total_bins);

        switch params.rate_type
            case 'binned'

                % Take median of pupil within bin
                for bin = 1:length(time_bins)
                    pupil_binned(bin) = median(pupil_trigger_data{t}((time_bins(bin) <= (pupil_trigger_times{t} - pupil_trig_time(t))) & ((pupil_trigger_times{t} - pupil_trig_time(t)) <= time_bins(bin)+params.bin_width)),'omitnan');
                end

            otherwise
                throw(makeME('exit','Error in %s: Did not make method for calculating firing rate for time-locked trigger pupil data',{mfilename}))
        end
    else
        pupil_binned = NaN(1,total_bins);
    end



    %%%%%%%% Neural stuff %%%%%%%%
    if ~isempty(neural_data)
        % Get exact times of events for neural data
        neural_start_time(t) = neural_data.trigger.time(t)-params.time_before;
        neural_trig_time(t)  = neural_data.trigger.time(t);
        neural_stop_time(t)  = neural_data.trigger.time(t)+params.time_after;

        % First, find spikes during the time lock
        neural_trigger_data{t} = neural_data.spike((neural_start_time(t) <= neural_data.spike) & ...
            (neural_data.spike <= neural_stop_time(t)))';
        neural_trigger_times{t} = neural_trigger_data{t} - neural_trig_time(t);

        % Second, get firing rate based on rate_type
        neural_binned = NaN(1,total_bins);

        switch params.rate_type
            case 'binned'

                % Count spikes per bin
                for bin = 1:length(time_bins)
                    neural_binned(bin) = sum( (time_bins(bin) <= (neural_trigger_data{t}-neural_trig_time(t))) & ((neural_trigger_data{t}-neural_trig_time(t)) <= time_bins(bin)+params.bin_width) ) / params.bin_width; % Convert to per second (aka Hz)
                end

            otherwise
                throw(makeME('exit','Error in %s: Did not make method for calculating firing rate for time-locked trigger neural data',{mfilename}))
        end
    else
        neural_binned = NaN(1,total_bins);
    end


    % Extract the data into cell arrays. Trials may not be the same size,
    % but that will be fixed later
    pupil2corr_trial{t} = pupil_binned;
    neural2corr_trial{t} = neural_binned;


    % Check to make sure there are enough neural bins
    if numel(neural2corr_trial{t}) ~= total_bins

        warning('Neural data to correlated did not have enough points. Setting this trial to NaNs!')
        neural_ind{t} = NaN(1,total_bins);
        pupil_ind{t} = NaN(1,total_bins);
        neural_time{t} = NaN(1,total_bins);
        pupil_time{t} = NaN(1,total_bins);
        pupil2corr_trial{t} = NaN(1,total_bins);
        neural2corr_trial{t} = NaN(1,total_bins);

    end

end % Loop through triggers



%%%%% CHOP OFF the end of the trials so all trials are the same size
if ~all(cellfun(@length, pupil2corr_trial) == cellfun(@length, neural2corr_trial))
    throw(makeME('ignore','Error in %s: Pupil and neural data are not the same size',{mfilename}))
end


pupil = struct(...
    'type',         params.pupil_type,...
    'rate_type',    params.rate_type,...
    'raw',          struct(... % times, index, raw data, ect
    'time',struct('start',pupil_start_time,'trig',pupil_trig_time,'stop',pupil_stop_time,'time',{pupil_trigger_times}),...
    'data',{pupil_trigger_data}), ...
    'data',         cell2mat(pupil2corr_trial)); % 10 x 16 double


neural = struct(...
    'type',         params.neural_type,...
    'rate_type',    params.rate_type,...
    'raw',          struct(... % times, index, raw data, ect
    'time',         struct('start',neural_start_time,'trig',neural_trig_time,'stop',neural_stop_time,'time',{neural_trigger_times}),...
    'data',         {neural_trigger_data}), ...
    'data',         cell2mat(neural2corr_trial)); % 10 x 16 double


trigger_data = struct( ...
    'pupil',         pupil, ...
    'neural',        neural, ...
    'trigger_num',   trigger_num,...
    'trigger_total', trigger_total,...
    'time_before',   params.time_before,...
    'bin_width',     params.bin_width,...
    'time_after',    params.time_after,...
    'time_bins',     time_bins,...
    'total_bins',    total_bins);


disp(['Done loading data for ',strjoin({mouse_id, date_str, file_num, exper_cond_field},', '),'.'])


end % End load data