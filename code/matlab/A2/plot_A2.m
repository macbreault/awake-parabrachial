% plot_A2.m
%
% Script that plots responses from A2
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 01-13-2023

clear
close all
clc

% p = genpath('/Users/Jessesmith 1/Desktop/facemap-main/JESSE FIG FILES/Library');
% addpath(p);
% set(0,'DefaultFigureWindowStyle','docked') %docks the figures onto the Matlab window

step_plot = @(x) reshape([x;x],1,2*numel(x));

save_it = 1;

%% Load data

% ~~~~~~~~~~~ Set mouse information ~~~~~~~~~~~
date_str = '11-19-20';   % Date of recording
%date_str = '11-11-20';   % Date of recording
%date_str = '12-2-20';   % Date of recording




rate_type = 'binned'; % Use 'binned' or 'sliding' window for firing rate
bin_width = 0.05; % seconds (MUST BE LARGER THAN 1/FRAME RATE=0.03333)
neural_type = 'raw'; % 'raw' or 'filter'

time_before = 5; % Seconds before trigger
time_after = 10; % Seconds after trigger


% ~~~~~~~~~~~ Other variables ~~~~~~~~~~~
% Figure configurations
fig_params = struct('Position',[0 0 8 3.5], 'Units','Inches', 'WindowStyle','normal', 'FontSize', 10);



% ~~~~~~~~~~~ Whose computer dis is? ~~~~~~~~~~~

paths = struct('dropbox',[],'neural',[]);

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

addpath(genpath(fullfile(paths.dropbox,'code'))); % Add libraries to directory
plot_path = fullfile(paths.dropbox,'plots','A2', date_str);
if ~exist(plot_path,'dir') && save_it; mkdir(plot_path); end


% ~~~~~~~~~~~ Get pinch folder information ~~~~~~~~~~~
a2_folder = fullfile(paths.dropbox,'data/A2',date_str,'Sorted');

% Find folders inside Sorted
listings = dir(a2_folder);
folder_ind = [listings.isdir] & contains({listings.name},'File');
listings(~folder_ind) = [];



%% Loop through Files

folders = {listings.name};
for folder = folders


    % Skip if File does not have a MATLAB folder
    folder_path = fullfile(a2_folder, folder{:}, 'MATLAB');
    if ~exist(folder_path,'dir')
        continue
    end

    % Check if event, spike, and pinch files exist
    if not (exist(fullfile(folder_path,'pinch.mat'),'file') && exist(fullfile(folder_path,'event.mat'),'file') && exist(fullfile(folder_path,'spike.mat'),'file'))
        continue
    end


    %% Process spike data

    % ~~~~~~~~~~~ Load spike data ~~~~~~~~~~~
    spike_struct = load(fullfile(folder_path,'spike.mat'));

    % Fix files
    spike_name = fieldnames(spike_struct);

    if numel(spike_name) == 1
        spike_data = spike_struct.(spike_name{:}); % Time each spike is occuring
    else
        throw(makeME('ignore','Error in %s: Found more than 1 spike file',{mfilename}))
    end

    % Process
    % Round based on smpRate
    smpRate = 2.5e-05; % 40kHz -> 40000 Hz -> 1/40000 Rate
    spike_data = round(spike_data,6);
    neural_data.spike = spike_data;

    % Calculate firing rate (and assigns neural_data.time and neural_data.raw)
    switch rate_type
        case 'binned'

            % Use binned count
            time_end = ceil(spike_data(end)); % Time of last spike
            neural_data.time = (0:bin_width:time_end)';
            neural_data.raw = NaN(size(neural_data.time)); % Initalize vector that keeps track of how many spikes are in each 1 second bin
            for t = 1:length(neural_data.time)
                neural_data.raw(t) = sum(  (neural_data.time(t) <= spike_data) & (spike_data < (neural_data.time(t)+bin_width)) ) / bin_width; % Convert to per second (aka Hz)
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
            k = rectwin(find(neural_data.time == bin_width)-1); % Put bin width here!
            firing_rate = conv(spike_train,k);
            firing_rate = firing_rate(ceil(length(k)/2):end-floor(length(k)/2));

            neural_data.raw = firing_rate ./ bin_width;

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


    %% Process Local field potentials (pinch)

    % ~~~~~~~~~~~ Load pinch data ~~~~~~~~~~~
    pinch_struct = load(fullfile(folder_path,'pinch.mat'));

    % Fix files
    pinch_name = fieldnames(pinch_struct);

    if numel(pinch_name) == 1
        mech_data = pinch_struct.(pinch_name{:}); % Time each spike is occuring
    else
        throw(makeME('ignore','Error in %s: Found more than 1 pinch file',{mfilename}))
    end

    %pinch_data_raw = resample(mech_data(:,1),mech_data(:,2),10);
    %pinch_data_time = linspace(mech_data(1,2), mech_data(end,2), numel(pinch_data_raw_1));

    pinch_data = struct('pinch',mech_data, 'time', mech_data(:,2), 'raw', mech_data(:,1));


    %% Process data into events

    event_struct = load(fullfile(folder_path,'event.mat'));

    % Fix files
    event_name = fieldnames(event_struct);

    if numel(event_name) == 1
        event_data = event_struct.(event_name{:}); % Time each spike is occuring
    else
        throw(makeME('ignore','Error in %s: Found more than 1 event file',{mfilename}))
    end

    time_bins = -1*time_before : bin_width : time_after;
    total_bins = numel(time_bins); % Check that each trigger has this many bins in it

    trigger_total = numel(event_data);
    trigger_num = 1:numel(trigger_total);

    trigger_loop = 1:trigger_total; % Triggers to actually loop through

    % Initialize



    % Run
    for t = trigger_loop

        % ~~~~~~~~~~~~~~~~ Pinch stuff ~~~~~~~~~~~~~~~~
        if ~isempty(pinch_data)
            % Get exact times of events for pinch data
            pinch_start_time(t) = event_data(t,1) - time_before;
            pinch_trig_time(t)  = event_data(t,1);
            pinch_stop_time(t)  = event_data(t,1) + time_after;

            % Extract pinch data during trigger
            pinch_trigger_ind = (pinch_start_time(t) <= pinch_data.time) & (pinch_data.time <= pinch_stop_time(t));
            pinch_trigger_times{t} = pinch_data.time(pinch_trigger_ind);
            pinch_trigger_data{t} = pinch_data.raw(pinch_trigger_ind);

            % Get average pinch size per bin (using similar method as neural data)
            pinch_binned = NaN(1,total_bins);

            switch rate_type
                case 'binned'

                    % Take median of pinch within bin
                    for bin = 1:length(time_bins)
                        pinch_binned(bin) = median(pinch_trigger_data{t}((time_bins(bin) <= (pinch_trigger_times{t} - pinch_trig_time(t))) & ((pinch_trigger_times{t} - pinch_trig_time(t)) <= time_bins(bin)+bin_width)),'omitnan');
                    end

                otherwise
                    throw(makeME('exit','Error in %s: Did not make method for calculating firing rate for time-locked trigger pinch data',{mfilename}))
            end
        else
            pinch_binned = NaN(1,total_bins);


        end


        % ~~~~~~~~~~~~~~~~ Neural ~~~~~~~~~~~~~~~~
        if ~isempty(neural_data)
            % Get exact times of events for neural data
            neural_start_time(t) = event_data(t,1)-time_before;
            neural_trig_time(t)  = event_data(t,1);
            neural_stop_time(t)  = event_data(t,1)+time_after;

            % First, find spikes during the time lock
            neural_trigger_data{t} = neural_data.spike((neural_start_time(t) <= neural_data.spike) & ...
                (neural_data.spike <= neural_stop_time(t)))';
            neural_trigger_times{t} = neural_trigger_data{t} - neural_trig_time(t);

            % Second, get firing rate based on rate_type
            neural_binned = NaN(1,total_bins);

            switch rate_type
                case 'binned'

                    % Count spikes per bin
                    for bin = 1:length(time_bins)
                        neural_binned(bin) = sum( (time_bins(bin) <= (neural_trigger_data{t}-neural_trig_time(t))) & ((neural_trigger_data{t}-neural_trig_time(t)) <= time_bins(bin)+bin_width) ) / bin_width; % Convert to per second (aka Hz)
                    end

                otherwise
                    throw(makeME('exit','Error in %s: Did not make method for calculating firing rate for time-locked trigger neural data',{mfilename}))
            end
        else
            neural_binned = NaN(1,total_bins);
        end


        % Extract the data into cell arrays. Trials may not be the same size,
        % but that will be fixed later
        pinch2corr_trial{t} = pinch_binned;
        neural2corr_trial{t} = neural_binned;


        % Check to make sure there are enough neural bins
        if numel(neural2corr_trial{t}) ~= total_bins

            warning('Neural data to correlated did not have enough points. Setting this trial to NaNs!')
            neural_ind{t} = NaN(1,total_bins);
            pinch_ind{t} = NaN(1,total_bins);
            neural_time{t} = NaN(1,total_bins);
            pinch_time{t} = NaN(1,total_bins);
            pinch2corr_trial{t} = NaN(1,total_bins);
            neural2corr_trial{t} = NaN(1,total_bins);

        end

    end % Loop through triggers


    %%%%% CHOP OFF the end of the trials so all trials are the same size
    if ~all(cellfun(@length, pinch2corr_trial) == cellfun(@length, neural2corr_trial))
        throw(makeME('ignore','Error in %s: pinch and neural data are not the same size',{mfilename}))
    end



    pinch = struct(...
        'rate_type',    rate_type,...
        'raw',          struct(... % times, index, raw data, ect
        'time',         struct('start',pinch_start_time,'trig',pinch_trig_time,'stop',pinch_stop_time,'time',{pinch_trigger_times}),...
        'data',         {pinch_trigger_data}), ...
        'data',         cell2mat(pinch2corr_trial')); % 10 x 16 double


    neural = struct(...
        'type',         neural_type,...
        'rate_type',    rate_type,...
        'raw',          struct(... % times, index, raw data, ect
        'time',         struct('start',neural_start_time,'trig',neural_trig_time,'stop',neural_stop_time,'time',{neural_trigger_times}),...
        'data',         {neural_trigger_data}), ...
        'data',         cell2mat(neural2corr_trial')); % 10 x 16 double


    trigger_data = struct( ...
        'pinch',         pinch, ...
        'neural',        neural, ...
        'trigger_num',   trigger_num,...
        'trigger_total', trigger_total,...
        'time_before',   time_before,...
        'bin_width',     bin_width,...
        'time_after',    time_after,...
        'time_bins',     time_bins,...
        'total_bins',    total_bins);



    %% Plot punch per day

    fig = figure(1);
    clf(fig)
    fig.WindowStyle = fig_params.WindowStyle;
    fig.Units = fig_params.Units;
    fig.Position = fig_params.Position;

    layout = tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
    nexttile(layout)

    % Plot SPIKES (left)
    yyaxis left
    hold on
    plot_time = step_plot(trigger_data.time_bins);
    plot_data = step_plot(mean(trigger_data.neural.data,1,'omitnan'));
    plot_CI = step_plot(calculateCR(trigger_data.neural.data')');
    H = shadedErrorBar(plot_time(2:end), plot_data(1:end-1), plot_CI(1:end-1),'lineProps',{'r'});
    H.mainLine.DisplayName = 'Neural';
    hold off


    % Plot PINCH (right)
    yyaxis right
    hold on
    plot_time = step_plot(trigger_data.time_bins);
    plot_data = step_plot(mean(trigger_data.pinch.data,1,'omitnan'));
    plot_CI = step_plot(calculateCR(trigger_data.pinch.data')');
    H = shadedErrorBar(plot_time(2:end), plot_data(1:end-1), plot_CI(1:end-1),'lineProps',{'b'});
    H.mainLine.DisplayName = 'Pinch';
    hold off


    xline(0,'Tag','trigger')

    % Format
    xline(0,'Tag','trigger')
    box off
    xlim([-time_before, time_after])
    xticks(unique([-time_before:1:0,0:1:time_after]))

    % Label
    title(['Number of triggers = ',num2str(trigger_data.trigger_total)])
    xtickangle(0)

    xlabel('Time since trigger (s)')

    yyaxis left
    set(gca,'YColor','r')
    ylabel('Spike rate (Hz)')


    yyaxis right
    set(gca,'YColor','b')
    ylabel('Pinch (mV)')


    % ~~~~~~~~~~ Uniformly format plots ~~~~~~~~~~

    yyaxis left
    uniformFormat(fig)
    yyaxis right
    uniformFormat(fig)


    % ~~~~~~~~~~ Save ~~~~~~~~~~
    if save_it
        if ~exist(fullfile(plot_path,folder{:}),'dir'); mkdir(fullfile(plot_path,folder{:})); end
        savefig(fig,fullfile(plot_path,folder{:},'1_average.fig'))
        exportgraphics(fig,fullfile(plot_path,folder{:},'1_average.pdf'),'Resolution',300,'ContentType','vector')
    end




end
