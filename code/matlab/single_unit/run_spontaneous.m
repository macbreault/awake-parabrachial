% run_spontaneous.m
%
% Script that plots spontanious activity vs. pupil size based on 5 minutes before first trigger in LASER files
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 08-17-2022
% Modified: 08-08-2022 - Added threshold that removes firing rate below thresh from being correlated


clear
clc
set(0,'DefaultFigureWindowStyle','docked') %docks the figures onto the Matlab window


CI = @(x) 1.96 * std(x,[],1,'omitnan') / sqrt(size(x,2)); % Function
step_plot = @(x) reshape([x;x],1,2*numel(x));


%% Set data

% ~~~~~~~~~~~ Set mouse information ~~~~~~~~~~~
mouse_id = 'C7';                        % Mouse ID
date_str = '7-1-2022';                   % Date of recording
file_num = 'File1';                     % File number ('File1')


% Other variables to change
run_all = 0;  % Boolean as to whether to run all conditions (1) or not (0)
save_it = 0; % Whether to save plots (1) or not (0)

remove_blink = 1; % Remove blink (1) or not (0)
rate_type = 'binned'; % Use 'binned' or 'sliding' window for firing rate
bin_width = 1; % seconds (MUST BE LARGER THAN 1/FRAME RATE=0.03333)
norm_it = 1; % Normalize pupil size by largest size (1) or not (0)

pupil_type = 'raw'; % 'raw' or 'filter'
neural_type = 'raw'; % 'raw' or 'filter'
T = 2; % Trigger number to plot in example

thresh = 10; % Filter for neural data
%thresh = '10%'; % Filter for neural data
% Either use filter on Hz (ex. 0.1, or 10) OR could use '10%' cut-off (ignore neural data under 10%)




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
time_after = 0; % Seconds after trigger
time_before = 300; % Seconds before trigger

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


        % Get threshold info
        if isa(thresh,'char')
            thresh_prct = str2double(erase(thresh,'%'));
            thresh_val = prctile(trigger_data.neural.data(:), thresh_prct);
        else
            thresh_val = thresh;
        end


        %%%%%% Prepare for plotting %%%%%%

        % Create experiment title for all plots
        exper_title = {['Mouse ',num2str(mouse_id),', ',date_str,', ',file_num, ', ', exper_cond];...
            ['rate type = ',rate_type,', bin width = ',num2str(bin_width),' (s), remove blink=',num2str(remove_blink),', thresh=',num2str(thresh),' Hz']};

        params_folder = strjoin({...
                    [rate_type,'_',num2str(bin_width)],...
                    ['cutoff_',num2str(thresh)]
                    },'|');

        plot_path  = fullfile(paths.dropbox,'plots',findPlotFolder,'spontaneous', params_folder, mouse_id, date_str, file_num);
        if ~exist(plot_path,'dir') && save_it; mkdir(plot_path); end




        %% Plot raw data

        run('plot_raw_data.m')



        %% Plot spontaneous activity (time locked to 5 minutes before first trigger)
        
        T = 1; % Trigger number

        % ~~~~~~~~~ Check data ~~~~~~~~~

        % Check that trigger 1 exists
        if ~any(pupil_data.trigger.num == T) || ~any(neural_data.trigger.num == T)
            warning('Error in %s: Trigger 1 was removed, cannot accurately get first 5 minutes of data. Skip!',mfilename)
            continue
        end

        % Find time of first trigger (in seconds)
        first_trigger_time_pupil = pupil_data.trigger.times(pupil_data.trigger.num == T,1);
        first_trigger_time_neural = neural_data.trigger.time(neural_data.trigger.num == T);

        % Check that there is at least 5 minutes (300 seconds) before first trigger
        if (first_trigger_time_pupil < 300) || (first_trigger_time_neural < 300)
            warning('Error in %s: Less than 5 minutes of data before trigger 1. Skip!',mfilename)
            continue
        end
        


        % ~~~~~~~~~ Plot ~~~~~~~~~
        % Inialize figure
        fig = figure(2);
        clf(fig)
        fig.WindowStyle = fig_params.WindowStyle;
        fig.Units = fig_params.Units;
        fig.Position = fig_params.Position;

        layout = tiledlayout(2,1,'TileSpacing','normal','Padding','tight');

        title(layout, exper_title)

        % Plotting constants
        binned_time = trigger_data.time_bins;



        % ~~~~~~~~~~ Plot RAW DATA ~~~~~~~~~~ 

        nexttile(layout)
        hold on

        %%%%% Plot neural %%%%%
        yyaxis left
        plot_time = step_plot(binned_time);
        plot_data = step_plot(trigger_data.neural.data(T,:));
        plot(plot_time(2:end), plot_data(1:end-1),'-','Color','r','DisplayName','Neural','LineWidth',1)
        axis tight

        % Plot spike
        plot(trigger_data.neural.raw.time.time{T}, ones(size(trigger_data.neural.raw.time.time{T})),'|','MarkerSize',2.5,'Color','k','DisplayName','Spikes');
        
        % Plot line for threshold pupil
        hold on
        yline(thresh_val,'Tag','thresh')
        hold off

        %%%%% Plot pupil %%%%%
        yyaxis right
        plot_time = step_plot(binned_time);
        plot_data = step_plot(trigger_data.pupil.data(T,:));
        plot(plot_time(2:end), plot_data(1:end-1),'-','Color','b','DisplayName','Pupil','LineWidth',0.5)
        axis tight


        % Format
        box on


        % Label
        xtickangle(0)
        xlabel(['Time before trigger ',num2str(T), ' (s)'])

        yyaxis right
        set(gca,'YColor','b')
        if norm_it
            ylabel('Norm pupil area (pixels^2)')
        else
            ylabel('Pupil area (pixels^2)')
        end

        yyaxis left
        set(gca,'YColor','r')
        ylabel('Spike rate (Hz)')
        title(['Spontaneous activity before trigger ',num2str(T)])




        % ~~~~~~~~~~ Plot CORRELATION DATA ~~~~~~~~~~ 
        nexttile(layout)

        thresh_ind = trigger_data.neural.data(T,:)' >= thresh_val; % Indices to plot and correlate, whose neural data is at least about the threshold only

        % Data
        [r,p] = corr(trigger_data.pupil.data(T,thresh_ind)', trigger_data.neural.data(T,thresh_ind)', 'type', 'Pearson','rows','pairwise');

        % Plot
        plot(trigger_data.pupil.data(T,thresh_ind), trigger_data.neural.data(T,thresh_ind),'o','DisplayName',['r=',num2str(r)],'MarkerSize',4,'MarkerFaceColor','none','MarkerEdgeColor',[0.5 0.5 0.5])
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
        if norm_it
            xlabel(['Norm \bf{',lower(trigger_data.pupil.type(1:end)),'} \rm{pupil area (pixels^2)}'])
        else
            xlabel(['\bf{',upper(trigger_data.pupil.type(1)),lower(trigger_data.pupil.type(2:end)),'} \rm{pupil area (pixels^2)}'])
        end
        ylabel(['\bf{',upper(trigger_data.neural.type(1)),lower(trigger_data.neural.type(2:end)),'} \rm{spiking rate (Hz)}'])
        legend('location','eastoutside')


        % Save
        uniformFormat(fig)
        if save_it
            savefig(fig,fullfile(plot_path,'2_spont.fig'))
            exportgraphics(fig,fullfile(plot_path,'2_spont.pdf'),'Resolution',300,'ContentType','vector')
        end



        disp(['             Done running ',date_str,', ',file_num, ', ',exper_cond])

    end % End loop through conditions

    disp(' ')

end % End loop through mice

disp(' ')