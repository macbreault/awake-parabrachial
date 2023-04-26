% plot_corr_conditions.m
%
% Script that plots STUFF using NEW conditioning experiment
%
% Uses info from video files and arduino for timing of tone
%
% Gets correlations of pupil size and firing rate immediately before,
% during and after the stimulation
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 06-14-2022 (v0)
% Modified: 08-04-2022 (v1) - Changed way data is loaded to use a function (loadData)
% Modified: 08-18-2022 (v2) - Changed so it now uses binned version of pupil and neural
% Modified: 09-15-2022 - Replaced anonymous function of CI with defined CI function, which uses a different method for calculating


clear
clc
% p = genpath('/Users/Jessesmith 1/Desktop/facemap-main/JESSE FIG FILES/Library');
% addpath(p);
set(0,'DefaultFigureWindowStyle','docked') %docks the figures onto the Matlab window

CI = @(x) 1.96 * std(x,[],1,'omitnan') / sqrt(size(x,2)); % Function
step_plot = @(x) reshape([x;x],1,2*numel(x));



%% Set data

% ~~~~~~~~~~~ Set mouse information ~~~~~~~~~~~
mouse_id = 'C3';              % Mouse ID
date_str = '9-6-22';         % Date of recording
file_num = 'File1';           % File number ('File1')
exper_cond = 'Extinction';         % What experiment were you recording
%{
- Tone Alone        Expect nothing
- Conditioning
- Laser
- Test
- Recondition
- Retest
- Extinction
- Reextinction
%}


% Other variables to change
run_all = 0;  % Boolean as to whether to run all conditions (1) or not (0)
save_it = 0; % Whether to save plots (1) or not (0)

remove_blink = 1; % Remove blink (1) or not (0)
rate_type = 'binned'; % Use 'binned' or 'sliding' window for firing rate
bin_width = 0.1; % seconds (MUST BE LARGER THAN 1/FRAME RATE=0.03333)norm_it = 1; % Normalize pupil size by largest size (1) or not (0)
norm_it = 1; % Normalize pupil size by largest size (1) or not (0)

pupil_type = 'raw'; % 'raw' or 'filter'
neural_type = 'raw'; % 'raw' or 'filter'
T = 2; % Trigger number to plot in example

time_before = 5; % Seconds
time_after = 10; % Seconds



% ~~~~~~~~~~~~~~~~~~~~~~ DO NOT CHANGE BELOW ~~~~~~~~~~~~~~~~~~~~~~

% Figure configurations
fig_params = struct('Position',[0 0 8 5], 'Units','Inches', 'WindowStyle','normal', 'FontSize', 10);
plot_cross = 0; % Plot cross-corr results (1) or not (0)


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

% Find path to neural data by first checking SEAGATE, then checking RED, then assigning DROPBOX
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

        plot_path  = fullfile(paths.dropbox,'plots',findPlotFolder,'conditioning', params_folder, mouse_id, date_str, file_num, exper_cond);
        if ~exist(plot_path,'dir') && save_it; mkdir(plot_path); end



        %% Plot raw data

        run('plot_raw_data.m')



        %% Plot split data

        run('plot_split_data.m')



        %% Plot average data

        run('plot_average_data.m')



        %% Plot correlation - NOT AVERAGED

        % Start figure
        fig = figure(2);
        clf(fig)
        fig.WindowStyle = fig_params.WindowStyle;
        fig.Units = fig_params.Units;
        fig.Position = fig_params.Position;

        layout = tiledlayout(2,1,'TileSpacing','compact','Padding','none');

        title(layout, exper_title)

        [r,p] = corr(trigger_data.pupil.data(:), trigger_data.neural.data(:), 'type', 'Pearson','rows','pairwise');


        % ~~~~~~~~~~~~~~~~~~~ Plot scatter plot of correlation between pupil and spikes ~~~~~~~~~~~~~~~~~~~
        nexttile(layout)
        hold on
        plot(trigger_data.pupil.data(:), trigger_data.neural.data(:),'o','DisplayName',['r=',num2str(r)],'MarkerSize',4,'MarkerFaceColor','none','MarkerEdgeColor',[0.5 0.5 0.5])
        ls = lsline; % Make least square line
        ls.DisplayName = ['p=',num2str(p)];
        if p <= 0.05
            ls.Color = 'm'; % If significant, then lsline turns red
            ls.LineWidth = 2;
        else
            ls.Color = 'k';
        end

        %%%%% Plot data used for trigger example on correlation plot %%%%%
        plot(trigger_data.pupil.data(T,:),trigger_data.neural.data(T,:),'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName',['Trigger ',num2str(T)])
        hold off

        % Format
        axis square
        box on

        % Label
        legend('location','eastoutside')
        xlabel(['\bf{',upper(trigger_data.pupil.type(1)),lower(trigger_data.pupil.type(2:end)),'} \rm{pupil area (pixels^2)}'])
        ylabel(['\bf{',upper(trigger_data.neural.type(1)),lower(trigger_data.neural.type(2:end)),'} \rm{spiking rate (Hz)}'])
        title([num2str(trigger_data.trigger_total),' trials NOT averaged'])




        % ~~~~~~~~~~~~~~~~~~~ Plot example of raw data used for correlation ~~~~~~~~~~~~~~~~~~~
        nexttile(layout)

        % Shift time so spike and pupil line up
        binned_time = trigger_data.time_bins;

        % Plot
        hold on

        %%%%% Plot binned neural %%%%%
        yyaxis left
        plot_time = step_plot(binned_time);
        plot_data = step_plot(trigger_data.neural.data(T,:));
        plot(plot_time(2:end), plot_data(1:end-1),'-','Color','r','DisplayName','Binned neural','LineWidth',1.5)

        %%%%% Plot spiking data %%%%%
        yyaxis left
        plot(trigger_data.neural.raw.time.time{T}, ones(size(trigger_data.neural.raw.time.time{T})),'|','MarkerSize',2.5,'Color','k','DisplayName','Spikes');

        %%%%% Plot raw pupil %%%%%
        yyaxis right
        plot(trigger_data.pupil.raw.time.time{T} - trigger_data.pupil.raw.time.trig(T),   trigger_data.pupil.raw.data{T},'-','Color','b','DisplayName','Raw pupil','LineWidth',0.5)

        %%%%% Plot binned pupil %%%%%
        yyaxis right
        plot_time = step_plot(binned_time);
        plot_data = step_plot(trigger_data.pupil.data(T,:));
        plot(plot_time(2:end),   plot_data(1:end-1), '-','Color','b','DisplayName','Binned pupil','LineWidth',2)

        xline(0,'LineWidth',2.5,'Color',[0.4940 0.1840 0.5560],'DisplayName','Trigger')

        hold off

        % Format
        box on
        xlim([-trigger_data.time_before, trigger_data.time_after])
        xticks(unique([-trigger_data.time_before:1:0,0:1:trigger_data.time_after]))

        % Label
        legend('location','eastoutside')
        title(['Trigger ',num2str(trigger_data.trigger_num(T))])
        xlabel('Time since trigger (s)')

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

        set(findobj(fig,'-property','FontSize'),'FontSize',fig_params.FontSize)


        if save_it
            savefig(fig,fullfile(plot_path,'2_corr_raw.fig'))
            exportgraphics(fig,fullfile(plot_path,'2_corr_raw.pdf'),'Resolution',300,'ContentType','vector')
        end



        %% Plot correlation - AVERAGED

        % Start figure
        fig = figure(3);
        clf(fig)
        fig.WindowStyle = fig_params.WindowStyle;
        fig.Units = fig_params.Units;
        fig.Position = fig_params.Position;

        layout = tiledlayout(2,1,'TileSpacing','compact','Padding','none');

        title(layout, exper_title)

        plot_pupil = mean(trigger_data.pupil.data,1,'omitnan');
        %CI_pupil = CI(trigger_data.pupil.data);
        CI_pupil = calculateCR(trigger_data.pupil.data')';
        plot_neural = mean(trigger_data.neural.data,1,'omitnan');
        %CI_neural = CI(trigger_data.neural.data);
        CI_neural = calculateCR(trigger_data.neural.data')';

        [r,p] = corr(plot_pupil', plot_neural', 'Type','Pearson','rows','pairwise');

        % Plot SCATTER DATA
        nexttile(layout)
        hold on
        plot(plot_pupil, plot_neural,'o','DisplayName',['r=',num2str(r)],'MarkerSize',4,'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerFaceColor','none')
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
        legend('location','eastoutside')
        xlabel(['\bf{',upper(trigger_data.pupil.type(1)),lower(trigger_data.pupil.type(2:end)),'} \rm{pupil area (pixels^2)}'])
        ylabel(['\bf{',upper(trigger_data.neural.type(1)),lower(trigger_data.neural.type(2:end)),'} \rm{spiking rate (Hz)}'])
        title([num2str(trigger_data.trigger_total), ' trials averaged'])



        % Plot AVERAGED data
        nexttile(layout)
        binned_time = trigger_data.time_bins;

        yyaxis left
        plot_time = step_plot(binned_time);
        plot_data = step_plot(plot_neural);
        plot_CI = step_plot(CI_neural);
        H = shadedErrorBar(plot_time(2:end), plot_data(1:end-1), plot_CI(1:end-1),'lineProps',{'r'});
        delete(H.edge)
        set(gca,'YColor','r')
        ylabel('Spike rate (Hz)')

        yyaxis right
        plot_time = step_plot(binned_time);
        plot_data = step_plot(plot_pupil);
        plot_CI = step_plot(CI_pupil);
        H = shadedErrorBar(plot_time(2:end), plot_data(1:end-1), plot_CI(1:end-1),'lineProps',{'b'});
        delete(H.edge)
        set(gca,'YColor','b')
        if norm_it
            ylabel('Norm pupil area (pixels^2)')
        else
            ylabel('Pupil area (pixels^2)')
        end

        box on
        xlim([-trigger_data.time_before, trigger_data.time_after])
        xline(0,'LineWidth',2.5,'Color',[0.4940 0.1840 0.5560])
        xticks(unique([-trigger_data.time_before:1:0,0:1:trigger_data.time_after]))
        xlabel('Time since trigger (s)')


        if save_it
            savefig(fig,fullfile(plot_path,'3_corr_avg.fig'))
            exportgraphics(fig,fullfile(plot_path,'3_corr_avg.pdf'),'Resolution',300,'ContentType','vector')
        end



        %% Plot cross-correlation - NOT AVERAGED

        if plot_cross

            error('Modify code to add cross-correlation as was done in run_laser_corr.m')
            % Start figure
            fig = figure(4);
            clf(fig)
            fig.WindowStyle = fig_params.WindowStyle;
            fig.Units = fig_params.Units;
            fig.Position = fig_params.Position;

            layout = tiledlayout(2,2,'TileSpacing','compact','Padding','none');

            title(layout, exper_title)

            x_pupil = trigger_data.pupil.data(:);
            x_neural = trigger_data.neural.data(:);
            NumLags = 10/bin_width;

            [R,lag] = crosscorr(x_pupil, x_neural, NumLags); % Maximum number of lag = 10 seconds  NumLags = NumLags
            [R_max,ind_max] = max(abs(R));
            lag_max = lag(ind_max);

            if lag_max > 0
                pupil_ind = 1:(length(x_pupil)-lag_max);
                neural_ind = (lag_max+1):length(x_neural);
            elseif lag_max == 0
                pupil_ind = 1:length(x_pupil);
                neural_ind = 1:length(x_neural);
            elseif lag_max < 0
                pupil_ind = (-1*lag_max+1):length(x_pupil);
                neural_ind = 1:(length(x_pupil)-lag_max*-1);
            end
            [r,p] = corr(x_pupil(pupil_ind), x_neural(neural_ind), 'Type','Pearson','rows','pairwise');



            % Plot CROSS-CORRELATION
            nexttile(layout)
            crosscorr(x_pupil, x_neural, NumLags); % NumLags = NumLags
            %legend(['Lag (s) = ',num2str(lag_max*bin_width)])
            xlabel('Lag (s)')
            obj = get(gca,'Children');
            arrayfun(@(i) set(obj(i),'XData',get(obj(i),'XData')*bin_width), 1:numel(obj))
            title('')
            hold on
            x = xline(lag_max*bin_width, 'DisplayName',['Lag (s) = ',num2str(lag_max*bin_width)],'Color',[0.4940 0.1840 0.5560],'LineWidth',2.5);
            legend(x)
            hold off


            % Plot SCATTER DATA
            nexttile(layout)
            hold on
            plot(x_pupil(pupil_ind), x_neural(neural_ind),'o','DisplayName',['r=',num2str(r)],'MarkerSize',4,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor','none')
            ls = lsline; % Make least square line
            ls.DisplayName = ['p=',num2str(p)];
            if p<=0.05
                ls.Color = 'm'; % If significant, then lsline turns red
                ls.LineWidth = 2;
            else
                ls.Color = 'k';
            end

            % Format
            axis square
            box on

            % Label
            legend('location','eastoutside')
            xlabel(['\bf{',upper(trigger_data.pupil.type(1)),lower(trigger_data.pupil.type(2:end)),'} \rm{pupil area (pixels^2)}'])
            ylabel(['\bf{',upper(trigger_data.neural.type(1)),lower(trigger_data.neural.type(2:end)),'} \rm{spiking rate (Hz)}'])
            title([num2str(trigger_data.trigger_total), ' trials NOT averaged'])



            % Plot AVERAGED data
            nexttile(layout, [1,2])
            % Shift time so spike and pupil line up
            binned_time = linspace(-trigger_data.time_before, trigger_data.time_after, trigger_data.num_bins);
            full_time = linspace(-trigger_data.time_before, trigger_data.time_after, numel(trigger_data.pupil.time(T,1):trigger_data.pupil.time(T,2)));

            % Get trial ind
            if lag_max > 0
                pupil_trial_ind = 1:(length(binned_time)-lag_max);
                neural_trial_ind = (lag_max+1):length(binned_time);
            elseif lag_max == 0
                pupil_trial_ind = 1:length(binned_time);
                neural_trial_ind = 1:length(binned_time);
            elseif lag_max < 0
                pupil_trial_ind = (-1*lag_max+1):length(binned_time);
                neural_trial_ind = 1:(length(binned_time)-lag_max*-1);
            end
            [r,p] = corr(trigger_data.pupil.data(T,pupil_trial_ind)', trigger_data.neural.data(T,neural_trial_ind)', 'Type','Pearson','rows','pairwise');

            % Plot
            hold on
            %%%%% Plot raw neural %%%%%
            yyaxis left
            x = neural.raw(trigger_data.neural.time(T,1):trigger_data.neural.time(T,2));
            plot(((1:length(neural_trial_ind))-1)*bin_width, x(neural_trial_ind),'-','Color','r','DisplayName','Raw neural','LineWidth',0.5)
            x = neural.raw(trigger_data.neural.ind(T,:));
            plot(((1:length(neural_trial_ind))-1)*bin_width, x(neural_trial_ind), '*', 'Color','r','DisplayName','Raw neural ind','MarkerSize',4)

            %%%%% Plot filtet neural %%%%%
            yyaxis left
            x = neural.filter(trigger_data.neural.time(T,1):trigger_data.neural.time(T,2));
            plot(((1:length(neural_trial_ind))-1)*bin_width, x(neural_trial_ind),'-','Color','m','DisplayName','Filter neural','LineWidth',0.5)
            x = neural.filter(trigger_data.neural.ind(T,:));
            plot(((1:length(neural_trial_ind))-1)*bin_width, x(neural_trial_ind), '*', 'Color','m','DisplayName','Filter neural ind','MarkerSize',4)

            %%%%% Plot raw pupil %%%%%
            yyaxis right
            x = pupil.raw(trigger_data.pupil.time(T,1):trigger_data.pupil.time(T,2));

            if lag_max > 0
                [~,max_lag_ind] = min(abs(full_time - binned_time(pupil_trial_ind(end))));
                full_time_ind = 1:max_lag_ind;
            elseif lag_max == 0
                full_time_ind = 1:length(full_time);
            elseif lag_max < 0
                [~,max_lag_ind] = min(abs(full_time - binned_time(pupil_trial_ind(1))));
                full_time_ind = max_lag_ind:length(full_time);
            end

            plot(full_time(full_time_ind) - full_time(full_time_ind(1)),   x(full_time_ind), '-','Color','b','DisplayName','Raw pupil','LineWidth',0.5)
            x = pupil.raw(trigger_data.pupil.ind(T,:));
            plot(((1:length(pupil_trial_ind))-1)*bin_width, x(pupil_trial_ind), '*', 'Color','b','DisplayName','Raw pupil ind','MarkerSize',4)

            %%%%% Plot filter pupil %%%%%
            x = pupil.filter(trigger_data.pupil.time(T,1):trigger_data.pupil.time(T,2));
            plot(full_time(full_time_ind) - full_time(full_time_ind(1)),   x(full_time_ind),'-','Color','c','DisplayName','Filter pupil','LineWidth',0.5)
            x = pupil.filter(trigger_data.pupil.ind(T,:));
            plot(((1:length(pupil_trial_ind))-1)*bin_width, x(pupil_trial_ind), '*', 'Color','c','DisplayName','Filter pupil ind','MarkerSize',4)
            hold off

            % Format
            box on
            %axis
            %xlim([-trigger_data.time_before, trigger_data.time_after])
            %xticks(unique([-trigger_data.time_before,0:10:trigger_data.time_after]))

            % Label
            legend('location','eastoutside')
            title(['Trigger ',num2str(trigger_data.trigger_num(T))])
            xlabel('Time (s)')

            yyaxis right
            set(gca,'YColor','b')
            ylabel('Pupil area (pixels^2)')

            yyaxis left
            set(gca,'YColor','r')
            ylabel('Spike rate (Hz)')

            set(findobj(fig,'-property','FontSize'),'FontSize',fig_params.FontSize)


            if save_it
                savefig(fig,fullfile(plot_path,'4_cross_corr_raw.fig'))
                exportgraphics(fig,fullfile(plot_path,'4_cross_corr_raw.pdf'),'Resolution',300, 'ContentType','vector')
            end
        end




        %% Plot cross-correlation - AVERAGED

        if plot_cross

            % Start figure
            fig = figure(5);
            clf(fig)
            fig.WindowStyle = fig_params.WindowStyle;
            fig.Units = fig_params.Units;
            fig.Position = fig_params.Position;

            layout = tiledlayout(2,2,'TileSpacing','compact','Padding','none');

            title(layout, exper_title)

            x_pupil = mean(trigger_data.pupil.data,1,'omitnan');
            x_neural = mean(trigger_data.neural.data,1,'omitnan');
            NumLags = 10/bin_width;

            [R,lag] = crosscorr(x_pupil, x_neural, NumLags); % Shift second one over!!! NumLags = NumLags
            [R_max,ind_max] = max(abs(R));
            lag_max = lag(ind_max);

            if lag_max > 0
                pupil_ind = 1:(length(x_pupil)-lag_max);
                neural_ind = (lag_max+1):length(x_neural);
            elseif lag_max == 0
                pupil_ind = 1:length(x_pupil);
                neural_ind = 1:length(x_neural);
            elseif lag_max < 0
                pupil_ind = (-1*lag_max+1):length(x_pupil);
                neural_ind = 1:(length(x_pupil)-lag_max*-1);
            end
            [r,p] = corr(x_pupil(pupil_ind)', x_neural(neural_ind)', 'Type','Pearson','rows','pairwise');



            % Plot CROSS-CORRELATION
            nexttile(layout)
            crosscorr(x_pupil, x_neural, NumLags); % NumLags = NumLags
            %legend(['Lag (s) = ',num2str(lag_max*bin_width)])
            xlabel('Lag (s)')
            obj = get(gca,'Children');
            arrayfun(@(i) set(obj(i),'XData',get(obj(i),'XData')*bin_width), 1:numel(obj))
            title('')
            hold on
            x = xline(lag_max*bin_width, 'DisplayName',['Lag (s) = ',num2str(lag_max*bin_width)],'Color',[0.4940 0.1840 0.5560],'LineWidth',2.5);
            legend(x)
            hold off


            % Plot SCATTER DATA
            nexttile(layout)
            hold on
            plot(x_pupil(pupil_ind), x_neural(neural_ind),'o','DisplayName',['r=',num2str(r)],'MarkerSize',4,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor','none')
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
            legend('location','eastoutside')
            xlabel(['\bf{',upper(trigger_data.pupil.type(1)),lower(trigger_data.pupil.type(2:end)),'} \rm{pupil area (pixels^2)}'])
            ylabel(['\bf{',upper(trigger_data.neural.type(1)),lower(trigger_data.neural.type(2:end)),'} \rm{spiking rate (Hz)}'])
            title([num2str(trigger_data.trigger_total), ' trials averaged'])



            % Plot AVERAGED data
            nexttile(layout, [1,2])
            % Shift time so spike and pupil line up
            yyaxis left
            H = shadedErrorBar(((1:length(neural_ind))-1)*bin_width, mean(trigger_data.neural.data(:,neural_ind),1,'omitnan'), CI(trigger_data.neural.data(:,neural_ind)),'lineProps',{'r'});
            delete(H.edge)
            set(gca,'YColor','r')
            ylabel('Spike rate (Hz)')

            yyaxis right
            H = shadedErrorBar(((1:length(pupil_ind))-1)*bin_width, mean(trigger_data.pupil.data(:,pupil_ind),1,'omitnan'), CI(trigger_data.pupil.data(:,pupil_ind)),'lineProps',{'b'});
            delete(H.edge)
            set(gca,'YColor','b')
            ylabel('Pupil area (pixels^2)')

            box on
            %xlim([-trigger_data.time_before, trigger_data.time_after])
            %xline(0,'k','LineWidth',2.5)
            xlabel('Time (s)')


            % Label
            title(['Trigger ',num2str(trigger_data.trigger_num(T))])
            xlabel('Time (s)')

            yyaxis right
            set(gca,'YColor','b')
            ylabel('Pupil area (pixels^2)')

            yyaxis left
            set(gca,'YColor','r')
            ylabel('Spike rate (Hz)')

            set(findobj(fig,'-property','FontSize'),'FontSize',fig_params.FontSize)



            if save_it
                savefig(fig,fullfile(plot_path,'5_cross_corr_avg.fig'))
                exportgraphics(fig,fullfile(plot_path,'5_cross_corr_avg.pdf'),'Resolution',300, 'ContentType','vector')
            end
        end


        %% Run final stuff before ending loop

        disp(['             Done running ',date_str,', ',file_num, ', ',exper_cond])

    end % End loop through conditions

    disp(' ')

end % End loop through mice

disp(' ')


