% plot_figure_9.m
%
% Script to makes figure 9 for publication, about single units EVOKED data
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Modified: 10-21-2022

%clear
clc
set(0,'DefaultFigureWindowStyle','normal') %docks the figures onto the Matlab window

save_it = 1; % Whether to save plots (1) or not (0)



%% Set data


% ~~~~~~~~~~~ Set mouse information ~~~~~~~~~~~

% from run_evoked.m
params = struct('mouse_id',  'C6',...
    'date_str',  '7-27-2022',...
    'file_num', 'File1', ...
    'remove_blink', 1,...
    'rate_type', 'binned',...
    'bin_width', 0.1,...
    'norm_it', 1,...
    'pupil_type','raw',...
    'neural_type','raw',...
    'time_before',5, ...
    'time_after',5);





% ~~~~~~~~~~~~~~~~~~~~~~ DO NOT CHANGE BELOW ~~~~~~~~~~~~~~~~~~~~~~
% Figure configurations
fig_params = struct('Position',[0 0 12 9], 'Units','centimeters', 'WindowStyle','normal', 'FontSize', 10);
% 18 cm = max




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

plot_path  = fullfile(paths.dropbox,'plots',findPlotFolder,'1-draft');
if ~exist(plot_path,'dir') && save_it; mkdir(plot_path); end



%% Load data

exper_conds = {'test','tone alone','retest','extinction'}; % All conditions that are needed to complete this plot

%~~~~~~~~~~~~ Load table data ~~~~~~~~~~~~

condition_file = fullfile(paths.dropbox,'tables','Conditioning Recordings.xlsx');
condition_table = readtable(condition_file, 'Sheet', params.mouse_id,'PreserveVariableNames',true); %'VariableNamingRule','preserve')


% Find correct row
datetime_str = datetime(params.date_str,'InputFormat','M-d-yy');
date_inds = condition_table.Date == datetime_str;
file_inds = strcmp(condition_table.File,params.file_num);
condition_inds = any(cell2mat(cellfun(@(exper_cond)strcmpi(condition_table.Condition, exper_cond), exper_conds, 'un',0)),2);
exper_inds = date_inds & file_inds & condition_inds;
condition_table(~exper_inds,:) = [];

if min(size(condition_table)) ~= 4
    error('Did not find all conditions')
end

%~~~~~~~~~~~~ Load data ~~~~~~~~~~~~
date_str = strjoin(cellfun(@(num) num2str(str2double(num)),{datestr(unique(condition_table.Date),'mm'),datestr(unique(condition_table.Date),'dd'),datestr(unique(condition_table.Date),'yy')},'un',0),'-');
file_num = cell2mat(unique(condition_table.File));

data_params = struct('remove_blink',params.remove_blink,'rate_type',params.rate_type,'bin_width',params.bin_width,'norm_it',params.norm_it,'pupil_type',params.pupil_type,'neural_type',params.neural_type,'time_before',params.time_before,'time_after',params.time_after);

% % Initialize structure that hold onto data
% exper_data = struct();
% 
% for exper_cond = exper_conds
%     exper_cond_field = replace(exper_cond{:},' ','_');
%     exper_data.(exper_cond_field) = struct('pupil_data',struct(),'neural_data',struct(),'trigger_data',struct());
%     [exper_data.(exper_cond_field).pupil_data, exper_data.(exper_cond_field).neural_data, exper_data.(exper_cond_field).trigger_data] = loadData(params.mouse_id, date_str, file_num, exper_cond_field, paths, data_params);
% end



%% Initialize figure

fig = figure(9);
clf(fig)
fig.WindowStyle = fig_params.WindowStyle;
fig.Units = fig_params.Units;
fig.Position = fig_params.Position;

layout_main = tiledlayout(2,2,'TileSpacing','compact','Padding','none');



%% Plot A

plot_evoked(exper_data, 'tone_alone', 'test', 'neural', layout_main, 'A')



%% Plot B

plot_evoked(exper_data, 'tone_alone', 'test', 'pupil', layout_main, 'B')



%% Plot C

plot_evoked(exper_data, 'extinction', 'retest', 'neural', layout_main, 'C')



%% Plot D

plot_evoked(exper_data, 'extinction', 'retest', 'pupil', layout_main, 'D')



%% Save

sameLimits(findobj(fig,'Type','Axes'), 'x')
sameLimits(findobj(fig,'Type','Axes','UserData','neural'),'y')
sameLimits(findobj(fig,'Type','Axes','UserData','pupil'),'y')
uniformFormat_pub(fig)

if save_it
    exportgraphics(fig,fullfile(plot_path,['figure_',num2str(fig.Number),'.pdf']),'Resolution',300,'ContentType','vector')
end












%% ~~~~~~~~~~~~~~~~~~~~~~~~~~ Functions for figures ~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Plot evoked trials
function [] = plot_evoked(data, cond1, cond2, data_type, layout, tag)
% cond1 = black
% cond2 = color

step_plot = @(x) reshape([x;x],1,2*numel(x));

color = struct('neural','r','pupil','b');
legend_label = struct('tone_alone','Pre Conditioning',...
                      'test', 'Post Conditioning 1',...
                      'extinction', 'Extinction',...
                      'retest', 'Post Conditioning 2');

%~~~~~~~~~~~~ Prepare data for plotting ~~~~~~~~~~~~

% Collect data for each trigger for data.(cond1)
unbinned1 = struct('time',   data.(cond1).trigger_data.time_bins,...
                  'neural', NaN(data.(cond1).trigger_data.trigger_total, numel(data.(cond1).trigger_data.time_bins)),...
                  'pupil',  NaN(data.(cond1).trigger_data.trigger_total, numel(data.(cond1).trigger_data.time_bins)));
binned1 = struct('time',     step_plot(unbinned1.time),...
                'neural', NaN(data.(cond1).trigger_data.trigger_total, numel(step_plot(unbinned1.time))),...
                'pupil',  NaN(data.(cond1).trigger_data.trigger_total, numel(step_plot(unbinned1.time))));

for t = 1:data.(cond1).trigger_data.trigger_total

    % Collect neural data
    unbinned1.neural(t,:) = data.(cond1).trigger_data.neural.data(t,:);
    binned1.neural(t,:) = step_plot(unbinned1.neural(t,:));

    % Collect pupil data
    unbinned1.pupil(t,:) = data.(cond1).trigger_data.pupil.data(t,:);
    binned1.pupil(t,:) = step_plot(unbinned1.pupil(t,:));

end

data2_before = unbinned1.(data_type)(:,unbinned1.time < 0);
thresh1 = max(calculateCI(data2_before(:)'));



% Collect data for each trigger for data.(cond2)
unbinned2 = struct('time',   data.(cond2).trigger_data.time_bins,...
                  'neural', NaN(data.(cond2).trigger_data.trigger_total, numel(data.(cond2).trigger_data.time_bins)),...
                  'pupil',  NaN(data.(cond2).trigger_data.trigger_total, numel(data.(cond2).trigger_data.time_bins)));
binned2 = struct('time',     step_plot(unbinned2.time),...
                'neural', NaN(data.(cond2).trigger_data.trigger_total, numel(step_plot(unbinned2.time))),...
                'pupil',  NaN(data.(cond2).trigger_data.trigger_total, numel(step_plot(unbinned2.time))));

for t = 1:data.(cond2).trigger_data.trigger_total

    % Collect neural data
    unbinned2.neural(t,:) = data.(cond2).trigger_data.neural.data(t,:);
    binned2.neural(t,:) = step_plot(unbinned2.neural(t,:));

    % Collect pupil data
    unbinned2.pupil(t,:) = data.(cond2).trigger_data.pupil.data(t,:);
    binned2.pupil(t,:) = step_plot(unbinned2.pupil(t,:));

end

data2_before = unbinned2.(data_type)(:,unbinned2.time < 0);
thresh2 = max(calculateCI(data2_before(:)'));



%~~~~~~~~~~~~ Plot ~~~~~~~~~~~~
nexttile(layout)

hold on

% Plot cond1
plot_time = binned1.time;
plot_data = mean(binned1.(data_type),1,'omitnan');
plot_CR = calculateCR(binned1.(data_type)')';
H = shadedErrorBar(plot_time(2:end), plot_data(:,1:end-1), plot_CR(:,1:end-1),'lineProps',{'k'});
H.mainLine.DisplayName = legend_label.(cond1);
delete(H.edge)

% Plot cond2
plot_time = binned2.time;
plot_data = mean(binned2.(data_type),1,'omitnan');
plot_CR = calculateCR(binned2.(data_type)')';
H = shadedErrorBar(plot_time(2:end), plot_data(:,1:end-1), plot_CR(:,1:end-1),'lineProps',{color.(data_type)});
H.mainLine.DisplayName = legend_label.(cond2);
delete(H.edge)

% Plot other stuff
xline(0,'Tag','trigger')
yline(thresh1,'Color','k','Alpha',1,'LineStyle',':')
yline(thresh2,'Color',color.(data_type),'Alpha',1,'LineStyle',':')

hold off

%~~~~~~~~~~~~ Format ~~~~~~~~~~~~
axis tight

%~~~~~~~~~~~~ Label ~~~~~~~~~~~~
xlabel('Time (s)')
switch data_type
    case 'neural'
        ylabel('Firing Rate (Hz)')
    case 'pupil'
        ylabel('Normalized pupil area')
end
set(gca,'Tag',tag)
set(gca,'UserData',data_type)

lg = legend('Location','northwest');

end

