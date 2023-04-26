% plot_figure_8.m
%
% Script to makes figure 8 for publication, about single units
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Modified: 10-20-2022

%clear
clc
set(0,'DefaultFigureWindowStyle','normal') %docks the figures onto the Matlab window

save_it = 1; % Whether to save plots (1) or not (0)

step_plot = @(x) reshape([x;x],1,2*numel(x));


%% Set data


params=struct('ab',struct(),'c',struct(),'d',struct(),'e',struct());

% ~~~~~~~~~~~ Set mouse information ~~~~~~~~~~~

% from run_sponteneous.m
params.ab = struct('mouse_id',  'C7',...
    'date_str',  '7-1-2022',...
    'file_num', 'File1', ...
    'exper_cond', 'Laser', ...
    'remove_blink', 1,...
    'rate_type', 'binned',...
    'bin_width', 1,...
    'norm_it', 1,...
    'pupil_type','raw',...
    'neural_type','raw',...
    'thresh',10,...
    'time_before',300, ...
    'time_after',0);


% from run_laser_corr.m
params.c = struct('mouse_id',  'C2',...
    'date_str',  '7-26-2022',...
    'exper_cond', 'Laser', ...
    'file_num', 'File1', ...
    'remove_blink', 1,...
    'rate_type', 'binned',...
    'bin_width', 0.1,...
    'norm_it', 1,...
    'pupil_type','raw',...
    'neural_type','raw',...
    'time_before',5, ...
    'time_after',10);


params.d.column = 'Rmax';
params.e.column = 'Lag';





% ~~~~~~~~~~~~~~~~~~~~~~ DO NOT CHANGE BELOW ~~~~~~~~~~~~~~~~~~~~~~
% Figure configurations
fig_params = struct('Position',[0 0 18 9], 'Units','centimeters', 'WindowStyle','normal', 'FontSize', 10);
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



%% Initialize figure

fig = figure(8);
clf(fig)
fig.WindowStyle = fig_params.WindowStyle;
fig.Units = fig_params.Units;
fig.Position = fig_params.Position;

layout_main = tiledlayout(2,3,'TileSpacing','normal','Padding','none');


% Prepare figure A
tile_a = nexttile(layout_main,[1,2]); axis off;
layout_a = tiledlayout(layout_main, 1, 1,'TileSpacing','normal','Padding','none');
layout_a.Layout.Tile = 1;
layout_a.Layout.TileSpan = [1,2];

% Prepare figure B
tile_b = nexttile(layout_main,1); axis off;
layout_b = tiledlayout(layout_main, 1, 1,'TileSpacing','normal','Padding','none');
layout_b.Layout.Tile = 3;

% Prepare figure C
tile_c = nexttile(layout_main,[1,2]); axis off;
layout_c = tiledlayout(layout_main, 1, 1,'TileSpacing','normal','Padding','none');
layout_c.Layout.Tile = 4;
layout_c.Layout.TileSpan = [1,2];

% Prepare figure D+E
tile_de = nexttile(layout_main,1); axis off;
layout_de = tiledlayout(layout_main, 1, 2,'TileSpacing','normal','Padding','none');
layout_de.Layout.Tile = 6;


% Initialize for later
% if save_it
%     data = struct('ab',struct(),'c',struct());
% end



%% A

%~~~~~~~~~~~~ Load data ~~~~~~~~~~~~

if ~any(strcmp(fieldnames(data),'ab')) || (numel(fieldnames(data.ab)) == 0)
    data.ab = loadConditionData(params.ab, paths);
end


%~~~~~~~~~~~~ Prepare data ~~~~~~~~~~~~

% Collect data for each trigger
unbinned = struct('time',   sort(abs(data.ab.trigger_data.time_bins)),...
                  'neural', data.ab.trigger_data.neural.data(1,:),...
                  'pupil',  data.ab.trigger_data.pupil.data(1,:));
binned = struct('time',     step_plot(unbinned.time),...
                'neural',   step_plot(unbinned.neural),...
                'pupil',    step_plot(unbinned.pupil));


%~~~~~~~~~~~~ Plot ~~~~~~~~~~~~
nexttile(layout_a)
set(gca,'UserData','keep_color')

hold on

%%%%% Plot neural %%%%%
yyaxis left
plot(binned.time(2:end), binned.neural(1:end-1),'-','Color','r')


%%%%% Plot pupil %%%%%
yyaxis right
plot(binned.time(2:end), binned.pupil(1:end-1),'-','Color','b')


hold off


%~~~~~~~~~~~~ Format ~~~~~~~~~~~~
xtickangle(0)


%~~~~~~~~~~~~ Label ~~~~~~~~~~~~
xlabel('Time (s)')

yyaxis left
set(gca,'YColor','r')
ylabel('Spike rate (Hz)')
y_lim = get(gca,'YLim');
y_tick = get(gca,'YTick');

yyaxis right
set(gca,'YColor','b')
if params.ab.norm_it
    ylabel('Normalized pupil area')
else
    ylabel('Pupil area (pixels^2)')
end
x_lim = get(gca,'YLim');
x_tick = get(gca,'YTick');

title('Spontaneous activity')
set(gca,'Tag','A')



%% B

%~~~~~~~~~~~~ Load data ~~~~~~~~~~~~

if ~any(strcmp(fieldnames(data),'ab')) || (numel(fieldnames(data.ab)) == 0)
    data.ab = loadConditionData(params.ab, paths);
end


%~~~~~~~~~~~~ Prepare data ~~~~~~~~~~~~

% Collect data for each trigger
unbinned = struct('time',   data.ab.trigger_data.time_bins,...
                  'neural', data.ab.trigger_data.neural.data(1,:),...
                  'pupil',  data.ab.trigger_data.pupil.data(1,:));
binned = struct('time',     step_plot(unbinned.time),...
                'neural',   step_plot(unbinned.neural),...
                'pupil',    step_plot(unbinned.pupil));

x = unbinned.pupil;
y = unbinned.neural;

[r,p] = corr(x',y','Type','Pearson','Rows','Pairwise');



%~~~~~~~~~~~~ Plot ~~~~~~~~~~~~
nexttile(layout_b)

plot(x,y,'oy')

ls = lsline;


%~~~~~~~~~~~~ Format ~~~~~~~~~~~~
axis square

set(gca,{'XLim','XTick','XColor'},{x_lim,x_tick,'b'})
set(gca,{'YLim','YTick','YColor'},{y_lim,y_tick,'r'})


%~~~~~~~~~~~~ Label ~~~~~~~~~~~~

ylabel('Firing Rate (Hz)')

if params.ab.norm_it
    xlabel('Normalized pupil area')
else
    xlabel('Pupil area (pixels^2)')
end

% Add correlation to plot
txt = text(0.05,0.95, {['r=',num2str(round(r,3))];['p=',num2str(p,'%.2e')]}, 'HorizontalAlignment','left','VerticalAlignment','top','Units','normalized');


set(gca,'Tag','B')
set(gca,'UserData','keep_color')



%% C

%~~~~~~~~~~~~ Load data ~~~~~~~~~~~~

if ~any(strcmp(fieldnames(data),'c')) || (numel(fieldnames(data.c)) == 0)
    data.c = loadConditionData(params.c, paths);
end


%~~~~~~~~~~~~ Prepare data ~~~~~~~~~~~~

% Collect data for each trigger
unbinned = struct('time',   data.c.trigger_data.time_bins,...
                  'neural', NaN(data.c.trigger_data.trigger_total, numel(data.c.trigger_data.time_bins)),...
                  'pupil',  NaN(data.c.trigger_data.trigger_total, numel(data.c.trigger_data.time_bins)));
binned = struct('time',     step_plot(unbinned.time),...
                'neural', NaN(data.c.trigger_data.trigger_total, numel(step_plot(unbinned.time))),...
                'pupil',  NaN(data.c.trigger_data.trigger_total, numel(step_plot(unbinned.time))));

for t = 1:data.c.trigger_data.trigger_total

    % Collect neural data
    unbinned.neural(t,:) = data.c.trigger_data.neural.data(t,:);
    binned.neural(t,:) = step_plot(unbinned.neural(t,:));

    % Collect pupil data
    unbinned.pupil(t,:) = data.c.trigger_data.pupil.data(t,:);
    binned.pupil(t,:) = step_plot(unbinned.pupil(t,:));

end


%~~~~~~~~~~~~ Plot ~~~~~~~~~~~~
nexttile(layout_c)
set(gca,'UserData','keep_color')

hold on
%%%%% Plot binned neural %%%%%
yyaxis left
plot_mean = mean(binned.neural,1,'omitnan');
plot_CR = calculateCR(binned.neural')';
shadedErrorBar(binned.time(2:end), plot_mean(:,1:end-1), plot_CR(:,1:end-1),'lineProps',{'r'})

%%%%% Plot binned pupil %%%%%
yyaxis right
plot_mean = mean(binned.pupil,1,'omitnan');
plot_CR = calculateCR(binned.pupil')';
shadedErrorBar(binned.time(2:end), plot_mean(:,1:end-1), plot_CR(:,1:end-1),'lineProps',{'b'})

xline(0,'Tag','trigger')
hold off


%~~~~~~~~~~~~ Format ~~~~~~~~~~~~
box off
xlim([-data.c.trigger_data.time_before, data.c.trigger_data.time_after])
xticks(unique([-data.c.trigger_data.time_before:5:0,0:5:data.c.trigger_data.time_after]))
xtickangle(0)


%~~~~~~~~~~~~ Label ~~~~~~~~~~~~
xlabel('Time (s)')

yyaxis left
set(gca,'YColor','r')
ylabel('Firing Rate (Hz)')

yyaxis right
set(gca,'YColor','b')
if params.c.norm_it
    ylabel('Normalized pupil area')
else
    ylabel('Pupil area (pixels^2)')
end

title('Heat evoked')
set(gca,'Tag','C')



%% D

%~~~~~~~~~~~~ Load data ~~~~~~~~~~~~
table_file = fullfile(paths.dropbox,'tables','laser_corr_values.xlsx');
corr_table = readtable(table_file,'ReadRowNames',1,'PreserveVariableNames',1);

%~~~~~~~~~~~~ Prepare data ~~~~~~~~~~~~

%plot_ind = ~contains(corr_table{:,'Significant?'},'No');
plot_data = corr_table{:,params.d.column};


%~~~~~~~~~~~~ Plot ~~~~~~~~~~~~
nexttile(layout_de)

setNotBoxPlotTags(notBoxPlot(plot_data,'jitter',0.5));
yline(0,'Tag','thresh')


%~~~~~~~~~~~~ Format ~~~~~~~~~~~~
xlim([0,2])
xticks([])
set(gca,'XColor','none')


%~~~~~~~~~~~~ Label ~~~~~~~~~~~~
ylabel('Correlation value')
set(gca,'Tag','D')



%% D

%~~~~~~~~~~~~ Load data ~~~~~~~~~~~~
table_file = fullfile(paths.dropbox,'tables','laser_corr_values.xlsx');
corr_table = readtable(table_file,'ReadRowNames',1,'PreserveVariableNames',1);

%~~~~~~~~~~~~ Prepare data ~~~~~~~~~~~~

%plot_ind = ~contains(corr_table{:,'Significant?'},'No');
plot_data = corr_table{:,params.e.column};


%~~~~~~~~~~~~ Plot ~~~~~~~~~~~~
nexttile(layout_de)

setNotBoxPlotTags(notBoxPlot(plot_data,'jitter',0.5));
yline(0,'Tag','thresh')


%~~~~~~~~~~~~ Format ~~~~~~~~~~~~
xlim([0,2])
xticks([])
set(gca,'XColor','none')


%~~~~~~~~~~~~ Label ~~~~~~~~~~~~
ylabel('Lag (s)')
set(gca,'Tag','E')



%% Save

uniformFormat_pub(fig)

if save_it
    exportgraphics(fig,fullfile(plot_path,['figure_',num2str(fig.Number),'.pdf']),'Resolution',300,'ContentType','vector')
end












%% ~~~~~~~~~~~~~~~~~~~~~~~~~~ Functions for figures ~~~~~~~~~~~~~~~~~~~~~~~~~~


%% Load condition data
function data_struct = loadConditionData(params, paths)

%~~~~~~~~~~~~ Load table data ~~~~~~~~~~~~

condition_file = fullfile(paths.dropbox,'tables','Conditioning Recordings.xlsx');
condition_table = readtable(condition_file, 'Sheet', params.mouse_id,'PreserveVariableNames',true); %'VariableNamingRule','preserve')

% Find correct row
datetime_str = datetime(params.date_str,'InputFormat','M-d-yy');
date_inds = condition_table.Date == datetime_str;
file_inds = strcmp(condition_table.File,params.file_num);
condition_inds = strcmpi(condition_table.Condition, params.exper_cond);
exper_inds = date_inds & file_inds & condition_inds;
condition_table(~exper_inds,:) = [];

if min(size(condition_table)) ~= 1
    error('Found 0 or more than 1 condition in table')
end

%~~~~~~~~~~~~ Load data ~~~~~~~~~~~~
date_str = strjoin(cellfun(@(num) num2str(str2double(num)),{datestr(condition_table.Date,'mm'),datestr(condition_table.Date,'dd'),datestr(condition_table.Date,'yy')},'un',0),'-');
file_num = condition_table.File{1};
exper_cond = condition_table.Condition{1};

data_params = struct('remove_blink',params.remove_blink,'rate_type',params.rate_type,'bin_width',params.bin_width,'norm_it',params.norm_it,'pupil_type',params.pupil_type,'neural_type',params.neural_type,'time_before',params.time_before,'time_after',params.time_after);

[pupil_data, neural_data, trigger_data] = loadData(params.mouse_id, date_str, file_num, exper_cond, paths, data_params);

data_struct.pupil_data = pupil_data;
data_struct.neural_data = neural_data;
data_struct.trigger_data = trigger_data;

end


%% Sets tags for boxplot for formatting later
function H = setNotBoxPlotTags(H)

    for field = reshape(fieldnames(H),1,[])
        if isgraphics(H.(field{:}))
            H.(field{:}).Tag = field{:};
        end
    end

end

