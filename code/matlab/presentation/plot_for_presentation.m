% plot_for_presentation.m
%
% Script to loads and plots data for presentation
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Modified: 09-23-2022

clear
clc


%% Set data

% ~~~~~~~~~~~ Set mouse information ~~~~~~~~~~~
mouse_id = 'M3';


% Other variables
save_it = 0;


%%%%%% Set paths %%%%%%
paths = struct('dropbox',[]);

if contains(pwd,'/mac/')
    paths.dropbox = '/Users/mac/Dropbox (MIT)/Jesse';
else
    paths.dropbox = '/Users/Jessesmith 1/Dropbox/Jesse';
end

if ~exist(paths.dropbox,'dir'); error('Error in %s: Dropbox path not found', mfilename); end



%% Load data


behavior = readtable(fullfile(paths.dropbox,'tables','Photometry Behavior.xlsx'),'ReadRowNames',1,'PreserveVariableNames',1);
traces = readtable(fullfile(paths.dropbox,'tables','Photometry Traces.xlsx'),'Sheet','Sheet2','ReadRowNames',1,'PreserveVariableNames',1);


%%%%%% Set paths %%%%%%
paths = struct('dropbox',[]);

if contains(pwd,'/mac/')
    paths.dropbox = '/Users/mac/Dropbox (MIT)/Jesse';
else
    paths.dropbox = '/Users/Jessesmith 1/Dropbox/Jesse';
end 

addpath(genpath(fullfile(paths.dropbox,'code'))); % Add libraries to directory

plot_path  = fullfile(paths.dropbox,'plots',findPlotFolder,date,'other');
if ~exist(plot_path,'dir') && save_it; mkdir(plot_path); end


% Figure configurations
fig_params = struct('Position',[0 0 8 4], 'Units','Inches', 'WindowStyle','normal', 'FontSize', 10);



%% Plot

% Prepare figure
fig = figure(1);
clf(fig)
fig.WindowStyle = fig_params.WindowStyle;
fig.Units = fig_params.Units;
fig.Position = fig_params.Position;

layout = tiledlayout(1,2,'TileSpacing','normal','Padding','compact');
title(layout, {mouse_id,''})



%%%%% Plot trace %%%%%


% Get data
plot_time = cellfun(@str2num, traces.Properties.RowNames);
plot_data = traces{:,mouse_id};
plot_CI = calculateCI(plot_data(0 > plot_time)'); % Calculate CI using [-1,-2,-3]


% Plot
nexttile(layout)
plot(plot_time, plot_data,'g-o','MarkerFaceColor','g')
yline(plot_CI,'Tag','thresh')
xline(0,'Tag','trigger')


% Format
axis padded
box off
%xlim([min(plot_time),max(plot_time)])




% Label
title('Population CGRP Activity')
xlabel('Days Post CFA')
ylabel('Duration (s)')



%%%%% Plot behavior %%%%%

% Get data for plotting
plot_time = cellfun(@str2num, behavior.Properties.RowNames);
plot_data = behavior{:,mouse_id};

% Plot
nexttile(layout)
plot(plot_time, plot_data,'k-o','MarkerFaceColor','k')


% Format
axis padded
box off


% Label
title('Mechanical Withdrawal Threshold')
ylabel({'Mechanical Withdraw','Threshold (g)'})
xlabel('Days Post CFA')





% Save
uniformFormat(fig)
if save_it
    exportgraphics(fig,fullfile(plot_path,[mouse_id,'.pdf']),'Resolution',300,'ContentType','vector')
    savefig(fig,fullfile(plot_path,[mouse_id,'.fig']))
end


