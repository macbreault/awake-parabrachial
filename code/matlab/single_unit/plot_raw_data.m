% plot_raw_data.m
%
% Script that plots raw video and neural data for comparison
%
% Among other things, it mainly uses pupil_data and neural_data that are outputs of loadData
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 08-05-2022
% Modified: 08-18-2022 - v1 - Plots raw neural data as steps and plots spikes

% Inialize figure
fig = figure(1000);
clf(fig)
fig.WindowStyle = fig_params.WindowStyle;
fig.Units = fig_params.Units;
fig.Position = fig_params.Position;
fig.Position(3) = fig.Position(3)*2;

layout = tiledlayout(2,1,'TileSpacing','compact','Padding','none');

title(layout, exper_title)



% ~~~~~~~~~~ Plot 1: Pupil ~~~~~~~~~~

% Plot
nexttile(layout)
hold on
h1 = plot(pupil_data.time, pupil_data.raw,'-','Color','b','DisplayName','Raw pupil','LineWidth',0.5);
h2 = plot(pupil_data.time, pupil_data.filter,'--','Color','c','DisplayName','Filter pupil','LineWidth',0.5);
h3 = arrayfun(@(t) xline(pupil_data.trigger.times(t,1),'Tag','trigger'),1:pupil_data.trigger.total);
hold off

% Format
axis tight
box on
xticks(sort([0, round(pupil_data.trigger.times(:,1))']))

% Label
xlabel('Time relative to starting video (s)')
if norm_it
    ylabel('Norm pupil area (pixels^2)')
else
    ylabel('Pupil area (pixels^2)')
end
legend([h1,h2],'location','eastoutside')



% ~~~~~~~~~~ Plot neural data (hz) ~~~~~~~~~~
nexttile(layout)

plot_time = reshape([neural_data.time';neural_data.time'],1,2*numel(neural_data.time));
plot_neural = reshape([neural_data.raw';neural_data.raw'],1,2*numel(neural_data.raw));

% Plot
hold on
h1 = plot(plot_time(2:end), plot_neural(1:end-1),'-','Color','r','DisplayName','Raw neural','LineWidth',0.5);
h2 = plot(neural_data.time + (params.bin_width/2), neural_data.filter,'--','Color','m','DisplayName','Filter neural','LineWidth',0.5);
if any(ismember(fieldnames(neural_data),'spike'))
    h3 = plot(neural_data.spike, ones(size(neural_data.spike)),'|','MarkerSize',2.5,'Color','k','DisplayName','Spikes');
else
    h3 = [];
end
h4 = arrayfun(@(t) xline(neural_data.trigger.time(t),'Tag','trigger'),1:neural_data.trigger.total);
hold off

% Format
axis tight
box on
xticks(sort([0, round(neural_data.trigger.time)']))

% Label
xlabel('Time relative to starting neural (s)')
ylabel('Spike rate (Hz)')
legend([h1,h2,h3],'location','eastoutside')


% ~~~~~~~~~~ Uniformly format plots ~~~~~~~~~~
uniformFormat(fig)


if save_it
    savefig(fig,fullfile(plot_path,'0_raw.fig'))
    exportgraphics(fig,fullfile(plot_path,'0_raw.pdf'),'Resolution',300,'ContentType','vector')
end




% %% Plotting as steps
% 
% clf
% plot_time = reshape([neural_data.time';neural_data.time'],1,2*numel(neural_data.time));
% plot_neural = reshape([neural_data.raw';neural_data.raw'],1,2*numel(neural_data.raw));
% 
% plot_time = plot_time(1:100);
% plot_neural = plot_neural(1:100);
% 
% plot(plot_time(2:end),plot_neural(1:end-1),'-')


% %% Plotting spikes
% 
% clf
% hold on
% stem(neural_data.spike, ones(size(neural_data.spike)),'Marker','none','Color','r')
% h3 = arrayfun(@(t) xline(neural_data.trigger.time(t),'LineWidth',2,'DisplayName','Trigger neural start','Color',[0.5 0.5 0.5]),1:neural_data.trigger.total);