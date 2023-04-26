% plot_split_data.m
%
% Script that plots raw video and neural data split into triggers for comparison
%
% Among other things, it mainly uses trigger_data that is an output of loadData
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 08-05-2022
% Modified: 08-18-2022 - v1 - Plots raw neural data as steps and plots spikes

step_plot = @(x) reshape([x;x],1,2*numel(x));


% Inialize figure
fig = figure(10);
clf(fig)
fig.WindowStyle = fig_params.WindowStyle;
fig.Units = fig_params.Units;
fig.Position = fig_params.Position;
fig.Position(3) = fig.Position(3)*1.5;

% Get optimal number of subplots
layout = tiledlayout(1,4,'TileSpacing','compact','Padding','compact');



%% Collect data for each trigger

unbinned = struct('time',   trigger_data.time_bins,...
                  'neural', NaN(trigger_data.trigger_total, numel(trigger_data.time_bins)),...
                  'pupil',  NaN(trigger_data.trigger_total, numel(trigger_data.time_bins)));
binned = struct('time',     step_plot(unbinned.time),...
                'neural', NaN(trigger_data.trigger_total, numel(step_plot(unbinned.time))),...
                'pupil',  NaN(trigger_data.trigger_total, numel(step_plot(unbinned.time))));


for t = 1:trigger_data.trigger_total

    % Collect neural data
    unbinned.neural(t,:) = trigger_data.neural.data(t,:);
    binned.neural(t,:) = step_plot(unbinned.neural(t,:));

    % Collect pupil data
    unbinned.pupil(t,:) = trigger_data.pupil.data(t,:);
    binned.pupil(t,:) = step_plot(unbinned.pupil(t,:));

end



%% Plot average

nexttile(layout,[1,3])

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

% Format
box off
xlim([-trigger_data.time_before, trigger_data.time_after])
xticks(unique([-trigger_data.time_before:1:0,0:1:trigger_data.time_after]))

% Label
title(['Number of triggers = ',num2str(trigger_data.trigger_total)])
xtickangle(0)

xlabel('Time since trigger (s)')

yyaxis left
set(gca,'YColor','r')
ylabel('Spike rate (Hz)')


yyaxis right
set(gca,'YColor','b')
if norm_it
    ylabel('Normalized pupil area (pixels^2)')
else
    ylabel('Pupil area (pixels^2)')
end



% ~~~~~~~~~~ Uniformly format plots ~~~~~~~~~~

yyaxis left
uniformFormat(fig)
yyaxis right
uniformFormat(fig)



%% Plot correlation

% Get data
time_ind = (0 < unbinned.time) & (unbinned.time <= time_after); % BOUNDARY BETWEEN 0 AND time_after SECONDS
mean_neural = mean(unbinned.neural(:,time_ind), 1,'omitnan')';
mean_pupil = mean(unbinned.pupil(:,time_ind),  1,'omitnan')';
[r,p] = corr(mean_neural,mean_pupil,'type','Pearson','rows','pairwise');

% Plot
nexttile(layout)
scatter(mean_pupil,mean_neural, 'o','MarkerEdgeColor','k','MarkerFaceColor','none','LineWidth',1,'MarkerFaceAlpha',0.25,'DisplayName',['r=',num2str(r)])
ls = lsline;
ls.DisplayName = ['p=',num2str(p)];

% Format
axis square padded
box on

% Label
title(['Correlated between (0,',num2str(time_after),'] s'])
if norm_it
    xlabel('Normalized pupil area (pixels^2)')
else
    xlabel('Pupil area (pixels^2)')
end
ylabel('Spike rate (Hz)')
lg = legend('location','southoutside','UserData','ignore_text');


% Edit for significant correlations
if p <= 0.05
    ls.Color = 'm';
    lg.FontWeight = 'bold';
else
    ls.Color = [0.5 0.5 0.5];
    lg.FontWeight = 'normal';
end



%% Add cross-correlation value to title

time_ind = (0 < unbinned.time) & (unbinned.time <= time_after);
mean_neural = mean(trigger_data.neural.data(:,time_ind), 1,'omitnan')';
mean_pupil = mean(trigger_data.pupil.data(:,time_ind),  1,'omitnan')';
[r_max, lag_max, p_max] = calculateCrossCorr(mean_neural,mean_pupil,bin_width);
lag_max_time = lag_max*bin_width;

title(layout, {exper_title; ...
    ['\bfr_{max} = ',num2str(round(r_max,3)),', lag_{max} = ',num2str(round(lag_max_time,3)),' s']}) % Add to title





%% Save it

uniformFormat(fig)
if save_it
    savefig(fig,fullfile(plot_path,'1_average.fig'))
    exportgraphics(fig,fullfile(plot_path,'1_average.pdf'),'Resolution',300,'ContentType','vector')
end
