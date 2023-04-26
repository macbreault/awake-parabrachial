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
fig = figure(1);
clf(fig)
fig.WindowStyle = fig_params.WindowStyle;
fig.Units = fig_params.Units;
fig.Position = fig_params.Position;
fig.Position(3) = fig.Position(3)*2;

% Get optimal number of subplots
P = numSubplots(trigger_data.trigger_total);

layout = tiledlayout(P(1),P(2),'TileSpacing','compact','Padding','none');



%% Collect data for plotting

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



%% Plot for each trigger

cross_corr = struct('r',NaN(trigger_data.trigger_total,1), 'lag',NaN(trigger_data.trigger_total,1));

for t = 1:trigger_data.trigger_total

    nexttile(layout)

    % Shift time so spike and pupil line up

    hold on
    %%%%% Plot binned neural %%%%%
    yyaxis left
    plot(binned.time(2:end), binned.neural(t,1:end-1),'-','Color','r','DisplayName','Binned neural','LineWidth',1.5)

    %%%%% Plot spiking data %%%%%
    yyaxis left
    plot(trigger_data.neural.raw.time.time{t}, ones(size(trigger_data.neural.raw.time.time{t})),'|','MarkerSize',2.5,'Color','k','DisplayName','Spikes');
        
    %%%%% Plot raw pupil %%%%%
    yyaxis right
    plot(trigger_data.pupil.raw.time.time{t} - trigger_data.pupil.raw.time.trig(t),   trigger_data.pupil.raw.data{t},'-','Color','b','DisplayName','Raw pupil','LineWidth',0.5)

    %%%%% Plot binned pupil %%%%%
    yyaxis right
    plot(binned.time(2:end),   binned.pupil(t,1:end-1), '-','Color','b','DisplayName','Binned pupil','LineWidth',2)

    xline(0,'Tag','trigger')

    hold off

    % Format
    box off
    xlim([-trigger_data.time_before, trigger_data.time_after])
    xticks(unique([-trigger_data.time_before:1:0,0:1:trigger_data.time_after]))

    % Label
    title(['Trigger ',num2str(trigger_data.trigger_num(t))])
    xtickangle(0)

    [y,x] = ind2sub(flip(P),t);
    %title(['(y=',num2str(y),', x=',num2str(x),')'])

    if x == P(1)
        xlabel('Time since trigger (s)')
    else
        xticklabels('')
    end

    yyaxis right
    set(gca,'YColor','b')
    if y == P(2)
        if norm_it
            ylabel('Normalized pupil area (pixels^2)')
        else
            ylabel('Pupil area (pixels^2)')
        end
    else
        %yticklabels('')
    end

    yyaxis left
    set(gca,'YColor','r')
    if y == 1
        ylabel('Spike rate (Hz)')
    else
        %yticklabels('')
    end


    % Calculate cross-correlation
    time_ind = (0 < unbinned.time) & (unbinned.time <= time_after);
    x = unbinned.neural(t,time_ind);
    y = unbinned.pupil(t,time_ind);
    [r_max, lag_max, p_max] = calculateCrossCorr(x,y,bin_width);
    cross_corr.r(t) = r_max;
    cross_corr.lag(t) = lag_max*bin_width;

end

% Add cross-correlation to title
R_max = mean(cross_corr.r,'omitnan');
lag_max = mean(cross_corr.lag,'omitnan');
title(layout, {exper_title; ...
    ['\bfAverage r_{max} = ',num2str(round(R_max,3)),', Average lag_{max} = ',num2str(round(lag_max,3)),' s']}) % Add to title



% ~~~~~~~~~~ Uniformly format plots ~~~~~~~~~~
AX = findobj(fig,'Type','Axes');

% Set left and right axes to have same limits
for side = {'left','right'}
    arrayfun(@(i) yyaxis(AX(i),side{:}),1:numel(AX))
    sameLimits(AX)
end

uniformFormat(fig)

if save_it
    savefig(fig,fullfile(plot_path,'1_split.fig'))
    exportgraphics(fig,fullfile(plot_path,'1_split.pdf'),'Resolution',300,'ContentType','vector')
end
