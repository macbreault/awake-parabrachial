% plot_stuff.m
%
% Script that plots STUFF
% Gets correlations of pupil size and firing rate immediately before,
% during and after the stimulation
%
% Jesse Smith
% Created: 03-10-2022

error('OUTDATED')

clear
clc

p = genpath('/Users/Jessesmith 1/Desktop/facemap-main/JESSE FIG FILES/Library');
addpath(p);
set(0,'DefaultFigureWindowStyle','docked') %docks the figures onto the Matlab window
%% Load data

%%%%~~~~~NEED TO KNOW~~~~~
%Mouse ID
%Date
%File name
%video laser time points in seconds



% ~~~~~~~~~~~ Set mouse information ~~~~~~~~~~~
%%%%% CHANGE THIS STUFF ONLY %%%%%
mouse_id = 14;          % Mouse ID
date_str = '2-2-22';    % Date of recording
file_num = 'File1';     % File number

% ~~~~~~~~~~~ Load video laser trigger CHANGE THIS TOO ~~~~~~~~~~~
run('laser_video_triggers.m') % ADD LASER VIDEO TRIGGER INTO THIS SCRIPT FOR EACH CASE


% ~~~~~~~~~~~ Load pupil data ~~~~~~~~~~~
pupil_file = fullfile('/Volumes/Untitled/Recordings',['M',num2str(mouse_id)],date_str,file_num,['M',num2str(mouse_id),' ',date_str,' ',file_num,'_proc.mat']);
load(pupil_file,'proc') % pupil_file,'proc','-struct',pupil)

%~~~~~~~~~~~~ Convert pupil data from frames to seconds~~~~~~~~~~~~
pupil_data = double(proc.pupil(1).area);

%~~~~~~~~~~~~ Filter pupil data using moving average ~~~~~~~~~~~~
pupil_filter = smoothdata(pupil_data,'movmean',100); % https://www.mathworks.com/help/matlab/ref/smoothdata.html#bvhejau-method


%~~~~~~~~~~~~ Convert frames to seconds for pupil~~~~~~~~~
%pupil_time =  (1:proc.nframes)./25;
%plot(pupil_time,pupil_data,'-')

% ~~~~~~~~~~~ Load spike data ~~~~~~~~~~~
neural_file = fullfile('/Volumes/Seagate Backup Plus Drive/Electrophysiology/Awake Recordings/Mouse Recordings/Neurotar',...
    ['Mouse ',num2str(mouse_id)],date_str,'Sorted',file_num,'MATLAB','spike.mat');
neural_struct = load(neural_file);

% Extract spike data (regardless of the name of the spike file)
neural_name = fieldnames(neural_struct);

if numel(neural_name) == 1
    spike_data = neural_struct.(neural_name{:});
else
    error('Error in %s: Found more than 1 spike file. Check!',mfilename)
end


% ~~~~~~~~~~~ Calculate spike rate ~~~~~~~~~~~
bin_width = 1; % Seconds %%%%%%%%%%%%%%%%% CHANGE IF YOU WANT? %%%%%%%%%%%%%%%%%
T = ceil(spike_data(end)); % Time of last spike
rate_time = 0:bin_width:T;
neural_data = NaN(1,length(rate_time)); % Initalize vector that keeps track of how many spikes are in each 1 second bin
for t = 1:length(rate_time)
    neural_data(t) = sum(  (rate_time(t) <= spike_data) & (spike_data < (rate_time(t)+bin_width)) ) / bin_width; % Convert to per second (aka Hz)
end


% ~~~~~~~~~~~ Filter spike data ~~~~~~~~~~~
neural_filter = smoothdata(neural_data,'movmean',5);




% ~~~~~~~~~~~ Load neural laser trigger ~~~~~~~~~~~
laser_neural_file = fullfile('/Volumes/Seagate Backup Plus Drive/Electrophysiology/Awake Recordings/Mouse Recordings/Neurotar', ...
    ['Mouse ' num2str(mouse_id)],date_str,'Sorted',file_num,'MATLAB','event.mat');
laser_neural_struct = load(laser_neural_file);

% Extract neural video data (regardless of the name of the event file)
laser_neural_name = fieldnames(laser_neural_struct);

if numel(laser_neural_name) == 1
    laser_neural_data = laser_neural_struct.(laser_neural_name{:});
else
    error('Error in %s: Found more than 1 laser neural file. Check!',mfilename)
end


%% Plot raw signals

% HYPOTHETICAL DATA (REMOVE LATER)
%time = 1:10;
%data = rand(1,10);


% Inialize figure
fig = figure(1);
clf
layout = tiledlayout(2,1,'TileSpacing','compact');

title(layout,['Mouse ',num2str(mouse_id),', ',date_str,', ',file_num])



% ~~~~~~~~~~ Plot 1: Pupil ~~~~~~~~~~
nexttile(layout)

% Plot
%~~~~~~~~~~~~ Convert frames to seconds for pupil~~~~~~~~~
pupil_time =  (0:(double(proc.nframes)-1))./ 25;
hold on
plot(pupil_time, pupil_filter,'--','Color','r','LineWidth',1,'DisplayName','Filter pupil')
plot(pupil_time, pupil_data,'-','Color','b','DisplayName','Raw pupil')
stem(laser_video, max(pupil_data).*ones(size(laser_video)),'k-','LineWidth',1,'MarkerEdgeColor','none','DisplayName','Laser video')
hold off

% Format
axis tight
box on
set(gca, 'Children',flipud(get(gca, 'Children'))) % Re-order so trigger is behind signal


% Label
xlabel('Time (s)')
ylabel('Pupil area (pixels^2)')
legend('location','bestoutside')



% ~~~~~~~~~~ Plot neural data (hz) ~~~~~~~~~~
nexttile(layout)

% % Check to make sure all spikes are counted for in spike_rate
% if not((sum(spike_rate)*bin_width) == numel(spike_data))
%     error('Error in %s: Number of spikes from rate does not match spike_data. Check',mfilename)
% end

% Plot
hold on
plot(rate_time, neural_filter,'--','Color','r','LineWidth',1,'DisplayName','Filter neural')
plot(rate_time, neural_data,'-','Color','b','DisplayName','Raw neural')
stem(laser_neural_data, max(neural_data).*ones(size(laser_neural_data)),'-k','LineWidth',1,'MarkerEdgeColor','none','DisplayName','Laser neural')
hold off

% Format
axis tight
box on
set(gca, 'Children',flipud(get(gca, 'Children'))) % Re-order so trigger is behind signal

% Label
xlabel('Time (s)')
ylabel('Spike rate (Hz)')
legend('location','bestoutside')


% ~~~~~~~~~~ Uniformly format plots ~~~~~~~~~~
AX = findobj(fig,'Type','Axes');
set(AX,'XLim',[min(cell2mat(get(AX,'XLim')),[],'all'), max(cell2mat(get(AX,'XLim')),[],'all')]) % Make X limits the same for both axes
set(findobj(fig,'-property','FontSize'),'FontSize',15)





%% Plot correlation plot

%%%% Things to change %%%%
pupil_type = 'raw'; % 'raw' or 'filter'
neural_type = 'raw'; % 'raw' or 'filter'
T = 5; % Trigger number to plot in example

time_before = 5; % In seconds
time_after = 30; % mean([diff(laser_video); diff(laser_neural_data)]); % In seconds




%%%% STOP, don't change below here! %%%%
% Choose which pupil data to use for correlating
pupil_2_corr = [];
switch pupil_type
    case 'raw'
        pupil_2_corr = pupil_data;
    case 'filter'
        pupil_2_corr = pupil_filter;
    otherwise
        error('Error in %s: Did not recognize type of pupil',mfilename)
end

neural_2_corr = [];
switch neural_type
    case 'raw'
        neural_2_corr = neural_data;
    case 'filter'
        neural_2_corr = neural_filter;
    otherwise
        error('Error in %s: Did not recognize type of neural',mfilename)
end

% Time lock using respective triggers

trigger_num = numel(laser_video);

pupil2corr = [];
neural2corr = [];
pupil2corr_trail = cell(trigger_num,1);
neural2corr_trial = cell(trigger_num,1);
pupil_start = NaN(trigger_num,1);
pupil_end = NaN(trigger_num,1);
neural_start = NaN(trigger_num,1);
neural_end = NaN(trigger_num,1);
neural_ind = cell(trigger_num,1);
pupil_ind = cell(trigger_num,1);

trigger_loop = 1:trigger_num; % Triggers to actually loop through

for t = trigger_loop
    
    % Find first and last index of trigger into pupil data
    [~,pupil_start(t)] = min(abs(pupil_time - (laser_video(t)-time_before)));
    [~,pupil_end(t)] = min(abs(pupil_time - (laser_video(t)+time_after)));
    
    % Find first and last index of trigger into neural data
    [~,neural_start(t)] = min(abs(rate_time - (laser_neural_data(t)-time_before)));
    [~,neural_end(t)] = min(abs(rate_time - (laser_neural_data(t)+time_after)));
    
    % Make index into pupil_data and neural_data
    neural_ind{t} = neural_start(t):neural_end(t);
    pupil_ind{t}  = round( linspace(pupil_start(t), pupil_end(t), numel(neural_ind{t})) );
    
    % Extract the data into cell arrays. Trials may not be the same size,
    % but that will be fixed later
    pupil2corr_trail{t} = pupil_2_corr( pupil_ind{t} );
    neural2corr_trial{t} = neural_2_corr( neural_ind{t} )';
    
end

%%%%% CHOP OFF the end of the trials so all trials are the same size
min_bins = min(cellfun(@length, [pupil2corr_trail; neural2corr_trial])); % Size all trials should be

% Initialize matrices
pupil2corr = NaN(min_bins, trigger_num);
neural2corr = NaN(min_bins, trigger_num);

% Fill it in!
for t = trigger_loop
    pupil2corr(:,t) = pupil2corr_trail{t}(1:min_bins);
    neural2corr(:,t) = neural2corr_trial{t}(1:min_bins);
end


% Correlate
[r,p] = corr(pupil2corr(:), neural2corr(:),'Type','Pearson');


% Start figure
fig = figure(2);
clf


% ~~~~~~~~~~~~~~~~~~~ Plot scatter plot of correlation between pupil and spikes ~~~~~~~~~~~~~~~~~~~
subplot(2,1,1)
hold on
plot(pupil2corr(:), neural2corr(:),'o','DisplayName',['r=',num2str(r)])
ls = lsline; % Make least square line
ls.DisplayName = ['p=',num2str(p)];
if p<=0.05
    ls.Color = 'r'; % If significant, then lsline turns red
    ls.LineWidth = 2;
else
    ls.Color = 'k';
end
%%%%% Plot data used for trigger example on correlation plot %%%%%
plot(pupil2corr(:,T),neural2corr(:,T),'o','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName',['Trigger ',num2str(T)])
hold off

% Format
axis square
box on

% Label
legend('location','eastoutside')
xlabel(['\bf{',upper(pupil_type(1)),lower(pupil_type(2:end)),'} \rm{pupil area (pixels^2)}'])
ylabel(['\bf{',upper(neural_type(1)),lower(neural_type(2:end)),'} \rm{spiking rate (Hz)}'])
title(['Mouse ',num2str(mouse_id),', ',date_str,', ',file_num])



% ~~~~~~~~~~~~~~~~~~~ Plot example of raw data used for correlation ~~~~~~~~~~~~~~~~~~~
subplot(2,1,2)

% Shift time so spike and pupil line up
binned_time = linspace(-time_before, time_after, min_bins);
full_time = linspace(-time_before, time_after, numel(pupil_start(T):pupil_end(T)));

% Plot
hold on

%%%%% Plot raw neural %%%%%
yyaxis left
plot(binned_time,neural_data(neural_start(T):neural_end(T)),'-','Color','r','DisplayName','Raw neural')
plot(binned_time, neural_data(neural_ind{T}), '*', 'Color','r','DisplayName','Raw neural ind')

%%%%% Plot filtet neural %%%%%
yyaxis left
plot(binned_time,neural_filter(neural_start(T):neural_end(T)),'-','Color','m','DisplayName','Filter neural')
plot(binned_time, neural_filter(neural_ind{T}), '*', 'Color','m','DisplayName','Filter neural ind')

%%%%% Plot raw pupil %%%%%
yyaxis right
plot(full_time,pupil_data(pupil_start(T):pupil_end(T)),'-','Color','b','DisplayName','Raw pupil')
plot(binned_time,pupil_data(pupil_ind{T}), '*', 'Color','b','DisplayName','Raw pupil ind')

%%%%% Plot filter pupil %%%%%
plot(full_time,pupil_filter(pupil_start(T):pupil_end(T)),'-','Color','c','DisplayName','Filter pupil')
plot(binned_time,pupil_filter(pupil_ind{T}), '*', 'Color','c','DisplayName','Filter pupil ind')


hold off

% Format
box on
xlim([-time_before, time_after])
xticks(unique([-time_before,0:10:time_after]))

% Label
legend('location','eastoutside')
title(['Trigger ',num2str(T)])
xlabel('Time since trigger [s]')

yyaxis right
set(gca,'YColor','b')
ylabel('Pupil area (pixels^2)')

yyaxis left
set(gca,'YColor','r')
ylabel('Spike rate (Hz)')

set(findobj(fig,'-property','FontSize'),'FontSize',15)


%% %%%%%%%%% CROSS CORRELATION%%%%%%%%%
fig = figure (3);

% Gets cross-corr value
[R,lag] = crosscorr(pupil2corr(:,T), neural2corr(:,T)); %https://www.mathworks.com/help/econ/crosscorr.html
[R_max,ind_max] = max(R);
lag_max = lag(ind_max);

% Plots correlation value
crosscorr(pupil2corr(:,T), neural2corr(:,T)); %https://www.mathworks.com/help/econ/crosscorr.html



% %plot(lag21,R,[t21 t21],[-0.5 1],'r:')
% stem(lag,R)
% xline(lag_max,'r:')
% text(lag_max+10,0.5,['\bf{Lag: ' int2str(lag_max),'}'])
%
% ylabel('Correlation value (R)')
xlabel({['Lag (index = ',num2str(bin_width),' sec)'],['Lag < 0 = Neural before pupil                                 '...
    'Lag > 0 = Pupil before neural']})
% axis tight
% title('Cross-Correlations')

% Add legend with lag value
legend(['\bf{Lag: ' int2str(lag_max),'}'])

set(findobj(fig,'-property','FontSize'),'FontSize',15)






%%%%% Plot data used to run cross-correlation %%%%%

fig = figure(4);
clf

time_concat = repmat(binned_time,1,numel(trigger_loop));

% Pupil data
subplot(2,1,1)
hold on
plot(pupil2corr(:))
arrayfun(@(x) xline(x,'k-','LineWidth',2),find(time_concat == 0))
arrayfun(@(x) xline(x,'k:'),find(time_concat == min(time_concat))-0.5)
hold off
axis tight
title('Pupil data')
xlabel(['Index = ',num2str(bin_width),' sec'])
ylabel('Pupil area (pixels^2)')


% Neural data
subplot(2,1,2)
hold on
plot(neural2corr(:))
arrayfun(@(x) xline(x,'k-','LineWidth',2),find(time_concat == 0))
arrayfun(@(x) xline(x,'k:'),find(time_concat == min(time_concat))-0.5)
hold off
axis tight
title('Neural data')
xlabel(['Index = ',num2str(bin_width),' sec'])
ylabel('Firing rate (Hz)')

set(findobj(fig,'-property','FontSize'),'FontSize',15)




%%%%% Plot data SHIFTED by MAX cross-correlation compared to LAG = 0%%%%%

fig = figure(5);
clf

time_concat = repmat(binned_time,1,trigger_num);
pupil2corr_concat = pupil2corr(:);
neural2corr_concat = neural2corr(:);


% ORIGINAL
subplot(2,1,1)
% Pupil data
yyaxis right
plot(pupil2corr(:),'b-')
hold off
axis tight
set(gca,'YColor','b')
xlabel(['Index = ',num2str(bin_width),' sec'])
ylabel('Pupil area (pixels^2)')

% Pupil data
yyaxis left
plot(neural2corr(:),'r-')
hold off
axis tight
xlabel(['Index = ',num2str(bin_width),' sec'])
set(gca,'YColor','r')
ylabel('Firing rate (Hz)')
title({['R = ',num2str(R(lag == 0))];'Lag = 0'})


% SHIFTED
subplot(2,1,2)
% Pupil data
yyaxis right
hold on
plot(1:numel(time_concat),pupil2corr_concat,'b-')
%arrayfun(@(x) xline(x,'k-','LineWidth',2),find(time_concat == 0))
%arrayfun(@(x) xline(x,'k:'),find(time_concat == min(time_concat))-0.5)
hold off
axis tight
xlabel(['Index = ',num2str(bin_width),' sec'])
ylabel('Pupil area (pixels^2)')
set(gca,'YColor','b')

% Neural data
yyaxis left
hold on
plot((1:numel(time_concat))-lag_max,neural2corr_concat,'r-')
%arrayfun(@(x) xline(x,'k-','LineWidth',2),find(time_concat == 0))
%arrayfun(@(x) xline(x,'k:'),find(time_concat == min(time_concat))-0.5)
hold off
axis tight
xlabel(['Index = ',num2str(bin_width),' sec'])
ylabel('Firing rate (Hz)')
set(gca,'YColor','r')

title({['Max R = ',num2str(R_max)];['Max lag index = ',num2str(lag_max)]})
set(findobj(fig,'-property','FontSize'),'FontSize',15)






%% Plot second derivative?

dpdt  = gradient(pupil_data)./gradient(pupil_time'); % First derivative
ddpdt = gradient(dpdt)./gradient(pupil_time'); % Second derivative

dndt  = gradient(neural_data)./gradient(rate_time); % First derivative
ddndt = gradient(dndt)./gradient(rate_time); % Second derivative


fig = figure(6);
subplot(3,1,1)
plot(pupil_time, pupil_data)
axis tight
ylabel('Pupil area (pixels^2)')

subplot(3,1,2)
plot(pupil_time, dpdt)
axis tight
ylabel({'First derivative of';'pupil area (pixels^2) per second'})

subplot(3,1,3)
plot(pupil_time, ddpdt)
axis tight
ylabel({'Second derivative of';'pupil area (pixels^2) per second^2'})
xlabel('Time (s)')

set(findobj(fig,'-property','FontSize'),'FontSize',15)

fig = figure(7);
subplot(3,1,1)
plot(rate_time, neural_data)
axis tight
ylabel('Firing rate (Hz)')

subplot(3,1,2)
plot(rate_time, dndt)
axis tight
ylabel({'First derivative of';'firing rate (Hz) per second'})

subplot(3,1,3)
plot(rate_time, ddndt)
axis tight
ylabel({'Second derivative of';'firing rate (Hz) per second^2'})
xlabel('Time (s)')

set(findobj(fig,'-property','FontSize'),'FontSize',15)

fig = figure(8);
clf
subplot(3,1,1)
plot(pupil_time, ddpdt, 'Color','b')
%stem(laser_video, max(pupil_data).*ones(size(laser_video)),'k-','LineWidth',1,'MarkerEdgeColor','none','DisplayName','Laser video')
axis tight
ylabel({'Second derivative of';'pupil area (pixels^2) per second^2'})
xlabel('Time (s)')

subplot(3,1,2)
plot(rate_time, ddndt,'Color','r')
%stem(laser_neural_data, max(neural_data).*ones(size(laser_neural_data)),'-k','LineWidth',1,'MarkerEdgeColor','none','DisplayName','Laser neural')
axis tight
ylabel({'Second derivative of';'firing rate (Hz) per second^2'})
xlabel('Time (s)')

subplot(3,1,3)
hold on
yyaxis left
plot(pupil_time, ddpdt, '-','Color','b','DisplayName','Pupil 2nd d')
axis tight
ylabel({'Second derivative of';'pupil area (pixels^2) per second^2'})
xlabel('Time (s)')
stem(laser_video, max(ddpdt).*ones(size(laser_video)),'k-','LineWidth',1,'MarkerEdgeColor','none','DisplayName','Laser video')
yyaxis right
plot(rate_time, ddndt, '-','Color','r','DisplayName','neural 2nd d')
stem(laser_neural_data, max(ddndt).*ones(size(laser_neural_data)),'-k','LineWidth',1,'MarkerEdgeColor','g','DisplayName','Laser neural')
axis tight
ylabel({'Second derivative of';'firing rate (Hz) per second^2'})

hold off

%%%
fig = figure (9);
clf
hold on
plot(neural2corr)
hold on
plot(mean(neural2corr,2),'k-','LineWidth',2)
title({['Firing rate']})

fig = figure (10);
clf
plot(pupil2corr)
hold on
plot(mean(pupil2corr,2),'b-','LineWidth',2)
title({['Pupil']})

fig = figure (11);
clf
hold on
yyaxis right
plot(mean(pupil2corr,2),'b-','LineWidth',2)
yyaxis left
plot(mean(neural2corr,2),'r-','LineWidth',2)
%crosscorr(mean(pupil2corr,2),mean(neural2corr,2))
%clf
%crosscorr(mean(pupil2corr,2),mean(neural2corr,2))
[R,lag] = crosscorr(mean(pupil2corr,2),mean(neural2corr,2));


max(abs(R));


max(abs(R)) == R;


lag(max(abs(R)) == R);


std_pupil = std(pupil2corr,[],2);


mean_pupil = mean(pupil2corr,2);

ci_pupil = [mean_pupil - 1.96*(std_pupil./sqrt(numel(pupil2corr))), mean_pupil + 1.96*(std_pupil./sqrt(numel(pupil2corr)))];
ci_pupil_plot = [ci_pupil(:,1) - ci_pupil(:,2)];

fig = figure (12);
clf
shadedErrorBar(binned_time,mean_pupil, ci_pupil_plot,'lineProps','-b')

std_neural = std(neural2corr, [],2);
mean_neural = mean(neural2corr,2);

ci_neural = [mean_neural - 1.96*(std_neural./sqrt(numel(neural2corr))), mean_neural + 1.96*(std_neural./sqrt(numel(neural2corr)))];
ci_neural_plot = [ci_neural(:,1) - ci_neural(:,2)];

fig = figure (13);
clf
shadedErrorBar(binned_time,mean_neural, ci_neural_plot,'lineProps','-r')

fig = figure (14);
clf
subplot (2,1,2)
hold on
yyaxis left
axis tight
shadedErrorBar(binned_time,mean_pupil, ci_pupil_plot,'lineProps','-b')
ylabel({'pupil area (pixels^2)'})
xlabel('Time (s)')
yyaxis right
shadedErrorBar(binned_time,mean_neural, ci_neural_plot,'lineProps','-r')
ylabel({'Firing rate (Hz)'})
title(['Mouse ',num2str(mouse_id),', ',date_str,', ',file_num,])
legend('location','bestoutside')
%title({['PSTH pupil and firing rate']})
ignore_this = xline(0,'k-','LineWidth',2);
set(get(get(ignore_this, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off')

%fig = figure (15);
%clf
subplot (2,1,1)
[r,p] = corr(mean_pupil, mean_neural,'Type','Pearson');
hold on
plot(mean_neural,mean_pupil, 'o','DisplayName',['r=',num2str(r)])
legend('location','eastoutside')
ls = lsline; % Make least square line
ls.DisplayName = ['p=',num2str(p)];
if p<=0.05
    ls.Color = 'r'; % If significant, then lsline turns red
    ls.LineWidth = 2;
else
    ls.Color = 'k';
end
xlabel(['\bf{',upper(pupil_type(1)),lower(pupil_type(2:end)),'} \rm{pupil area (pixels^2)}'])
ylabel(['\bf{',upper(neural_type(1)),lower(neural_type(2:end)),'} \rm{spiking rate (Hz)}'])
title(['Mouse ',num2str(mouse_id),', ',date_str,', ',file_num])

median_pupil = median(pupil2corr,2);
median_neural = median(neural2corr,2);

fig = figure (15);
clf
subplot (2,1,2)
hold on
yyaxis left
axis tight
shadedErrorBar(binned_time,median_pupil, ci_pupil_plot,'lineProps','-b')
ylabel({'pupil area (pixels^2)'})
xlabel('Time (s)')
yyaxis right
shadedErrorBar(binned_time,median_neural, ci_neural_plot,'lineProps','-r')
ylabel({'Firing rate (Hz)'})
title(['MEDIAN Mouse ',num2str(mouse_id),', ',date_str,', ',file_num,])
legend('location','bestoutside')
%title({['PSTH pupil and firing rate']})
ignore_this = xline(0,'k-','LineWidth',2);
set(get(get(ignore_this, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off')

subplot (2,1,1)
[r,p] = corr(median_pupil, median_neural,'Type','Pearson');
plot(median_neural,median_pupil, 'o','DisplayName',['r=',num2str(r)])
legend('location','eastoutside')
ls = lsline; % Make least square line
ls.DisplayName = ['p=',num2str(p)];
if p<=0.05
    ls.Color = 'r'; % If significant, then lsline turns red
    ls.LineWidth = 2;
else
    ls.Color = 'k';
end
xlabel(['\bf{',upper(pupil_type(1)),lower(pupil_type(2:end)),'} \rm{pupil area (pixels^2)}'])
ylabel(['\bf{',upper(neural_type(1)),lower(neural_type(2:end)),'} \rm{spiking rate (Hz)}'])
title(['Mouse ',num2str(mouse_id),', ',date_str,', ',file_num])

fig = figure (16);

shadedErrorBar(binned_time,mean_pupil, ci_pupil_plot,'lineProps','-r')





%% Make plot for presentation

% Set figure
fig = figure(17);
clf
%fig.WindowStyle = 'Normal';
fig.Units = 'inches';
fig.Position = [0,0,13,6];
layout = tiledlayout(1,3,'TileSpacing','normal','Padding','compact');


%~~~~~~~~ Plot shaded error ~~~~~~~~
nexttile(layout,[1,2])

%%%% Pupil %%%%
% Plot
yyaxis left
h1 = shadedErrorBar(binned_time,mean_pupil, ci_pupil_plot,'lineProps','-b');

% Format
set(gca,'YColor','b')
set(gca,'TickDir','out')
delete(h1.edge)
x_lim = get(gca,'YLim');

% Label
ylabel({'Pupil area (pixels^2)'})



%%%% Pupil %%%%
% Plot
yyaxis right
hold on
h2 = shadedErrorBar(binned_time,mean_neural, ci_neural_plot,'lineProps','-r');
ln = plot(binned_time(end)+[-10,-5], 1.25*mean(mean_neural).*ones(1,2),'k-');
hold off

% Format
set(gca,'XLim',[binned_time(1), binned_time(end)])
set(gca,'YColor','r')
set(gca,'TickDir','out')
delete(h2.edge)
y_lim = get(gca,'YLim');
set(gca,{'XTick','XTickLabels','XColor'},{[],'','none'})

% Label
ylabel('Firing rate (Hz)','Rotation',270,'VerticalAlignment','baseline')
%xlabel('Time (s)')
text(binned_time(end)+-10, 1.2*mean(mean_neural),'5 sec')
legend('Average pupil area','Average firing rate','Location','northeast')





%~~~~~~~~ Correlation ~~~~~~~~
nexttile(layout,[1,1])

% Plot
[r,p] = corr(mean_pupil,mean_neural,'Type','Pearson');
scatter(mean_pupil,mean_neural, 'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.5)

% Format
set(gca,'TickDir','out')
axis square
box on
set(gca,'XLim',x_lim)
set(gca,'YLim',y_lim)

% Label
ylabel('Firing rate (Hz)')
xlabel({'Pupil area (pixels^2)'})
corr_str = {['$r=',num2str(round(r,2,'significant')),'$'];...
    ['$p=',replace(num2str(round(p,3,'significant')),'e',' \times 10^{'),'}$']};
text(gca,x_lim(1)*1.1,y_lim(end)*0.95,corr_str,'interpret','latex','VerticalAlignment','top')

% Add line 
hold on
ls = lsline;
ls.Color = [0.7,0.7,0.7];
hold off



%~~~~~~~~ Uniform formatting ~~~~~~~~
set(findall(fig,'-property','FontSize'),'FontSize',14)
set(findall(fig,'-property','FontName'),'FontName','Arial')
set(findobj(fig,'-property','LineWidth'),'LineWidth',1)
ls.LineWidth = 2;
ln.LineWidth = 5;



% hold on
% yyaxis left
% axis tight
% shadedErrorBar(binned_time,mean_pupil, ci_pupil_plot,'lineProps','-b')
% ylabel({'pupil area (pixels^2)'})
% xlabel('Time (s)')
% yyaxis right
% shadedErrorBar(binned_time,mean_neural, ci_neural_plot,'lineProps','-r')
% ylabel({'Firing rate (Hz)'})
% title(['Mouse ',num2str(mouse_id),', ',date_str,', ',file_num,])
% legend('location','bestoutside')
% %title({['PSTH pupil and firing rate']})
% ignore_this = xline(0,'k-','LineWidth',2);
% set(get(get(ignore_this, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off')


% [r,p] = corr(mean_pupil, mean_neural,'Type','Pearson');
% hold on
% plot(mean_neural,mean_pupil, 'o','DisplayName',['r=',num2str(r)])
% legend('location','eastoutside')
% ls = lsline; % Make least square line
% ls.DisplayName = ['p=',num2str(p)];
% if p<=0.05
%     ls.Color = 'r'; % If significant, then lsline turns red
%     ls.LineWidth = 2;
% else
%     ls.Color = 'k';
% end
% xlabel(['\bf{',upper(pupil_type(1)),lower(pupil_type(2:end)),'} \rm{pupil area (pixels^2)}'])
% ylabel(['\bf{',upper(neural_type(1)),lower(neural_type(2:end)),'} \rm{spiking rate (Hz)}'])
% title(['Mouse ',num2str(mouse_id),', ',date_str,', ',file_num])
%% Make crosscorrelation for each stimulation
% 
% % Set figure
% fig = figure(18);
% clf


