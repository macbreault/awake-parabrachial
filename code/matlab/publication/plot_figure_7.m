% plot_figure_7.m
%
% Script to makes figure 7 for publication, about photometry
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Modified: 10-17-2022

clear
clc
set(0,'DefaultFigureWindowStyle','normal') %docks the figures onto the Matlab window

save_it = 1; % Whether to save plots (1) or not (0)
overlay_it = 0; % Make overlaed version (1) or not (0)



%% Set data


params=struct('a',struct(),'b',struct(),'c',struct(),'d',struct());

% ~~~~~~~~~~~ Set mouse information ~~~~~~~~~~~

params.a = struct('mouse_id',  'F3',...
    'date_str',  '08242022',...
    'cond_type', 'pre', ...
    'cond_day',  0, ...
    'stim_type', 'heat', ...
    'type',      'zscore', ...
    'pre_time', 5,...
    'post_time', 30);

params.b = struct('mouse_id',  'M3',...
    'date_str',  '08182022',...
    'cond_type', 'pre', ...
    'cond_day',  0, ...
    'stim_type', 'heat', ...
    'type',      'zscore', ...
    'pre_time', 5,...
    'post_time', 30);

params.c = struct('mouse_id',  'F3',...
    'date_str',  '09012022',...
    'cond_type', 'cfa', ...
    'cond_day',  4, ...
    'stim_type', 'heat', ...
    'type',      'zscore', ...
    'pre_time', 5,...
    'post_time', 30);

params.d = struct('mouse_id',  'M3',...
    'date_str',  '09012022',...
    'cond_type', 'cfa', ...
    'cond_day',  4, ...
    'stim_type', 'heat', ...
    'type',      'zscore', ...
    'pre_time', 5,...
    'post_time', 30);



params.ef = struct('trace_type','Z Score Heat Duration');



min_response_time_on  = 0.5; % seconds
min_response_time_off = 0.5; % seconds





% ~~~~~~~~~~~~~~~~~~~~~~ DO NOT CHANGE BELOW ~~~~~~~~~~~~~~~~~~~~~~
% Figure configurations
fig_params = struct('Position',[0 0 18 12], 'Units','centimeters', 'WindowStyle','normal', 'FontSize', 10);
% 18 cm = max




%%%%%% Set paths %%%%%%
paths = struct('dropbox',[],'neural',[]);

if contains(pwd,'/mac/')
    paths.dropbox = '/Users/mac/Dropbox (MIT)/Jesse';
else
    paths.dropbox = '/Users/Jessesmith 1/Dropbox/Jesse';
end

red_driver_path = '/Volumes/T7';
seagate_driver_path = '/Volumes/Seagate Backup Plus Drive';

% Find path to neural data by first checking SEAGATE, then checking RED, then assigning DROPBOX
paths.neural = fullfile(seagate_driver_path, 'Photometry');
if ~exist(paths.neural,'dir') % Check seagate
    paths.neural = replace(paths.neural,seagate_driver_path,red_driver_path);
    if ~exist(paths.neural,'dir')% Check red drive
        paths.neural = fullfile(paths.dropbox,'data','photometry'); % Else, assign dropbox
    end
end


% Check that both paths exist
if ~exist(paths.dropbox,'dir'); error('Error in %s: Dropbox path not found', mfilename); end
if ~exist(paths.neural,'dir'); error('Error in %s: Neural hard drive not found', mfilename); end

addpath(genpath(fullfile(paths.dropbox,'code'))); % Add libraries to directory

plot_path  = fullfile(paths.dropbox,'plots',findPlotFolder,'1-draft');
if ~exist(plot_path,'dir') && save_it; mkdir(plot_path); end



%% Initialize figure

fig = figure(7);
clf(fig)
fig.WindowStyle = fig_params.WindowStyle;
fig.Units = fig_params.Units;
fig.Position = fig_params.Position;

layout_main = tiledlayout(4,1,'TileSpacing','compact','Padding','none');


% Prepare figure A+B
tile_ab = nexttile(layout_main); axis off;
layout_ab = tiledlayout(layout_main, 1, 2,'TileSpacing','tight','Padding','none');
layout_ab.Layout.Tile = 1;

% Prepare figure C+D
tile_cd = nexttile(layout_main); axis off;
layout_cd = tiledlayout(layout_main, 1, 2,'TileSpacing','tight','Padding','none');
layout_cd.Layout.Tile = 2;

% Prepare figure E+F
if overlay_it
    tile_ef = nexttile(layout_main); axis off;
    layout_ef = tiledlayout(layout_main, 1, 2,'TileSpacing','tight','Padding','none');
    layout_ef.Layout.Tile = 3;
else

    % Prepare figure E+F
    tile_ef = nexttile(layout_main,[2,1]); axis off;
    layout_ef = tiledlayout(layout_main, 1, 2,'TileSpacing','tight','Padding','none');
    layout_ef.Layout.Tile = 3;
    layout_ef.Layout.TileSpan = [2,1];

    % Prepare figure E
    tile_e = nexttile(layout_ef); axis off;
    layout_e = tiledlayout(layout_ef, 5, 1,'TileSpacing','tight','Padding','none','TileIndexing','columnmajor');
    layout_e.Layout.Tile = 1;
    %layout_e.Layout.TileSpan = [1,2];

    % Prepare figure C_2
    tile_f = nexttile(layout_ef); axis off;
    layout_f = tiledlayout(layout_ef, 5, 1,'TileSpacing','tight','Padding','none','TileIndexing','columnmajor');
    layout_f.Layout.Tile = 2;
    %layout_f.Layout.TileSpan = [1,2];
end



%% A

plot_recording(params.a, layout_ab, paths, 'a')



%% B

plot_recording(params.b, layout_ab, paths, 'b')



%% C

plot_recording(params.c, layout_cd, paths, 'c')



%% D

plot_recording(params.d, layout_cd, paths, 'd')



%% E

if overlay_it
    plot_traces_overlay(params.ef, {'F5','F4','F3','F2','F1'}, layout_ef, paths, 'e')
else
    plot_traces(params.ef, {'F5','F4','F3','F2','F1'}, layout_e, paths, 'e')
end



%% F

if overlay_it
    plot_traces_overlay(params.ef, {'M5','M4','M3','M2','M1'}, layout_ef, paths, 'f')
else
    plot_traces(params.ef, {'M5','M4','M3','M2','M1'}, layout_f, paths, 'f')
end



%% Save

sameLimits(findobj(fig,'Type','Axes','UserData','same_c'),'x')
%sameLimits(findobj(fig,'Type','Axes','UserData','same_c'),'y')
%sameLimits(findobj(fig,'Type','Axes','UserData','same_ab'))
uniformFormat_pub(fig)

if save_it
    exportgraphics(fig,fullfile(plot_path,['figure_',num2str(fig.Number),'.pdf']),'Resolution',300,'ContentType','vector')
end
















%% ~~~~~~~~~~~~~~~~~~~~~~~~~~ Functions for figures ~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Plot recording
function [] = plot_recording(params, layout, paths, tag)

min_response_time_on  = 0.5; % seconds
min_response_time_off = 0.5; % seconds
heat_duration = 5; % seconds

%~~~~~~~~~~~~ Load data ~~~~~~~~~~~~
condition_file = fullfile(paths.dropbox,'tables','Photometry Recordings.xlsx');
condition_table = readtable(condition_file, 'Sheet', params.mouse_id,'PreserveVariableNames',true); %'VariableNamingRule','preserve')

% Filter out the conditions that do not match
datetime_str = datetime(params.date_str,'InputFormat','MMddyyyy');
date_inds = condition_table.Date == datetime_str;

cond_inds = strcmpi(condition_table.Condition, params.cond_type);

if ~strcmpi(params.cond_type,'pre')
    day_inds = condition_table.('Days Post CFA') == params.cond_day;
else
    day_inds = isnan(condition_table.('Days Post CFA')) | (condition_table.('Days Post CFA') == 0);
end

stim_ind = contains(condition_table.('Stimulation type'), params.stim_type,'IgnoreCase',1);

condition_table(~(date_inds & cond_inds & day_inds & stim_ind),:) = [];

if min(size(condition_table)) ~= 1
    error('Found 0 or more than 1 condition in table')
end

% Load data
mouse_id = params.mouse_id;
type = params.type;
date_str = datestr(condition_table.Date,'mmddyyyy');
cond_type = condition_table.Condition{1};
if strcmpi(cond_type,'pre')
    cond_day = '0';
else
    cond_day = num2str(condition_table.('Days Post CFA'));
end
stim_type = strrep(lower(condition_table.('Stimulation type'){1}),' ','-');
stim_type(1) = upper(stim_type(1));

% Pick file name based on type of data to process
switch type
    case 'zscore'
        file_name = '_ZScorePETH';
    case 'delta'
        file_name = '_DeltaFPETH';
end

% Load .fig of Z-Score data
folder_name = strjoin({['#',mouse_id],cond_type,cond_day,stim_type},'-');
data_path = fullfile(paths.neural, date_str, folder_name, 'Figures');
listings = dir(data_path);
file_ind = reshape(find(contains({listings.name}, file_name) & contains({listings.name},'.fig')),1,[]);
if ~exist(data_path,'dir') || ~any(file_ind)
    warning('Could not find Figure folder or .fig file. Skip!')
end
if numel(file_ind) == 0
    error('Listing not found')
elseif numel(file_ind) > 1
    file_ind = file_ind(contains({listings(file_ind).name}, num2str(params.post_time)) & ...
                        contains({listings(file_ind).name}, num2str(params.pre_time)));
    if numel(file_ind) ~= 1
        error('Too many listing files')
    end
end
listing = listings(file_ind);

% Load figure as invisible
data_fig = openfig(fullfile(listing.folder, listing.name), 'invisible');

% Photometry data (as image object from data_fig)
img_obj = findobj(data_fig,'Type','Image');

if numel(img_obj) ~= 1
    figure(data_fig) % Make figure visible because of error
    error('Error in %s: Found more than one image plotted in figure that contains data. Check!',mfilename)
end

data = img_obj.CData; % trials x time bins

% Extract y-label from axes that has average plot on it
lines = findobj(data_fig,'Type','Line');
y_label = lines(1).Parent.YLabel.String;

% Trial information
total_trials = size(data,1);
good_trials = true(total_trials,1);
remove_triggers = condition_table.('Triggers removed?'); % Get triggers to ignore
if ~isnan(remove_triggers)
    good_trials(remove_triggers(total_trials >= remove_triggers)) = 0;
end

% Time information
pre_time = img_obj.Parent.XLim(1); % -5
post_time = img_obj.Parent.XLim(2); % 10
total_bins = size(data,2);
time_bins  = linspace(pre_time, post_time, total_bins);
%time_bins = img_obj.XData;

% Close figure
close(data_fig)



%~~~~~~~~~~~~ Prepare data ~~~~~~~~~~~~
plot_time = time_bins;
plot_data = mean(data(good_trials,:),1,'omitnan');
plot_CI =  calculateCR(data(good_trials,:)');

% Calculate duration of significant response(s)

duration = NaN;

ind_before = find(time_bins < 0);
ind_after  = find(time_bins > 0);

data_before = plot_data(ind_before);
data_after  = plot_data(ind_after);

% 1. Use upper bound of CI from baseline as threshold
thresh = max(calculateCI(data_before));

% 2. Refine responses, so consecutive bins after threshold above/below threshold to be considered on/off
min_response_bin_on = round(min_response_time_on / unique(round(diff(time_bins),5))); % Number of bins corresponding to min_response_time
min_response_bin_off = round(min_response_time_off / unique(round(diff(time_bins),5))); % Number of bins corresponding to min_response_time

% Onset
onset = thresh <= data_after;
labeledVector = bwlabel(onset);
measurements = regionprops(labeledVector, data_after, 'Area', 'PixelValues','PixelIdxList');
measurements([measurements.Area] < min_response_bin_on) = []; % Remove indices that are too small
onset_inds = arrayfun(@(ind) min(measurements(ind).PixelIdxList), 1:numel(measurements));


% Offset
offset = thresh > data_after;
labeledVector = bwlabel(offset);
measurements = regionprops(labeledVector, data_after, 'Area', 'PixelValues','PixelIdxList');
measurements([measurements.Area] < min_response_bin_off) = []; % Remove indices that are too small
offset_inds = arrayfun(@(ind) min(measurements(ind).PixelIdxList), 1:numel(measurements));


% Get duration for every onset...
duration = zeros(size(onset_inds));
plot_duration = cell(size(onset_inds));
for j = 1:numel(onset_inds)

    onset_time = plot_time(ind_after(onset_inds(j)));

    % Find offset index that matches onset index
    offset_ind = min(offset_inds(onset_inds(j) < offset_inds));
    if isempty(offset_ind)
        offset_ind = numel(ind_after);
    end
    offset_time = plot_time(ind_after(offset_ind));

    duration(j) = offset_time - onset_time;
    plot_duration{j} = ind_after(onset_inds(j)):ind_after(offset_ind);

end

% Check that there is no overlap between durations
overlap = NaN(numel(duration));
for j1 = 1:numel(duration)
    for j2 = 1:numel(duration)
        overlap(j1,j2) = all(ismember(plot_duration{j1},plot_duration{j2}));
    end
end
actual_overlap = triu(overlap,1) | tril(overlap,-1);
[row,col] = ind2sub(size(actual_overlap), find(actual_overlap));
for j = 1:numel(row)
    % Keep the index with the largest duration
    if duration(row(j)) > duration(col(j))
        duration(col(j)) = NaN;
        plot_duration{col(j)} = {};
    elseif duration(row(j)) < duration(col(j))
        duration(row(j)) = NaN;
        plot_duration{row(j)} = {};
    end
end

% Remove empty
duration(isnan(duration)) = [];
plot_duration(cellfun(@isempty,plot_duration)) = [];

% Calculate AREA under significant response(s)
AUC = NaN(size(duration));
plot_AUC = cell(size(duration));

for j = 1:numel(AUC)

    % Calculate AUC between these curves:
    x  = plot_time(plot_duration{j});
    y1 = plot_data(plot_duration{j}); % Response
    y2 = repmat(thresh ,size(x));     % Threshold

    % 1. Find all points where these curves intersect (not including ends)
    [x0,y0] = intersections(x,y1,x,y2);

    % 2. Add end points
    x_intersect = unique([x(1); x0; x(end)]);

    % 3. Find y1 value for each intersect
    x_with_intersects =  [x, x0'];
    y1_with_intersects = [y1,y0'];
    [x_with_intersects_sorted, order] = sort(x_with_intersects);
    y1_with_intersects_sorted = y1_with_intersects(order);

    % 4. Integrate between each point of intersection
    AUC_temp = NaN(numel(x_intersect)-1,1);
    for jj = 1:(numel(x_intersect)-1)
        mask = (x_intersect(jj) <= x_with_intersects_sorted) & (x_with_intersects_sorted <= x_intersect(jj+1));
        AUC_temp(jj) = trapz(x_with_intersects_sorted(mask), abs(y1_with_intersects_sorted(mask)-thresh));
    end

    % 5. Sum all areas together
    AUC(j) = sum(AUC_temp); % trapz(x,abs(y1-thresh)) % Can't use this one because of areas that go below threshold

    % Make mask of AUC for plotting
    plot_AUC{j}(1,:) = [x, fliplr(x)]; % x-value
    plot_AUC{j}(2,:) = [y1, fliplr(y2)]; % y-value

end


%~~~~~~~~~~~~ Plot ~~~~~~~~~~~~
nexttile(layout)
H = shadedErrorBar(plot_time, plot_data, plot_CI);
delete(H.edge)
hold on
arrayfun(@(j) fill(plot_AUC{j}(1,:), plot_AUC{j}(2,:), 'g', 'EdgeColor','none','FaceAlpha',0.25,'DisplayName', ['AUC ',num2str(j),' = ',num2str(AUC(j))]), 1:numel(AUC))
arrayfun(@(j) plot(plot_time(plot_duration{j}), plot_data(plot_duration{j}),'-','Color','g','DisplayName', ['Duration ',num2str(j),' = ',num2str(duration(j)),' s']), 1:numel(duration))
xline(0,'Tag','trigger')
yline(thresh,'Tag','thresh')

% Add rectangle for heat
heat_pos = [0, min(ylim), heat_duration, sum(abs(ylim))]; %[x y w h]
rectangle('Position', heat_pos)

hold off

%~~~~~~~~~~~~ Format ~~~~~~~~~~~~
axis tight

%~~~~~~~~~~~~ Label ~~~~~~~~~~~~
set(gca,'Tag',tag)
xlabel('Time (s)')
if strcmp(params.cond_type,'pre')
    title({'Pre CFA',params.mouse_id})
elseif strcmp(params.cond_type,'cfa')
    title({[num2str(params.cond_day), ' days Post CFA'],params.mouse_id})
end
ylabel('Z Score Heat')
set(gca,'UserData','same_ab')

end % End plot recording





%% Sets tags for boxplot for formatting later
function H = setNotBoxPlotTags(H)

    for field = reshape(fieldnames(H),1,[])
        if isgraphics(H.(field{:}))
            H.(field{:}).Tag = field{:};
        end
    end

end





%% Plot traces
function [] = plot_traces(params, mice, layout, paths, tag)

%~~~~~~~~~~~~ Prepare data ~~~~~~~~~~~~
behavior_file = fullfile(paths.dropbox,'tables','Photometry Behavior.xlsx');
trace_file = fullfile(paths.dropbox,'tables','Photometry Traces.xlsx');

behavior_all = readtable(behavior_file,'ReadRowNames',1,'PreserveVariableNames',1);
trace = readtable(trace_file,'Sheet',params.trace_type,'ReadRowNames',1,'PreserveVariableNames',1);

mouse_ids = flip(sort(behavior_all.Properties.VariableNames));
mouse_ids = cellfun(@(mouse) mouse_ids(ismember(mouse_ids,mouse)), mice);

for mouse = 1:numel(mouse_ids)

    mouse_id = mouse_ids{mouse};

    % Behavior
    mouse_ind = strcmp(behavior_all.Properties.VariableNames,mouse_id);
    plot_time = cellfun(@str2num, behavior_all(:,mouse_ind).Properties.RowNames);
    plot_data = behavior_all{:,mouse_ind};
    pain_thresh = min(plot_data(plot_time <= 0)); %%%%%%%%%%%% CHANGE!!!
    behave_table = [plot_time, plot_data];
    pain_ind = plot_data == pain_thresh;
    behave_table(pain_ind,:) = []; % Remove days without pain

    % Traces
    mouse_ind = strcmp(trace.Properties.VariableNames,mouse_id);
    plot_time = cellfun(@str2num, trace(:,mouse_ind).Properties.RowNames);
    plot_data = trace{:,mouse_id};
    plot_CI = calculateCI(plot_data(0 > plot_time)'); % Calculate CI using [-1,-2,-3]

    trace_table = [plot_time, plot_data];

    % Remove NaN
    bad_ind = any(isnan(behave_table),2);
    behave_table(bad_ind,:) = [];

    bad_ind = any(isnan(trace_table),2);
    trace_table(bad_ind,:) = [];

    %~~~~~~~~~~~~ Plot ~~~~~~~~~~~~
    nexttile(layout)
    hold on
    plot(trace_table(:,1), trace_table(:,2), '-','Color','k','Tag','no_legend') % Plot line traces
    plot(trace_table(:,1), trace_table(:,2), 'o', 'MarkerFaceColor','k','MarkerEdgeColor','none','DisplayName','Normal pain','MarkerSize',5) % Plot pain days
    
    % Overlay days of pain
    trace_2_behave = ismember(trace_table(:,1), behave_table(:,1));
    plot(trace_table(trace_2_behave,1), trace_table(trace_2_behave,2), 'o', 'MarkerFaceColor','y','MarkerEdgeColor','none','DisplayName','Hyperalgesia','MarkerSize',5) % Plot pain days
    xline(0,'Tag','trigger')
    yline(plot_CI,'Tag','thresh')
    hold off


    %~~~~~~~~~~~~ Format ~~~~~~~~~~~~
    axis tight
    box off

    %~~~~~~~~~~~~ Label ~~~~~~~~~~~~
    if mouse == 1
        if      contains(mouse_id, 'F','IgnoreCase',1)
            title({'Females',mouse_id})
        elseif  contains(mouse_id, 'M','IgnoreCase',1)
            title({'Males',mouse_id})
            legend('Location','northeast')
        end
        set(gca,'Tag',tag)
    else
        title(mouse_id)
    end

    if (mouse == ceil(numel(mice)/2)) && ~isempty(tag)
        if contains(params.trace_type,'duration','IgnoreCase',1)
            ylabel([params.trace_type, ' (s)'])
        else
            ylabel(params.trace_type)
        end
    end

    if mouse ~= numel(mouse_ids)
        set(get(gca,'XAxis'),'Color','none')
        xticks([])
    else
        xlabel('Days Post CFA')
    end

    set(gca,'UserData','same_c')
 
end

end



%% Plot traces and overlays
function [] = plot_traces_overlay(params, mice, layout, paths, tag)

%~~~~~~~~~~~~ Prepare data ~~~~~~~~~~~~
behavior_file = fullfile(paths.dropbox,'tables','Photometry Behavior.xlsx');
trace_file = fullfile(paths.dropbox,'tables','Photometry Traces.xlsx');

behavior_all = readtable(behavior_file,'ReadRowNames',1,'PreserveVariableNames',1);
trace = readtable(trace_file,'Sheet',params.trace_type,'ReadRowNames',1,'PreserveVariableNames',1);

mouse_ids = flip(sort(behavior_all.Properties.VariableNames));
mouse_ids = cellfun(@(mouse) mouse_ids(ismember(mouse_ids,mouse)), mice);

nexttile(layout)

for mouse = 1:numel(mouse_ids)

    mouse_id = mouse_ids{mouse};

    % Behavior
    mouse_ind = strcmp(behavior_all.Properties.VariableNames,mouse_id);
    plot_time = cellfun(@str2num, behavior_all(:,mouse_ind).Properties.RowNames);
    plot_data = behavior_all{:,mouse_ind};
    pain_thresh = min(plot_data(plot_time <= 0)); %%%%%%%%%%%% CHANGE!!!
    behave_table = [plot_time, plot_data];
    pain_ind = plot_data == pain_thresh;
    behave_table(pain_ind,:) = []; % Remove days without pain

    % Traces
    mouse_ind = strcmp(trace.Properties.VariableNames,mouse_id);
    plot_time = cellfun(@str2num, trace(:,mouse_ind).Properties.RowNames);
    plot_data = trace{:,mouse_id};
    plot_CI = calculateCI(plot_data(0 > plot_time)'); % Calculate CI using [-1,-2,-3]

    trace_table = [plot_time, plot_data];

    % Remove NaN
    bad_ind = any(isnan(behave_table),2);
    behave_table(bad_ind,:) = [];

    bad_ind = any(isnan(trace_table),2);
    trace_table(bad_ind,:) = [];

    %~~~~~~~~~~~~ Plot ~~~~~~~~~~~~
    hold on
    plot(trace_table(:,1), trace_table(:,2), '-','Color','k','Tag','no_legend') % Plot line traces
    plot(trace_table(:,1), trace_table(:,2), 'o', 'MarkerFaceColor','k','MarkerEdgeColor','none','DisplayName','Normal pain','MarkerSize',5) % Plot pain days
    
    % Overlay days of pain
    trace_2_behave = ismember(trace_table(:,1), behave_table(:,1));
    plot(trace_table(trace_2_behave,1), trace_table(trace_2_behave,2), 'o', 'MarkerFaceColor','y','MarkerEdgeColor','none','DisplayName','Hyperalgesia','MarkerSize',5) % Plot pain days
    xline(0,'Tag','trigger')
    yline(max(plot_CI),'Tag','thresh')
    hold off


    %~~~~~~~~~~~~ Format ~~~~~~~~~~~~
    axis tight
    box off

    %~~~~~~~~~~~~ Label ~~~~~~~~~~~~
%     if mouse == 1
        if      contains(mouse_id, 'F','IgnoreCase',1)
            title('Females')
        elseif  contains(mouse_id, 'M','IgnoreCase',1)
            title('Males')
%             legend('Location','northeast')
        end
        set(gca,'Tag',tag)
%     else
%         title(mouse_id)
%     end

%     if (mouse == ceil(numel(mice)/2)) && ~isempty(tag)
        if contains(params.trace_type,'duration','IgnoreCase',1)
            ylabel([params.trace_type, ' (s)'])
        else
            ylabel(params.trace_type)
        end
%     end

%     if mouse ~= numel(mouse_ids)
%         set(get(gca,'XAxis'),'Color','none')
%         xticks([])
%     else
        xlabel('Days Post CFA')
%     end

    set(gca,'UserData','same_c')
 
end

end