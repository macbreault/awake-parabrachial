% plot_recordings_v1.m
%
% Script to loads and processes photometry data for Jesse
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 09-12-2022 - v0
% Modified: 10-05-2022 - v1 - Now, instead of relying on Excel for time before and after and time bin, it now relys on .fig that preprocessing saves to import data

clear
clc
close all
set(0,'DefaultFigureWindowStyle','docked') %docks the figures onto the Matlab window

warning('TODO: 1. Combine overlapping significant response durations, 2. Add duration and AUC to table')



%% Set data

% ~~~~~~~~~~~ Set mouse information ~~~~~~~~~~~
mouse_id = 'F4';                     % Mouse ID
date_str = '02062023';               % Date of recording
cond_type = 'cfa';                   % Pain type ('pre','cfa','cci')
cond_day  = 57;                      % Day since pain type (0 for 'pre')
stim_type = 'heat';                  % Type of stimulation applied ('heat' or 'tail')



% Other variables to change
run_all = 0;  % Boolean as to whether to run all conditions (1) or not (0)
save_it = 1; % Whether to save plots (1) or not (0)


% Data to run
types = {'zscore','delta'};


min_response_time_on  = 0.5; % seconds
min_response_time_off = 0.5; % seconds





% ~~~~~~~~~~~~~~~~~~~~~~ DO NOT CHANGE BELOW ~~~~~~~~~~~~~~~~~~~~~~
% Figure configurations
fig_params = struct('Position',[0 0 8 5], 'Units','Inches', 'WindowStyle','normal', 'FontSize', 10);




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




%% Loop through all conditions

condition_file = fullfile(paths.dropbox,'tables','Photometry Recordings.xlsx');
mouse_ids = sheetnames(condition_file);

if ~run_all
    mouse_ids = mouse_ids(ismember(mouse_ids,mouse_id));
end

% Loop through all mice...
for mouse = 1:numel(mouse_ids)

    mouse_id = mouse_ids{mouse};

    disp(['Running mouse ',mouse_id,'...'])

    %%%%%% Import their table %%%%%%

    opts = detectImportOptions(condition_file, 'Sheet', mouse_id,'PreserveVariableNames',true);
    opts.VariableTypes(strcmp(opts.VariableNames,'Triggers removed?')) = {'char'};
    condition_table = readtable(condition_file, opts);
    %condition_table = readtable(condition_file, 'Sheet', mouse_id,'PreserveVariableNames',true); %'VariableNamingRule','preserve')

    % Filter out the conditions that do not match
    if ~run_all

        % ~~~~~~~ Check which conditions should be run (if not empty) ~~~~~~~
        date_inds = true(size(condition_table,1),1);
        cond_inds = true(size(condition_table,1),1);
        day_inds = true(size(condition_table,1),1);
        stim_ind = true(size(condition_table,1),1);

        % Date
        if ~isempty(date_str)
            datetime_str = datetime(date_str,'InputFormat','MMddyyyy');
            date_inds = condition_table.Date == datetime_str;
        end

        % Condition
        if ~isempty(cond_inds)
            cond_inds = strcmpi(condition_table.Condition,cond_type);
        end

        % Pain day
        if ~isempty(day_inds)
            if ~strcmpi(cond_type,'pre')
                day_inds = condition_table.('Days Post CFA') == cond_day;
            else
                day_inds = isnan(condition_table.('Days Post CFA')) | (condition_table.('Days Post CFA') == 0);
            end
        end

        % Stimulation type
        if ~isempty(stim_ind)
            stim_ind = contains(condition_table.('Stimulation type'), stim_type,'IgnoreCase',1);
        end


        % Remove rows from condition_table that do not match conditions
        exper_inds = date_inds & cond_inds & day_inds & stim_ind;
        condition_table(~exper_inds,:) = [];

    end % Only keep conditions requested



    %% Loop through all conditions

    for i = 1:size(condition_table,1)

        % Extract conditions (if running all files in condition table) to make it match folders/files
        try
            date_str = datestr(condition_table(i,:).Date,'mmddyyyy');
            cond_type = condition_table(i,:).Condition{1};

            if strcmpi(cond_type,'pre')
                cond_day = '0';
            else
                cond_day = num2str(condition_table(i,:).('Days Post CFA'));
            end

            stim_type = strrep(lower(condition_table(i,:).('Stimulation type'){1}),' ','-');
            stim_type(1) = upper(stim_type(1));

            if strcmp(stim_type,'pinch')
                warning('RENAME!')
            end
        catch
            warning('Invalid date in excel table. Skip!')
            continue
        end

        disp(' ')
        disp(['             Running ',strjoin({date_str,cond_type,cond_day,stim_type},', '),'...'])



        %% Plot data for Z-Score and Delta
        for type = reshape(types,1,[])

            % Pick file name based on type of data to process
            switch type{:}
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
                continue
            end



            %% Extract data for every .fig that meets the condition
            for k = file_ind

                % Load figure as invisible
                data_fig = openfig(fullfile(listings(k).folder, listings(k).name), 'invisible');

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
                remove_triggers = condition_table(i,:).('Triggers removed?'); % Get triggers to ignore

                if iscell(remove_triggers) && (numel(remove_triggers) == 1)% Is a cell of size 1
                    if (numel(remove_triggers) == 1)
                        remove_triggers = cell2mat(remove_triggers);
                    else
                        error('Error in %s: Triggers to removed are not formatted for matlab. MSB must add new condition to this line', mfilename)
                    end
                end
                
                if isa(remove_triggers,'char')
                    remove_triggers = cellfun(@str2double, split(remove_triggers,','))';
                end
               
                if isa(remove_triggers,'double') % If a double
                    good_trials(remove_triggers(total_trials >= remove_triggers)) = 0;
                else
                    error('Error in %s: Triggers to removed are not formatted for matlab. MSB must add new condition to this line', mfilename)
                end

                % Time information
                pre_time = img_obj.Parent.XLim(1); % -5
                post_time = img_obj.Parent.XLim(2); % 10
                total_bins = size(data,2);
                time_bins  = linspace(pre_time, post_time, total_bins);
                %time_bins = img_obj.XData;


                % Close figure
                close(data_fig)



                %%%%%% Prepare for plotting %%%%%%

                % Create experiment title for all plots
                fig_title = ['Mouse ',strjoin({mouse_id,date_str,cond_type,cond_day,stim_type},', ')];

                extra = erase(cell2mat(regexp(listings(k).name,[file_name,'.*fig'],'match')),{file_name,'.fig'});
                params_folder = strjoin({...
                    ['time_(',num2str(pre_time),',',num2str(post_time),')',extra],...
                    },'|');

                plot_path  = fullfile(paths.dropbox,'plots',findPlotFolder,'recording',date_str, folder_name, params_folder);
                if ~exist(plot_path,'dir') && save_it; mkdir(plot_path); end



                %% Plot trials - SPLIT

                % Prepare figure
                fig = figure(sub2ind([numel(types),2],1,find(strcmp(type,types))));
                clf(fig)
                fig.WindowStyle = fig_params.WindowStyle;
                fig.Units = fig_params.Units;
                fig.Position = fig_params.Position;

                num = numSubplots(sum(good_trials));
                layout = tiledlayout(num(1),num(2),'TileSpacing','compact','Padding','none');
                title(layout, fig_title)

                % Loop through all trials
                for t = 1:total_trials

                    if any(~isnan(remove_triggers) & (t == remove_triggers))
                        continue
                    end

                    nexttile(layout)

                    % Plot
                    hold on
                    plot(time_bins, data(t,:))
                    xline(0,'Tag','trigger')
                    hold off

                    % Format
                    axis tight

                    % Label
                    title(['Trial ',num2str(t)])

                end

                xlabel(layout,'Time since stimulation [s]')
                ylabel(layout,y_label)

                % ~~~~~~~~~~~ Uniform formatting ~~~~~~~~~~~
                %sameLimits(fig,'y')
                uniformFormat(fig)

                % Save
                if save_it

                    % Save figure
                    exportgraphics(fig,fullfile(plot_path,[type{:},'_1_trials.pdf']),'Resolution',300,'ContentType','vector')
                    savefig(fig,fullfile(plot_path,[type{:},'_1_trials.fig']))

                end



                %% Plot trials - AVERAGED

                % Get data for plotting
                plot_time = time_bins;
                plot_data = mean(data(good_trials,:),1,'omitnan');
                plot_CI =  calculateCR(data(good_trials,:)');

                %         test_data = (1:3)'.*repmat(1:10,3,1);
                %         plot_time = 1:10;
                %         plot_data = mean(test_data,1,'omitnan');
                %         plot_CI =  calculateCR(test_data');


                % ~~~ Calculate duration of significant response(s) ~~~

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

                %{
        % Check to make sure durations do not overlap. If so, then combine
        if numel(plot_duration) > 1
            % Do something
            old_plot_duration = plot_duration;
            plot_duration = {};
            plot_duration{1} = old_plot_duration{1};
            new_plot_duration = {};
            
            for j = 2:numel(old_plot_duration)
                if any(intersect(plot_duration{j-1}, old_plot_duration{j}))
                    plot_duration{j-1} = unique([plot_duration{j-1}, old_plot_duration{j}]);
                end
            end
        end
                %}



                % ~~~ Calculate AREA under significant response(s) ~~~
                AUC = NaN(size(duration));
                plot_AUC = cell(size(duration));

                for j = 1:numel(AUC)

                    % Calculate AUC between these curves:
                    x  = plot_time(plot_duration{j});
                    y1 = plot_data(plot_duration{j}); % Response
                    y2 = repmat(thresh ,size(x));     % Threshold

%                     % Test cases
%                     thresh = 0.5;
%                     x = [0,1];
%                     y1 = [0,1];
%                     y2 = repmat(thresh,size(x));

                    % 1. Find all points where these curves intersect (not including ends)
                    [x0,y0] = intersections(x,y1,x,y2);

                    % 2. Add end points
                    x_intersect = unique([x(1); x0; x(end)]);

                    % 3. Find y1 value for each intersect
                    x_with_intersects =  [x, x0'];
                    y1_with_intersects = [y1,y0'];
                    [x_with_intersects_sorted, order] = sort(x_with_intersects);
                    y1_with_intersects_sorted = y1_with_intersects(order);
                    %{
                    figure
                    clf
                    hold on
                    plot(x,y1,'o-')
                    plot(x,y2)
                    plot(x_intersect,repmat(thresh,size(x_intersect)),'*')
                    plot(x_with_intersects_sorted,y1_with_intersects_sorted,'x')
                    %}
                    
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





                % ~~~ Figure ~~~

                % Prepare figure
                fig = figure(sub2ind([numel(types),2],2,find(strcmp(type,types))));
                clf(fig)
                fig.WindowStyle = fig_params.WindowStyle;
                fig.Units = fig_params.Units;
                fig.Position = fig_params.Position;

                sig_color = [0.3020    0.6863    0.2902];

                layout = tiledlayout(1,1,'TileSpacing','normal','Padding','normal');
                title(layout, {fig_title; ['Minimum response time (onset, offset) = (',num2str(min_response_time_on),', ',num2str(min_response_time_off),') s']})

                % Plot
                nexttile(layout)
                H = shadedErrorBar(plot_time, plot_data, plot_CI);
                delete(H.edge)
                hold on
                arrayfun(@(j) fill(plot_AUC{j}(1,:), plot_AUC{j}(2,:), sig_color,'EdgeColor','none','FaceAlpha',0.25,'DisplayName', ['AUC ',num2str(j),' = ',num2str(AUC(j))]), 1:numel(AUC))
                arrayfun(@(j) plot(plot_time(plot_duration{j}), plot_data(plot_duration{j}),'.-','Color',sig_color,'DisplayName', ['Duration ',num2str(j),' = ',num2str(duration(j)),' s']), 1:numel(duration))
                xline(0,'Tag','trigger')
                yline(thresh,'Tag','thresh')
                hold off

                % Format
                axis tight

                % Label
                title(['Trials averaged (',num2str(sum(good_trials)),')'])
                xlabel('Time since stimulation [s]')
                ylabel(y_label)
                H.mainLine.DisplayName = 'Average';
                H.mainLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
                if ~isempty(AUC)
                    legend('location','southoutside','Orientation','horizontal','NumColumns',numel(AUC))
                end


                % ~~~~~~~~~~~ Uniform formatting ~~~~~~~~~~~
                uniformFormat(fig)

                % Save
                if save_it

                    % Save plot
                    exportgraphics(fig,fullfile(plot_path,[type{:},'_2_average.pdf']),'Resolution',300,'ContentType','vector')
                    savefig(fig,fullfile(plot_path,[type{:},'_2_average.fig']))

                    % Save duration into condition table
                    warning('Save duration to spreadsheet')
                    %condition_table.("Duration of Significatn Response") =

                end

            end % Loop through files


        end % Loop through plot types


        disp(['             Done running ',strjoin({date_str,cond_type,cond_day,stim_type},', '),'.'])



    end % End loop through all conditions

    if size(condition_table,1) == 0
        disp(['No data found for mouse ',mouse_id,' and/or conditions. Skip!'])
    else
        disp(['Done running mouse ',mouse_id,'.'])
    end

    disp(' ')

end % End loop through mice



