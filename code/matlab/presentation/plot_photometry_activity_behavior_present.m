% plot_photometry_table_present.m
%
% Script that plots behavioral and activity results from photometry experiments for each animal
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Modified: 10-10-2022

clear
clc


%% Set data

% ~~~~~~~~~~~ Set mouse information ~~~~~~~~~~~



% Other variables
save_it = 1;




% ~~~~~~~~~~~~~~~~~~~~~~ DO NOT CHANGE BELOW ~~~~~~~~~~~~~~~~~~~~~~

%%%%%% Set paths %%%%%%
paths = struct('dropbox',[]);

if contains(pwd,'/mac/')
    paths.dropbox = '/Users/mac/Dropbox (MIT)/Jesse';
else
    paths.dropbox = '/Users/Jessesmith 1/Dropbox/Jesse';
end

addpath(genpath(fullfile(paths.dropbox,'code'))); % Add libraries to directory

plot_path  = fullfile(paths.dropbox,'plots',findPlotFolder,date,'photometry','activity_behavior');
if ~exist(plot_path,'dir') && save_it; mkdir(plot_path); end


% Figure configurations
fig_params = struct('Position',[0 0 17.8 8.9], 'Units','centimeters', 'WindowStyle','normal', 'FontSize', 10);




%% Load data

% Get file names
behavior_file = fullfile(paths.dropbox,'tables','Photometry Behavior.xlsx');
trace_file = fullfile(paths.dropbox,'tables','Photometry Traces.xlsx');

% Get types of traces to plot
trace_types = sheetnames(trace_file);
trace_types = trace_types(~strcmp(sheetnames(trace_file),'TEST'));

if ~save_it
    trace_types = trace_types(1:2,:);
end

% Load data
behavior_all = readtable(behavior_file,'ReadRowNames',1,'PreserveVariableNames',1);
traces_all = struct();
for i = 1:numel(trace_types)
    traces_all.(strrep(trace_types(i,:),' ','_')) = readtable(trace_file,'Sheet',trace_types(i,:),'ReadRowNames',1,'PreserveVariableNames',1);
end



%% For each type of condition...

for i = 1:numel(trace_types)

    %% Initialize figure

    fig = figure(1);
    clf(fig)
    fig.WindowStyle = fig_params.WindowStyle;
    fig.Units = fig_params.Units;
    fig.Position = fig_params.Position;

    layout_main = tiledlayout(1,3,'TileSpacing','normal','Padding','none');



    %% ~~~~~~~~~~ A. ~~~~~~~~~~

    % ~~~~~~~~~  Get trace data ~~~~~~~~~
    mouse_ids = behavior_all.Properties.VariableNames;
    trace_type = strrep(trace_types{i},' ','_');
    trace = traces_all.(trace_type);

    % Get pre-CFA
    pre_days = trace.Properties.RowNames(cellfun(@str2num, trace.Properties.RowNames) <= 0);

    % ONE POINT PER ANIMAL
    %males = mean(trace{pre_days, contains(mouse_ids,'M')},1,'omitnan');
    %females = mean(trace{pre_days, contains(mouse_ids,'F')},1,'omitnan');

    % ALL POINTS PER ANIMAL
    males   = trace{pre_days, contains(mouse_ids,'M')}(:);
    females = trace{pre_days, contains(mouse_ids,'F')}(:);


    % ~~~~~~~~~ Plot ~~~~~~~~~ 
    ax = nexttile(layout_main,1);
    hold on

    % Males
    H_male = notBoxPlot(males, 1 * ones(size(males)));
    delete(H_male.sdPtch)
    H_male.mu.Color = 'k';
    H_male.semPtch.EdgeColor = 'none';
    H_male.semPtch.FaceColor = [0.75 0.75 0.75];
    if ~isempty(H_male.data)
        H_male.data.Color = 'k';
        H_male.data.MarkerEdgeColor = 'k';
        H_male.data.MarkerFaceColor = 'none';
        H_male.data.LineWidth = 1;
        H_male.data.MarkerSize = 3;
    end


    % Females
    H_female = notBoxPlot(females, 2 * ones(size(females)));
    delete(H_female.sdPtch)
    H_female.mu.Color = 'k';
    H_female.semPtch.EdgeColor = 'none';
    H_female.semPtch.FaceColor = [0.75 0.75 0.75];
    if ~isempty(H_female.data)
        H_female.data.Color = 'k';
        H_female.data.MarkerEdgeColor = 'k';
        H_female.data.MarkerFaceColor = 'none';
        H_female.data.LineWidth = 1;
        H_female.data.MarkerSize = 3;
    end


    % ~~~~~~~~~ Format ~~~~~~~~~
    hold off
    axis padded
    box off
    xticks([1,2])


    % ~~~~~~~~~ Label ~~~~~~~~~ 
    ylabel(trace_types{i})
    xticklabels({'Males','Females'})
    xtickangle(0)
    title('Pre CFA')
    ax.Tag = 'A';




    %% ~~~~~~~~~~ B. ~~~~~~~~~~

    % CHOOSE which mice to plot
    %mouse_ids = sort({'F2','M1','F3','M3'});
    mouse_ids = flip(sort(behavior_all.Properties.VariableNames));

    % Prepare new layout
    tile_b = nexttile(layout_main,[1,2]); axis off;
    layout_b = tiledlayout(layout_main,numel(mouse_ids)/2,2,'TileSpacing','compact','Padding','none','TileIndexing','columnmajor');
    layout_b.Layout.Tile = 2;
    layout_b.Layout.TileSpan = [1 2];


    % Loop through mice to plot
    for mouse = 1:numel(mouse_ids)

        mouse_id = mouse_ids{mouse};

        %~~~~~~ Get data for plotting ~~~~~~

        % Behavior
        mouse_ind = strcmp(behavior_all.Properties.VariableNames,mouse_id);
        plot_time = cellfun(@str2num, behavior_all(:,mouse_ind).Properties.RowNames);
        plot_data = behavior_all{:,mouse_ind};
        pain_thresh = min(plot_data(plot_time <= 0)); %%%%%%%%%%%% CHANGE!!!
        warning('Add proper formual for mechanical threshold')
        behave_table = [plot_time, plot_data];
        pain_ind = plot_data < pain_thresh;
        behave_table(~pain_ind,:) = []; % Remove days without pain

        % Traces
        trace_type = strrep(trace_types{i},' ','_');
        trace = traces_all.(trace_type);
        mouse_ind = strcmp(trace.Properties.VariableNames,mouse_id);
        plot_time = cellfun(@str2num, trace(:,mouse_ind).Properties.RowNames);
        plot_data = trace{:,mouse_id};
        plot_CI = calculateCI(plot_data(0 > plot_time)'); % Calculate CI using [-1,-2,-3]

        trace_table = [plot_time, plot_data];
        

        % Plot
        ax = nexttile(layout_b);
        hold on
        plot(trace_table(:,1), trace_table(:,2), '-','Color','k','Tag','no_legend') % Plot line traces
        plot(trace_table(:,1), trace_table(:,2), 'o', 'MarkerFaceColor','k','MarkerEdgeColor','none','DisplayName','No pain') % Plot pain days
        
        % Overlay days of pain
        trace_2_behave = ismember(trace_table(:,1), behave_table(:,1));
        plot(trace_table(trace_2_behave,1), trace_table(trace_2_behave,2), 'o', 'MarkerFaceColor','y','MarkerEdgeColor','none','DisplayName','Pain') % Plot pain days
        
        if ~any(isnan(plot_CI))
            yline(plot_CI,'Tag','thresh')
        end
        xline(0,'Tag','trigger')
        hold off


        % Format
        axis padded
        box off

        % Label
        title(mouse_id)
        if any((mouse == [1,((numel(mouse_ids)/2)+1)]))
            if contains(mouse_id,'F')
                title({'Female',mouse_id})
            elseif contains(mouse_id,'M')
                title({'Male',mouse_id})
            end
        end
        if mouse == 1
            ax.Tag = 'B';
        end
        if mouse == ceil((numel(mouse_ids)/2) + (numel(mouse_ids)/2/2))
            lg = legend('Orientation','vertical','Box','off');
            lg.Location = 'eastoutside';
            %lg.Layout.Tile = 'east';
        end


    end % Loop through mice


    % Format
    AX = findobj(layout_b,'Type','Axes');
    % sameLimits(AX)
    
    % Label
    ylabel(layout_b,trace_types{i})
    xlabel(layout_b,'Days Post CFA')




    %% Uniform plotting

    uniformFormat_pub(fig)

    % Save
    if save_it
        exportgraphics(fig,fullfile(plot_path,[trace_types{i},'.pdf']),'Resolution',300,'ContentType','vector')
        savefig(fig,fullfile(plot_path,[trace_types{i},'.fig']))
    end

end % Loop through trace types




