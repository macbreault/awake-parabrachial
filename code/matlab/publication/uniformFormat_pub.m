function [] = uniformFormat_pub(fig)
% function [] = uniformFormat_pub(fig)
%
% Function with figure as an input that makes the plot uniformly formatted for PUBLICATION QUALITY
%
% Inputs:
%   - fig = figure handle (optional) If not provided, uses gcf
%
% Outputs:
%   - None, just a nice looking figure!
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 10-05-2022


%% Check argument

if isempty(whos('fig'))
    fig = gcf;
end



%% Intialize formating variables

% Set colors
colors = struct('black', [0,0,0],...
                'white', [1,1,1],...
                'grey',  [0.5,0.5,0.5],...
                'gray',  [0.5,0.5,0.5],...
                'purple',[0.5961    0.3059    0.6392],... % From cbrewer
                'red',   [0.8941    0.1020    0.1098],... % From cbrewer
                'blue',  [0.2157    0.4941    0.7216],... % From cbrewer
                'brown', [0.6510    0.3373    0.1569],... % From cbrewer
                'green', [0.3020    0.6863    0.2902]... % From cbrewer
                ); 


% Set variables
font_name = 'Arial';
font_size = 7;
font_weight = 'normal';
font_angle = 'normal';
font_color = colors.black;

legend_font_size = 5;
legend_font_weight = 'normal';

line_width = 0.75;
constant_line_width = 1;

marker_sz = 3;
marker_line_width = 1;

axes_line_width = 1;
axes_line_color = colors.black;
axes_back_color = 'none';

tick_dir = 'out'; % 'in' | 'out' | 'both' | 'none'
tick_length = 0.15; % cm

box_plot_line_width = 1;


font_info = struct('panel', struct('FontName','Arial','FontSize',10,'FontWeight','bold','FontAngle','normal','Color',colors.black));



%% Format AXs

AX = findobj(fig,'Type','Axes');

set(AX,'LineWidth', axes_line_width)
set(AX,'Color',     axes_back_color)
set(AX,'Box',       'off')

set(findobj(AX,'flat','-not','XColor','none','-not','UserData','keep_color'),'XColor',axes_line_color)
set(findobj(AX,'flat','-not','YColor','none','-not','UserData','keep_color'),'YColor',axes_line_color)
set(findobj(AX,'flat','-not','ZColor','none','-not','UserData','keep_color'),'ZColor',axes_line_color)


%% Format ticks

set(AX,'TickDir',tick_dir)

% Set tick length for each axes
for i = 1:numel(AX)
    

    % 1. Change units of AX and get largest size
    old_units = AX(i).Units;
    AX(i).Units = 'centimeters';


    % 2. Find appropriate size
    tick_sz = tick_length / max(AX(i).Position(3:4));

    % 3. Set size
    set(AX(i),'TickLength', [tick_sz, 1])
    set(AX(i),'Units',old_units)

end
%set(AX,'TickLength',[0.05, 0.05])

% warning('Set tick length')



%% LINES (not markers?)

% Format lines that:
format_lines = findobj(fig);

% 1. Have Line Width as a property
format_lines = findobj(format_lines,'flat','-property','LineWidth');

% 2. Are not tagged as edge (for nondirected graphs)
format_lines = findobj(format_lines,'flat','-not','Tag','edge');

% 3. Are not text
format_lines = findobj(format_lines,'flat','-not','Type','Text');

% 4. If they have a marker property, it is set to 'none'. Otherwise, keep objects without marker as a property
format_lines = [findobj(format_lines,'flat','-not','-property','Marker');...
                findobj(format_lines,'flat','Marker','none')];

% 5. They are not Axes
format_lines = findobj(format_lines,'flat','-not','Type','Axes');

% 5. And they are not an lsline
format_lines = findobj(format_lines,'flat','-not','Tag','lsline');

% End. Format
set(format_lines,'LineWidth',line_width)



%% Markers

format_markers = findall(fig,'Marker','o');
set(format_markers, 'MarkerSize', marker_sz)
%set(format_markers, 'LineWidth', marker_line_width)



%% Text (general)

text_ignore = 'stat_text';

set(findall(fig,'-property','FontName','-not','Tag',text_ignore),'FontName', font_name)
set(findall(fig,'-property','FontAngle','-not','Tag',text_ignore),'FontAngle', font_angle)
set(findall(fig,'-property','FontWeight','-not','Tag',text_ignore),'FontWeight', font_weight)
set(findall(fig,'-property','FontSize','-not','Tag',text_ignore),'FontSize', font_size)

set(findall(fig,'Type','Text'),'Color',colors.black)



%% Format Legends

lgs = findall(fig,'Type','Legend');
set(lgs,'FontSize', legend_font_size)
set(lgs,'FontWeight', legend_font_weight)
set(lgs,'Box','off')

% Make icon in legend smaller
for i = 1:numel(lgs)
    entries = lgs(i).EntryContainer.Children;
    for entry = entries(:)'
        T = entry.Icon.Transform;
        T.Matrix(1) = T.Matrix(1) / 2;
    end
end

% Remove legend on things tagged 'no_legend'
remove_legend = findall(fig,'Tag','no_legend');
arrayfun(@(i) set(get(get(remove_legend(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off'), 1:numel(remove_legend))



%% Sublabels (found as non-empty AX.Tag)

arrayfun(@(i) addSublabel(AX(i),AX(i).Tag, font_info), find(~cellfun(@isempty,{AX.Tag})));



%% Format specific TAGS

% ~~~~~~~~~~~~~~~~~~~~~~~ shadedErrorBar ~~~~~~~~~~~~~~~~~~~~~~~
if ~isempty(findobj(fig,'Tag','shadedErrorBar_mainLine'))

    % shadedErrorBar_mainLine
    mainLines = findall(fig,'Tag','shadedErrorBar_mainLine');

    % shadedErrorBar_edge
    edges = findall(fig,'Tag','shadedErrorBar_edge');
    delete(edges)

    % shadedErrorBar_patch
    patches = findall(fig,'Tag','shadedErrorBar_patch');

end


% ~~~~~~~~~~~~~~~~~~~~~~~ notBoxPlot ~~~~~~~~~~~~~~~~~~~~~~~
if ~isempty(findobj(fig,'Tag','mu'))
    
    % sdPtch
    delete(findobj(fig,'Tag','sdPtch'))

    % semPtch
    semPtch = findobj(fig,'Tag','semPtch');
    for i = 1:numel(semPtch)
        center = mean(semPtch(i).XData(1:2));
        hold(semPtch(i).Parent,'on')
        plot(semPtch(i).Parent, center+((semPtch(i).XData(1:2)-center)*0.5), semPtch(i).YData(1:2),'-','Tag','mu')
        plot(semPtch(i).Parent, center+((semPtch(i).XData(3:4)-center)*0.5), semPtch(i).YData(3:4),'-','Tag','mu')
        plot(semPtch(i).Parent, [mean(semPtch(i).XData(1:2)); mean(semPtch(i).XData(3:4))], unique(semPtch(i).YData),'-','Tag','mu')
        hold(semPtch(i).Parent,'off')
    end
    delete(findobj(fig,'Tag','semPtch'))
    
    % mu
    set(findobj(fig,'Tag','mu'),{'Color','LineWidth'}, {colors.black, box_plot_line_width})
    cellfun(@(h) set(h.LegendInformation,'IconDisplayStyle','off'),get(findobj(fig,'Tag','mu'),'Annotation'),'un',0)
    
    % data
     set(findobj(fig,'Tag','data'),{'MarkerFaceColor','MarkerEdgeColor','MarkerSize'}, {'none', colors.grey, 5})
     data = findobj(fig,'Tag','data');
     arrayfun(@(i) uistack(data(i),'bottom'), 1:numel(data))

end



% ~~~~~~~~~~~~~~~~~~~~~~~ ConstantLine ~~~~~~~~~~~~~~~~~~~~~~~
constant_line = findobj(fig,'Type','ConstantLine');

if ~isempty(constant_line)

    for i = 1:numel(constant_line)

        switch constant_line(i).Tag
            case 'trigger'
                set(constant_line(i),...
                    {'Color', 'LineStyle','LineWidth','Alpha'},...
                    {colors.black, ':', 0.5, 1})

            case 'thresh'
                set(constant_line(i),...
                   {'Color', 'LineStyle','LineWidth','Alpha'},...
                   {colors.brown, ':', 0.5, 1})

            otherwise

        end

        set(get(get(constant_line(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')

    end
    
end



% ~~~~~~~~~~~~~~~~~~~~~~~ lsline ~~~~~~~~~~~~~~~~~~~~~~~

ls_line = findobj(fig,'Tag','lsline');

for i = 1:numel(ls_line)
    
    x_lim = ls_line(i).Parent.XLim;
    y_lim = ls_line(i).Parent.YLim;

    % Extend lsline (assuming axes limits are correct)
    p = polyfit(ls_line(i).XData, ls_line(i).YData, 1);
    XData_new = x_lim;
    YData_new = p(1).*XData_new + p(2);
    
    set(ls_line(i),...
                    {'Color', 'LineStyle','Marker','XData','YData'},...
                    {colors.black,'-','none',XData_new,YData_new})

    set(ls_line(i).Parent,{'XLim','YLim'},{x_lim,y_lim})
    
    %set(get(get(ls_line(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')

    % Send to back
    uistack(ls_line(i),'top')

end




% ~~~~~~~~~~~~~~~~~~~~~~~ rectangle ~~~~~~~~~~~~~~~~~~~~~~~

rect = findobj(fig,'Type','rectangle');

for i = 1:numel(rect)

    % 1. Change y and height to fit new axes size where: rect(i).Position = [x y w h];
    y_lim_new = rect(i).Parent.YLim;
    rect(i).Position(2) = y_lim_new(1);
    rect(i).Position(4) = sum(abs(y_lim_new));

    % 2. Format
    rect(i).FaceColor = [1.0000    1.0000    0.8980];%colors.brown;
    rect(i).EdgeColor = 'none';


    % 3. Send to back
    uistack(rect(i),'bottom')

end



%% Change COLORS

color_properties = {'Color','MarkerFaceColor','MarkerEdgeColor','FaceColor','EdgeColor','YColor','XColor'};

% Red
cellfun(@(property) set(findall(fig, property, 'r'), property, colors.red), color_properties)

% Blue
cellfun(@(property) set(findall(fig, property, 'b'), property, colors.blue), color_properties)

% Green
cellfun(@(property) set(findall(fig, property, 'g'), property, colors.green), color_properties)

% Purple
cellfun(@(property) set(findall(fig, property, 'y'), property, colors.purple), color_properties)



end % end uniformFormat

