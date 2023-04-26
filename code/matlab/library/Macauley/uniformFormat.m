function [] = uniformFormat(fig)
% function [] = uniformFormat(fig)
%
% Function with figure as an input that makes the plot uniformly formatted
%
% Inputs:
%   - fig = figure handle (optional) If not provided, uses gcf
%
% Outputs:
%   - None, just a nice looking figure!
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 11-18-2021

%{

% For testing

clf
x = 1:10;
y = rand(10,10);
shadedErrorBar(x,y,{@median,@std});
hold on
yline(0.5,'Tag','thresh')
xline(5,'Tag','0')
set(gca,'Tag','(a)')
fig = gcf;
uniformFormat(fig)

%}


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
font_size = 12;
font_weight = 'normal';
font_angle = 'normal';
font_color = colors.black;

legend_font_size = 10;

line_width = 1.5;
constant_line_width = 1.5;

marker_sz = 5;

axes_line_width = 1;
axes_line_color = colors.black;
axes_back_color = 'none';
tick_dir = 'out'; % 'in' | 'out' | 'both' | 'none'


font_info = struct('panel', struct('FontName','Arial','FontSize',8,'FontWeight','bold','FontAngle','normal','Color',colors.black));



%% Format AXs

AX = findobj(fig,'Type','Axes');

set(AX,'LineWidth',axes_line_width)
set(AX,'Color',axes_back_color)


% % DO NOT color axes that have the color 'none'
% AX_X = findobj(AX,'flat','Visible','on','-not','XColor','none');
% AX_Y = findobj(AX,'flat','Visible','on','-not','YColor','none');
% AX_Z = findobj(AX,'flat','Visible','on','-not','ZColor','none');
% 
% 
% set(AX_X,'XColor',axes_line_color)
% set(AX_Y,'YColor',axes_line_color)
% set(AX_Z,'ZColor',axes_line_color)



%% Format ticks

set(AX,'TickDir',tick_dir)

%warning('Set tick length')



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

% % Format markers that:
% format_markers = findobj(fig);
% 
% % 1. Have Marker as a property
% format_markers = findobj(format_markers,'flat','-property','Marker');
% 
% % 2. Are not tagged as edge (for nondirected graphs)
% % format_markers = findobj(format_markers,'flat','-not','Tag','edge');
% 
% % 3. If they have a marker property, it is NOT set to 'none'
% format_markers = findobj(format_markers,'flat','-not','Marker','none');
% 
% % End. Format
% set(format_markers,'MarkerEdgeColor','none')



%% Text (general)

text_ignore = 'ignore_text';

set(findall(fig,'-property','FontName','-not','UserData',text_ignore),'FontName', font_name)
set(findall(fig,'-property','FontAngle','-not','UserData',text_ignore),'FontAngle', font_angle)
set(findall(fig,'-property','FontWeight','-not','UserData',text_ignore),'FontWeight', font_weight)
set(findall(fig,'-property','FontSize','-not','UserData',text_ignore),'FontSize', font_size)



%% Format Legends

lgs = findall(fig,'Type','Legend');
set(lgs,'FontSize', legend_font_size)

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
if ~isempty(findobj(fig,'Tag','sdPtch'))
    
    % sdPtch
    set(findobj(fig,'Tag','sdPtch'),{'FaceColor','EdgeColor'},{ colors.grey,'none'})
    cellfun(@(h) set(h.LegendInformation,'IconDisplayStyle','off'),get(findobj(fig,'Tag','sdPtch'),'Annotation'),'un',0)
    
    % semPtch
    delete(findobj(fig,'Tag','semPtch'))
    
    % mu
    set(findobj(fig,'Tag','mu'),'Color', colors.black)
    cellfun(@(h) set(h.LegendInformation,'IconDisplayStyle','off'),get(findobj(fig,'Tag','mu'),'Annotation'),'un',0)
    
    % data
    %delete(findobj(fig,'Tag','data'))
end



% ~~~~~~~~~~~~~~~~~~~~~~~ ConstantLine ~~~~~~~~~~~~~~~~~~~~~~~
constant_line = findobj(fig,'Type','ConstantLine');

if ~isempty(constant_line)

    for i = 1:numel(constant_line)

        switch constant_line(i).Tag
            case 'trigger'
                set(constant_line(i),...
                    {'Color', 'LineStyle','LineWidth','Alpha'},...
                    {colors.black, ':', constant_line_width, 1})

            case 'thresh'
                set(constant_line(i),...
                   {'Color', 'LineStyle','LineWidth','Alpha'},...
                   {colors.brown, ':', constant_line_width, 1})

            otherwise

        end

        set(get(get(constant_line(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')

    end
    
    % Send to back
%     arrayfun(@(i) uistack(constant_line(i),'bottom'), 1:numel(constant_line))
    
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
                    {'LineStyle','Marker','XData','YData'},...
                    {'-','none',XData_new,YData_new})

    set(ls_line(i).Parent,{'XLim','YLim'},{x_lim,y_lim})
    
    %set(get(get(ls_line(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')

    % Send to back
    uistack(ls_line(i),'top')

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

