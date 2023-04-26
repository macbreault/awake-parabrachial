function txt = addSublabel(ax, str, font_info)
% function txt = addSublabel(ax, str, font_info)
%
%
% Function that adds a sublabel to the upper left corner of the given axis
%
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 06-15-2021

%{
annotation = figure property
text = axes property
        txt = text(ax,min(ax.XLim),min(ax.YLim),str);
%}


%% Get figure property

fig = ancestor(ax,'Figure');
fig_units = fig.Units;



%% Add

set(ax,'Units',fig_units)
txt = annotation(fig,'textbox','String',str,'Units',fig_units);



%% Format
field_names = fieldnames(font_info.panel)';
for field = field_names
    txt.(field{:}) = font_info.panel.(field{:});
end

txt.BackgroundColor = 'none';
txt.EdgeColor = 'none';
txt.Margin = 0;



%% Position witin container            % [left bottom width height]

is_tile = strcmp(class(ax.Parent.Parent),'matlab.graphics.layout.TiledChartLayout'); % Boolean if object is in a TILE

% If AX is NOT in a nested Tilelayout, do the following (aka container is Figure)
if ~is_tile
    
    set(ax,'Units',fig_units)
    set(txt,'Units',fig_units)
    
    %warning('Use TightInset instead of OuterPosition?')
    %warning('Also, check that parent is figure?')
    
    ax_outerposition = ax.OuterPosition;
    pause(0.01)
    
    txt.Position = [ax_outerposition(1),...
        ax_outerposition(2) + ax_outerposition(4) - txt.Position(4),...
        txt.Position(3),...
        txt.Position(4)];
    
    
% aka (aka container is TiledLayout)    
elseif is_tile

    % Check to see if tile is TILE in TILE
    is_tile_in_tile = strcmp(class(ax.Parent.Parent.Parent),'matlab.graphics.layout.TiledChartLayout'); % Boolean if object is in a TILE
    

    if ~is_tile_in_tile

        parent = ax.Parent;
        parent.Units = fig_units;

        % Convert axes position relative to TiledLayout parent to position relative to Figure
        % [left bottom width height]
        ax_outerposition_parent = ax.OuterPosition;
        pause(0.01)

        % Left
        ax_outerposition_figure(1) = ax_outerposition_parent(1) + parent.OuterPosition(1);
        % Bottom
        ax_outerposition_figure(2) = ax_outerposition_parent(2) + parent.OuterPosition(2);
        % Width
        ax_outerposition_figure(3) = ax_outerposition_parent(3);
        % Height
        ax_outerposition_figure(4) = ax_outerposition_parent(4);


        txt.Position = [ax_outerposition_figure(1),...
            ax_outerposition_figure(2) + ax_outerposition_figure(4) - txt.Position(4),...
            txt.Position(3),...
            txt.Position(4)];

    else

        parent = ax.Parent;
        parent.Units = fig_units;

        grandparent = parent.Parent;
        grandparent.Units = fig_units;

        txt.Position = parent.OuterPosition + [grandparent.OuterPosition(1), grandparent.OuterPosition(2), 0, 0];


    end

    
end


   
end % End addSublabel