function [] = sameLimits(obj,varargin)
% function [] = sameLimits(fig,varargin)
%
% Function that sets the x and/or y limits the same to all the axes in the figure
%
% Input:
%       - obj: < 1 x 1 > figure graphical object OR multiple axes object (optional) Default is gcf
%       - ARGUMENTS:
%           - 'x' or 'X' - < 1 x 1 > string to only set the x-axis the same
%           - 'y' or 'Y' - < 1 x 1 > string to only set the y-axis the same
%           - 'z' or 'Z' - < 1 x 1 > string to only set the z-axis the same
%
% Output:
%       - A figure
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 10-04-2018

%% Initialize arguments

Lim = {'XLim','YLim','ZLim'}; % Default limits to change
type = [];

% For only 1 input
if (nargin < 2)
    
    % If not inputs, then set to default font size for each device
    if nargin == 0
        
        obj = gcf;
        
        
    % Input was 'x' or 'y'
    elseif ischar(obj)
        
        Lim = validatestring(obj,Lim);
        obj = gcf;
        
    % Input was an figure
    elseif isa(obj,'matlab.ui.Figure')
        type = 'figure';
        
        % Input was an axes
    elseif isa(obj,'matlab.graphics.axis.Axes')
        type = 'axes';
        
    end
    
end

% For multiple inputs
if ~(isa(obj,'matlab.ui.Figure') || isa(obj,'matlab.graphics.axis.Axes')) || ~isempty(varargin)
    
    inputs = {obj,varargin{:}};
    
    % Assign object type
    try
        obj = inputs{cellfun(@(in) isa(in,'matlab.ui.Figure'),inputs)};
        type = 'figure';
        
    catch
        
        try 
            obj = inputs{cellfun(@(in) isa(in,'matlab.graphics.axis.Axes'),inputs)};
            type = 'axes';
            
        catch
            obj = gcf;
        end
        
    end
    
    % Assign Limits
    Lim = {validatestring(inputs{cellfun(@(in) ischar(in),inputs)},Lim)};
    
    
    % Assign varargins to correct fieldname
    varargin = inputs(cellfun(@(in) isa(in,'cell') & ~isempty(in),inputs));

    
end

% Extract arguments
getArgs(varargin);



%% Set font size

% Object handle
if ~isempty(whos('obj'))
    
    if strcmp(type,'figure')
        handles = findobj(obj,'Type','Axes');
    elseif strcmp(type,'axes')
        handles = obj;
    end
    
else
    handles = obj;
end

AX = findobj(handles,'Type','Axes');

if isempty(AX) || (numel(AX) == 1)
    return
end

if strcmp(type,'axes') && numel(Lim) == 3
    % Set YLim and XLim the same
    limits = cellfun(@(lim) get(AX,lim), Lim,'un',0);
    %limits = vertcat(limits{:});
    %cellfun(@(lim) set(AX,lim,[min(limits(:,1)),max(limits(:,2))]),Lim)
    arrayfun(@(l) set(AX,Lim{l},[min([limits{l}{:}]),max([limits{l}{:}])]),1:numel(Lim))
    
else
    cellfun(@(lim) set(AX,lim,[min(min(vertcat(AX.(lim)))), max(max(vertcat(AX.(lim))))]),Lim)
end

% Fix ticks (but skip if the ticks are the same across axes)
if any(contains(Lim,'X')) && ~isempty([AX.XTick]) && ~isequal(AX.XTick) && ~any(cellfun(@isempty,{AX.XTick}))
    set(AX,'XTick',linspace(min(cellfun(@min,{AX.XTick})),max(cellfun(@max,{AX.XTick})),max(cellfun(@numel,{AX.XTick}))))
end

if any(contains(Lim,'Y')) && ~isempty([AX.YTick]) && ~isequal(AX.YTick) && ~any(cellfun(@isempty,{AX.YTick}))
    set(AX,'YTick',linspace(min(cellfun(@min,{AX.YTick})),max(cellfun(@max,{AX.YTick})),max(cellfun(@numel,{AX.YTick}))))
end

if any(contains(Lim,'Z')) && ~isempty([AX.ZTick]) && ~isequal(AX.ZTick) && ~any(cellfun(@isempty,{AX.ZTick}))
    set(AX,'ZTick',linspace(min(cellfun(@min,{AX.ZTick})),max(cellfun(@max,{AX.ZTick})),max(cellfun(@numel,{AX.ZTick}))))
end

% Set same ticks if X and Y are the SAME
if any(contains(Lim,'X')) && any(contains(Lim,'Y')) && ~isempty([AX.XTick]) && ~isempty([AX.YTick]) &&...
        all(min(vertcat(AX.XTick),[],2) == min(vertcat(AX.YTick),[],2)) && ...
        all(max(vertcat(AX.XTick),[],2) == max(vertcat(AX.YTick),[],2)) && ...
        ~any(cellfun(@isempty,{AX.XTick})) && ~any(cellfun(@isempty,{AX.YTick}))
    
    % Set same axes based on which axes has the largest/smallest tick
    [min_tick,min_lim] = min([min([AX.XTick]),min([AX.YTick])]);
    [max_tick,max_lim] = max([max([AX.XTick]),max([AX.YTick])]);
    numel_tick = [max([AX.XTick]),max([AX.YTick])];
    
    % Save YTick
    if (min_lim == 2) && (max_lim == 2)
        
        new_ticks = AX.YTick;
        
    % Save XTick   
    elseif (min_lim == 1) && (max_lim == 1)
        
        new_ticks = AX.XTick;
     
    % Make new Tick
    else
        
        new_ticks = linspace(min_tick,max_tick,max(numel_tick));
        
    end
    
    % Set new ticks
    arrayfun(@(ax) set(AX(ax),'XTick',new_ticks),1:numel(AX))
    arrayfun(@(ax) set(AX(ax),'YTick',new_ticks),1:numel(AX))
    
end



end % end sameLimits