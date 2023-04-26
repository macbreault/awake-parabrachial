function arg = getArgs(varargin)
% function varargout = getArgs(varargin)
%
% Function that extracts the arguments and returns arguments to caller function as output and assignin
% Assumes format of varargin as: {'variable_name', variable_value}
% Which would output of: variable_name = variable_value
%
%
% Inputs:
%       - varargin: < 1 x V > cell array with arugments with format of {'variable_name',variable value} OR structure array such that varargin{:}.name = value
%
% Output:
%       - arg: < 1 x V > structure array with names and values of each argument
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 08-07-2018
% Modified: 01-03-2020 - Added ability to use structure format for aguments such as varargin.name = value


%% Initialize variables

argin = varargin{:};

% Check format of argin
if isstruct(argin)
    argin = cellfun(@(name,value) {name,value}, fieldnames(argin), struct2cell(argin), 'un',0);
end

arg = struct('name',[],'value',[]);



%% Check arguments

if numel(argin) > 0
   for v = 1:length(argin)
       
       % If the first argument is a string, then assign a variable name to that string to the second
       % argument
       if ischar(argin{v}{1})
           
           arg(v).name = argin{v}{1};
           arg(v).value = argin{v}{2};
           
           assignin('caller',argin{v}{1}, argin{v}{2})
           
       end
       
   end
end


end % end getArg