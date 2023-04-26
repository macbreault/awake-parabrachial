function [r_max, lag_max, p_max] = calculateCrossCorr(x, y, bin_width)
% function [r_max, lag_max, p_max] = calculateCrossCorr(x, y, bin_width)
%
% Function that caluculates the cross correlation between x and y using a maximum lag of 3 seconds
%
% Input:
%       - x: Data x
%       - y: Data y
%       - bin_width: Time resolution of data
%
% Output:
%       - r_max: Correlation of shifted data
%       - lag_max: Lag value used to shift x over
%       - p_max: P-value of correlation of shifted data
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 10-15-2022
% MATLAB R2021b


%% Initialize

r_max = NaN;
lag_max = NaN;
p_max = NaN;


% Cut-offs
min_seconds = 3; % in seconds. Minium number of data needed + Maximum number of shifting
max_lag = min_seconds / bin_width; % 3 seconds


% Reshape data
x = reshape(x,[],1);
y = reshape(y,[],1);


% Remove NaNs
if any(isnan(x)) || any(isnan(y))

    nan_ind = isnan(x) | isnan(y);

    warning('Found NaN while cross correlating. Removed %s indices out of %s',num2str(sum(nan_ind)),num2str(length(nan_ind)))
    
    % Remove indices
    x(nan_ind) = [];
    y(nan_ind) = [];

end


% Check that x and y are the same size
if length(x) ~= length(y)
    error('Error in %s: X and Y are not the same length',mfilename)
end


% Check that there is enough data to correlate (using max_lag)
if max_lag > (length(x) - 1)
    warning('Not enough data left to cross-correlate. Skip!')
    return
end




%% Calculate correlation

% Calculate cross-corr
[r,lags] = crosscorr(x,y,'NumLags',max_lag);
[~,i_max] = max(abs(r));
lag_max = lags(i_max); % How much to push x over

% Shift data to maximize correlation
if lag_max > 0
    x_new = nan(1,length(x)+lag_max);
    x_new((1:numel(x))+lag_max) = x;
    y_new = nan(1,length(x)+lag_max);
    y_new(1:numel(y)) = y;
elseif lag_max < 0
    x_new = nan(1,length(x)+abs(lag_max));
    x_new(1:numel(x)) = x;
    y_new = nan(1,length(x)+abs(lag_max));
    y_new((1:numel(y))+abs(lag_max)) = y;
elseif lag_max == 0
    x_new = x;
    y_new = y;
end
[r_max, p_max] = corr(reshape(x_new,[],1),reshape(y_new,[],1),'Type','Pearson','Rows','Pairwise');



%{
        % Lag>0
        x = [0 1 0 0];
        y = [0 0 1 0];
        % Lag<0
        x = [0 0 1 0];
        y = [0 1 0 0];

        max_lag = ceil(numel(x)*.1);
        [r,lags] = crosscorr(x,y,'NumLags',max_lag);
        [r_max,i_max] = max(r);
        lag_max = lags(i_max); % How much to push x over
        
        % Shift data to maximize correlation
        if lag_max > 0
            x_new = nan(1,length(x)+lag_max);
            x_new((1:numel(x))+lag_max) = x;
            y_new = nan(1,length(x)+lag_max);
            y_new(1:numel(y)) = y;
        elseif lag_max < 0
            x_new = nan(1,length(x)+abs(lag_max));
            x_new(1:numel(x)) = x;
            y_new = nan(1,length(x)+abs(lag_max));
            y_new((1:numel(y))+abs(lag_max)) = y;
        end
        corr(x_new',y_new','rows','pairwise')
%}





end % End main