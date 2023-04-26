function CR = calculateCR(x, test_it)
% function CR = calculateCR(x, test_it)
%
% Function that caluculates the RANGE of 95% confidence interval of the vector/matrix of x
%
% Input:
%       - x: < N x M > double of values of N samples across M trials
%       - test_it: < 1 x 1 > boolean whether to test CI function (1) or not (0)
%
% Output:
%       - CR: < N x 1 > double of confidenct interval of N samples
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 09-15-2022
% MATLAB R2021b


%% Validate attributes

CI_type = 'usual'; % Can use 'std' for 2 std away or 'student' to use student t-test or 'usual' for way that CI anonymous function used

% X
validateattributes(x,{'double'},{'2d'})
if isempty(x)
    x = NaN;
end


% PROBABILITIES
alpha = 0.95;
conf_level_value = ([1,-1].*(1-alpha)/2) + [0,1]; % Should be [0.025, 0.975]
z_score = 1.96;
std_factor = 2; % How many standard deviations away


% Check test cases to make sure CI si running correctly
if exist('test_it','var') && test_it
    test = getTestCases;
    for i = 1:numel(test)
        %if ~isequal(test(i).CI, calculateCI(test(i).x,0))
        if ~all(round(abs(test(i).CR - calculateCR(test(i).x,0))) == 0,'all')
            error('Error in %s: CI does not match test cases. Check!',mfilename)
        end
    end
end



%% Calculate CR based on if x is a vector or matrix

if isvector(x)

    switch CI_type
        case 'usual'

            sample_size = sum(~isnan(x));

            CR = 1.96 * (std(x,'omitnan') ./ sqrt(sample_size));

        case 'student'

            sample_size = sum(~isnan(x));

            % Standard Error
            SEM = std(x,'omitnan')/sqrt(sample_size);

            % T-score (student t-distribution)
            ts = tinv(conf_level_value, sample_size-1);

            % CR
            CR = max(ts)*SEM;

        case 'std'

            CR = std_factor * std(x,'omitnan');

        otherwise
            error('Error in %s: Method to calculate confidence range undefined',mfilename)
    end


elseif ismatrix(x)

    switch CI_type

        case 'usual'

            sample_size = sum(~isnan(x),2);

            CR = 1.96 * (std(x,[],2,'omitnan') ./ sqrt(sample_size));

        case 'student'

            sample_size = sum(~isnan(x),2);

            % Standard Error
            SEM = std(x,[],2,'omitnan') ./ sqrt(sample_size);

            % T-score (student t-distribution)
            ts = cell2mat(arrayfun(@(n) tinv(conf_level_value, sample_size(n)-1), (1:numel(sample_size))', 'un',0));

            % CR
            CR = max(ts,[],2).*SEM;

        case 'std'

            CR = std_factor * std(x,[],2,'omitnan');

        otherwise
            error('Error in %s: Method to calculate confidence range undefined',mfilename)
    end

else
    error('Error in %s: x must be a vector or matrix',mfilename)
end


% Check CR before returning
validateattributes(CR,{'double'},{'size',[NaN,1]})





%% Function that gets test cases
function test = getTestCases

    error('These are only true when CI us being calculated using student t-test')
    test = struct('x',[],'CR',[]);
    
    % < 1 x 10 >
    test(1).x = (1:10);
    test(1).CR = 2.1659;
    
    % < 10 x 1 >
    test(2).x = (1:10)';
    test(2).CR = 2.1659;
    
    % < 10 x 3 >
    test(3).x = (1:3).*repmat((1:10)',1,3);
    test(3).CR =   [2.48413771175033
                    4.96827542350066
                    7.45241313525099
                    9.93655084700132
                    12.4206885587516
                    14.9048262705020
                    17.3889639822523
                    19.8731016940026
                    22.3572394057530
                    24.8413771175033];
    
    % < 3 x 10 >
    test(4).x = (1:3)'.*repmat(1:10,3,1); 
    test(4).CR =   [2.1659
                    4.3317
                    6.4976];

end % End nested function



end % End main