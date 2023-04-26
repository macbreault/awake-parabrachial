function CI = calculateCI(x, test_it)
% function CI = calculateCI(x, test_it)
%
% Function that caluculates confidence interval of the vector/matrix of x
% Based on being 2 STD
%
% Input:
%       - x: < N x M > double of values of N samples across M trials
%       - test_it: < 1 x 1 > boolean whether to test CI function (1) or not (0)
%
% Output:
%       - CI: < N x 2 > double of confidenct interval of N samples
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 09-15-2022
% MATLAB R2021b


%% Validate attributes

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
        if ~all(round(abs(test(i).CI - calculateCI(test(i).x,0))) == 0,'all')
            error('Error in %s: CI does not match test cases. Check!',mfilename)
        end
    end
end



%% Calculate CI based on if x is a vector or matrix

% Calculate confidence range
CR = calculateCR(x);

% Calculate CI
if isvector(x)

    CI = mean(x,2,'omitnan') + [-1,1] * CR;

    % Check by plotting
    %{
    figure
    clf
    hold on
    histogram(x)
    xline(mean(x,'omitnan'),'k--','LineWidth',2,'DisplayName','mean')
    xline(mean(x,'omitnan')+std(x,'omitnan'),'r-','LineWidth',2,'DisplayName','+1std')
    xline(mean(x,'omitnan')-std(x,'omitnan'),'r-','LineWidth',2,'DisplayName','-1std')
    xline(CI(1), '-b','LineWidth',2,'DisplayName','Lower CI')
    xline(CI(2), '-b','LineWidth',2,'DisplayName','Upper CI')
    hold off
    ylabel('Count')
    xlabel('x')
    legend('location','southoutside')
    %}


elseif ismatrix(x)

    % CI
    CI = mean(x,2,'omitnan') + [-1,1] .* repmat(CR,1,2);

else
    error('Error in %s: x must be a vector or matrix',mfilename)
end


% Check CI before returning
validateattributes(CI,{'double'},{'size',[NaN,2]})




%% Function that gets test cases
function test = getTestCases
    
    error('These tests only work for student t-distribution')
    test = struct('x',[],'CI',[]);
    
    % < 1 x 10 >
    test(1).x = (1:10);
    test(1).CI = [3.3341    7.6659];
    
    % < 10 x 1 >
    test(2).x = (1:10)';
    test(2).CI = [3.3341    7.6659];
    
    % < 10 x 3 >
    test(3).x = (1:3).*repmat((1:10)',1,3);
    test(3).CI =   [-0.484137711750330	4.48413771175033
                    -0.968275423500660	8.96827542350066
                    -1.45241313525099	13.4524131352510
                    -1.93655084700132	17.9365508470013
                    -2.42068855875165	22.4206885587517
                    -2.90482627050198	26.9048262705020
                    -3.38896398225231	31.3889639822523
                    -3.87310169400264	35.8731016940026
                    -4.35723940575297	40.3572394057530
                    -4.84137711750330	44.8413771175033];
    
    % < 3 x 10 >
    test(4).x = (1:3)'.*repmat(1:10,3,1); 
    test(4).CI =   [3.3341    7.6659
                    6.6683   15.3317
                    10.0024   22.9976];

end % End nested function



end % End main