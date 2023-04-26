% run_all.m
%
% Script that runs all the plotting functions
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 08-18-2022


try
    run('run_conditions_corr.m')
end

try
    run('run_evoked.m')
end

try
    run('run_laser_corr.m')
end

try
    run('run_spontaneous.m')
end

try
    run('run_evoked_neural.m')
end

try
    run('run_evoked_pupil.m')
end