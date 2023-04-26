function [] = playSound()
% function [] = playSound()
%
% Function that makes a sound for times when code is running in the background and I'm doing something else
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 04-01-2021


Data = load('handel.mat');
sound(Data.y, Data.Fs)

% beep


end % end playSound