% check_blink.m
%
% Script that finds blinking cut-off parameter
%
% MATLAB R2021b
% Copyright (c) 2023, Macauley Breault  breault@mit.edu
% Created: 07-02-2022


blink = double(proc.blink.area);
blink_mean = mean(blink);
blink_std = std(blink);

thresh = blink_mean - 3*blink_std;
blink_ind = (blink < thresh);


% Get data
noblink = blink;
noblink(blink_ind) = NaN;
parea_blink = double(proc.pupil.area);
parea_noblink = parea_blink;
parea_noblink(blink_ind) = NaN;



%% Plot distribution of blink to find cut-off
figure(1)
clf

subplot(2,3,[1,4])
hold on
histogram(blink)
xline(blink_mean,'r','LineWidth',2)
xline(thresh,'r','LineWidth',2)
xlabel('Pixels (blinking)')
ylabel('Count')
axis square
hold off




%% Plot on actual data
ax1 = subplot(2,3,[2,3]);
hold on
plot(ax1, blink,'-b')
plot(ax1, noblink,'-r')
yline(thresh,'--k','LineWidth',2)
xlabel('Frame number')
ylabel('Blink (pixels)')
hold off

ax2 = subplot(2,3,[5,6]);
hold on
plot(ax2, parea_blink,'b')
plot(ax2, parea_noblink,'r')
xlabel('Frame number')
ylabel('Pupil (pixels)')
hold off

linkaxes([ax1,ax2],'x')


