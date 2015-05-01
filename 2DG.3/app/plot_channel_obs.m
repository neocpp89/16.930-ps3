clear all; close all; clc;

load('channel_obstacle.mat', 'meshes');
load('channel_obstacle.mat', 'solutions');
load('channel_obstacle.mat', 'p');

for i=1:numel(p)
    h = figure;
    set(h, 'units', 'inches', 'position', [1 1 4 4])
    set(h, 'PaperUnits','centimeters');
    set(h, 'Units','centimeters');
    pos=get(h,'Position');
    set(h, 'PaperSize', [pos(3) pos(4)]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
    subplot(2,1,1);
    title(sprintf('Mach number at t=100 for order %d', p(i)));
    scaplot(meshes{i}, eulereval(solutions{i}, 'M', 1.4), [], [], 1)
    subplot(2,1,2);
    title(sprintf('Pressure at t=100 for order %d', p(i)));
    scaplot(meshes{i}, eulereval(solutions{i}, 'p', 1.4), [], [], 1)
    print(sprintf('../../report/channel_obs_%d.pdf', p(i)), '-dpdf');
end
