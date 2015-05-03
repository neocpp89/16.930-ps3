clear all; close all; clc;

load('wave_problem.mat', 'meshes');
load('wave_problem.mat', 'solutions');
load('wave_problem.mat', 'incident');
load('wave_problem.mat', 'pl');

for i=1:numel(pl)
    u = squeeze(solutions{i}(:,3,:));
    ue = squeeze(incident{i}(:,3,:));
    h = figure;
    set(h, 'units', 'inches', 'position', [1 1 4 4])
    set(h, 'PaperUnits','centimeters');
    set(h, 'Units','centimeters');
    pos=get(h,'Position');
    set(h, 'PaperSize', [pos(3) pos(4)]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
    subplot(2,1,1);
    title(sprintf('Scattered Field for order %d', pl(i)));
    colormap('jet');
    scaplot(meshes{i}, u-ue, [-1.2, 1.2], [], 0);
    subplot(2,1,2);
    title(sprintf('|\\eta_s|  at (R = e) for order %d', pl(i)));
    [fidx, polar] = find_nodes_at_radius(meshes{i}, exp(1), 1e-4);
    theta = squeeze(polar(:,2,:));
    for x=1:size(fidx,2)
        if (any(fidx(:,x)) == 0)
            continue;
        end
        hold on;
        tt=theta(:,x);
        uu=u(:,x);
        ee=ue(:,x);
        plot(tt(fidx(:,x)), abs(uu(fidx(:,x)) - ee(fidx(:,x))), '-b.');
    end
    hold off;
    xlim([-pi, pi]);
    print(sprintf('../../report/wave_%d.pdf', pl(i)), '-dpdf');
end
