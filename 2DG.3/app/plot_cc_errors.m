% Plot errors in convection convergence study.
clear all; close all; clc;



load('cc_nodistort.mat');
lx = log(1./x);
le = log(errs);

% only take last 3 points for fits
r=3:5;
for pp=1:numel(p)
    pc = polyfit(lx(r), le(pp,r), 1);
    fprintf('Order %d Rate %g\n', p(pp), pc(1));
end
ll = cell(numel(p), 1);
for pp=1:numel(p)
    ll{pp} = sprintf('Order %d', p(pp));
end
h = figure;
set(h, 'units', 'inches', 'position', [1 1 4 4])
set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos=get(h,'Position');
set(h, 'PaperSize', [pos(3) pos(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
loglog(1./x, errs');
legend(ll, 'location', 'SouthEast');
xlabel('1/(m, n) [average element length]');
ylabel('L_2 error');
title('Error in Convection Equation');
print('../../report/cc_err_nodistort.pdf', '-dpdf');

clear all;
load('cc_distort.mat');
lx = log(1./x);
le = log(errs);

% only take last 3 points for fits
r=3:5;
for pp=1:numel(p)
    pc = polyfit(lx(r), le(pp,r), 1);
    fprintf('Order %d Rate %g\n', p(pp), pc(1));
end
ll = cell(numel(p), 1);
for pp=1:numel(p)
    ll{pp} = sprintf('Order %d', p(pp));
end
h = figure;
set(h, 'units', 'inches', 'position', [1 1 4 4])
set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos=get(h,'Position');
set(h, 'PaperSize', [pos(3) pos(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
loglog(1./x, errs');
legend(ll, 'location', 'SouthEast');
xlabel('1/(m, n) [average element length]');
ylabel('L_2 error');
title('Error in Convection Equation (distorted)');
print('../../report/cc_err_distort.pdf', '-dpdf');