%%% PSU + Hopkins: Find u*

clear; close all; clc;

means_path = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means';
caze       = 'WT4_WVC_AG0';

means = load(fullfile(means_path, strcat(caze, '_MEANS.mat')));
means = means.output;

x = unique(means.X) * 1E-3;
y = flipud(unique(means.Y)) * 1E-3;
u = mean(means.ensemble.u, 2, 'omitnan');

nan_mask = ~isnan(u);
y = y(nan_mask);
u = u(nan_mask);

% Crop near wall
mask = y > 0.02;
y_mask = y(mask);
u_mask = u(mask);

% Crop freestream
mask = y_mask < 0.11;
y_mask = y_mask(mask);
u_mask = u_mask(mask);



% Fit u*
Kappa = 0.4;
Nu = 1.48E-5;

% Fitting u* and c+
% options = optimset('MaxIter', 1E5, 'Display', 'off');
% LL_fun  = @(a) a(1) * ((1/Kappa) * log(y_mask * a(1) / Nu) + a(2));
% Fit_fun = @(a) sqrt(sum((LL_fun(a) - u_mask).^2, 'omitnan'));
% params  = fminsearch(Fit_fun, [0.5, 5], options);
% u_star  = params(1);
% C = params(2);
% 
% fprintf("u* = %1.3f m/s\nC = %2.3f\n\n", u_star, C)

% Fitting u* only
C = 5;
options = optimset('MaxIter', 1E5, 'Display', 'off');
LL_fun  = @(a) a * ((1/Kappa) * log(y_mask * a / Nu) + C);
Fit_fun = @(a) sqrt(sum((LL_fun(a) - u_mask).^2, 'omitnan'));
params  = fminsearch(Fit_fun, 0.5, options);
u_star  = params(1);

fprintf("u* = %1.3f m/s\nC = %2.3f\n\n", u_star, C)



% Plot
lw = 2;
FS = 18;
LL_fun_plot = @(a,b) a * ((1/Kappa) * log(b * a / Nu) + C);
y_plot = linspace(min(y), max(y), 100);

ax = figure('Name', 'u* Fit');
hold on
plot(y, u, 'LineWidth', lw, 'color', 'black')
plot(y_plot, LL_fun_plot(u_star, y_plot), 'LineWidth', lw, 'LineStyle', '--', 'color', 'red')
xline(min(y_mask), 'linewidth', lw)
xline(max(y_mask), 'linewidth', lw)
hold off
set(gca, 'XScale', 'Log')
xlabel('$y$ [m]', 'Interpreter', 'latex', 'FontSize', FS)
ylabel('$u$ [m/s]', 'Interpreter', 'latex', 'FontSize', FS)
tit = strcat(caze, sprintf("\nu* = %1.3f m/s, C = %2.3f", u_star, C));
title(tit, 'interpreter', 'none', 'FontSize', 12)
grid on







