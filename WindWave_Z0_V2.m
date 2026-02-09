%%% PSU + Hopkins: Find Z_0

clear; close all; clc;

means_path = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means';
caze       = 'WT4_WVD_AG0';

means = load(fullfile(means_path, strcat(caze, '_MEANS.mat')));
means = means.output;

x  = unique(means.X) * 1E-3;
y  = flipud(unique(means.Y)) * 1E-3;
u  = mean(means.ensemble.u, 2, 'omitnan');
uv = mean(means.ensemble.uv, 2, 'omitnan');

kinematic_viscosity = 1.48E-5;
u_star = sqrt(max(-1 * uv));
u_plus = u / u_star;
y_plus = (y * u_star) / kinematic_viscosity;

mask = y_plus > 0;
y_plus = y_plus(mask);
u_plus = u_plus(mask);


upper_limit = length(y_plus)-1;
options     = optimset('MaxIter', 1E5, 'Display', 'off');

for limit = upper_limit:-1:1
    disp(limit)
    LL_fun  = @(a) a(1) * log(y_plus(end - limit:end)) + a(2);
    Fit_fun = @(a) sqrt(sum((LL_fun(a) - u_plus(end - limit:end)).^2, 'omitnan'));
    params  = fminsearch(Fit_fun, [1,1], options);
    z0      = params(2);
    if sign(z0) == 1
        break
    end
end



disp(z0)

LL_fun_plot = @(a,b,c) a * log(c) + b;
y_plus_plot = 0:max(y_plus(end - limit:end));

lw = 2;
figure()
hold on
plot(y_plus, u_plus, 'LineWidth', lw)
plot(y_plus_plot, LL_fun_plot(params(1),params(2),y_plus_plot), 'LineWidth', lw, 'LineStyle', '--')
hold off
set(gca, 'XScale', 'Log')
xlabel('$y^+$', 'Interpreter', 'latex')
ylabel('$u^+$', 'Interpreter', 'latex')
xlim([0,max(y_plus)])
ylim([0,30])
title(caze, 'interpreter', 'none')
grid on







