%%% Curvilinear gradients and integrals

addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
means_path = fullfile(project_path, 'means');
curvilinear_path = fullfile(project_path, 'curvilinear');

% Case
caze = 'WT8_WVA_AG0';

% Load data
curvilinear = load(fullfile(curvilinear_path, strcat(caze, '_CURVILINEAR.mat')));
curvilinear = curvilinear.output;

% Wave Parameters
wave_type       = caze(strfind(caze, 'WV') + 2);
wave_parameters = readcell("Offshore_Waves.xlsx");
wavelength      = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 2};
amplitude       = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 3};

clear means_path curvilinear_path project_path wave_parameters wave_type

%% Test taking gradients

phase = 2;

% X (xi) in mm
X = curvilinear.X;
x = X(1,:);
dx = mean(diff(X(1,:)));
vertical_lines = curvilinear.phase(phase).vertical_lines;

% Y (zeta) in mm
Y = curvilinear.Y;
y = Y(:,1);
dy = mean(diff(Y(:,1)));
horizontal_lines = curvilinear.phase(phase).horizontal_lines;

u = curvilinear.phase(phase).u;
v = curvilinear.phase(phase).v;


[dxi_dx, dxi_dy] = gradient(vertical_lines, 1, 1);
[dzeta_dx, dzeta_dy] = gradient(horizontal_lines, 1, 1);

% Should only be using dxi_dx and dzeta_dy
% Units are in mm

levels = 10;
figure()
tiledlayout(2,2)

h(1) = nexttile;
hold on
contourf(vertical_lines, horizontal_lines, dxi_dx, levels)
plot(x, curvilinear.phase(phase).wave_profile, 'color', 'red', 'linewidth', 3)
for i = 1:10:150
    plot(x, horizontal_lines(i,:), 'color', 'black', 'linewidth', 2)
end
hold off
axis equal
colorbar
title("$d\xi / dx$", 'Interpreter', 'latex')

h(2) = nexttile;
hold on
contourf(vertical_lines, horizontal_lines, dxi_dy, levels)
plot(X(1,:), curvilinear.phase(phase).wave_profile, 'color', 'red', 'linewidth', 3)
for i = 1:10:150
    plot(x, horizontal_lines(i,:), 'color', 'black', 'linewidth', 2)
end
hold off
axis equal
colorbar
title("$d\xi / dy$", 'Interpreter', 'latex')

h(3) = nexttile;
hold on
contourf(vertical_lines, horizontal_lines, dzeta_dx, levels)
plot(X(1,:), curvilinear.phase(phase).wave_profile, 'color', 'red', 'linewidth', 3)
for i = 1:10:171
    plot(vertical_lines(:,i), y + curvilinear.phase(phase).wave_profile(i), 'color', 'black', 'linewidth', 2)
end
hold off
axis equal
colorbar
title("$d\zeta / dx$", 'Interpreter', 'latex')

h(4) = nexttile;
hold on
contourf(vertical_lines, horizontal_lines, dzeta_dy, levels)
plot(X(1,:), curvilinear.phase(phase).wave_profile, 'color', 'red', 'linewidth', 3)
for i = 1:10:171
    plot(vertical_lines(:,i), y + curvilinear.phase(phase).wave_profile(i), 'color', 'black', 'linewidth', 2)
end
hold off
axis equal
colorbar
title("$d\zeta / dy$", 'Interpreter', 'latex')

linkaxes(h, 'xy')
ylim([-30, 200])

%%

buffer = 1;
u(vertical_lines < 1) = nan;
u(vertical_lines > 235) = nan;
u(horizontal_lines > 201) = nan;
u(horizontal_lines < curvilinear.phase(phase).max_wave_profile + buffer) = nan;


figure()
hold on
contourf(vertical_lines, horizontal_lines, u, 20)
plot(x, curvilinear.phase(phase).max_wave_profile, 'color', 'red')
plot(x, curvilinear.phase(phase).max_wave_profile + buffer, 'color', 'blue')
plot(x, curvilinear.phase(phase).wave_profile, 'color', 'black')
% plot(x, curvilinear.phase(phase).wave_profile + 11, 'color', 'black', 'linestyle', '--')
hold off
axis equal
colorbar

%%

[du_dxi, du_dzeta] = gradient(u, 1, 1);

% Convert dxi/dzeta into meters
du_dxi = du_dxi ./ (dxi_dx * 1E-3);
du_dzeta = du_dzeta ./ (dzeta_dy * 1E-3);
buffer = 3;

du_dxi(horizontal_lines < curvilinear.phase(phase).max_wave_profile + buffer) = nan;
du_dzeta(horizontal_lines < curvilinear.phase(phase).max_wave_profile + buffer) = nan;


figure()
tiledlayout(1,2)

nexttile
hold on
contourf(vertical_lines, horizontal_lines, du_dxi, 50, 'linestyle', 'none')
plot(x, curvilinear.phase(phase).wave_profile)
plot(x, curvilinear.phase(phase).max_wave_profile + buffer)
hold off
axis equal
colorbar
title("$du / d\xi$", "interpreter", "latex")

nexttile
hold on
contourf(vertical_lines, horizontal_lines, du_dzeta, 50, 'linestyle', 'none')
plot(x, curvilinear.phase(phase).wave_profile)
plot(x, curvilinear.phase(phase).max_wave_profile + buffer)
hold off
axis equal
colorbar
title("$du / d\zeta$", "interpreter", "latex")


%% Integration

% From ChatGPT
u_zeta_cumulative = nan(size(u));
for i = 1:size(u, 2)
    zeta_column = horizontal_lines(:,i) * 1E-3;
    data_column = u(:,i);
    nan_mask = ~isnan(zeta_column) & ~isnan(data_column);
    
    if sum(nan_mask) > 1
        u_zeta_cumulative(nan_mask, i) = flipud(cumtrapz(flipud(zeta_column(nan_mask)), flipud(data_column(nan_mask))));
    end
end

figure()
nexttile
hold on
contourf(vertical_lines, horizontal_lines, u_zeta_cumulative, 50, 'linestyle', 'none')
plot(x, curvilinear.phase(phase).wave_profile)
plot(x, curvilinear.phase(phase).max_wave_profile + buffer)
hold off
axis equal
colorbar

%% Made into a function

test = curvilinearIntegrate(u, horizontal_lines);

figure()
nexttile
hold on
contourf(vertical_lines, horizontal_lines, test, 50, 'linestyle', 'none')
plot(x, curvilinear.phase(phase).wave_profile)
plot(x, curvilinear.phase(phase).max_wave_profile + buffer)
hold off
axis equal
colorbar


%% Make it a function

function output = curvilinearIntegrate(quantity, horizontal_lines)
    output = nan(size(quantity));
    for i = 1:size(quantity, 2)
        zeta_column = horizontal_lines(:,i) * 1E-3;
        data_column = quantity(:,i);
        nan_mask = ~isnan(zeta_column) & ~isnan(data_column);
        
        if sum(nan_mask) > 1
            output(nan_mask, i) = flipud(cumtrapz(flipud(zeta_column(nan_mask)), flipud(data_column(nan_mask))));
        end
    end
end
