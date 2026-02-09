%%% Checking if we should crop near the surface for ensemble averages with
%%% waves

clear; close all; clc;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions/');
mat_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";

caze = "WT4_WVA_AG0";

%% cropped data

cropped_data = load(fullfile(mat_path, 'means', strcat(caze, "_MEANS.mat")));
cropped_data = cropped_data.output;

%% uncropped data

uncropped_data = load(fullfile(mat_path, 'data', strcat(caze, "_DATA.mat")));
uncropped_data = uncropped_data.output;

%% compute uncropped mean

X = uncropped_data.X;
Y = uncropped_data.Y;

x = X(1,:);
y = Y(:,1);

uncropped_U = mean(uncropped_data.U, 3, 'omitnan');
uncropped_V = mean(uncropped_data.V, 3, 'omitnan');

uncropped_U(Y < 0) = nan;
uncropped_V(Y < 0) = nan;

%% compare contours and profiles

figure()
tiledlayout(1,2)

nexttile
contourf(X, Y, uncropped_U, 500, 'LineStyle', 'none')
axis equal
colorbar()
clim([0, 3])
title('Uncropped')

nexttile
contourf(X, Y, cropped_data.ensemble.u, 500, 'LineStyle', 'none')
axis equal
colorbar()
clim([0, 3])
title('Uncropped')


idx = 86;
lw = 3;
figure()
hold on
plot(uncropped_U(:, idx), y, 'LineWidth', 3)
plot(cropped_data.ensemble.u(:, idx), y, 'LineWidth', 3)
hold off

%% compute uncropped stresses

uncropped_U_flucs = uncropped_data.U - mean(uncropped_data.U, 3, 'omitnan');
uncropped_V_flucs = uncropped_data.V - mean(uncropped_data.V, 3, 'omitnan');

uncropped_UU = mean(uncropped_U_flucs.^2, 3, 'omitnan');
uncropped_VV = mean(uncropped_V_flucs.^2, 3, 'omitnan');
uncropped_UV = mean(uncropped_U_flucs .* uncropped_V_flucs, 3, 'omitnan');

%% compare stresses

figure()
tiledlayout(1,2)

nexttile
contourf(X, Y, uncropped_UV, 500, 'LineStyle', 'none')
axis equal
colorbar()
clim([-0.03, 0])
title('Uncropped')

nexttile
contourf(X, Y, cropped_data.ensemble.uv, 500, 'LineStyle', 'none')
axis equal
colorbar()
clim([-0.03, 0])
title('Uncropped')


idx = 86;
lw = 3;
figure()
hold on
plot(uncropped_UV(:, idx), y, 'LineWidth', 3)
plot(cropped_data.ensemble.uv(:, idx), y, 'LineWidth', 3)
hold off

%% plot phase average profiles

figure()
hold on
for i = 1:4
    plot(cropped_data.phase(i).u(:, 86), y)
end
hold off