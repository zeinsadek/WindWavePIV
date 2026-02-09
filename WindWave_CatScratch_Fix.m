%% Data was taken with a bad camera, trying to fix the cat scratch in the middle of the frame
% Zein Sadek

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions/');
wave_parameters = readcell('Offshore_Waves.xlsx');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data
experiment_name = 'WT4_WVA_AG0';

% Experiment Specifics
tunnel_freq = experiment_name(strfind(experiment_name, 'WT') + 2);
wave_type   = experiment_name(strfind(experiment_name, 'WV') + 2);

% Save paths
results_path  = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/';
mtlb_file     = strcat(results_path, 'data'   , '/', experiment_name, '_DATA.mat');
mean_file     = strcat(results_path, 'means'  , '/', experiment_name, '_MEANS.mat');
phase_file    = strcat(results_path, 'phase'  , '/', experiment_name, '_PHASE.mat');
BL_file       = strcat(results_path, 'boundary_layer/', experiment_name, '_BL.mat');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

means = load(mean_file); 
means = means.output;

% Coordinates
Y = means.Y;
y = Y(:,1);
X = means.X;
x = X(1,:);

dy = ((y(3) - y(2)) * 10^(-3)); % [m]
dx = ((x(3) - x(2)) * 10^(-3)); % [m]

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRY DIFFERENT FILTERING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phase = 3;
component = 'u';
u = means.phase(phase).(component);
max_wave_profile = means.phase(phase).max_wave_profile;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARTIFICIAL VERTICAL SHIFT: SINGLE X-PROFILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

row = 10;
slice = u(row, :);
LHS_value = slice(86);
RHS_value = slice(87);
difference = RHS_value - LHS_value;

slice_corrected = slice;
slice_corrected(87:end) = slice_corrected(87:end) - difference;

figure()
hold on
plot(x, slice, 'color', 'black')
scatter(x, slice, 10, 'filled', 'MarkerFaceColor', 'black')

plot(x, slice_corrected, 'color', 'red')
scatter(x, slice_corrected, 10, 'filled', 'MarkerFaceColor', 'red')

xline(x(86))
xline(x(87))
hold off

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARTIFICIAL VERTICAL SHIFT: ARRAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LHS_values = u(:,86);
RHS_values = u(:,87);
differences = RHS_values - LHS_values;
differences(isnan(differences)) = 0;

u_corrected = u;
u_corrected(:,87:end) = u_corrected(:,87:end) - differences;

% Contours
close all
figure()
tiledlayout(1,2)
nexttile
hold on
contourf(X, Y, u, 10, 'LineStyle', 'none')
plot(x, max_wave_profile)
hold off
xline(x(86))
xline(x(87))
axis equal
title('Original data')

nexttile
hold on
contourf(X, Y, u_corrected, 10, 'LineStyle', 'none')
plot(x, max_wave_profile)
hold off
xline(x(86))
xline(x(87))
axis equal
title('Shifted data')

% Profiles
figure()
hold on
for i = 1:20:length(u)
    plot(x, u(i,:) + 0.01 * i, 'color', 'black')
    plot(x, u_corrected(i,:) + 0.01 * i, 'color', 'red')
end
hold off
xline(x(86))
xline(x(87))
title('Profiles of shifted data')
xlim([0, 235])

% Gradients
[dudx, ~] = gradient(u, dx);
[du_shiftdx, ~] = gradient(u_corrected, dx);

figure()
tiledlayout(1,2)
nexttile 
contourf(X, Y, dudx, 100, 'linestyle', 'none')
axis equal
title('Orignal data')

nexttile 
contourf(X, Y, du_shiftdx, 100, 'linestyle', 'none')
axis equal
title('Shifted data')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LINEAR INTERPOLATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Where to interpolate between
[~, left] = min(abs(x - 115));
[~, right] = min(abs(x - 120));

% Interpolate across strip
u_lin_interp = u_corrected;
u_lin_interp(:, left:right) = nan;
u_lin_interp = fillmissing(u_lin_interp, 'makima', 2);
u_lin_interp(Y < max_wave_profile) = nan;

% Contours
figure()
tiledlayout(1,2)
nexttile
contourf(X, Y, u, 10, 'linestyle', 'none')
axis equal
title('Original data')

nexttile
contourf(X, Y, u_lin_interp, 10, 'linestyle', 'none')
axis equal
title('Shifted + Interpolated data')

% Profiles
figure()
hold on
for i = 1:10:length(u)
    plot(x, u(i,:) + 0.01 * i, 'color', 'black')
    plot(x, u_lin_interp(i,:) + 0.01 * i, 'color', 'red')
end
hold off
xline(x(left))
xline(x(right))
xlim([0, 235])

% Gradients
[du_shift_interp_dx, ~] = gradient(u_lin_interp, dx);

figure()
tiledlayout(1,2)
nexttile 
contourf(X, Y, dudx, 50, 'linestyle', 'none')
axis equal
title('Original data')

nexttile 
contourf(X, Y, du_shift_interp_dx, 50, 'linestyle', 'none')
axis equal
title('Shifted + Interpolated data')





