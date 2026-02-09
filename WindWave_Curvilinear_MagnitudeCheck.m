%%% Checking agnitude between curvilinear and cartesian coordinates

clc; close all; clear
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')

project_path = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV';
caze = 'WT4_WVC_AG0';

ensemble = load(fullfile(project_path, 'means', strcat(caze, '_MEANS.mat')));
ensemble = ensemble.output;

curvilinear = load(fullfile(project_path, 'curvilinear', strcat(caze, '_CURVILINEAR.mat')));
curvilinear = curvilinear.output;

clear project_path

%% Plotting magnitude

levels = 30;
phase = 4;

% Compute magnitudes
ensemble_magnitude = sqrt(ensemble.phase(phase).u.^2 + ensemble.phase(phase).v.^2);
curvilinear_magnitude = sqrt(curvilinear.phase(phase).u.^2 + curvilinear.phase(phase).v.^2);

% Plot magnitudes
close all;
figure('color', 'white')
tiledlayout(1,2)
sgtitle(sprintf('%s: Phase %1.0f', caze, phase), 'interpreter', 'none')

nexttile
contourf(ensemble.X, ensemble.Y, ensemble_magnitude, levels, 'LineStyle', 'none')
axis equal
xlim([-20, 250])
ylim([-20, 200])
colorbar()
title('Cartesian Magnitude')

nexttile
contourf(curvilinear.phase(phase).vertical_lines, curvilinear.phase(phase).horizontal_lines, curvilinear_magnitude, levels, 'LineStyle', 'none')
axis equal
xlim([-20, 250])
ylim([-20, 200])
colorbar()
title('Curvilinear Magnitude')

% Match color ranges
ax = findall(gcf, 'Type', 'axes');
clims = cell2mat(get(ax, 'CLim'));
global_clim = [min(clims(:,1)), max(clims(:,2))];
set(ax, 'CLim', global_clim);


% Plot difference between magnitudes
% figure('color', 'white')
% contourf(curvilinear.phase(phase).vertical_lines, curvilinear.phase(phase).horizontal_lines, ensemble_magnitude - curvilinear_magnitude, levels, 'LineStyle', 'none')
% axis equal
% xlim([0, 200])
% ylim([0, 200])
% colorbar()
% colormap('coolwarm')
% title('Ensemble Magnitude - Curvilinear Magnitude')
% clim([-1, 1])

clc;
fprintf('Average difference between Ensemble and Curvilinear: %1.4f [m/s]\n\n', mean(ensemble_magnitude - curvilinear_magnitude, 'all', 'omitnan'))

clear levels phase ensemble_magnitude curvilinear_magnitude ax clims global_clim



%% Compare vertical velocity

levels = 50;
phase = 1;

% Plot magnitudes
close all;
figure('color', 'white')
tiledlayout(1,2)
sgtitle(sprintf('%s: Phase %1.0f', caze, phase), 'interpreter', 'none')

nexttile
contourf(ensemble.X, ensemble.Y, ensemble.phase(phase).v, levels, 'LineStyle', 'none')
axis equal
xlim([-20, 250])
ylim([-20, 200])
colorbar()
colormap('coolwarm')
title('Cartesian Vertical Velocity')

nexttile
contourf(curvilinear.phase(phase).vertical_lines, curvilinear.phase(phase).horizontal_lines, -curvilinear.phase(phase).v, levels, 'LineStyle', 'none')
axis equal
xlim([-20, 250])
ylim([-20, 200])
colorbar()
colormap('coolwarm')
title('Curvilinear Vertical Velocity')

% Match color ranges
ax = findall(gcf, 'Type', 'axes');
clims = cell2mat(get(ax, 'CLim'));
global_clim = [min(clims(:,1)), max(clims(:,2))];
set(ax, 'CLim', global_clim);

clear levels phase ax clims global_clim