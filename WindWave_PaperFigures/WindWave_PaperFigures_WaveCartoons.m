% Cartoon figures for hysteresis figure

clc; clear; close all; 
save_folder = '/Users/zeinsadek/Desktop/WindWaveCartoons';

%% Generate the top surface wave image

% Wave details
x = 0:0.01:5;
phase_shift = 0.5;
lambda = 1;
amplitude = 0.5;
wave = amplitude * cos(2*pi*(1/lambda)*(x - phase_shift)) + (amplitude);


% Plot details
linewidth = 1.5;
wave_surface_color = 'black';
below_surface_alpha = 0.25;

% Figure
figure('color', 'white', 'units', 'centimeters', 'Position', [4,4,10,2]);
tiledlayout(1,1,'padding','tight')
ax = nexttile;
hold on
plot(x, wave, 'color', wave_surface_color, 'linewidth', linewidth)

% Shaded region
hFill = patch( ...
[x, fliplr(x)], ...
[wave, -5 * ones(size(wave))], ...
wave_surface_color, ...
'FaceAlpha', below_surface_alpha, ...       
'EdgeColor', 'none', ...     
'HandleVisibility', 'off'); 
uistack(hFill, 'bottom')

hold off
xlim([0, max(x)])
ylim([-1,  5])
set(gca,'XColor', 'none','YColor','none')


% Save figure
pause(3)
exportgraphics(ax, fullfile(save_folder, 'WindWave_LongWaveCartoon.pdf'), 'Resolution', 600)
close all


%% Generate phase pictures

% Wave details
x = 0:0.01:5;
phase_shift = 1;
lambda = 1;
amplitude = 1;
wave = amplitude * cos(2*pi*(1/lambda)*(x - phase_shift)) + (amplitude);


% Plot details
linewidth = 1.5;
wave_surface_color = 'black';
below_surface_alpha = 0.25;

% Figure
figure('color', 'white', 'units', 'centimeters', 'Position', [4,4,3,3]);
tiledlayout(1,1,'padding','tight')
ax = nexttile;
hold on
plot(x, wave, 'color', wave_surface_color, 'linewidth', linewidth)

% Shaded region
hFill = patch( ...
[x, fliplr(x)], ...
[wave, -5 * ones(size(wave))], ...
wave_surface_color, ...
'FaceAlpha', below_surface_alpha, ...       
'EdgeColor', 'none', ...     
'HandleVisibility', 'off'); 
uistack(hFill, 'bottom')

hold off
axis square
xlim([0, lambda/2])
ylim([-1,  5])
set(gca,'XColor', 'none','YColor','none')
% set(gca, 'box', 'on')


% Save figure
pause(3)
exportgraphics(ax, fullfile(save_folder, 'WindWave_SmallWaveCartoon.pdf'), 'Resolution', 600)
close all




