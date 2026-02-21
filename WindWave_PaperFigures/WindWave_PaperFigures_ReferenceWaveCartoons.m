%% Cartoon wave at different phases for figures

clc; close all; clear

x = -0.25:0.01:0.25;

% Phase 1: 0
% Phase 2: 0.25
% Phase 3: 0.5
% Phase 4: 0.75

phase = 0.25;
wave = cos(2 * pi * (x - phase)) + 2;

fig = figure('color', 'white');
hold on
area(x, wave)
plot(x, wave, 'linewidth', 10, 'color', 'black')
hold off
ylim([-1.5 + 2, 1.5 + 3])
xlim([min(x), max(x)])

ax = gca;
ax.XAxis.TickLength = [0 0];
ax.YAxis.TickLength = [0 0];
set(gca,'XTick',[], 'YTick', [])
box off
axis off 


% figure_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/WaveReferences';
% figure_name = strcat('WaveReference_Phase_2.png');
% 
% fprintf('Saving Figure\n\n')
% exportgraphics(fig, fullfile(figure_folder, figure_name), 'Resolution', 200)
% close all; clc; fprintf('Figure Saved!\n\n')

