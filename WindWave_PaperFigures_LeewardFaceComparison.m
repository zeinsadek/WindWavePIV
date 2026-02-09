%%% Wind + Wave paper curvilinear coordinates figure for demonstration

% WindWave paper figures: Curvilinear Phase Averages

clc; clear; close all
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions/Inpaint_nans/Inpaint_nans')
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps/slanCM')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
figure_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/pdf_test6';
curvilinear_path = fullfile(project_path, 'curvilinear_new');
cartesian_path = fullfile(project_path, 'means');

% Cases
wind_speed = 'WT6';
wave = 'C';
caze = strcat(wind_speed, '_WV', wave, '_AG0');

if ismember(wind_speed(end), {'4'})
    u_inf = 2.4181;
elseif ismember(wind_speed(end), {'6'})
    u_inf = 3.8709;
elseif ismember(wind_speed(end), {'8'})
    u_inf = 5.4289;
end

% Load data
% instantaneous = load(fullfile(data_path, strcat(caze, '_DATA.mat')));
cartesian = load(fullfile(cartesian_path, strcat(caze, '_MEANS.mat')));
cartesian = cartesian.output;
curvilinear = load(fullfile(curvilinear_path, strcat(caze, '_CURVILINEAR.mat')));
curvilinear = curvilinear.output;


clc; fprintf('All %s cases loaded\n', wind_speed)
clear caze tmp w no_wave_caze cartesian_path curvilinear_path data_path



%% Make a sexy ass plot

% Plot settings
phase = 4;
component = 'u';

% Cartesian Coordinates
X = cartesian.X;
x = X(1,:);
Y = cartesian.Y;
y = Y(:,1);

% Curvilinear Coordinates
vertical_lines = curvilinear.phase(phase).vertical_lines;
horizontal_lines = curvilinear.phase(phase).horizontal_lines;

% Load data
max_wave = cartesian.phase(phase).max_wave_profile;
reference_wave = cartesian.phase(phase).reference_wave;

% Plot cartesian grid on top, but crop out waves
X_grid = X;
Y_grid = Y;
[rows, columns] = size(X);


% Plot phase avg between cartesian and curvilinear side-by-side
gridLineWidth = 1;
linewidth = 3;
gridStep = 10;
levels = 10;
line_trans = 0.1;

% Quiver options
qstep = 8;
qlinewidth = 1;

% Near-wave cropping
vetical_crop_shift = 2;

% Y and color limits
ymin = -20;
ymax = 100;
cmin = 0.4;
cmax = 1;


% Plotting
figure('color', 'white')
tiledlayout(2,2)

% Cartesian
h(1) = nexttile(1);
title('Cartesian')
hold on

% Load velocities
cartesian_u = cartesian.phase(phase).u;
cartesian_v = cartesian.phase(phase).v;
cartesian_u_norm = cartesian_u / u_inf;
cartesian_v_norm = cartesian_v / u_inf;

% Manually crop to remove max_wave profile noise
cartesian_u_norm(Y <= reference_wave.' + vetical_crop_shift) = nan;
cartesian_v_norm(Y <= reference_wave.' + vetical_crop_shift) = nan;

% Plot contour
contourf(X, Y, cartesian_u_norm, ...
         levels, 'linestyle', 'none')
plot(x, reference_wave, 'color', 'black', 'linewidth', linewidth)

% Quivers
quiver(X(1:qstep:end, 1:qstep:end), Y(1:qstep:end, 1:qstep:end), ...
       cartesian_u_norm(1:qstep:end, 1:qstep:end), ...
       cartesian_v_norm(1:qstep:end, 1:qstep:end), ...
       'color', 'black', 'linewidth', qlinewidth)


% Vertical lines
X_grid(Y <= reference_wave.' + vetical_crop_shift) = nan;
Y_grid(Y <= reference_wave.' + vetical_crop_shift) = nan;
for i = 1:gridStep:rows
    P = plot(X_grid(:,i), y, 'color', 'black', 'linewidth', gridLineWidth);
    P.Color(4) = line_trans;
end


% Horizontal lines
for i = 1:gridStep:rows
    P = plot(x, Y_grid(i,:), 'color', 'black', 'linewidth', gridLineWidth);
    P.Color(4) = line_trans;
end

% Sample line
idx = round(length(x)/2);
plot(X_grid(:,idx), y, 'color', 'black', 'linewidth', 2);

hold off
axis equal

colorbar()
colormap((slanCM('Reds')))
clim([cmin, cmax])



% Curvilinear
h(2) = nexttile(3);
title('Curvilinear')
hold on

% Load velocities
curvilinear_u = curvilinear.phase(phase).u;
curvilinear_v = curvilinear.phase(phase).v;
curvilinear_u_norm = curvilinear_u / u_inf;
curvilinear_v_norm = curvilinear_v / u_inf;

% Manually crop to remove max_wave profile noise
curvilinear_u_norm(Y <= vetical_crop_shift) = nan;
curvilinear_v_norm(Y <= vetical_crop_shift) = nan;

% Plot contour
contourf(vertical_lines, horizontal_lines, curvilinear_u_norm, ...
         levels, 'linestyle', 'none')
plot(x, reference_wave, 'color', 'black', 'linewidth', linewidth)

% Quivers
quiver(vertical_lines(1:qstep:end, 1:qstep:end), horizontal_lines(1:qstep:end, 1:qstep:end), ...
       curvilinear_u_norm(1:qstep:end, 1:qstep:end), ...
       curvilinear_v_norm(1:qstep:end, 1:qstep:end), ...
       'color', 'black', 'linewidth', qlinewidth)

% Sample line
idx = round(length(x)/2);
plot(vertical_lines(:,idx), y + 1.05 * reference_wave(idx), 'color', 'black', 'linewidth', 2);

% Vertical lines
for i = 1:gridStep:rows
    P = plot(vertical_lines(:,i), y + 1.05 * reference_wave(i), 'color', 'black', 'linewidth', gridLineWidth);
    P.Color(4) = line_trans;
end
% Horizontal lines
for i = 1:gridStep:rows
    P = plot(x, horizontal_lines(i,:), 'color', 'black', 'linewidth', gridLineWidth);
    P.Color(4) = line_trans;
end
hold off
axis equal
colorbar()
clim([cmin, cmax])


linkaxes(h, 'xy')
ylim([ymin, ymax])



% Compare velocity profiles
h(3) = nexttile([2,1]);
idx = round(length(x)/2);
hold on
plot(cartesian_u_norm(:, idx), y, 'color', 'black', 'linewidth', 2, 'displayname', 'Cartesian')
plot(curvilinear_u_norm(:, idx), horizontal_lines(:, idx), 'color', 'red', 'linewidth', 2, 'displayname', 'Curvilinear')
hold off

legend







%% Compare uv stresses


wave_transparency = 0.25;


% Fontsizes
tickFontSize = 8;
labelFontSize = 10;
colorbar_fontsize = 10;
titleFontSize = 12;

% Plot phase avg between cartesian and curvilinear side-by-side
gridLineWidth = 1;
linewidth = 2;
gridStep = 8;
line_trans = 0.2;
levels = 20;



% Near-wave cropping
vetical_crop_shift = 0;

% Y and color limits
ymin = -20;
ymax = 120;
cmax = 7E-3;
cmin = 0;


% Plotting
clc; close all
figure('color', 'white', 'units', 'centimeters', 'Position', [10,10,13,5])
t = tiledlayout(1, 2, 'padding', 'compact', 'TileSpacing', 'tight');

%%% Cartesian
h(1) = nexttile(1);
set(h(1), 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
title("-$\langle u' v' \rangle \mathbin{/} u_{\infty}^2$", 'interpreter', 'latex', 'fontsize', titleFontSize)
hold on

% Load velocities
cartesian_uv = cartesian.phase(phase).uv;
cartesian_uv_norm = cartesian_uv / (u_inf^2);

% Manually crop to remove max_wave profile noise
cartesian_uv_norm(Y <= reference_wave.' + vetical_crop_shift) = nan;
cartesian_uv_norm(Y <= reference_wave.' + vetical_crop_shift) = nan;

% Plot contour
contourf(X, Y, -cartesian_uv_norm, ...
         levels, 'linestyle', 'none')

% Plot wave
plot(x, reference_wave, 'color', 'black', 'linewidth', linewidth)

% Shaded region below wave
hFill = patch( ...
[x, fliplr(x)], ...
[reference_wave.', -100 * ones(size(reference_wave.'))], ...
'k', ...
'FaceAlpha', wave_transparency, ...      
'EdgeColor', 'none', ...    
'HandleVisibility', 'off'); 
uistack(hFill, 'bottom')

% Vertical lines
X_grid(Y <= reference_wave.' + vetical_crop_shift) = nan;
Y_grid(Y <= reference_wave.' + vetical_crop_shift) = nan;
for i = 1:gridStep:rows
    P = plot(X_grid(:,i), y, 'color', 'black', 'linewidth', gridLineWidth);
    P.Color(4) = line_trans;
end

% Horizontal lines
for i = 1:gridStep:rows
    P = plot(x, Y_grid(i,:), 'color', 'black', 'linewidth', gridLineWidth);
    P.Color(4) = line_trans;
end
hold off
xticks(0:100:200)
yticks(0:50:300)
axis equal

colormap((slanCM('RdPu')))
clim([cmin, cmax])




%%% Curvilinear
h(2) = nexttile(2);
set(h(2), 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
title("-$\langle u_{\xi}' u_{\zeta}' \rangle \mathbin{/} u_{\infty}^2$", 'interpreter', 'latex', 'fontsize', titleFontSize)
hold on

% Load velocities
curvilinear_uv = curvilinear.phase(phase).uv;
curvilinear_uv_norm = curvilinear_uv / (u_inf^2);

% Manually crop to remove max_wave profile noise
curvilinear_uv_norm(Y <= vetical_crop_shift) = nan;

% Plot contour
contourf(vertical_lines, horizontal_lines, -curvilinear_uv_norm, ...
         levels, 'linestyle', 'none')

% Plot wave
plot(x, reference_wave, 'color', 'black', 'linewidth', linewidth)

% Shaded region below wave
hFill = patch( ...
[x, fliplr(x)], ...
[reference_wave.', -100 * ones(size(reference_wave.'))], ...
'k', ...
'FaceAlpha', wave_transparency, ...      
'EdgeColor', 'none', ...    
'HandleVisibility', 'off'); 
uistack(hFill, 'bottom')

% Vertical lines
for i = 1:gridStep:rows
    P = plot(vertical_lines(:,i), y + 1.05 * reference_wave(i), 'color', 'black', 'linewidth', gridLineWidth);
    P.Color(4) = line_trans;
end
% Horizontal lines
for i = 1:gridStep:rows
    P = plot(x, horizontal_lines(i,:), 'color', 'black', 'linewidth', gridLineWidth);
    P.Color(4) = line_trans;
end
hold off
xticks(0:100:200)
yticks(0:50:300)
axis equal

colorbar()
colormap((slanCM('RdPu')))
clim([cmin, cmax])

C = colorbar;
C.Label.Interpreter = 'latex';
C.Label.FontSize = colorbar_fontsize;
C.TickLabelInterpreter = 'latex';
C.Ruler.Exponent = -3;
C.Ruler.TickLabelFormat = '%2.0f';

% Turn off the y-axis
ax = gca;
ax.YAxis.Visible = 'off';


linkaxes(h, 'xy')
ylim([ymin, ymax])
xlim([0, 230])
xlabel(t, '$x$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(t, '$y$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)



% Save figure
% pause(3)
% figure_name = sprintf('LogLaw_%s_WV%s_Phase_%1.0f_uv_Comparison.pdf', wind_speed, wave, phase);
% exportgraphics(t, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 600, 'ContentType', 'image');
% close all
% fprintf('Generated figure: %s\n\n', figure_name)






%% Compare uu stresses


% wave_transparency = 0.25;
% 
% 
% % Fontsizes
% tickFontSize = 14;
% labelFontSize = 16;
% colorbar_fontsize = 16;
% titleFontSize = 16;
% 
% % Plot phase avg between cartesian and curvilinear side-by-side
% gridLineWidth = 1;
% linewidth = 3;
% gridStep = 8;
% line_trans = 0.2;
% levels = 20;
% 
% % Near-wave cropping
% vetical_crop_shift = 0;
% 
% % Y and color limits
% ymin = -20;
% ymax = 120;
% cmax = 35E-3;
% cmin = 0;
% 
% % Plotting
% clc; close all
% figure('color', 'white', 'units', 'centimeters', 'Position', [2,4,24,8])
% t = tiledlayout(1, 2, 'padding', 'tight', 'TileSpacing', 'tight');
% 
% %%% Cartesian
% h(1) = nexttile(1);
% set(h(1), 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% title("-$\langle u' u' \rangle / u_{\infty}^2$", 'interpreter', 'latex', 'fontsize', titleFontSize)
% hold on
% 
% % Load velocities
% cartesian_uu = cartesian.phase(phase).uu;
% cartesian_uu_norm = cartesian_uu / (u_inf^2);
% 
% % Manually crop to remove max_wave profile noise
% cartesian_uu_norm(Y <= reference_wave.' + vetical_crop_shift) = nan;
% cartesian_uu_norm(Y <= reference_wave.' + vetical_crop_shift) = nan;
% 
% % Plot contour
% contourf(X, Y, cartesian_uu_norm, ...
%          levels, 'linestyle', 'none')
% 
% % Plot wave
% plot(x, reference_wave, 'color', 'black', 'linewidth', linewidth)
% 
% % Shaded region below wave
% hFill = patch( ...
% [x, fliplr(x)], ...
% [reference_wave.', -100 * ones(size(reference_wave.'))], ...
% 'k', ...
% 'FaceAlpha', wave_transparency, ...      
% 'EdgeColor', 'none', ...    
% 'HandleVisibility', 'off'); 
% uistack(hFill, 'bottom')
% 
% % Vertical lines
% X_grid(Y <= reference_wave.' + vetical_crop_shift) = nan;
% Y_grid(Y <= reference_wave.' + vetical_crop_shift) = nan;
% for i = 1:gridStep:rows
%     P = plot(X_grid(:,i), y, 'color', 'black', 'linewidth', gridLineWidth);
%     P.Color(4) = line_trans;
% end
% 
% % Horizontal lines
% for i = 1:gridStep:rows
%     P = plot(x, Y_grid(i,:), 'color', 'black', 'linewidth', gridLineWidth);
%     P.Color(4) = line_trans;
% end
% hold off
% xticks(0:100:200)
% yticks(0:50:300)
% axis equal
% 
% colormap((slanCM('Reds')))
% clim([cmin, cmax])
% 
% 
% 
% 
% %%% Curvilinear
% h(2) = nexttile(2);
% set(h(2), 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% title("-$\langle u_{\xi}' u_{\xi}' \rangle / u_{\infty}^2$", 'interpreter', 'latex', 'fontsize', titleFontSize)
% hold on
% 
% % Load velocities
% curvilinear_uu = curvilinear.phase(phase).uu;
% curvilinear_uu_norm = curvilinear_uu / (u_inf^2);
% 
% % Manually crop to remove max_wave profile noise
% curvilinear_uu_norm(Y <= vetical_crop_shift) = nan;
% 
% % Plot contour
% contourf(vertical_lines, horizontal_lines, curvilinear_uu_norm, ...
%          levels, 'linestyle', 'none')
% 
% % Plot wave
% plot(x, reference_wave, 'color', 'black', 'linewidth', linewidth)
% 
% % Shaded region below wave
% hFill = patch( ...
% [x, fliplr(x)], ...
% [reference_wave.', -100 * ones(size(reference_wave.'))], ...
% 'k', ...
% 'FaceAlpha', wave_transparency, ...      
% 'EdgeColor', 'none', ...    
% 'HandleVisibility', 'off'); 
% uistack(hFill, 'bottom')
% 
% % Vertical lines
% for i = 1:gridStep:rows
%     P = plot(vertical_lines(:,i), y + 1.05 * reference_wave(i), 'color', 'black', 'linewidth', gridLineWidth);
%     P.Color(4) = line_trans;
% end
% % Horizontal lines
% for i = 1:gridStep:rows
%     P = plot(x, horizontal_lines(i,:), 'color', 'black', 'linewidth', gridLineWidth);
%     P.Color(4) = line_trans;
% end
% hold off
% xticks(0:100:200)
% yticks(0:50:300)
% axis equal
% 
% colorbar()
% colormap((slanCM('Reds')))
% clim([cmin, cmax])
% 
% C = colorbar;
% C.Label.Interpreter = 'latex';
% C.Label.FontSize = colorbar_fontsize;
% C.TickLabelInterpreter = 'latex';
% C.Ruler.Exponent = -3;
% C.Ruler.TickLabelFormat = '%2.0f';
% 
% % Turn off the y-axis
% ax = gca;
% ax.YAxis.Visible = 'off';
% 
% 
% linkaxes(h, 'xy')
% ylim([ymin, ymax])
% xlim([0, 230])
% xlabel(t, '$x$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylabel(t, '$y$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)
% 
% 
% 
% % Save figure
% % pause(3)
% % figure_name = sprintf('LogLaw_%s_WV%s_Phase_%1.0f_uu_Comparison.pdf', wind_speed, wave, phase);
% % exportgraphics(t, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 600, 'ContentType', 'image');
% % close all
% % fprintf('Generated figure: %s\n\n', figure_name)






%% Compare vv stresses

% wave_transparency = 0.25;
% 
% 
% % Fontsizes
% tickFontSize = 14;
% labelFontSize = 16;
% legendFontSize = 12;
% annotationFontSize = 10;
% colorbar_fontsize = 16;
% titleFontSize = 16;
% sz = 30;
% 
% % Plot phase avg between cartesian and curvilinear side-by-side
% gridLineWidth = 1;
% linewidth = 3;
% gridStep = 8;
% levels = 20;
% line_trans = 0.2;
% 
% % Near-wave cropping
% vetical_crop_shift = 0;
% 
% % Y and color limits
% ymin = -20;
% ymax = 120;
% cmax = 11E-3;
% cmin = 0;
% 
% 
% % Plotting
% clc; close all
% figure('color', 'white', 'units', 'centimeters', 'Position', [2,4,24,8])
% t = tiledlayout(1, 2, 'padding', 'tight', 'TileSpacing', 'tight');
% 
% %%% Cartesian
% h(1) = nexttile(1);
% set(h(1), 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% title("-$\langle v' v' \rangle / u_{\infty}^2$", 'interpreter', 'latex', 'fontsize', titleFontSize)
% hold on
% 
% % Load velocities
% cartesian_vv = cartesian.phase(phase).vv;
% cartesian_vv_norm = cartesian_vv / (u_inf^2);
% 
% % Manually crop to remove max_wave profile noise
% cartesian_vv_norm(Y <= reference_wave.' + vetical_crop_shift) = nan;
% cartesian_vv_norm(Y <= reference_wave.' + vetical_crop_shift) = nan;
% 
% % Plot contour
% contourf(X, Y, cartesian_vv_norm, ...
%          levels, 'linestyle', 'none')
% 
% % Plot wave
% plot(x, reference_wave, 'color', 'black', 'linewidth', linewidth)
% 
% % Shaded region below wave
% hFill = patch( ...
% [x, fliplr(x)], ...
% [reference_wave.', -100 * ones(size(reference_wave.'))], ...
% 'k', ...
% 'FaceAlpha', wave_transparency, ...      
% 'EdgeColor', 'none', ...    
% 'HandleVisibility', 'off'); 
% uistack(hFill, 'bottom')
% 
% % Vertical lines
% X_grid(Y <= reference_wave.' + vetical_crop_shift) = nan;
% Y_grid(Y <= reference_wave.' + vetical_crop_shift) = nan;
% for i = 1:gridStep:rows
%     P = plot(X_grid(:,i), y, 'color', 'black', 'linewidth', gridLineWidth);
%     P.Color(4) = line_trans;
% end
% 
% % Horizontal lines
% for i = 1:gridStep:rows
%     P = plot(x, Y_grid(i,:), 'color', 'black', 'linewidth', gridLineWidth);
%     P.Color(4) = line_trans;
% end
% hold off
% xticks(0:100:200)
% yticks(0:50:300)
% axis equal
% 
% colormap((slanCM('Reds')))
% clim([cmin, cmax])
% 
% 
% 
% 
% %%% Curvilinear
% h(2) = nexttile(2);
% set(h(2), 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% title("-$\langle u_{\zeta}' u_{\zeta}' \rangle / u_{\infty}^2$", 'interpreter', 'latex', 'fontsize', titleFontSize)
% hold on
% 
% % Load velocities
% curvilinear_vv = curvilinear.phase(phase).vv;
% curvilinear_vv_norm = curvilinear_vv / (u_inf^2);
% 
% % Manually crop to remove max_wave profile noise
% curvilinear_vv_norm(Y <= vetical_crop_shift) = nan;
% 
% % Plot contour
% contourf(vertical_lines, horizontal_lines, curvilinear_vv_norm, ...
%          levels, 'linestyle', 'none')
% 
% % Plot wave
% plot(x, reference_wave, 'color', 'black', 'linewidth', linewidth)
% 
% % Shaded region below wave
% hFill = patch( ...
% [x, fliplr(x)], ...
% [reference_wave.', -100 * ones(size(reference_wave.'))], ...
% 'k', ...
% 'FaceAlpha', wave_transparency, ...      
% 'EdgeColor', 'none', ...    
% 'HandleVisibility', 'off'); 
% uistack(hFill, 'bottom')
% 
% % Vertical lines
% for i = 1:gridStep:rows
%     P = plot(vertical_lines(:,i), y + 1.05 * reference_wave(i), 'color', 'black', 'linewidth', gridLineWidth);
%     P.Color(4) = line_trans;
% end
% % Horizontal lines
% for i = 1:gridStep:rows
%     P = plot(x, horizontal_lines(i,:), 'color', 'black', 'linewidth', gridLineWidth);
%     P.Color(4) = line_trans;
% end
% hold off
% xticks(0:100:200)
% yticks(0:50:300)
% axis equal
% 
% colorbar()
% colormap((slanCM('Reds')))
% clim([cmin, cmax])
% 
% C = colorbar;
% C.Label.Interpreter = 'latex';
% C.Label.FontSize = colorbar_fontsize;
% C.TickLabelInterpreter = 'latex';
% C.Ruler.Exponent = -3;
% C.Ruler.TickLabelFormat = '%2.0f';
% 
% % Turn off the y-axis
% ax = gca;
% ax.YAxis.Visible = 'off';
% 
% 
% linkaxes(h, 'xy')
% ylim([ymin, ymax])
% xlim([0, 230])
% xlabel(t, '$x$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylabel(t, '$y$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)
% 
% 
% 
% % Save figure
% % pause(3)
% % figure_name = sprintf('LogLaw_%s_WV%s_Phase_%1.0f_vv_Comparison.pdf', wind_speed, wave, phase);
% % exportgraphics(t, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 600, 'ContentType', 'image');
% % close all
% % fprintf('Generated figure: %s\n\n', figure_name)



%% Generate figure with 6 tiles comparing stresses between cartesian and curvilinear


% Below wave color
wave_transparency = 0.25;

% Fontsizes
tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;
annotationFontSize = 10;
colorbar_fontsize = 10;
titleFontSize = 12;
sz = 30;

% Grid lines
gridLineWidth = 1;
linewidth = 2;
gridStep = 8;
levels = 20;
line_trans = 0.2;

% Near-wave cropping
vetical_crop_shift = 0;



% Order to plot stresses
stresses = {'uv', 'uu', 'vv'};
cmins = [0, 0, 0];
cmaxs = [9E-3, 35E-3, 11E-3];

colors = {'RdPu', 'OrRd', 'YlGnBu'};

% Stress titles
cartesian_stress_labels.('uu') = "$\langle u' u' \rangle \mathbin{/} u_{\infty}^2$";
cartesian_stress_labels.('vv') = "$\langle v' v' \rangle \mathbin{/} u_{\infty}^2$";
cartesian_stress_labels.('uv') = "-$\langle u' v' \rangle \mathbin{/} u_{\infty}^2$";

curvilinear_stress_labels.('uu') = "$\langle u_{\xi}' u_{\xi}' \rangle \mathbin{/} u_{\infty}^2$";
curvilinear_stress_labels.('vv') = "$\langle u_{\zeta}' u_{\zeta}' \rangle \mathbin{/} u_{\infty}^2$";
curvilinear_stress_labels.('uv') = "-$\langle u_{\xi}' u_{\zeta}' \rangle \mathbin{/} u_{\infty}^2$";

% Plot
clc; close all
figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,13])
t = tiledlayout(3, 2, 'TileSpacing', 'tight', 'padding', 'loose');

c = 1;
 % Loop through stresses
for s = 1:length(stresses)
    stress = stresses{s};
    cmin = cmins(s);
    cmax = cmaxs(s);

    if strcmp(stress, 'uv')
        mult = -1;
    else
        mult = 1;
    end

    for i = 1:2
        h(c) = nexttile;
        set(h(c), 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
        hold on

        % Plot wave
        plot(x, reference_wave, 'color', 'black', 'linewidth', linewidth)
        
        % Shaded region below wave
        hFill = patch( ...
        [x, fliplr(x)], ...
        [reference_wave.', -100 * ones(size(reference_wave.'))], ...
        'k', ...
        'FaceAlpha', wave_transparency, ...      
        'EdgeColor', 'none', ...    
        'HandleVisibility', 'off'); 
        uistack(hFill, 'bottom')

        % Cartesian
        if i == 1
            title(cartesian_stress_labels.(stress), 'interpreter', 'latex', 'fontsize', titleFontSize)
            % Load velocities
            cartesian_vv_norm = mult * cartesian.phase(phase).(stress) / (u_inf^2);
            % Manually crop to remove max_wave profile noise
            cartesian_vv_norm(Y <= reference_wave.' + vetical_crop_shift) = nan;
            cartesian_vv_norm(Y <= reference_wave.' + vetical_crop_shift) = nan;
            % Plot contour
            contourf(X, Y, cartesian_vv_norm, ...
                     levels, 'linestyle', 'none')
            clim([cmin, cmax])

            % Vertical lines
            X_grid(Y <= reference_wave.' + vetical_crop_shift) = nan;
            Y_grid(Y <= reference_wave.' + vetical_crop_shift) = nan;
            for j = 1:gridStep:rows
                P = plot(X_grid(:,j), y, 'color', 'black', 'linewidth', gridLineWidth);
                P.Color(4) = line_trans;
            end
            
            % Horizontal lines
            for j = 1:gridStep:rows
                P = plot(x, Y_grid(j,:), 'color', 'black', 'linewidth', gridLineWidth);
                P.Color(4) = line_trans;
            end
            % colormap((slanCM(colors{s})))

        elseif i == 2
            title(curvilinear_stress_labels.(stress), 'interpreter', 'latex', 'fontsize', titleFontSize)
            % Load velocities
            curvilinear_tmp_norm = mult * curvilinear.phase(phase).(stress) / (u_inf^2);
            % Manually crop to remove max_wave profile noise
            curvilinear_tmp_norm(Y <= vetical_crop_shift) = nan;
            % Plot contour
            contourf(vertical_lines, horizontal_lines, curvilinear_tmp_norm, ...
                     levels, 'linestyle', 'none')
            clim([cmin, cmax])

            % Vertical lines
            for j = 1:gridStep:rows
                P = plot(vertical_lines(:,j), y + 1.05 * reference_wave(j), 'color', 'black', 'linewidth', gridLineWidth);
                P.Color(4) = line_trans;
            end
            % Horizontal lines
            for j = 1:gridStep:rows
                P = plot(x, horizontal_lines(j,:), 'color', 'black', 'linewidth', gridLineWidth);
                P.Color(4) = line_trans;
            end
            % Turn off the y-axis
            ax = gca;
            ax.YAxis.Visible = 'off';

            % Colorbar
            C = colorbar;
            C.Label.Interpreter = 'latex';
            C.Label.FontSize = colorbar_fontsize;
            C.TickLabelInterpreter = 'latex';
            C.Ruler.Exponent = -3;
            C.Ruler.TickLabelFormat = '%2.0f';
        end

        colormap(h(c), (slanCM(colors{s})))
        hold off
        axis equal

        % if c < 5
            % % Turn off the y-axis
            % ax = gca;
            % ax.XAxis.Visible = 'off';
        % end

        c = c + 1;
    end
end


linkaxes(h, 'xy')
ylim([ymin, ymax])
yticks(0:50:200)
xlim([0, max(x)])
xlabel(t, '$x$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(t, '$y$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)

% Add labels
addPanelLabels(h, {'a', 'b', 'c', 'd', 'e', 'f'}, 'FontSize', 10, 'Offset', [-0.08, 1.18])


% Save figure
% pause(3)
% figure_name = sprintf('LogLaw_%s_WV%s_Phase_%1.0f_AllStressComparisons.pdf', wind_speed, wave, phase);
% exportgraphics(t, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 600, 'ContentType', 'image');
% close all
% fprintf('Generated figure: %s\n\n', figure_name)











%% compare cartesian normal stress anisotropy




% Fontsizes
tickFontSize = 14;
labelFontSize = 16;
colorbar_fontsize = 16;
titleFontSize = 16;

% Plot phase avg between cartesian and curvilinear side-by-side
gridLineWidth = 1;
linewidth = 3;
gridStep = 8;
line_trans = 0.2;
levels = 20;

% Near-wave cropping
vetical_crop_shift = 0;

% Y and color limits
ymin = -20;
ymax = 120;
cmax = 35E-3;
cmin = 0;

% Plotting
clc; close all
figure('color', 'white', 'units', 'centimeters', 'Position', [2,4,24,8])
t = tiledlayout(1, 4, 'padding', 'tight', 'TileSpacing', 'tight');

for phase = 1:4
    %%% Cartesian
    h(phase) = nexttile;
    hold on
    
    % Load velocities
    cartesian_uu = cartesian.phase(phase).uu;
    cartesian_vv = cartesian.phase(phase).vv;
    
    % Plot contour
    contourf(X, Y, cartesian_uu - cartesian_vv, ...
             levels, 'linestyle', 'none')

    hold off

    xticks(0:100:200)
    yticks(0:50:300)
    axis equal
    
    colormap((slanCM('Reds')))
    colorbar()
end










%% Functions

function addPanelLabels(ax, labels, varargin)
% addPanelLabels(ax, labels) adds (a),(b),... just OUTSIDE top-left of each axes.
% ax     : array of axes handles (e.g., from tiledlayout / findall)
% labels : cellstr like {'a','b','c'} or string array ["a" "b" "c"]
%
% Optional name-value:
% 'Offset'   : [dx dy] in normalized axes units (default [-0.10 1.02])
% 'FontSize' : default 12
% 'FontName' : default 'Times New Roman'

p = inputParser;
addParameter(p,'Offset',[-0.10 1.1]);
addParameter(p,'FontSize', 10);
addParameter(p,'FontName','Times New Roman');
parse(p,varargin{:});
off = p.Results.Offset;

labels = string(labels);
for k = 1:numel(ax)
    if ~isgraphics(ax(k),'axes'), continue; end

    % Plain parentheses + italic letter:
    % TeX interpreter: \it turns italic ON, \rm returns to roman.
    s = sprintf('(\\ita\\rm)');              % placeholder
    s = sprintf('(\\it%s\\rm)', labels(k));  % actual label

    text(ax(k), off(1), off(2), s, ...
        'Units','normalized', ...
        'Interpreter','tex', ...           % keeps italics control simple
        'FontName',p.Results.FontName, ...
        'FontSize',p.Results.FontSize, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'Clipping','off');                 % critical: allow outside axes
    end
end


function hAnn = addPanelLabelsFixed(fig, ax, labels, varargin)
% Places (a),(b),... just outside top-left with FIXED physical offsets.
%
% fig    : figure handle (e.g., totalFigure)
% ax     : array of axes handles
% labels : {'a','b',...} or ["a","b",...]
%
% Name-value:
% 'OffsetPts' : [dx dy] in points (default [-10 6])  (left, up)
% 'FontSize'  : default 12
% 'FontName'  : default 'Times New Roman'

p = inputParser;
addParameter(p,'OffsetPts',[-10 6]);
addParameter(p,'FontSize',12);
addParameter(p,'FontName','Times New Roman');
parse(p,varargin{:});

labels = string(labels);
offPts = p.Results.OffsetPts;

ppi = get(0,'ScreenPixelsPerInch');
offPx = offPts/72 * ppi;   % points -> pixels

hAnn = gobjects(numel(ax),1);

for k = 1:numel(ax)
    if ~isgraphics(ax(k),'axes'), continue; end

    % Axes outer position in figure pixels
    op = getpixelposition(ax(k), true);  % [x y w h] in pixels relative to figure

    % Anchor point: outside top-left of axes
    x = op(1) + offPx(1);
    y = op(2) + op(4) + offPx(2);

    % TeX gives roman parentheses + italic letter
    str = sprintf('(\\it%s\\rm)', labels(k));

    hAnn(k) = annotation(fig,'textbox', ...
        'Units','pixels', ...
        'Position',[x y 40 20], ...   % small box; big enough for "(a)"
        'String',str, ...
        'Interpreter','tex', ...
        'FontName',p.Results.FontName, ...
        'FontSize',p.Results.FontSize, ...
        'LineStyle','none', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom');
    end
end
