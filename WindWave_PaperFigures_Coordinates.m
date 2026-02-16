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
fig_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/pdf_test7';
curvilinear_path = fullfile(project_path, 'curvilinear_new');
cartesian_path = fullfile(project_path, 'means');
data_path = fullfile(project_path, 'data');

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
instantaneous = load(fullfile(data_path, strcat(caze, '_DATA.mat')));
cartesian = load(fullfile(cartesian_path, strcat(caze, '_MEANS.mat')));
cartesian = cartesian.output;
curvilinear = load(fullfile(curvilinear_path, strcat(caze, '_CURVILINEAR.mat')));
curvilinear = curvilinear.output;


clc; fprintf('All %s cases loaded\n', wind_speed)
clear caze tmp w no_wave_caze cartesian_path curvilinear_path data_path

%% Make a sexy ass plot

% Plot settings
phase = 1;
image_index = 4;

% Cartesian Coordinates
X = cartesian.X;
x = X(1,:);
Y = cartesian.Y;
y = Y(:,1);

% Curvilinear Coordinates
vertical_lines = curvilinear.phase(phase).vertical_lines;
horizontal_lines = curvilinear.phase(phase).horizontal_lines;

% Load data
u_data = instantaneous.output.U(:,:,cartesian.phase(phase).idxs(image_index));
v_data = instantaneous.output.V(:,:,cartesian.phase(phase).idxs(image_index));
crop_wave = imresize(cartesian.phase(phase).waves(image_index,:), [1, length(x)]);
reference_wave = cartesian.phase(phase).reference_wave;

% Fill missing data
u_data = inpaint_nans(double(u_data));
u_data(Y < crop_wave) = nan;
v_data = inpaint_nans(double(v_data));
v_data(Y < crop_wave) = nan;

% Plot cartesian grid on top, but crop out waves
X_grid = X;
Y_grid = Y;

X_grid(Y < crop_wave) = nan;
Y_grid(Y < crop_wave) = nan;
[rows, columns] = size(X);

%% Interpolate data into curvilinear components

% Fill NaNs in original data just in case
u_filled = u_data;
v_filled = v_data;
u_filled(isnan(u_filled)) = 0;
v_filled(isnan(v_filled)) = 0;

% Interpolate
u_interp = interp2(X, Y, u_filled, vertical_lines, horizontal_lines, 'cubic', nan);
v_interp = interp2(X, Y, v_filled, vertical_lines, horizontal_lines, 'cubic', nan);
u_interp(u_interp == 0) = nan;
v_interp(v_interp == 0) = nan;

% Crop wave
u_interp(horizontal_lines < crop_wave) = nan;
v_interp(horizontal_lines < crop_wave) = nan;

% Angles
dx = mean(diff(x), 'all');
dy = mean(diff(y), 'all');
[dz_dx, dz_dy] = gradient(horizontal_lines, dx, dy);
dz_dx(horizontal_lines < reference_wave) = nan;
dz_dy(horizontal_lines < reference_wave) = nan;
angles = atan2(dz_dx, dz_dy);
angles(horizontal_lines < reference_wave) = nan;

% Compute unit vectors
% Tangent
xi_hat_x = cos(angles);    
xi_hat_y = sin(angles);

% Normal
zeta_hat_x = -sin(angles);   
zeta_hat_y =  cos(angles);

% Project Cartesian velocities onto curvilinear axes
u_tangent = u_interp .* xi_hat_x + v_interp .* xi_hat_y;
u_normal  = -(u_interp .* zeta_hat_x + v_interp .* zeta_hat_y);


%% Plotting


wave_transparency = 0.25;

% Plot details
linewidth = 1.5;
gridStep = 10;
levels = 200;

black_trans = 1;
white_trans = 0;
black_color = [0.6, 0.6, 0.6];
gridLineWidth = 0.5;
white_multiplier = 4;

tickFontSize = 8;
labelFontSize = 10;
colorbarFontSize = 10;

color_map = 'plasma';
color_lower = 0.4;
color_upper = 1;


%%% JFM MAX WIDTH: 13 CM

clc; close all;
totalFigure = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,11,5.5]);
t = tiledlayout(1,2, 'TileSpacing', 'compact', 'padding', 'tight');

%%% Cartesian
h(1) = nexttile;
ax = gca;
set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
hold on

% Velocity
contourf(X, Y, u_data / u_inf, levels, 'linestyle', 'none')

% % Vertical lines
% Layer 2: White halo lines
for i = 1:gridStep:rows
    P1 = plot(X_grid(:,i), y, 'color', 'white', 'linewidth', white_multiplier * gridLineWidth);
    P1.Color(4) = white_trans;
end
for i = 3:gridStep:rows
    P1 = plot(x, Y_grid(i,:), 'color', 'white', 'linewidth', white_multiplier * gridLineWidth);
    P1.Color(4) = white_trans;
end

% Layer 3: Black lines on top of white
for i = 1:gridStep:rows
    P2 = plot(X_grid(:,i), y, 'color', black_color, 'linewidth', gridLineWidth);
    P2.Color(4) = black_trans;
end
for i = 3:gridStep:rows
    P2 = plot(x, Y_grid(i,:), 'color', black_color, 'linewidth', gridLineWidth);
    P2.Color(4) = black_trans;
end


% Plot waves
% plot(x, crop_wave, 'color', 'red', 'linewidth', linewidth)
plot(x, reference_wave, 'color', 'black', 'linewidth', linewidth)
colormap(flipud(slanCM(color_map)))

% Shaded region below wave
hFill = patch( ...
[x, fliplr(x)], ...
[crop_wave, -100 * ones(size(crop_wave))], ...
'k', ...
'FaceAlpha', wave_transparency, ...      
'EdgeColor', 'none', ...    
'HandleVisibility', 'off'); 
uistack(hFill, 'bottom')
hold off

% Colorbar
clim([color_lower, color_upper])
axis equal
xlim([0, 235])
ylim([-25, 200])
xticks(0:100:230)
% yticks(-50:50:200)
yticks(0:100:200)
title('$u \mathbin{/} u_{\infty}$', 'interpreter', 'latex', 'fontsize', labelFontSize)








%%% Curvilinear
h(2) = nexttile;
ax = gca;
set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
ax.YTickLabel = [];
hold on

% Velocity
contourf(vertical_lines, horizontal_lines, u_tangent / u_inf, levels, 'linestyle', 'none')
% % Vertical lines
% Layer 2: White halo lines
for i = 1:gridStep:rows
    P1 = plot(vertical_lines(:,i), y + 1.05 * reference_wave(i), 'color', 'white', 'linewidth', white_multiplier * gridLineWidth);
    P1.Color(4) = white_trans;
end
for i = 3:gridStep:rows
    P1 = plot(x, horizontal_lines(i,:), 'color', 'white', 'linewidth', white_multiplier * gridLineWidth);
    P1.Color(4) = white_trans;
end

% Layer 3: Black lines on top of white
for i = 1:gridStep:rows
    P2 = plot(vertical_lines(:,i), y + 1.05 * reference_wave(i), 'color', black_color, 'linewidth', gridLineWidth);
    P2.Color(4) = black_trans;
end
for i = 3:gridStep:rows
    P2 = plot(x, horizontal_lines(i,:), 'color', black_color, 'linewidth', gridLineWidth);
    P2.Color(4) = black_trans;
end

% Plot waves
% plot(x, crop_wave, 'color', 'red', 'linewidth', linewidth)
plot(x, reference_wave, 'color', 'black', 'linewidth', linewidth)
colormap(flipud(slanCM(color_map)))

% Shaded region below wave
hFill = patch( ...
[x, fliplr(x)], ...
[crop_wave, -100 * ones(size(crop_wave))], ...
'k', ...
'FaceAlpha', wave_transparency, ...      
'EdgeColor', 'none', ...    
'HandleVisibility', 'off'); 
uistack(hFill, 'bottom')
hold off

% Colorbarxticks(0:100:230)
clim([color_lower, color_upper + 0])
c = colorbar;
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
c.FontSize = tickFontSize;
c.Label.FontSize = colorbarFontSize;
c.Ticks = color_lower:0.2:color_upper + 0.05;
title('$u_{\xi} \mathbin{/} u_{\infty}$', 'interpreter', 'latex', 'fontsize', labelFontSize)


axis equal
linkaxes(h, 'xy')
xlim([0, 234])
ylim([-25, 200])
xticks(0:100:230)
% yticks(-50:50:200)
yticks(0:100:200)
xlabel(t, '$x$ [mm]', 'Interpreter','latex', 'FontSize',labelFontSize)
ylabel(t, '$y$ [mm]', 'Interpreter','latex', 'FontSize',labelFontSize)

% Add (a) (b) labels
addPanelLabels(h, {'a', 'b'},  'Offset', [-0.18, 1.2])




% Save figure
% pause(3)
% fig_name = strcat(wind_speed, '_WV', wave, '_CurvilinearCoordinates.pdf');
% exportgraphics(totalFigure, fullfile(fig_folder, fig_name), 'Resolution', 600, 'ContentType', 'image')
% close all 





%% Function

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
