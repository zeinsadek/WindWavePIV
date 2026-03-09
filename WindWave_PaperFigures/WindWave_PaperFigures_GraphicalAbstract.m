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
figure_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/pdf_test7';
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
cartesian = load(fullfile(cartesian_path, strcat(caze, '_MEANS.mat')));
cartesian = cartesian.output;
curvilinear = load(fullfile(curvilinear_path, strcat(caze, '_CURVILINEAR.mat')));
curvilinear = curvilinear.output;


clc; fprintf('All %s cases loaded\n', wind_speed)
clear caze tmp w no_wave_caze cartesian_path curvilinear_path data_path



%% Graphical abstract of curvilinear shear stress

% Plot settings
phase = 4;

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

% Near-wave cropping
vetical_crop_shift = 0;

cmin = 1E-3;
cmax = 6E-3;
levels = 200;


% Plot phase avg between cartesian and curvilinear side-by-side
gridLineWidth = 1;
linewidth = 3;
gridStep = 10;
line_trans = 0.4;
linecolor = 0.5 .* ones(1,3);

% Axes limits
ymin = -18;
ymax = 175;
y_range = ymax - ymin;  % 205


clc; close all
fig = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,6,5]);
tiledlayout(1,1,'Padding','tight', 'TileSpacing', 'none')

ax = gca;
ax.Units = 'normalized';
ax.Position = [0 0 1 1];   % fill entire figure
ax.LooseInset = [0 0 0 0]; % remove MATLAB padding

hold on

% Load velocities
curvilinear_uv = -curvilinear.phase(phase).uv;
curvilinear_uv_norm = curvilinear_uv / (u_inf^2);

% Manually crop to remove max_wave profile noise
curvilinear_uv_norm(Y <= vetical_crop_shift) = nan;
% curvilinear_uv_norm(y <= sgolayfilt(max_wave, 2, 5)) = nan;
% curvilinear_uv_norm(y <= max_wave) = nan;

% Plot contour
contourf(vertical_lines, horizontal_lines, curvilinear_uv_norm, ...
         levels, 'linestyle', 'none')


% Vertical lines
for i = 2 + gridStep:gridStep:rows
    P = plot(vertical_lines(:,i), y + 1.05 * reference_wave(i), 'color', linecolor, 'linewidth', gridLineWidth);
    P.Color(4) = line_trans;
end
% Horizontal lines
for i = 1:gridStep:rows
    P = plot(x, horizontal_lines(i,:), 'color', linecolor, 'linewidth', gridLineWidth);
    P.Color(4) = line_trans;
end
Pw = plot(x, reference_wave, 'color', 'black', 'linewidth', linewidth);
uistack(Pw, 'top')

axis equal
clim([cmin, cmax])

cmap = slanCM('RdPu');
n_blend = 50;  % number of entries to blend to white
for i = 1:n_blend
    frac = (i-1) / (n_blend-1);  % 0 at first entry, 1 at n_blend
    cmap(i,:) = (1-frac)*[1 1 1] + frac*cmap(n_blend,:);
end
colormap(cmap)


% Crop x inward by ~5 mm on each side
xmin = -0.5;
xmax = 235;
x_range = xmax - xmin;  % 220

% Compute y range to preserve 1.2:1 aspect ratio (figure is 6x5 cm)
y_range = x_range / (6/5);  % 183.33

% Anchor y range: set the bottom, compute the top
ymin = -20;
ymax = ymin + y_range;

axis equal
xlim([xmin, xmax])
ylim([ymin, ymax])

patch([x(:); flipud(x(:))], [reference_wave(:); ymin*ones(size(reference_wave(:)))], 'k', ...
      'FaceAlpha', 0.25, ...
      'EdgeColor', 'none', ...
      'handlevisibility', 'off');

hold off

set(gca,'XColor', 'none','YColor','none')

fig = gcf;
fig.Units = 'centimeters';
fig.Position = [10 10 6 5];   % [left bottom width height]


% Test save
% pause(3)
% exportgraphics(fig, fullfile(figure_folder, sprintf('WindWave_JFM_GraphicalAbstract_%s_WV%s_Phase_%1.0f.jpg', wind_speed, wave, phase)), 'resolution', 600);
% close all

