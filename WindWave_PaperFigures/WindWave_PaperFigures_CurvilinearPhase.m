% WindWave paper figures: Curvilinear Phase Averages

clc; clear; close all
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
curvilinear_path = fullfile(project_path, 'curvilinear');

% Figure Folder
figure_folder = 'pdf_test4';

% Cases
wind_speed = 'WT6';
waves = {'A', 'B', 'C', 'D'};

% Approximate wavelengths in mm for labeling plots
wavelengths.A = '410';
wavelengths.B = '313';
wavelengths.C = '189';
wavelengths.D = '124';

% Wavelengths in mm
wavelength_values.A = 410.8;
wavelength_values.B = 313.3;
wavelength_values.C = 189.6;
wavelength_values.D = 124.3;

% Wave steepnesses 
steepnesses.A = '0.180';
steepnesses.B = '0.211';
steepnesses.C = '0.305';
steepnesses.D = '0.267';

if ismember(wind_speed(end), {'4'})
    u_inf = 2.4181;
elseif ismember(wind_speed(end), {'6'})
    u_inf = 3.8709;
elseif ismember(wind_speed(end), {'8'})
    u_inf = 5.4289;
end


% Load all wave cases
for w = 1:length(waves)
    wave = waves{w};
    caze = strcat(wind_speed, '_WV', wave, '_AG0');
    disp(caze)

    % Load data
    tmp = load(fullfile(curvilinear_path, strcat(caze, '_CURVILINEAR.mat')));
    tmp = tmp.output;

    % Save data
    data.(caze) = tmp;
end
cazes = fields(data);
clc; fprintf('All %s cases loaded\n', wind_speed)

clear caze tmp w wave no_wave_caze


%% Plot all phases single case

levels = 100;
linewidth = 2;

colorbar_fontsize = 16;
tickFontSize = 14;
labelFontSize = 16;
legendFontSize = 14;

wave = 'C';
component = 'u';
caze = strcat(wind_speed, '_WV', wave, '_AG0');
wave_length = 1;

phase_labels = {'$\varphi = 0$', '$\varphi = \lambda / 4$', ...
                '$\varphi = \lambda / 2$', '$\varphi = 3 \lambda / 4$'};


if strcmp(component, 'u')
    labelComponent = "u_{\xi}";
elseif strcmp(component, 'v')
    labelComponent = "u_{\zeta}";
elseif strcmp(component, 'uu')
    labelComponent = "u_{\xi}' u_{\xi}'";
elseif strcmp(component, 'vv')
    labelComponent = "u_{\zeta}' u_{\zeta}'";
elseif strcmp(component, 'uv')
    labelComponent = "u_{\xi}' u_{\zeta}'";
else
    labelComponent = component;
end

clc; close all;
totalFigure = figure('color', 'white', 'units', 'centimeters', 'position', [5,5,30,7.5]);
t = tiledlayout(1, 4, 'TileSpacing', 'tight', 'padding', 'tight');

% Normalize velocities vs stresses
if ismember(component, {'u', 'v'})
    norm = u_inf;
    colorbar_label = sprintf('$ \\langle %s \\rangle / u_{\\infty}$', labelComponent);
else
    norm = u_inf^2;
    colorbar_label = sprintf('$ \\langle %s \\rangle / u_{\\infty}^2$', labelComponent);
end

% Select appropriate colormap
if ismember(component, {'u', 'uu', 'vv'})
    cmap = 'parula';
else
    cmap = 'coolwarm';
end

% Loop through all phases
for phase = 1:4

    h(phase) = nexttile;

    % Make tick marks look fancy
    ax = gca;
    set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)

    tmp = data.(caze).phase(phase).(component) / norm;
    tmp(data.(caze).phase(phase).vertical_lines < 0) = nan;
    tmp(data.(caze).phase(phase).vertical_lines > 236) = nan;

    hold on
    % Plot contour
    contourf(data.(caze).phase(phase).vertical_lines / wave_length, ...
             data.(caze).phase(phase).horizontal_lines / wave_length, ...
             tmp, ...
             levels, 'linestyle', 'none')

    % Plot wave cropping
    plot(data.(caze).X(1,:) / wave_length, data.(caze).phase(phase).max_wave_profile / wave_length, 'linewidth', linewidth, 'color', 'red')

    % Plot Reference wave profile
    plot(data.(caze).X(1,:) / wave_length, data.(caze).phase(phase).wave_profile / wave_length, 'linewidth', linewidth, 'color', 'black')
    hold off

    axis equal
    ylim([-20 / wave_length, 200 / wave_length])
    title(phase_labels{phase}, 'interpreter', 'latex', 'FontSize', labelFontSize)

    if phase == 1
        ylabel('$y$ [mm]', 'interpreter', 'latex', 'FontSize', labelFontSize)
    end

    if phase ~= 1
        ax.YTickLabel = [];
    end

    if phase == 4
        colormap(cmap)
        C = colorbar;
        C.Label.String = colorbar_label;
        C.Label.Interpreter = 'latex';
        C.Label.FontSize = colorbar_fontsize;
        C.TickLabelInterpreter = 'latex';
    end

    clim([0.25, 1])
    clear c
end

% Sync color ranges
ax = findall(gcf, 'Type', 'axes');
clims = cell2mat(get(ax, 'CLim'));
global_clim = [min(clims(:,1)), max(clims(:,2))];
set(ax, 'CLim', global_clim);

% Add a shared x-axis label to the entire layout
xlabel(t, '$x$ [mm]', 'interpreter', 'latex', 'FontSize', labelFontSize);
linkaxes(h,'xy')

% Save figure
fig_folder = fullfile('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/', figure_folder, 'Phase');
fig_name = strcat(caze(1:end - 4), '_AllPhases_', component, '.pdf');
exportgraphics(totalFigure, fullfile(fig_folder, fig_name), 'Resolution', 600, 'ContentType', 'image')
close all 

clear ax C clims cmap colorbar_fontsize colorbar_label component global_clim h 
clear levels linewidth norm phase wave


%% Plot single phase all cases

wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};

levels = 100;
linewidth = 2;

colorbar_fontsize = 16;
tickFontSize = 14;
labelFontSize = 16;
legendFontSize = 14;


phase = 2;
component = 'u';

% Convert to curvilinear notation
if strcmp(component, 'u')
    labelComponent = "u_{\xi}";
elseif strcmp(component, 'v')
    labelComponent = "u_{\zeta}";
elseif strcmp(component, 'uu')
    labelComponent = "u_{\xi}' u_{\xi}'";
elseif strcmp(component, 'vv')
    labelComponent = "u_{\zeta}' u_{\zeta}'";
elseif strcmp(component, 'uv')
    labelComponent = "u_{\xi}' u_{\zeta}'";
else
    labelComponent = component;
end


clc; close all;
totalFigure = figure('color', 'white', 'units', 'centimeters', 'position', [5,5,33,15]);
tiledlayout(2, 4, 'TileSpacing', 'tight', 'padding', 'tight');

% Normalize velocities vs stresses
if ismember(component, {'u', 'v'})
    norm = u_inf;
    colorbar_label = sprintf('$\\langle %s \\rangle / u_{\\infty} $', labelComponent);
else
    norm = u_inf^2;
    colorbar_label = sprintf('$\\langle %s \\rangle / u_{\\infty}^2$', labelComponent);
end

% Select appropriate colormap
if ismember(component, {'u', 'uu', 'vv'})
    cmap = 'parula';
else
    cmap = 'coolwarm';
end

% Lopp through different waves
for c = 1:length(cazes)

    disp(c)
    caze = cazes{c};
    wave = caze(7);

    % wave_length = wavelength_values.(wave);
    wave_length = 1;

    % Make tick marks look fancy
    h(c) = nexttile;
    ax = gca;
    set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)

    tmp = data.(caze).phase(phase).(component) / norm;
    tmp(data.(caze).phase(phase).vertical_lines < 0) = nan;
    tmp(data.(caze).phase(phase).vertical_lines > 236) = nan;

    hold on

    % Plot contour
    contourf(data.(caze).phase(phase).vertical_lines / wave_length, ...
             data.(caze).phase(phase).horizontal_lines / wave_length, ...
             tmp, ...
             levels, 'linestyle', 'none')

    % Plot line of constant \xi
    index = round(size(tmp,2)/2);
    P = plot(data.(caze).phase(phase).vertical_lines(:, index) / wave_length, data.(caze).Y(:,1) / wave_length, ...
             'color', 'black', 'linestyle', ':', 'linewidth', 2, 'HandleVisibility', 'off');
    P.Color(4) = 0.25;

    % Plot wave crop
    plot(data.(caze).X(1,:) / wave_length, data.(caze).phase(phase).max_wave_profile / wave_length, 'linewidth', linewidth, 'color', 'red')

    % Plot reference wave
    plot(data.(caze).X(1,:) / wave_length, data.(caze).phase(phase).wave_profile / wave_length, 'linewidth', linewidth, 'color', 'black')

    hold off


    % Scale axes to wavelength but keep aspect ratio 1:1
    axis equal

    % absolute limits in [mm]
    xmin = -10;
    xmax = 250;
    ymin = -20;
    ymax = 200;

    x_range = (xmax) - (xmin);
    y_range = (ymax) - (ymin);

    aspect = x_range / y_range;

    xlim([xmin/wave_length, xmax/wave_length])
    ylim([ymin/wave_length, ymax/wave_length])
    xticks(0:50:200)
    yticks(0:50:200)

    % 
    % % Force identical inner plot boxes for every tile
    % h(c).PositionConstraint = "innerposition";
    % axis(h(c),'normal')
    % pbaspect(h(c), [aspect 1 1])
    
    % ticks (your logic is fine)
    % yticks(0:0.25:(200 / wave_length))
    % xticks(0:0.25:(250 / wave_length))


    % Titles
    if ismember(wave, waves)
        title_label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelengths.(wave), steepnesses.(wave));
    else
        title_label = sprintf('$\\lambda_{0}$');
    end
    title(title_label, 'interpreter', 'latex', 'fontsize', labelFontSize)

    if c == 1
        ylabel('$y$ [mm]', 'interpreter', 'latex', 'FontSize', labelFontSize)
    end

    if c == 2
        xlabel('     ', 'interpreter', 'latex', 'FontSize', labelFontSize)
    end

    if c ~= 1
        ax.YTickLabel = [];
    end

    if c == 4
        colormap(cmap)
        C = colorbar;
        C.Label.String = colorbar_label;
        C.Label.Interpreter = 'latex';
        C.Label.FontSize = colorbar_fontsize;
        C.TickLabelInterpreter = 'latex';
    end

    if strcmp(component, 'u')
        clim([0.2, 1])
    elseif strcmp(component, 'uv')
        clim([-10E-3, 0])
    elseif strcmp(component, 'vv')
        clim([0, 12E-3])
    end

    clear c
end

% Sync color ranges
ax = findall(gcf, 'Type', 'axes');
clims = cell2mat(get(ax, 'CLim'));
global_clim = [min(clims(:,1)), max(clims(:,2))];
set(ax, 'CLim', global_clim);

% Plot profiles
nexttile([1,4])

% Make tick marks look fancy
ax = gca;
set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)

hold on
for c = 1:length(cazes)

    caze = cazes{c};
    wave = caze(7);
    wave_length = 1;

    tmp = data.(caze).phase(phase).(component) / norm;
    tmp(data.(caze).phase(phase).vertical_lines < 0) = nan;
    tmp(data.(caze).phase(phase).vertical_lines > 236) = nan;

    index = round(size(tmp,2)/2);

    y_profile = data.(caze).phase(phase).horizontal_lines(:, index);
    profile = tmp(:, index);

    valid = ~isnan(y_profile) & ~isnan(profile);

    y_profile = y_profile(valid);
    profile = profile(valid);

    buffer = 2;
    y_profile = y_profile(1:end - buffer);
    profile = profile(1:end - buffer);

    % Titles
    if ismember(wave, waves)
        label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelengths.(wave), steepnesses.(wave));
    else
        label = sprintf('$\\lambda_{0}$');
    end

    % Color profile based on wave
    if strcmp(wave, 'A')
        color = wave_colors{1};
    elseif strcmp(wave, 'B')
        color = wave_colors{2};
    elseif strcmp(wave, 'C')
        color = wave_colors{3};
    elseif strcmp(wave, 'D')
        color = wave_colors{4};
    else 
        color = 'black';
    end

    % Plot a slice
    plot(y_profile / wave_length, profile, ...
         'linewidth', linewidth, 'displayname', label, 'color', color)

    % clear y_profile profile
end
hold off
xlabel('$\zeta$ [mm]', 'interpreter', 'latex', 'FontSize', labelFontSize)

% ylim([0.4, 1.1])
% yticks(0.4:0.1:1.1)
% xlim([0, 1.75])
% xticks(0:0.25:1.75)

ylabel(colorbar_label, 'interpreter', 'latex', 'FontSize', colorbar_fontsize);

if strcmp(component, 'u')
    loc = 'southeast';
elseif strcmp(component, 'vv')
    loc = 'northeast';
elseif strcmp(component, 'uv')
    loc = 'southeast';
end

legend('location', loc, 'interpreter', 'latex', 'box', 'off', 'fontsize', legendFontSize)


% Add a centered xlabel across the top row
annotation(totalFigure, 'textbox', [0.23, 0.46, 0.5, 0.04], ...
    'String', '$x$ [mm]', ...
    'Interpreter', 'latex', ...
    'FontSize', labelFontSize, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'EdgeColor', 'none');



% Save figure
% fig_folder = fullfile('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/', figure_folder, 'Phase');
% fig_name = strcat(wind_speed, '_Phase', num2str(phase), '_', component, '_WavelengthScaled.pdf');
% exportgraphics(totalFigure, fullfile(fig_folder, fig_name), 'Resolution', 600, 'ContentType', 'image')
% close all 


clear ax C caze clims cmap colorbar_fontsize colorbar_label component global_clim h 
clear levels linewidth norm phase wave h
