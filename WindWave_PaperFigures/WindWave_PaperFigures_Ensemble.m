% WindWave paper figures: Ensemble Averages

clc; clear; close all
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')
addpath('/Users/zeinsadek/Documents/MATLAB/MatlabFunctions/Inpaint_nans')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
means_path = fullfile(project_path, 'means');

% Cases
wind_speed = 'WT6';
waves = {'A', 'B', 'C', 'D'};

% Approximate wavelengths in mm for labeling plots
wavelengths.A = '410';
wavelengths.B = '313';
wavelengths.C = '189';
wavelengths.D = '124';

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

% Load No-Wave case
no_wave_caze = strcat(wind_speed, '_WV0_AGP');
disp(no_wave_caze)

% Load data
tmp = load(fullfile(means_path, strcat(no_wave_caze, '_MEANS.mat')));
tmp = tmp.output;

% Save data
data.(no_wave_caze) = tmp;

% Load all wave cases
for w = 1:length(waves)
    wave = waves{w};
    caze = strcat(wind_speed, '_WV', wave, '_AG0');
    disp(caze)

    % Load data
    tmp = load(fullfile(means_path, strcat(caze, '_MEANS.mat')));
    tmp = tmp.output;

    % Save data
    data.(caze) = tmp;
end
clc; fprintf('All %s cases loaded\n', wind_speed)

clear caze tmp w wave no_wave_caze


%% Repair cat-scratch

cazes = fields(data);
components = {'u', 'v', 'uu', 'vv', 'uv'};

clc;
for w = 1:length(cazes)
    caze = cazes{w};
    for c = 1:length(components)
        component = components{c};
    
        % Monitor
        fprintf("%s\n", component)
    
        % Load component
        tmp = data.(caze).ensemble.(component);
    
        % Shift RHS
        LHS_values = tmp(:,86);
        RHS_values = tmp(:,87);
        differences = RHS_values - LHS_values;
        differences(isnan(differences)) = 0;
        corrected = tmp;
        corrected(:,87:end) = corrected(:,87:end) - differences;
    
        % Interpolate across strip
        x = data.(caze).X(1,:);
        [~, left] = min(abs(x - 115));
        [~, right] = min(abs(x - 120));
        corrected(:, left:right) = nan;
        corrected = fillmissing(corrected, 'makima', 2);
        corrected(data.(caze).Y < data.(caze).ensemble.max_wave_profile) = nan;
        
        % Make sure original nans are in place
        corrected(isnan(tmp)) = nan;
    
    
        % Resave data
        % means_corrected.phase(phase).(component) = corrected;
        data.(caze).ensemble.(component) = corrected;
        clear corrected
    
        fprintf("\n")
    end
end
clc; fprintf('Cat scratch repaired on all components\n')

clear means components tmp LHS_values RHS_values differences left right
clear c component levels caze w x



%% Plot: contours with profiles below

wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};

levels = 100;
linewidth = 3;

colorbar_fontsize = 16;
tickFontSize = 8;
labelFontSize = 16;
legendFontSize = 10;

component = 'uv';

% Add tickmarks to the stresses
if strcmp(component, 'uu')
    labelComponent = "u'u'";
elseif strcmp(component, 'vv')
    labelComponent = "v'v'";
elseif strcmp(component, 'uv')
    labelComponent = "u'v'";
else
    labelComponent = component;
end


clc; close all;
totalFigure = figure('color', 'white', 'units', 'centimeters', 'position', [2,2,30,15]);
tiledlayout(2, length(cazes), 'TileSpacing', 'tight', 'padding', 'tight')

% Normalize velocities vs stresses
if ismember(component, {'u', 'v'})
    norm = u_inf;
    colorbar_label = sprintf('$\\overline{%s} / u_{\\infty}$', labelComponent);
else
    norm = u_inf^2;
    colorbar_label = sprintf('$\\overline{%s} / u_{\\infty}^2$', labelComponent);
end

% Select appropriate colormap
if ismember(component, {'u', 'uu', 'vv'})
    cmap = 'parula';
else
    cmap = 'coolwarm';
end

% Plot contours
for c = 1:length(cazes)
    caze = cazes{c};

    wave = caze(7);
    
    h(c) = nexttile;

    % Make tick marks look fancy
    ax = gca;
    set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)

    % Fill in slivers of missing data that looked annoying
    tmp = data.(caze).ensemble.(component) / norm;
    tmp = inpaint_nans(double(tmp));
    tmp(data.(caze).Y < data.(caze).ensemble.max_wave_profile) = nan;

    hold on
    contourf(data.(caze).X, data.(caze).Y, tmp, levels, 'linestyle', 'none')
    plot(data.(caze).X(1,:), data.(caze).ensemble.max_wave_profile, 'linewidth', linewidth, 'color', 'red')
    hold off

    axis equal
    ylim([-20, 200])
    yticks(-50:50:200)
    
    % Titles
    if ismember(wave, waves)
        title_label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelengths.(wave), steepnesses.(wave));
    else
        % title_label = sprintf('$\\lambda_{0}$');
        title_label = 'No Waves';
    end
    title(title_label, 'interpreter', 'latex', 'fontsize', labelFontSize)

    
    if c == 1
        ylabel('$y$ [mm]', 'interpreter', 'latex', 'FontSize', labelFontSize)
    end

    if c == 3
        xlabel('$x$ [mm]', 'interpreter', 'latex', 'FontSize', labelFontSize)
    end

    if c~= 1
        ax.YTickLabel = [];
    end

    if c == length(cazes)
        colormap(cmap)
        C = colorbar;
        C.Label.String = colorbar_label;
        C.Label.Interpreter = 'latex';
        C.Label.FontSize = colorbar_fontsize;
        C.TickLabelInterpreter = 'latex';
    end

    clear c ax tmp C caze wave
end

% Sync color ranges
ax = findall(gcf, 'Type', 'axes');
clims = cell2mat(get(ax, 'CLim'));
global_clim = [min(clims(:,1)), max(clims(:,2))];
set(ax, 'CLim', global_clim);

% Plot profiles
nexttile([1,5])

% Make tick marks look fancy
ax = gca;
set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)

hold on
for c = 1:length(cazes)
    caze = cazes{c};
    wave = caze(7);

    tmp = data.(caze).ensemble.(component) / norm;

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

    % Titles
    if ismember(wave, waves)
        label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelengths.(wave), steepnesses.(wave));
    else
        % label = sprintf('$\\lambda_{0}$');
        label = 'No Waves';
    end

    % Plot a slice
    plot(data.(caze).Y(:,1), tmp(:, round(size(tmp,2)/2)), ...
         'linewidth', linewidth, 'displayname', label, 'color', color)

    % Plot a streamwise average
    % tmp(data.(caze).Y <= max(data.(caze).ensemble.max_wave_profile)) = nan;
    % plot(data.(caze).Y(:,1), mean(tmp, 2, 'omitnan'), 'linewidth', linewidth, 'displayname', wave)

    clear c caze wave
end
hold off
xlabel('$y$ [mm]', 'interpreter', 'latex', 'FontSize', labelFontSize)
xlim([0, 200])
xticks(-50:50:200)
ylabel(colorbar_label, 'interpreter', 'latex', 'FontSize', colorbar_fontsize);

if strcmp(component, 'u')
    loc = 'southeast';
elseif strcmp(component, 'v')
    loc = 'northeast';
elseif strcmp(component, 'uv')
    loc = 'southeast';
end

legend('location', loc, 'interpreter', 'latex', 'box', 'off', 'fontsize', legendFontSize)

% xscale('log')

linkaxes(h,'xy')

% Save figure
fig_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/pdf_test2/Ensemble';
fig_name = strcat(wind_speed, '_Ensemble_', component, '_w_Profiles.pdf');
exportgraphics(totalFigure, fullfile(fig_folder, fig_name), 'Resolution', 600, 'ContentType', 'image')
close all 


    
     