%%% Cartesian vs Curvilinear Velocities and Stresses

addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
means_path = fullfile(project_path, 'means');
curvilinear_path = fullfile(project_path, 'curvilinear');

% Case
caze = 'WT8_WVD_AG0';

% Figures
figure_path = fullfile(project_path, "figures", caze, "cartesian_vs_curvilinear");

if ~exist(figure_path, 'dir')
    mkdir(figure_path)
end

% Free stream
tunnel_freq = caze(strfind(caze, 'WT') + 2);

if ismember(tunnel_freq, {'4'})
    u_inf = 2.4181;
elseif ismember(tunnel_freq, {'6'})
    u_inf = 3.8709;
elseif ismember(tunnel_freq, {'8'})
    u_inf = 5.4289;
end

fprintf("WT%s = %1.2f m/s Freestream\n\n", tunnel_freq, u_inf)

% Load data
cartesian = load(fullfile(means_path, strcat(caze, '_MEANS.mat')));
cartesian = cartesian.output;

curvilinear = load(fullfile(curvilinear_path, strcat(caze, '_CURVILINEAR.mat')));
curvilinear = curvilinear.output;

% Wave Parameters
wave_type       = caze(strfind(caze, 'WV') + 2);
wave_parameters = readcell("Offshore_Waves.xlsx");
wavelength      = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 2};
amplitude       = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 3};

clear means_path curvilinear_path project_path wave_parameters wave_type

%% Plots

components = {'u', 'v', 'uu', 'vv', 'uv'};

% Grid lines
lw = 2;
spacing = 5;
scale = 0.4;
quiver_color = 'black';


for c = 1:length(components)

    % Component
    component = components{c};
    if ismember(component, {'v'})
        map = 'coolwarm';
    else
        map = 'parula';
    end

    % How to normalize (velocity vs stress)
    if ismember(component, {'u', 'v'})
        normalization = u_inf;
    elseif ismember(component, {'uu', 'vv', 'uv'})
        normalization = u_inf^2;
    end
    
    % Colorbar limits and normalize
    cmin = min([cartesian.phase(:).(component), curvilinear.phase(:).(component)], [], 'all') / normalization;
    cmax = max([cartesian.phase(:).(component), curvilinear.phase(:).(component)], [], 'all') / normalization;
    
    close all
    ax = figure('position', [0, 0, 1920, 1080]);
    tiledlayout(2,4)
    sgtitle(sprintf("%s: Phase Averaged %s Cartesian vs Curvilinear", caze, component), 'interpreter', 'none')
    
    % Cartesian
    for p = 1:4
        h(p) = nexttile;
        hold on
        contourf(cartesian.X, cartesian.Y, cartesian.phase(p).(component) / normalization, 100, 'LineStyle', 'none')
        plot(cartesian.X(1,:), cartesian.phase(p).max_wave_profile, 'color', 'red', 'linewidth', lw)
        plot(cartesian.X(1,:), cartesian.phase(p).reference_wave, 'color', 'black', 'linewidth', lw)
    
        % Grid lines
        X = cartesian.X;
        Y = cartesian.Y;
        X(Y < cartesian.phase(p).reference_wave.') = nan;
        Y(Y < cartesian.phase(p).reference_wave.') = nan;
    
        for i = 1:spacing:length(X)
            P = plot(X(i,:), Y(i,:), 'color', 'black');
            P.Color(4) = 0.25;
        end
    
        for i = 1:spacing:length(Y)
            P = plot(X(:,i), Y(:,i), 'color', 'black');
            P.Color(4) = 0.25;
        end

        quiver(cartesian.X(1:spacing:end, 1:spacing:end), ...
               cartesian.Y(1:spacing:end, 1:spacing:end), ...
               cartesian.phase(p).u(1:spacing:end, 1:spacing:end) / normalization, ...
               cartesian.phase(p).v(1:spacing:end, 1:spacing:end) / normalization, scale, 'color', quiver_color)
    
        hold off
        axis equal
        xlim([0, 235])
        ylim([-30, 205])
        colorbar
        colormap(map)
        clim([cmin, cmax])
        title(sprintf("Cartesian Phase %1.0f", p))
    end
    
    % Curvilinear
    for p = 1:4
        h(p + 4) = nexttile;
        hold on
        contourf(curvilinear.phase(p).vertical_lines, curvilinear.phase(p).horizontal_lines, curvilinear.phase(p).(component) / normalization, 100, 'LineStyle', 'none')
        plot(cartesian.X(1,:), cartesian.phase(p).max_wave_profile, 'color', 'red', 'linewidth', lw)
        plot(cartesian.X(1,:), cartesian.phase(p).reference_wave, 'color', 'black', 'linewidth', lw)
    
        % Grid lines
        [row, col] = size(curvilinear.phase(p).vertical_lines);
        for i = 1:spacing:row
            P = plot(curvilinear.X(1,:), curvilinear.phase(p).horizontal_lines(i,:), 'color', 'black');
            P.Color(4) = 0.25;
        end
    
        for i = 1:spacing:col
            P = plot(curvilinear.phase(p).vertical_lines(:,i), curvilinear.Y(:,1) + curvilinear.phase(p).wave_profile(i), 'color', 'black');
            P.Color(4) = 0.25;
        end

       quiver(curvilinear.phase(p).vertical_lines(1:spacing:end, 1:spacing:end), ...
              curvilinear.phase(p).horizontal_lines(1:spacing:end, 1:spacing:end), ...
              curvilinear.phase(p).u(1:spacing:end, 1:spacing:end) / normalization, ...
              curvilinear.phase(p).v(1:spacing:end, 1:spacing:end) / normalization, scale, 'color', quiver_color)
    
        hold off
        axis equal
        xlim([0, 235])
        ylim([-30, 205])
        colorbar
        colormap(map)
        clim([cmin, cmax])
        title(sprintf("Curvilinear Phase %1.0f", p))
    end    
    linkaxes(h, 'xy')
    pause(2);

    % Save figure
    filename = strcat(caze, "_CartesianVsCurvilinear_", component, "_Normalized.png");
    exportgraphics(ax, fullfile(figure_path, filename), 'resolution', 300);
end

close all
clear h i p row col X Y P component spacing lw ax c cmax cmin map filename


