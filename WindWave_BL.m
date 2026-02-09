%% Computing Ensemble/Phase Average Boundary Layer Parameters
%%% Compute Boundary Layer Parameters
% Zein Sadek
% 11/2023

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Offshore/Offshore_Functions/');
experiment_log  = readcell('Offshore_Inflow_Log.xlsx');
wave_parameters = readcell('Offshore_Waves.xlsx');

s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data
experiment_name = 'WT8_WV0_AGP';

% Experiment Specifics
tunnel_freq = experiment_name(strfind(experiment_name, 'WT') + 2);
wave_type   = experiment_name(strfind(experiment_name, 'WV') + 2);

% Paths
results_path = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/';
mtlb_file    = strcat(results_path, 'data'   , '/', experiment_name, '_DATA.mat');
mean_file    = strcat(results_path, 'means'  , '/', experiment_name, '_MEANS.mat');
phase_file   = strcat(results_path, 'phase'  , '/', experiment_name, '_PHASE.mat');

if wave_type ~= '0'
    wave_file    = strcat(results_path, 'wave'   , '/', experiment_name, '_WAVE.mat');
end

% Save for BL thicknesses
BL_file       = strcat(results_path, 'boundary_layer/', experiment_name, '_BL.mat');
figure_folder = strcat(results_path, 'figures/', experiment_name, '/boundary_layer/');

if ~exist(figure_folder, 'dir')
    mkdir(figure_folder)
end

% Cropping Parameters
top_bound_value   = 205;       % relative to Y centered at still water
left_bound_value  = -121;      % relative to X centered at DaVis default
right_bound_value = 115;       % relative to X centered at DaVis default

% Wave Parameters
if wave_type == '0'
    wave_length    = 0;
    wave_amplitude = 0;
    phase_offset   = [0, wave_length/4, wave_length/2, 3*wave_length/4];
else
    wave_length    = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 2};
    wave_amplitude = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 3};
    wave_frequency = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 4};
    phase_offset   = [0, wave_length/4, wave_length/2, 3*wave_length/4];
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_figures = false;
figure_resolution = 300;
marker_alpha = 0.5;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB DATA TO ENSEMBLE/PHASE MEANS AND WAVES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if wave_type ~= '0'
%     waves = load(wave_file);
%     waves = waves.output;
% end

means = load(mean_file); 
means = means.output;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOUNDARY LAYER THICKNESS (Ensemble)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coordinates
Y = means.Y;
y = flipud(unique(Y));
X = means.X;
x = unique(X);
range = abs(left_bound_value) + abs(right_bound_value);
output.x = x;


% Selct center of image
x_middle_value    = range/2;
y_top_value       = (3/4) * top_bound_value;
[~, x_middle_idx] = min(abs(x - x_middle_value));
[~, y_top_idx]    = min(abs(y - y_top_value));

% Visualize where wer are sampling freestream velocity
box_size    = 20;
copy        = means.ensemble.u;
copy(Y < 0) = nan;
copy(y_top_idx - box_size:y_top_idx, x_middle_idx - box_size/2: x_middle_idx + box_size/2) = nan;

% Plot freestream sample location
ax = figure('Name', 'Freestream Sample');
hold on
contourf(X, Y, copy, 50, 'LineStyle', 'none')
yline(0, 'Color', 'black', 'LineWidth', 2)
plot(x, means.ensemble.max_wave_profile, 'Color', 'red', 'LineWidth', 3)
hold off

axis equal
xlim([0, range])
ylim([-20, top_bound_value])

C = colorbar();
C.Label.String = '${u}$ [m/s]';
C.Label.Interpreter = 'Latex';
xlabel('x [mm]', 'Interpreter', 'Latex')
ylabel('z [mm]', 'Interpreter', 'Latex')

% Average over sampled region
u_inf = mean(means.ensemble.u(y_top_idx - box_size:y_top_idx, x_middle_idx - box_size/2: x_middle_idx + box_size/2), 'all', 'omitnan');
title(strcat(experiment_name, {' '}, 'Freestream Velocity:', {' '}, num2str(round(u_inf, 2)), ' m/s'), 'Interpreter', 'none');
fprintf('Freestream Velocity: %3.3f m/s \n\n', u_inf)

% Save freestream
output.u_inf = u_inf;

% Save Figure
if save_figures == 1
    fprintf('Saving Figure...\n')
    figure_name = strcat(experiment_name, '_ens_avg_u_inf.png');
    exportgraphics(ax, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
    fprintf('Figure Saved\n\n')
    close all; clc;
end


%%% Find Boundary Layer Edge
% Define boundary layer edge velocity
smooth_kernel  = 60;
BL_crop        = means.ensemble.u;
BL_crop_smooth = juliaan_smooth(BL_crop, smooth_kernel);

% Save kernel
output.smooth_kernel = smooth_kernel;

%%% Visualize kernel size
% figure('Name', 'Kernel')
% hold on
% contourf(X, Y, means.ensemble.u, 500, 'LineStyle', 'none')
% axis equal
% for i = 1:smooth_kernel:length(x)
%    xline(x(i))
%    yline(y(i))
% end
% hold off
% xlim([0, range])
% ylim([-20, top_bound_value])
% title(strcat('Smoothing Kernel Size:', {' '}, num2str(smooth_kernel), ' x ', num2str(smooth_kernel), ' pixels'))

% Remove values outside boundary layer
BL_percent = 0.94;
BL_edge    = BL_percent * u_inf;

BL_crop(BL_crop_smooth > BL_edge)        = nan;
BL_crop_smooth(BL_crop_smooth > BL_edge) = nan;

% BL Thickness
[~, ens_I] = max(BL_crop_smooth, [], 1, 'omitnan');
ens_BL     = Y(ens_I);

%%% First and last points suck
ens_BL(1) = nan;
ens_BL(end) = nan;

% Save ensemble BL
output.ensemble.thickness = ens_BL;

%%% Plot Ensemble Avg Boundary Layer
ax = figure('Name', 'Ensemble BL');
hold on
contourf(X, Y, means.ensemble.u, 500, 'Linestyle', 'none')
yline(0, 'Color', 'black', 'LineWidth', 2)
plot(x, means.ensemble.max_wave_profile, 'Color', 'red', 'LineWidth', 3)
plot(x, ens_BL, 'Color', 'black', 'LineWidth', 2)
hold off

axis equal
xlim([0, range])
ylim([-20, top_bound_value])

C = colorbar();
C.Label.String = '${u}$ [m/s]';
C.Label.Interpreter = 'Latex';
xlabel('x [mm]', 'Interpreter', 'Latex')
ylabel('z [mm]', 'Interpreter', 'Latex')
title(strcat(experiment_name, ': Ens Avg Boundary Layer'), 'Interpreter', 'none')

% Save Figure
if save_figures  == 1
    fprintf('Saving Figure...\n')
    figure_name = strcat(experiment_name, '_ens_avg_BL.png');
    exportgraphics(ax, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
    fprintf('Figure Saved\n\n')
    close all; clc;
end


%%% Stacked boundary layer profiles
C = colormap(jet(length(x)));
ax = figure('Name', 'Ensemble Profiles');
hold on
for i = 1:length(x)
    u         = means.ensemble.u(:, i);
    mask      = ~isnan(u);
    u         = u(mask);
    y_mask    = y(mask);
    integrand = 1 - (u / u_inf);

    scatter(u, y_mask, 'MarkerFaceColor', C(i,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.1)
end

hold off
ylim([-10, top_bound_value])
xlabel('u [m/s]', 'Interpreter', 'Latex')
ylabel('z [mm]', 'Interpreter', 'Latex')
title(strcat(experiment_name, ': Ens Avg Velocity Profiles'), 'Interpreter', 'none')

% Save Figure
if save_figures  == 1
    fprintf('Saving Figure...\n')
    figure_name = strcat(experiment_name, '_ens_avg_BL_profiles.png');
    exportgraphics(ax, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
    fprintf('Figure Saved\n\n')
    close all; clc;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLACEMENT AND MOMENTUM THICKNESS (Ensemble)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% trapz(coordinates, data)

ens_DT = zeros(1, length(x));
ens_MT = zeros(1, length(x));

for i = 1:length(x)

    u         = BL_crop(:, i);
    mask      = ~isnan(u);
    u         = u(mask);
    y_mask    = y(mask);

    DT_integrand = 1 - (u / u_inf);
    MT_integrand = (u / u_inf) .* (1 - (u / u_inf));
    
    ens_DT(1,i)  = trapz(flipud(y_mask), flipud(DT_integrand));
    ens_MT(1,i)  = trapz(flipud(y_mask), flipud(MT_integrand));

end

%%%% First and last points suck
ens_DT(1,1)   = nan;
ens_DT(1,end) = nan;

ens_MT(1,1)   = nan;
ens_MT(1,end) = nan;


%%% Add offset
% ens_DT = ens_DT + means.ensemble.max_wave_profile;
% ens_MT = ens_MT + means.ensemble.max_wave_profile;

% Save displacement and momentum thicknesses
output.ensemble.displacement = ens_DT;
output.ensemble.momentum     = ens_MT;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENSEMBLE AVERAGE PROFILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colors = colormap(jet(length(x)));

ax = figure('Name', 'Integrands');
sgtitle('Ens Avg Integrands')
for i = 1:length(x)
    
    % u         = BL_crop(:, i);
    u         = means.ensemble.u(:, i);
    mask      = ~isnan(u);
    u         = u(mask);
    y_mask    = y(mask);

    marker_size = 10;
    marker_trans = 0.25;
    title_size = 6;
    xlabel_size = 12;

    % Velocity
    subplot(1,3,1)
    hold on
    scatter(u, y_mask, marker_size, 'markerfacecolor', colors(i,:), 'MarkerFaceAlpha', marker_trans, 'MarkerEdgeColor', 'none')
    xlabel('${u}$', 'Interpreter', 'Latex', 'FontSize', xlabel_size)
    ylabel('z [mm]', 'Interpreter', 'latex')
    title('Velocity Profiles', 'FontSize', title_size)
    ylim([0, top_bound_value])

    % Displacement
    subplot(1,3,2)
    hold on
    scatter(1 - (u/u_inf), y_mask, marker_size, 'markerfacecolor', colors(i,:), 'MarkerFaceAlpha', marker_trans, 'MarkerEdgeColor', 'none')
    xlabel('${1 - \frac{u}{u_{\infty}}}$', 'Interpreter', 'Latex', 'FontSize', xlabel_size)
    title('Displacement Thickness', 'FontSize', title_size)
    ylim([0, top_bound_value])
    
    % Momentum
    subplot(1,3,3)
    hold on
    scatter((u/u_inf) .* (1 - (u/u_inf)), y_mask, marker_size, 'markerfacecolor', colors(i,:), 'MarkerFaceAlpha', marker_trans, 'MarkerEdgeColor', 'none')
    xlabel('$\frac{u}{u_{\infty}}\left(1 - \frac{u}{u_{\infty}}\right)$', 'Interpreter', 'Latex', 'FontSize', xlabel_size)
    title('Momentum Thickness', 'FontSize', title_size)
    ylim([0, top_bound_value])
    
end
hold off
hold off
hold off

% Save figure
if save_figures == 1
    fprintf('Saving Figure...\n')
    figure_name = strcat(experiment_name, '_ens_avg_integrands.png');
    exportgraphics(ax, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
    fprintf('Figure Saved\n\n')
    close all; clc;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENSEMBLE AVERAGE PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% BL Thicknesses Contour
ax = figure('Name', 'Ens BL Lengths Contour');
hold on
contourf(X, Y, means.ensemble.u, 500, 'LineStyle', 'none', 'HandleVisibility', 'off')
plot(x, means.ensemble.max_wave_profile, 'Color', 'red', 'LineWidth', 3, 'HandleVisibility', 'off')
yline(0, 'Color', 'black', 'LineWidth', 3, 'HandleVisibility', 'off')
plot(x, ens_BL, 'Color', 'black', 'LineWidth', 2, 'DisplayName', '\delta')
plot(x, ens_DT, 'Color', 'black', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', '\delta^{*}')
plot(x, ens_MT, 'Color', 'black', 'LineStyle', '-.', 'LineWidth', 2, 'DisplayName', '\theta')


C = colorbar();
C.Label.String = '${u}$ [m/s]';
C.Label.Interpreter = 'Latex';

hold off
axis equal
xlim([0, range])
ylim([-10, top_bound_value])
title(strcat(experiment_name, ': Ens Avg BL Parameters'), 'Interpreter', 'none')
xlabel('x [mm]', 'Interpreter', 'Latex')
ylabel('x [mm]', 'Interpreter', 'Latex')
legend()

% Save figure
if save_figures  == 1
    fprintf('Saving Figure...\n')
    figure_name = strcat(experiment_name, '_ens_avg_BL_contour.png');
    exportgraphics(ax, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
    fprintf('Figure Saved\n\n')
    close all; clc;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENSEMBLE AVERAGE SHAPE FACTOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ens_H = ens_DT ./ ens_MT;
ens_H(1) = nan;
ens_H(end) = nan;

ax = figure('Name', 'Ens BL Shape Factor');
hold on
scatter(x, ens_H, 15, 'markerfacecolor', 'red', 'markeredgecolor', 'none')
plot(x, ens_H, 'linewidth', 2, 'color', 'red')
hold off

xlim([0, range])
xlabel('x [mm]', 'Interpreter', 'latex')
ylabel('H', 'Interpreter', 'latex')
title(strcat(experiment_name, {' '}, 'Ens Avg Shape Factor'), 'Interpreter', 'none')

% Save figure
if save_figures  == 1
    fprintf('Saving Figure...\n')
    figure_name = strcat(experiment_name, '_ens_avg_shape_factor.png');
    exportgraphics(ax, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
    fprintf('Figure Saved\n\n')
    close all; clc;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE AVERAGE PARAMETERS/PROFILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C = colormap(jet(length(x)));

if wave_type ~= '0'
    for j = 1:4
    
        % Plot
        ax = figure('Name', strcat('BL Profile Phase ', num2str(j)));
        hold on
        for i = 1:length(x)
            scatter(means.phase(j).u(:,i), y, 'MarkerFaceColor', C(i,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', marker_alpha)
        end
        hold off
    
        xlabel('u [m/s]', 'Interpreter', 'Latex')
        ylabel('z [mm]', 'Interpreter', 'Latex')
        ylim([-20 , top_bound_value])
        title(strcat(experiment_name, ': Phase', {' '}, num2str(j), ' Avg BL Profiles'), 'Interpreter', 'none')
    
        % Save figure
        if save_figures  == 1
            fprintf('Saving Figure...\n')
            figure_name = strcat(experiment_name, '_phase_avg_', num2str(j), '_profiles.png');
            exportgraphics(ax, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
            fprintf('Figure Saved\n\n')
            close all; clc;
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE AVERAGE CONTOURS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phase_BL = zeros(4, length(x));
phase_DT = zeros(4, length(x));
phase_MT = zeros(4, length(x));

if wave_type ~= '0'
    for i = 1:4
        
        % Wave
        ref_wave  = means.phase(i).reference_wave.';
        crop_wave = means.phase(i).max_wave_profile;
    
        % Define boundary layer edge velocity
        BL_crop        = means.phase(i).u;
        BL_crop_smooth = juliaan_smooth(BL_crop, smooth_kernel);
        
        % Remove values above boundary layer
        BL_crop(BL_crop_smooth > BL_edge)        = nan;
        BL_crop_smooth(BL_crop_smooth > BL_edge) = nan;
        
        % BL Thickness
        [~, phase_I] = max(BL_crop_smooth, [], 1, 'omitnan');
        BL           = Y(phase_I);
    
        % Remove ends
        BL(1)   = nan;
        BL(end) = nan;
        phase_BL(i,:) = BL;
        
        % Save phase avg BL
        output.phase(i).thickness = phase_BL(i,:);
        
        % Displacement & Momentum Thickness
        for j = 1:length(x)
            
            u            = BL_crop(:,j);
            mask         = ~isnan(u);
            u            = u(mask);
            y_mask       = y(mask);
            
            DT_integrand = 1 - (u / u_inf);
            MT_integrand = (u / u_inf) .* (1 - (u / u_inf));
    
            phase_DT(i,j) = trapz(flipud(y_mask), flipud(DT_integrand));
            phase_MT(i,j) = trapz(flipud(y_mask), flipud(MT_integrand));
       
        end
    
        %%% First and last points suck
        phase_DT(i,1)   = nan;
        phase_DT(i,end) = nan;
        
        phase_MT(i,1)   = nan;
        phase_MT(i,end) = nan;
    
        %%% Offset
        % phase_DT(i,:) = phase_DT(i,:) + crop_wave;
        % phase_MT(i,:) = phase_MT(i,:) + crop_wave;
        
        % Save
        output.phase(i).displacement = phase_DT(i,:);
        output.phase(i).momentum     = phase_MT(i,:);
        
        % Plot
        ax = figure('Name', strcat('BL Contour Phase ', num2str(i)));
        title(strcat('Phase', {' '}, num2str(i)))
        
        hold on
        contourf(X, Y, means.phase(i).u, 500, 'Linestyle', 'none', 'HandleVisibility', 'off')
        plot(x, ref_wave, 'Color', 'black', 'LineWidth', 3, 'HandleVisibility', 'off')
        plot(x, crop_wave, 'Color', 'red', 'LineWidth', 3, 'HandleVisibility', 'off')
    
        plot(x, phase_BL(i,:), 'Color', 'black', 'LineWidth', 2, 'DisplayName', '\delta')
        plot(x, phase_DT(i,:), 'Color', 'black', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', '\delta^{*}')
        plot(x, phase_MT(i,:), 'Color', 'black', 'LineStyle', '-.', 'LineWidth', 2, 'DisplayName', '\theta')
        hold off
    
        C = colorbar();
        C.Label.String = 'u [m/s]';
        C.Label.Interpreter = 'Latex';
    
        axis equal
        xlim([0+3, range-3])
        ylim([-20, top_bound_value - 3])
    
        xlabel('x [mm]', 'Interpreter', 'Latex')
        ylabel('z [mm]', 'Interpreter', 'Latex')
        legend()
        title(strcat(experiment_name, ' Phase', {' '}, num2str(i)), 'Interpreter', 'none')
        
        % Save figure
        if save_figures == 1
            fprintf('Saving Figure...\n')
            figure_name = strcat(experiment_name, '_phase_avg_', num2str(i), '_contours.png');
            exportgraphics(ax, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
            fprintf('Figure Saved\n\n')
            close all; clc;
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE AVERAGE INTEGRANDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if wave_type ~= '0'
    for j = 1:4
        ax = figure('Name', strcat('Phase', num2str(j), ' Integrands'));
        title_name = strcat('Phase', {' '}, num2str(j), ' Avg Integrands');
        sgtitle(title_name{1})
    
        for i = 1:length(x)
    
            u         = means.phase(j).u(:, i);
            mask      = ~isnan(u);
            u         = u(mask);
            y_mask    = y(mask);
        
            marker_size = 10;
            marker_trans = 0.25;
            title_size = 6;
            xlabel_size = 12;
        
            % Velcity
            subplot(1,3,1)
            hold on
            scatter(u, y_mask, marker_size, 'markerfacecolor', colors(i,:), 'MarkerFaceAlpha', marker_trans, 'MarkerEdgeColor', 'none')
            xlabel('${u}$', 'Interpreter', 'Latex', 'FontSize', xlabel_size)
            ylabel('z [mm]', 'Interpreter', 'latex')
            title('Velocity Profiles', 'FontSize', title_size)
            ylim([0, top_bound_value])
            
            % Displacement
            subplot(1,3,2)
            hold on
            scatter(1 - (u/u_inf), y_mask, marker_size, 'markerfacecolor', colors(i,:), 'MarkerFaceAlpha', marker_trans, 'MarkerEdgeColor', 'none')
            xlabel('${1 - \frac{u}{u_{\infty}}}$', 'Interpreter', 'Latex', 'FontSize', xlabel_size)
            title('Displacement Thickness', 'FontSize', title_size)
            ylim([0, top_bound_value])
            
            % Momentum
            subplot(1,3,3)
            hold on
            scatter((u/u_inf) .* (1 - (u/u_inf)), y_mask, marker_size, 'markerfacecolor', colors(i,:), 'MarkerFaceAlpha', marker_trans, 'MarkerEdgeColor', 'none')
            xlabel('$\frac{u}{u_{\infty}}\left(1 - \frac{u}{u_{\infty}}\right)$', 'Interpreter', 'Latex', 'FontSize', xlabel_size)
            title('Momentum Thickness', 'FontSize', title_size)
            ylim([0, top_bound_value])
            
        end
        hold off
        hold off
        hold off
        
        % Save figure
        if save_figures  == 1
            fprintf('Saving Figure...\n')
            figure_name = strcat(experiment_name, '_phase_avg_', num2str(j), '_integrands.png');
            exportgraphics(ax, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
            fprintf('Figure Saved\n\n')
            close all; clc;
        end
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE AVERAGE PROFILES PEAKS VS TROUGHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for j = 1:4
%     % Wave
%     ref_wave = means.phase(j).reference_wave.';
% 
%     % Find peaks/troughs
%     [x_peaks_troughs_value, peaks_troughs] = findpeaks(abs(ref));
% 
%     figure('Name', strcat('BL Integrand Phase ', num2str(j)));
%     for i = 1:length(peaks_troughs)
% 
%         u = means.phase(j).u(:, peaks_troughs(i));
% 
%         DT_integrand = 1 - (u / u_inf);
%         MT_integrand = (u / u_inf) .* (1 - (u / u_inf));
% 
%         hold on
%         hold on
%         if sign(ref(peaks_troughs(i))) == 1
%             subplot(1,2,1)
% 
%             %plot(u, y_mask, 'Color', colors(i,:), 'LineWidth', 2)
%             hold on
%             %area(DT_integrand, y_mask, min(y_mask), 'FaceColor', 'green')
%             plot(MT_integrand, y, 'LineWidth', 4)
%             hold off
%             xlabel('${1 - \frac{u}{u_{\infty}}}$', 'Interpreter', 'Latex', 'FontSize', 18)
%             ylabel('y [mm]', 'Interpreter', 'Latex')
%             ylim([-10, 80])
%             xlim([0, 0.6])
%             yline(wave_amplitude, 'Color', 'red', 'LineWidth', 2)
%             title('Peaks')
% 
%         elseif sign(ref(peaks_troughs(i))) == -1
%             subplot(1,2,2)
%             %plot(u, y_mask, 'Color', colors(i,:), 'LineWidth', 2)
%             hold on
%             %area(DT_integrand, y_mask, min(y_mask), 'FaceColor', 'green')
%             plot(MT_integrand, y, 'LineWidth', 4)
%             hold off
%             xlabel('${1 - \frac{u}{u_{\infty}}}$', 'Interpreter', 'Latex', 'FontSize', 18)
%             title(strcat('Phase', {' '}, num2str(j)))
%             ylim([-10, 80])
%             xlim([0, 0.6])
%             yline(-wave_amplitude, 'Color', 'red', 'LineWidth', 2)
%             title('Troughs')
%         end
%     end
% 
%     hold off
%     hold off
%     sgtitle(strcat(experiment_name, ' Phase', {' '},  num2str(j), ': Normalized Velocity Deficit'), 'Interpreter', 'none')
% 
% %     figure_name = strcat(experiment_name, '_phaseavg_', num2str(j), '_DT_integrand_filled.png');
% %     exportgraphics(ax, strcat(figure_folder, '/', figure_name))
% end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE AVERAGE BOUNDARY LAYER AND DISPLACEMENT/MOMENTUM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if wave_type ~= '0'

    %%% Boundary Layer Thickness
    ax1 = figure('Name', 'Phase Avg Boundary Layer Thickness');
    sgtitle('Phase Average Boundary Layer Thickness')
    hold on
    plot(x, ens_BL, 'Color', 'black', 'Linestyle', '--', 'LineWidth', 2, 'DisplayName', 'Ensemble')
    for i = 1:4
        ref_wave  = means.phase(i).reference_wave.';
        crop_wave = means.phase(i).max_wave_profile;
    
        name = strcat('Phase', {' '}, num2str(i));
        plot(x, phase_BL(i,:), 'LineWidth', 2, 'DisplayName', name{1})
    end
    hold off
    axis equal
    xlabel('x [mm]', 'Interpreter', 'Latex')
    ylabel('Boundary Layer Thickness $(\delta)$ [mm]', 'Interpreter', 'Latex')
    xlim([0, range])
    legend('Location', 'northoutside', 'Orientation', 'horizontal', 'fontsize', 6)
    
    % Save Figure
    if save_figures == 1
        fprintf('Saving Figure...\n')
        figure_name = strcat(experiment_name, '_phase_avg_BL_lines.png');
        exportgraphics(ax1, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
        fprintf('Figure Saved\n\n')
        close all; clc;
    end
    


    %%% Displacement Thickness
    ax2 = figure('Name', 'Phase Avg Displacement Thickness');
    sgtitle('Phase Average Displacement Thickness')
    hold on
    plot(x, ens_DT, 'Color', 'black', 'Linestyle', '--', 'LineWidth', 2, 'DisplayName', 'Ensemble')
    % plot(x, ens_DT ./ ens_BL, 'Color', 'black', 'Linestyle', '--', 'LineWidth', 2, 'DisplayName', 'Ensemble')
    for i = 1:4
        ref_wave  = means.phase(i).reference_wave.';
        crop_wave = means.phase(i).max_wave_profile;
    
        name = strcat('Phase', {' '}, num2str(i));
        plot(x, phase_DT(i,:), 'LineWidth', 2, 'DisplayName', name{1})
        % plot(x, phase_DT(i,:) ./ phase_BL(i,:), 'LineWidth', 2, 'DisplayName', name{1})
    end
    hold off
    
    xlabel('x [mm]', 'Interpreter', 'Latex')
    ylabel('$\delta^{*}$ [mm]', 'Interpreter', 'Latex')
    xlim([0, range])
    legend('Location', 'northoutside', 'Orientation', 'horizontal', 'fontsize', 6)
    
    % Save Figure
    if save_figures == 1
        fprintf('Saving Figure...\n')
        figure_name = strcat(experiment_name, '_phase_avg_DT_lines.png');
        exportgraphics(ax2, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
        fprintf('Figure Saved\n\n')
        close all; clc;
    end


    
    %%% Momentum Thickness
    ax3 = figure('Name', 'Phase Avg Momentum Thickness');
    sgtitle('Phase Average Momentum Thickness')
    hold on
    plot(x, ens_MT, 'Color', 'black', 'Linestyle', '--', 'LineWidth', 2, 'DisplayName', 'Ensemble')
    % plot(x, ens_MT ./ ens_BL, 'Color', 'black', 'Linestyle', '--', 'LineWidth', 2, 'DisplayName', 'Ensemble')
    for i = 1:4
        ref_wave  = means.phase(i).reference_wave.';
        crop_wave = means.phase(i).max_wave_profile;
    
        name = strcat('Phase', {' '}, num2str(i));
        plot(x, phase_MT(i,:), 'LineWidth', 2, 'DisplayName', name{1})
        % plot(x, phase_MT(i,:) ./ phase_BL(i,:), 'LineWidth', 2, 'DisplayName', name{1})
    end
    hold off
    
    xlabel('x [mm]', 'Interpreter', 'Latex')
    ylabel('$\theta$ [mm]', 'Interpreter', 'Latex')
    xlim([0, range])
    legend('Location', 'northoutside', 'Orientation', 'horizontal', 'fontsize', 6)
    
    % Save Figure
    if save_figures == 1
        fprintf('Saving Figure...\n')
        figure_name = strcat(experiment_name, '_phase_avg_MT_lines.png');
        exportgraphics(ax3, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
        fprintf('Figure Saved\n\n')
        close all; clc;
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE AVERAGE DISPLACEMENT/MOMENTUM NORMALIZED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if wave_type ~= '0'
    
    %%% Displacement Thickness
    ax2 = figure('Name', 'Normalized Phase Avg Displacement Thickness');
    sgtitle('Normalized Phase Average Displacement Thickness')
    hold on
    plot(x, ens_DT ./ ens_BL, 'Color', 'black', 'Linestyle', '--', 'LineWidth', 2, 'DisplayName', 'Ensemble')
    for i = 1:4
        ref_wave  = means.phase(i).reference_wave.';
        crop_wave = means.phase(i).max_wave_profile;
    
        name = strcat('Phase', {' '}, num2str(i)); 
        plot(x, phase_DT(i,:) ./ phase_BL(i,:), 'LineWidth', 2, 'DisplayName', name{1})
    end
    hold off
    
    xlabel('x [mm]', 'Interpreter', 'Latex')
    ylabel('$\delta^{*} / \delta$', 'Interpreter', 'Latex')
    xlim([0, range])
    legend('Location', 'northoutside', 'Orientation', 'horizontal', 'fontsize', 6)
    
    % Save Figure
    if save_figures == 1
        fprintf('Saving Figure...\n')
        figure_name = strcat(experiment_name, '_phase_avg_DT_lines_norm.png');
        exportgraphics(ax2, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
        fprintf('Figure Saved\n\n')
        close all; clc;
    end

    
    %%% Momentum Thickness
    ax3 = figure('Name', 'Normalized Phase Avg Momentum Thickness');
    sgtitle('Normalized Phase Average Momentum Thickness')
    hold on
    plot(x, ens_MT ./ ens_BL, 'Color', 'black', 'Linestyle', '--', 'LineWidth', 2, 'DisplayName', 'Ensemble')
    for i = 1:4
        ref_wave  = means.phase(i).reference_wave.';
        crop_wave = means.phase(i).max_wave_profile;
    
        name = strcat('Phase', {' '}, num2str(i));
        plot(x, phase_MT(i,:) ./ phase_BL(i,:), 'LineWidth', 2, 'DisplayName', name{1})
    end
    hold off
    
    xlabel('x [mm]', 'Interpreter', 'Latex')
    ylabel('$\theta / \delta$', 'Interpreter', 'Latex')
    xlim([0, range])
    legend('Location', 'northoutside', 'Orientation', 'horizontal', 'fontsize', 6)
    
    % Save Figure
    if save_figures == 1
        fprintf('Saving Figure...\n')
        figure_name = strcat(experiment_name, '_phase_avg_MT_lines_norm.png');
        exportgraphics(ax3, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
        fprintf('Figure Saved\n\n')
        close all; clc;
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE AVERAGE SHAPE FACTOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

new_cos_fit = @(b, v) b(1) * cos(2 * pi * (v - (range/2) - b(3)) / wave_length) + b(2);

if wave_type ~= '0'
    figure('Name', 'Phase Avg Shape Factor');
    t = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(t, strcat(experiment_name, ': Phase Average Shape Factor'), 'Interpreter', 'none')
    
    for i = 1:4
        
        h(i)    = nexttile;
        phase_H = phase_DT(i,:) ./ phase_MT(i,:);
        ref     = new_cos_fit([(max(phase_H) - min(phase_H)) / 2, mean(phase_H, 'omitnan'), phase_offset(i)], x);
        
        hold on
        % Reference wave position
        p          = plot(x, ref, 'Color', 'red', 'LineWidth', 2,'DisplayName', 'Wave');
        p.Color(4) = 0.25;
        
        % Ensemble shape factor
        ens_H_plot = plot(x,ens_H, 'Color', 'black', 'LineWidth', 2, 'Displayname', 'Ensemble');
        ens_H_plot.Color(4) = 0.25;
    
        % Phase avg shape factor
        plot(x, phase_H, 'Color', 'blue', 'LineWidth', 2, 'DisplayName', 'Shape Factor')
        hold off
    
        xlabel('x [mm]', 'interpreter', 'latex')
        ylabel('H', 'Interpreter', 'Latex')
        title(strcat('Phase', {' '}, num2str(i)));
        xlim([0, range])
      
    end
    
    linkaxes(h, 'xy')
    lh = legend({'Ensemble', 'Wave', 'Shape Factor'}, 'Interpreter', 'none', 'Location', 'NorthOutside', 'Orientation', 'Horizontal');
    lh.Layout.Tile = 'North';
    
    % Save Figure
    if save_figures == 1
        fprintf('Saving Figure...\n')
        figure_name = strcat(experiment_name, '_phase_avg_shape_factor.png');
        exportgraphics(t, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
        fprintf('Figure Saved\n\n')
        close all; clc;
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORRELATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colors = {'red', 'yellow', 'green', 'blue'};

ax = figure('Name', 'Castillo Correlation');
sgtitle(strcat(experiment_name, {' '}, 'Castillo Correlation'), 'Interpreter', 'none')
hold on
if wave_type ~= '0'
    for i = 1:4
        H      = (phase_DT(i,:)) ./ (phase_MT(i,:));
        DT_BLT = (phase_DT(i,:)) ./ (phase_BL(i,:));
        scatter(DT_BLT, H, 'filled', 'MarkerFaceColor', colors{i}, 'MarkerFaceAlpha', marker_alpha)
    end

    scatter(ens_DT ./ ens_BL, ens_H, 'filled', 'MarkerFaceColor', 'black', 'MarkerFaceAlpha', marker_alpha)
    hold off
    legend({'Phase 1', 'Phase 2', 'Phase 3', 'Phase 4', 'Ensemble'}, 'Location', 'northoutside', 'orientation', 'horizontal', 'fontsize', 6)

else
    scatter(ens_DT ./ ens_BL, ens_H, 'filled', 'MarkerFaceColor', 'black', 'MarkerFaceAlpha', marker_alpha)
    hold off
    legend({'Ensemble'}, 'Location', 'northoutside', 'orientation', 'horizontal', 'fontsize', 6)
end

xlabel('$\delta^* / \delta$', 'Interpreter', 'latex')
ylabel('H', 'Interpreter', 'latex')

% Save Figure
if save_figures == 1
    fprintf('Saving Figure...\n')
    figure_name = strcat(experiment_name, '_castillo_correlation.png');
    exportgraphics(ax, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
    fprintf('Figure Saved\n\n')
    close all; clc;
end


%% Save 

fprintf('Saving BL mat file...\n')
save(BL_file, 'output');
fprintf('BL mat file saved\n\n')







