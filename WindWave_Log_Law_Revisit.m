%% 2D Offshore PIV Processing: Ondrej Fercak, Zein Sadek, 1/2023
%%% Compute Log Law 
% Zein Sadek
% 11/2023

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions/');
% experiment_log  = readcell('Offshore_Inflow_Log.xlsx');
wave_parameters = readcell('Offshore_Waves.xlsx');

% s = settings;
% s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data
experiment_name = 'WT4_WV0_AGP';

% Experiment Specifics
tunnel_freq = experiment_name(strfind(experiment_name, 'WT') + 2);
wave_type   = experiment_name(strfind(experiment_name, 'WV') + 2);

% Save paths
results_path  = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/';
mtlb_file     = strcat(results_path, 'data'   , '/', experiment_name, '_DATA.mat');
mean_file     = strcat(results_path, 'means'  , '/', experiment_name, '_MEANS.mat');
phase_file    = strcat(results_path, 'phase'  , '/', experiment_name, '_PHASE.mat');
BL_file       = strcat(results_path, 'boundary_layer/', experiment_name, '_BL.mat');

if wave_type ~= '0'
    wave_file = strcat(results_path, 'wave'   , '/', experiment_name, '_WAVE.mat');
end

% Save for Log Law Values
LL_file       = strcat(results_path, 'log_law', '/', experiment_name, '_LL.mat');
figure_folder = strcat(results_path, 'figures/', experiment_name, '/log_law/');

if ~exist(figure_folder, 'dir')
    mkdir(figure_folder)
end

% Filtering Parameters
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
% IMPORT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if wave_type ~= '0'
%     waves = load(wave_file);
%     waves = waves.output;
% end

means = load(mean_file); 
means = means.output;

BL_data = load(BL_file);
BL_data = BL_data.output;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PICK LOCATION FOR SLICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coordinates
Y = means.Y;
y = Y(:,1);
X = means.X;
x = X(1,:);
range = abs(left_bound_value) + abs(right_bound_value);

% Selct center of phase avg waves
profile_x_value  = range/2;
[~, profile_idx] = min(abs(x - profile_x_value));
profile_buffer   = 20;
u_profile        = means.ensemble.u(:, profile_idx - profile_buffer:profile_idx + profile_buffer);
u_profile        = mean(u_profile, 2, 'omitnan');
u_profile(y < 0) = nan;

% Save Values
output.profile_idx        = profile_idx;
output.profile_x_value    = x(profile_idx);
output.ensemble.u_profile = u_profile;

figure()
plot(u_profile, y)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTEGRATED BOUNDARY LAYER EQUATION TERMS (Ensemble)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
nu = 1.48E-5;                   % [m^s / s]
dy = ((y(3) - y(2)) * 10^(-3)); % [m]
dx = ((x(3) - x(2)) * 10^(-3)); % [m]

%%% Viscous Shear Stress Term
dudy   = gradient(u_profile, dy);
term_1 = nu * dudy;

% Save Values
output.ensemble.IBL.visc = term_1;

% Reynolds Shear Stress
uv_profile = -1 * means.ensemble.uv(:, profile_idx - profile_buffer:profile_idx + profile_buffer);
uv_profile = mean(uv_profile, 2, 'omitnan');
uv_profile(y < 0) = nan;
term_2     = uv_profile;

% Save Values
output.ensemble.IBL.stress = term_2;

figure()
plot(term_2, y)

%%

%%% Mean Momentum Flux
% Compute integrad
[du2dx, ~] = gradient(means.ensemble.u.^2, dx);
du2dx_profile = du2dx(:, profile_idx - profile_buffer:profile_idx + profile_buffer);
du2dx_profile = mean(du2dx_profile, 2, 'omitnan');

% Create mask to remove nans
du2dx_profile_mask = ~isnan(du2dx_profile);
du2dx_profile      = du2dx_profile(du2dx_profile_mask);
y_mask             = y(du2dx_profile_mask) * 1E-3;

% Compute integral
int_du2dx_dy = -1 * flipud(cumtrapz(flipud(y_mask), flipud(du2dx_profile)));

% Second half of term: w/ continuity substitution
% umean_vmean = means.ensemble.u .* means.ensemble.v;
% umean_vmean = umean_vmean(:, profile_idx - profile_buffer:profile_idx + profile_buffer);
% umean_vmean = mean(umean_vmean, 2, 'omitnan');
% umean_vmean = -1 * umean_vmean;

% Second half of term: as shown
[dudx, ~] = gradient(means.ensemble.u, dx);
dudx_profile = dudx(:, profile_idx - profile_buffer:profile_idx + profile_buffer);
dudx_profile = mean(dudx_profile, 2, 'omitnan');

% Create mask to remove nans
dudx_profile_mask = ~isnan(dudx_profile);
dudx_profile      = dudx_profile(dudx_profile_mask);
y_mask            = y(dudx_profile_mask) * 1E-3;

% Compute integral
int_dudx_dy = flipud(cumtrapz(flipud(y_mask), flipud(dudx_profile)));

% Save term
% term_3 = int_du2dx_dy + umean_vmean(du2dx_profile_mask);
term_3 = int_du2dx_dy + (int_dudx_dy .* u_profile(dudx_profile_mask));
output.ensemble.IBL.mean_mom_y    = y_mask;
output.ensemble.IBL.mean_mom_mask = dudx_profile_mask;
output.ensemble.IBL.mean_mom      = term_3;

figure()
hold on
plot(term_3, y_mask)
plot(term_2, y * 1E-3)
hold on

%%
%%% Turbulent Momentum Flux
% Compute integrad
[duudx, ~] = gradient(means.ensemble.uu, dx);

% Pick Slice
duudx_profile = duudx(:, profile_idx - profile_buffer:profile_idx + profile_buffer);
duudx_profile = mean(duudx_profile, 2, 'omitnan');

% Create mask to remove nans
duudx_profile_mask = ~isnan(duudx_profile);
duudx_profile      = duudx_profile(duudx_profile_mask);
y_mask             = y(duudx_profile_mask) * 1E-3;

% Compute integrals
int_duudx_dy = -1 * flipud(cumtrapz(flipud(y_mask), flipud(duudx_profile)));

% Save
term_4 = int_duudx_dy;
output.ensemble.IBL.turb_mom_y    = y_mask;
output.ensemble.IBL.turb_mom_mask = duudx_profile_mask;
output.ensemble.IBL.turb_mom      = term_4;

figure()
hold on
plot(term_1, y * 1E-3)
plot(term_2, y * 1E-3)
plot(term_3, y_mask)
plot(term_4, y_mask)
hold off

%%
%%% Pressure Correction
% Compute integrad
[dvvdx, ~] = gradient(means.ensemble.vv, dx);

% Pick Slice
dvvdx_profile = dvvdx(:, profile_idx - profile_buffer:profile_idx + profile_buffer);
dvvdx_profile = mean(dvvdx_profile, 2, 'omitnan');

% Create mask to remove nans
dvvdx_profile_mask = ~isnan(dvvdx_profile);
dvvdx_profile      = dvvdx_profile(dvvdx_profile_mask);
y_mask             = y(dvvdx_profile_mask) * 1E-3;

% Compute integrals
int_dvvdx_dy = flipud(cumtrapz(flipud(y_mask), flipud(dvvdx_profile)));

% Save
term_5 = int_dvvdx_dy;
output.ensemble.IBL.pressure_y    = y_mask;
output.ensemble.IBL.pressure_mask = dvvdx_profile_mask;
output.ensemble.IBL.pressure      = term_5;

figure()
hold on
plot(term_1, y * 1E-3)
plot(term_2, y * 1E-3)
plot(term_3, y_mask)
plot(term_4, y_mask)
plot(term_5, y_mask)
hold off
ylim([0, 0.05])

%%
%%% Pressure Gradient
u_inf = BL_data.u_inf;
u_masked = imgaussfilt(means.ensemble.u(10,:), 5);
% u_masked(Y < BL_data.ensemble.thickness) = nan;
du_infdx = gradient(u_masked, dx);
du_infdx_avg = mean(du_infdx(:, profile_idx - profile_buffer:profile_idx + profile_buffer), 2);

term_6 = u_inf .* du_infdx_avg .* (y * 1E-3);
output.ensemble.IBL.pressure_grad = term_6;

figure()
plot(x, u_masked)
xlim([0, 200])
ylim([0, 1.1*u_inf])
yline(u_inf)

figure()
hold on
plot(x, du_infdx)
xline(x(profile_idx))
hold off

figure()
plot(term_6, y)

%% Comparing Terms

linewidth = 2;
marker_size = 10;
colors = turbo(linspace(0,6,1));

ax = figure('Name', 'IBL Terms');
title(strcat(experiment_name, ': IBL Terms'), 'Interpreter', 'none')

hold on
plot(term_1, y, 'color', colors(1,:), 'LineWidth', linewidth, 'DisplayName', 'Viscous Stress')
scatter(term_1, y, marker_size, 'markerfacecolor', colors(1,:),'MarkerEdgeColor', 'none', 'HandleVisibility', 'off')

plot(term_2, y, 'color', colors(2,:), 'LineWidth', linewidth, 'DisplayName', 'Reynolds Stress')
scatter(term_2, y, marker_size, 'markerfacecolor', colors(2,:),'MarkerEdgeColor', 'none', 'HandleVisibility', 'off')

plot(term_3, y(du2dx_profile_mask), 'color', colors(3,:), 'LineWidth', linewidth, 'DisplayName', 'Mean Momentum Flux')
scatter(term_3, y(du2dx_profile_mask), marker_size, 'markerfacecolor', colors(3,:),'MarkerEdgeColor', 'none', 'HandleVisibility', 'off')

plot(term_4, y(duudx_profile_mask), 'color', colors(4,:), 'LineWidth', linewidth, 'DisplayName', 'Turbuelent Flux')
scatter(term_4, y(duudx_profile_mask), marker_size, 'markerfacecolor', colors(4,:),'MarkerEdgeColor', 'none', 'HandleVisibility', 'off')

plot(term_5, y(dvvdx_profile_mask), 'color', colors(5,:), 'LineWidth', linewidth, 'DisplayName', 'Pressure Correction')
scatter(term_5, y(dvvdx_profile_mask), marker_size, 'markerfacecolor', colors(5,:),'MarkerEdgeColor', 'none', 'HandleVisibility', 'off')

% plot(term_6, y, 'LineWidth', 2, 'DisplayName', 'Pressure Gradient')
hold off

ylim([0, 80])
% set(gca,'yscale','log')
ylabel('y [mm]', 'Interpreter', 'latex')
xlabel('$m^2 / s^2$', 'Interpreter', 'latex')
legend()

if save_figures == 1
    fprintf('Saving Figure...\n')
    figure_name = strcat(experiment_name, '_IBL_terms.png');
    exportgraphics(ax, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
    fprintf('Figure Saved\n\n')
    close all; clc;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FRICTION VELOCITY [UV]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ens_u_star_profile           = sqrt(uv_profile);
[ens_u_star, ens_u_star_idx] = max(uv_profile);
ens_u_star                   = sqrt(ens_u_star);

% Save Values
output.ensemble.u_star.profile = ens_u_star_profile;
output.ensemble.u_star.u_star  = ens_u_star;
output.ensemble.u_star.idx     = ens_u_star_idx;

ax = figure('Name', 'u* profile');
title(strcat(experiment_name, {' '}, 'Friction Velocity:', {' '}, num2str(round(ens_u_star, 3)), ' m/s'), 'interpreter', 'none')
hold on
plot(ens_u_star_profile, y, 'Color', 'black')
yline(0, 'Color', 'red', 'LineWidth', 3);
scatter(ens_u_star, y(ens_u_star_idx), 'MarkerFaceColor', 'red')
hold off
ylim([-10, 210])
ylabel('y [mm]', 'Interpreter', 'Latex')
xlabel('$u^* = \sqrt{\frac{\tau_w}{\rho}}$', 'Interpreter', 'Latex', 'FontSize', 12)

if save_figures == 1
    fprintf('Saving Figure...\n')
    figure_name = strcat(experiment_name, '_ens_avg_u_star.png');
    exportgraphics(ax, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
    fprintf('Figure Saved\n\n')
    close all; clc;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOG LAW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ens_u_plus = u_profile / ens_u_star;
ens_y_plus = (y * (10^-3) * ens_u_star) / nu;

% Save Values
output.ensemble.u_plus = ens_u_plus;
output.ensemble.y_plus = ens_y_plus;

ax = figure('Name', 'Ensemble Log Law');
title(strcat(experiment_name, ': Ens Log Law'), 'Interpreter', 'none')
hold on
plot(ens_y_plus, ens_u_plus, 'color', 'red')
scatter(ens_y_plus, ens_u_plus, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'none')
hold off
set(gca,'xscale','log')
xlabel('$y^+$', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('$u^+$', 'Interpreter', 'Latex', 'FontSize', 18)
grid on

if save_figures == 1
    fprintf('Saving Figure...\n')
    figure_name = strcat(experiment_name, '_ens_avg_log_law.png');
    exportgraphics(ax, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
    fprintf('Figure Saved\n\n')
    close all; clc;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE AVERAGE FRICTION VELOCITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Phases
phase_offset = [0, wave_length/4, wave_length/2, 3*wave_length/4];

% Save values for each phase
phase_u_profiles = zeros(length(y), length(phase_offset));
phase_uv_profiles         = zeros(length(y), length(phase_offset));

phase_u_star_profile = zeros(length(y), length(phase_offset));
phase_u_star         = zeros(1, length(phase_offset));
phase_u_star_idx     = zeros(1, length(phase_offset));

phase_u_plus = zeros(length(y), length(phase_offset));
phase_y_plus = zeros(length(y), length(phase_offset));

phase_ens_u_plus = zeros(length(y), length(phase_offset));
phase_ens_y_plus = zeros(length(y), length(phase_offset));

if wave_type ~= '0'
    % Look through all phases
    for i = 1:length(phase_offset)
    
        %%% u Profiles
        phase_u_profile = means.phase(i).u(:, profile_idx - profile_buffer:profile_idx + profile_buffer);
        phase_u_profile = mean(phase_u_profile, 2, 'omitnan');
    
        % Save
        phase_u_profiles(:,i)     = phase_u_profile;
        output.phase(i).u_profile = phase_u_profile;
        
        
        %%% uv Profiles
        phase_uv_profile = -1 * means.phase(i).uv(:, profile_idx - profile_buffer:profile_idx + profile_buffer);
        phase_uv_profile = mean(phase_uv_profile, 2, 'omitnan');
           
        % Save Values
        phase_uv_profiles(:,i)     = phase_uv_profile;
        output.phase(i).uv_profile = phase_uv_profile;
        
        
        %%% Friction Velocity
        phase_u_star_profile(:,i) = sqrt(phase_uv_profile);
        [u_star, u_star_idx]  = max(phase_uv_profile);
        u_star                = sqrt(u_star);
    
        phase_u_star(1,i)     = u_star;
        phase_u_star_idx(1,i) = u_star_idx;
        
        % Save Values
        output.phase(i).u_star_profile = sqrt(phase_uv_profile);
        output.phase(i).u_star         = u_star;
        output.phase(i).u_star_idx     = u_star_idx;
        
        
        %%% Ensemble: Log Law
        phase_ens_u_plus(:,i) = phase_u_profile / ens_u_star;
        phase_ens_y_plus(:,i) = (y * (10^-3) * ens_u_star) / nu;
        
        % Phase: Log Law
        phase_u_plus(:,i) = phase_u_profile / u_star;
        phase_y_plus(:,i) = (y * (10^-3) * u_star) / nu;
        
        % Save Values
        output.phase(i).u_plus = phase_u_plus(:,i);
        output.phase(i).y_plus = phase_y_plus(:,i);
        
        output.phase(i).ens_u_plus = phase_ens_u_plus(:,i);
        output.phase(i).ens_y_plus = phase_ens_y_plus(:,i);
        
    end
end
    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE AVERAGE PROFILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if wave_type ~= '0'
    %%% Phase Average Velocity Profiles
    marker_size = 20;
    linewidth   = 2;
    colors = {'red', 'magenta', 'green', 'blue'};
    
    ax = figure('Name', 'Phase Avg Fixed Location Velocity Profiles');
    
    hold on
    plot(u_profile, y, 'Color', 'black', 'LineWidth', linewidth, 'DisplayName', 'Ensemble')
    scatter(u_profile, y, marker_size, 'markerfacecolor', 'black', 'markeredgecolor', 'none', 'HandleVisibility', 'off')
    
    for i = 1:4
        label = strcat('Phase', {' '}, num2str(i));
        plot(phase_u_profiles(:,i), y, 'color', colors{i}, 'LineWidth', linewidth, 'DisplayName', label{1})
        scatter(phase_u_profiles(:,i), y, marker_size, 'markerfacecolor', colors{i}, 'markeredgecolor', 'none', 'HandleVisibility', 'off')
    end
    hold off
    
    title(strcat(experiment_name, ': Phase Avg Velocity Profiles'), 'Interpreter', 'none')
    xlabel('u [m/s]', 'Interpreter', 'Latex')
    legend('Location', 'NorthWest')
    ylim([-2 * wave_amplitude, top_bound_value])
    
    % Save figure
    if save_figures == 1
        fprintf('Saving Figure...\n')
        figure_name = strcat(experiment_name, '_phase_avg_profiles.png');
        exportgraphics(ax, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
        fprintf('Figure Saved\n\n')
        close all; clc;
    end
end


%% Phase Average Stress Profiles

if wave_type ~= '0'
    marker_size = 10;
    linewidth   = 1;
    colors = {'red', 'magenta', 'green', 'blue'};
    
    ax = figure('Name', 'Phase Avg Fixed Location Stress Profiles');
    
    hold on
    plot(uv_profile, y, 'Color', 'black', 'LineWidth', linewidth, 'DisplayName', 'Ensemble')
    scatter(uv_profile, y, marker_size, 'markerfacecolor', 'black', 'markeredgecolor', 'none', 'HandleVisibility', 'off')
    
    for i = 1:4
        label = strcat('Phase', {' '}, num2str(i));
        plot(phase_uv_profiles(:,i), y, 'color', colors{i}, 'LineWidth', linewidth, 'DisplayName', label{1})
        scatter(phase_uv_profiles(:,i), y, marker_size, 'markerfacecolor', colors{i}, 'markeredgecolor', 'none', 'HandleVisibility', 'off')
    end
    hold off
    
    title(strcat(experiment_name, ': Phase Avg uv Profiles'), 'Interpreter', 'none')
    xlabel('$\bar{uv}$ [m/s]', 'Interpreter', 'Latex')
    legend('Location', 'NorthWest')
    ylim([-2 * wave_amplitude, top_bound_value])
    
    % Save figure
    if save_figures == 1
        fprintf('Saving Figure...\n')
        figure_name = strcat(experiment_name, '_phase_avg_uv_profiles.png');
        exportgraphics(ax, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
        fprintf('Figure Saved\n\n')
        close all; clc;
    end
end


%% Phase Avg Friction Velocity

if wave_type ~= '0'
    marker_size = 20;
    linewidth   = 2;
    colors = {'red', 'magenta', 'green', 'blue'};
    
    ax = figure('Name', 'Phase Avg Friction Velocity');
    
    hold on
    for i = 1:4
        label = strcat('Phase', {' '}, num2str(i));
        plot(phase_u_star_profile(:,i), y, 'Color', colors{i}, 'LineWidth', linewidth, 'DisplayName', label{1})
        scatter(phase_u_star_profile(:,i), y, marker_size, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'none', 'HandleVisibility', 'off')
    end
    
    plot(ens_u_star_profile, y, 'Color', 'black', 'LineWidth', linewidth, 'DisplayName', 'Ensemble')
    scatter(ens_u_star_profile, y, marker_size, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'none', 'HandleVisibility', 'off')
    
    hold off
    title(strcat(experiment_name, ': Phase Avg Friction Velocity'), 'Interpreter', 'none')
    ylim([-2*wave_amplitude, top_bound_value])
    xlabel('$u^*$ [m/s]', 'Interpreter', 'Latex')
    ylabel('y [mm]', 'Interpreter', 'Latex')
    legend()
    
    % Save figure
    if save_figures == 1
        fprintf('Saving Figure...\n')
        figure_name = strcat(experiment_name, '_phase_avg_friction_profiles.png');
        exportgraphics(ax, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
        fprintf('Figure Saved\n\n')
        close all; clc;
    end
end


%% Ensemble Avg Log Law

if wave_type ~= '0'
    marker_size = 10;
    linewidth   = 1;
    colors = {'red', 'magenta', 'green', 'blue'};
    
    ax = figure('Name', 'Ensemble Avg Log Law');
    title(strcat(experiment_name, ': Ensemble Avg Log Law'), 'Interpreter', 'none')
    hold on
    for i = 1:4
        label = strcat('Phase', {' '}, num2str(i));
        plot(phase_ens_y_plus(:,i), phase_ens_u_plus(:,i), 'Color', colors{i}, 'linewidth', linewidth,'HandleVisibility', 'off')
        scatter(phase_ens_y_plus(:,i), phase_ens_u_plus(:,i), marker_size, 'MarkerFaceColor', colors{i}, 'DisplayName', label{1}, 'MarkerEdgeColor', 'none')   
    end
    
    plot(ens_y_plus, ens_u_plus, 'Color', 'black', 'linewidth', linewidth, 'HandleVisibility', 'off')
    scatter(ens_y_plus, ens_u_plus, marker_size, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'none', 'DisplayName', 'Ensemble')
    hold off
    
    grid on
    legend('Location', 'NorthWest')
    set(gca,'xscale','log')
    ylabel('$u^+$', 'Interpreter', 'Latex')
    xlabel('$y^+$', 'Interpreter', 'Latex')
    
    % Save figure
    if save_figures == 1
        fprintf('Saving Figure...\n')
        figure_name = strcat(experiment_name, '_phase_avg_ens_avg_log_law.png');
        exportgraphics(ax, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
        fprintf('Figure Saved\n\n')
        close all; clc;
    end
end


%% Phase Avg Log Law

if wave_type ~= '0'
    marker_size = 10;
    linewidth   = 1;
    colors = {'red', 'magenta', 'green', 'blue'};
    
    ax = figure('Name', 'Phase Avg Log Law');
    title(strcat(experiment_name, ': Phase Avg Log Law'), 'Interpreter', 'none')
    hold on
    for i = 1:4
        label = strcat('Phase', {' '}, num2str(i));
        plot(phase_y_plus(:,i), phase_u_plus(:,i), 'Color', colors{i}, 'linewidth', linewidth, 'HandleVisibility', 'off')
        scatter(phase_y_plus(:,i), phase_u_plus(:,i), marker_size, 'MarkerFaceColor', colors{i}, 'DisplayName', label{1}, 'MarkerEdgeColor', 'none')   
    end
    
    plot(ens_y_plus, ens_u_plus, 'Color', 'black', 'linewidth', linewidth, 'HandleVisibility', 'off')
    scatter(ens_y_plus, ens_u_plus, marker_size, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'none', 'DisplayName', 'Ensemble')
    hold off
    
    grid on
    legend('Location', 'NorthWest')
    set(gca,'xscale','log')
    ylabel('$u^+$', 'Interpreter', 'Latex')
    xlabel('$y^+$', 'Interpreter', 'Latex')
    % xlim([5E0, 5E3])
    % ylim([0, 25])
    
    % Save figure
    if save_figures == 1
        fprintf('Saving Figure...\n')
        figure_name = strcat(experiment_name, '_phase_avg_log_law.png');
        exportgraphics(ax, strcat(figure_folder, '/', figure_name), 'Resolution', figure_resolution)
        fprintf('Figure Saved\n\n')
        close all; clc;
    end
end


%% Save Values to Mat File

fprintf('Saving LL mat file...\n')
save(LL_file, 'output')
fprintf('LL mat file saved\n\n')
clc;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLAUSER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
k = 0.4;
B = 5.0;

% Fitting
y_mask         = y(~isnan(u_profile));
u_profile_mask = u_profile(~isnan(u_profile));
clauser_fun    = @(Cf) ((1/k) * sqrt(Cf/2)) * log((y_mask * u_inf) / nu) + (((1/k) * sqrt(Cf/2)) * log(sqrt(Cf/2))) + B * sqrt(Cf/2);
fit_fun        = @(Cf) sum((clauser_fun(Cf) - (u_profile_mask / u_inf)).^2);
fitted_Cf      = fminsearch(fit_fun, 1E-4);
clauser_u_star = u_inf * sqrt((1/2 * fitted_Cf));

fprintf('Clauser Cf = %2.6f\n', fitted_Cf)
fprintf('Clauser u* = %3.3f\n', clauser_u_star)
fprintf('uv u* = %3.3f\n\n', ens_u_star)

% y0 = y(y > 0);
% clauser_fun_0 = @(Cf) ((1/k) * sqrt(Cf/2)) * log((y0 * u_inf) / nu) + (((1/k) * sqrt(Cf/2)) * log(sqrt(Cf/2))) + B * sqrt(Cf/2);
% inc = 1E-4;
% max_Cf = 0.01;
% Cfs = inc:inc:max_Cf;
% step = 20;
% Cf_colors = cool(linspace(0,length(Cfs),1));

linewidth = 2;

ax = figure('name', 'clauser lines');
title(strcat(experiment_name, ': Clauser Curves'), 'Interpreter', 'none')
hold on
plot(u_profile / u_inf, y, 'Color', 'black', 'LineStyle', '--', 'linewidth', linewidth, 'DisplayName', 'Data')
plot(clauser_fun(fitted_Cf), y_mask, 'Color', 'red', 'linewidth', linewidth, 'DisplayName', sprintf('Cf = %2.6f', fitted_Cf))
% for i = 1:step:length(Cfs)
%     Cf = Cfs(i);
%     plot(clauser_fun_0(Cf), y0, 'color', Cf_colors(i,:), 'linewidth', linewidth)
% end
hold off
ylim([0, top_bound_value])
xlabel('$u / u_{\infty}$', 'Interpreter', 'latex')
ylabel('y [mm]', 'Interpreter', 'latex')
legend('Location', 'northwest')



























