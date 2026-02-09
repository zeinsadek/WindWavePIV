% Plotting the hystersis observed between C_f and Re_{\theta} for the
% curvilinear velocity fields

addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps/slanCM')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
means_path = fullfile(project_path, 'means');
curvilinear_path = fullfile(project_path, 'curvilinear_new');
waves = {'A', 'B', 'C', 'D'};
wind_speeds = {'WT4', 'WT6', 'WT8'};

wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
wave_colors_rgb = [251/255, 54/255, 64/255;   % Red (A)
                   255/255, 195/255, 36/255;  % Yellow (B)
                   9/255, 129/255, 74/255;    % Green (C)
                   27/255, 231/255, 255/255]; % Cyan (D)

figure_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/pdf_test6';

% Approximate wavelengths in mm for labeling plots
wavelength_names.A = '410';
wavelength_names.B = '313';
wavelength_names.C = '189';
wavelength_names.D = '124';

steepnesses.A = '0.180';
steepnesses.B = '0.211';
steepnesses.C = '0.305';
steepnesses.D = '0.267';

% Wave parameters
wavelengths.A = 410.8 * 1E-3;
wavelengths.B = 313.3 * 1E-3;
wavelengths.C = 189.6 * 1E-3;
wavelengths.D = 124.3 * 1E-3;

amplitudes.A = 11.78 * 1E-3;
amplitudes.B = 10.53 * 1E-3;
amplitudes.C = 9.21 * 1E-3;
amplitudes.D = 5.29 * 1E-3;

frequencies.A = 1.96;
frequencies.B = 2.27;
frequencies.C = 3;
frequencies.D = 3.93;

% Wave phase speeds (c = f * lambda)
phase_speeds.A = frequencies.A * wavelengths.A;
phase_speeds.B = frequencies.B * wavelengths.B;
phase_speeds.C = frequencies.C * wavelengths.C;
phase_speeds.D = frequencies.D * wavelengths.D;

% Wave steepness (ak)
ak.A = 0.180;
ak.B = 0.211;
ak.C = 0.305;
ak.D = 0.267;

% Freestream velocities
freestreams.WT4 = 2.4181;
freestreams.WT6 = 3.8709;
freestreams.WT8 = 5.4289;

%% Curvilinear Friction Velocity

% Wind speed to consider
for s = 1:length(wind_speeds)

    % Loop through wind speeds
    wind_speed = wind_speeds{s};
 
    if ismember(wind_speed(end), {'4'})
        u_inf = 2.4181;
    elseif ismember(wind_speed(end), {'6'})
        u_inf = 3.8709;
    elseif ismember(wind_speed(end), {'8'})
        u_inf = 5.4289;
    end
    
    % Loop through waves
    for w = 1:length(waves)
    
        % Define case
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        disp(caze)
    
        % Load data
        curvilinear = load(fullfile(curvilinear_path, strcat(caze, '_CURVILINEAR.mat')));
        curvilinear = curvilinear.output;
    
        % Wave parameters
        wave_type       = caze(strfind(caze, 'WV') + 2);
        wave_parameters = readcell("Offshore_Waves.xlsx");
        wavelength      = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 2};
        amplitude       = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 3};

        wavelengths.(wave) = wavelength * 1E-3;
        amplitudes.(wave) = amplitude * 1E-3;
    
        for phase = 1:4
        
            % Get components
            uv = curvilinear.phase(phase).uv;
            % X (xi) in mm
            vertical_lines = curvilinear.phase(phase).vertical_lines * 1E-3;
            % Y (zeta) in mm
            horizontal_lines = curvilinear.phase(phase).horizontal_lines * 1E-3;
            % Waves
            wave_profile = curvilinear.phase(phase).wave_profile * 1E-3;
            max_wave_profile = curvilinear.phase(phase).max_wave_profile * 1E-3;

            %%% Save wave_profiles
            wave_profiles.(wave)(phase).profile = wave_profile;
            wave_profiles.(wave)(phase).x = curvilinear.X(1,:) * 1E-3;
            
            % Get profiles
            idx = round(length(uv)/2);
            uv_profile = uv(:, idx); 

            % Save profiles
            normalized_uv_profiles.(caze)(phase).uv = uv_profile / (u_inf^2);
            normalized_uv_profiles.(caze)(phase).y = horizontal_lines(:, idx) - wave_profile(idx);

        
            %%% Get Friction Velocity
            friction_velocity.(caze)(phase) = max(sqrt(-uv_profile), [], 'all', 'omitnan');
            skin_friction.(caze)(phase) = 2 * (friction_velocity.(caze)(phase) / u_inf)^2;

            %%% Get Skin Friciton profiles
            for i = 1:size(uv,2)
                column = uv(:,i);
                if all(isnan(column))
                    continue;
                end
        
                [M, I] = max(sqrt(-column), [], 'omitnan');
                cfs_tmp(i) = 2 * (M / u_inf).^2;
                x_tmp(i) = vertical_lines(I, i);
            end
        
        
            % Save values
            skin_friction_profiles.(caze)(phase).x = x_tmp;
            skin_friction_profiles.(caze)(phase).cf = cfs_tmp;
        
        end
    end
end

% Also include the no-wave cases
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    caze = strcat(wind_speed, '_WV0_AGP');

    % Loop through wind speeds
    wind_speed = wind_speeds{s};

    if ismember(wind_speed(end), {'4'})
        u_inf = 2.4181;
    elseif ismember(wind_speed(end), {'6'})
        u_inf = 3.8709;
    elseif ismember(wind_speed(end), {'8'})
        u_inf = 5.4289;
    end

    % Load data
    means = load(fullfile(means_path, strcat(caze, '_MEANS.mat')));
    means = means.output;

    % Get friction velocity and cf
    uv = means.ensemble.uv;

    % Get profiles
    idx = round(length(uv)/2);
    uv_profile = uv(:, idx); 

    % Save profiles
    normalized_uv_profiles.(caze).uv = uv_profile / (u_inf^2);
    normalized_uv_profiles.(caze).y = means.Y(:,1) * 1E-3;

    %%% Get Friction Velocity
    friction_velocity.(caze) = max(sqrt(-uv_profile), [], 'all', 'omitnan');
    skin_friction.(caze) = 2 * (friction_velocity.(caze) / u_inf)^2;
end

clc;


%% Boundary layer, integral parameters

% Boundary layer detection parameters
boundary_layer_percent = 0.96;
left_mask = 4E-3;
right_mask = 234E-3;

% Edge killing values
% OG was 0.008
threshold.thickness = 0.005;
threshold.displacement = 0.005;
threshold.momentum = 0.008;

% Data filtering and smoothing
smooth_kernel = 9;
gradient_tolerance = 1;
island_length = 5;
sgolay_length = 35;
n = 9;


% Loop through all waves and wind speeds
wind_speeds = {'WT4', 'WT6', 'WT8'};
waves = {'A', 'B', 'C', 'D'};

for s = 1:length(wind_speeds)

    % Loop through wind speeds
    wind_speed = wind_speeds{s};
    
    if ismember(wind_speed(end), {'4'})
        u_inf = 2.4181;
    elseif ismember(wind_speed(end), {'6'})
        u_inf = 3.8709;
    elseif ismember(wind_speed(end), {'8'})
        u_inf = 5.4289;
    end


    for w = 1:length(waves)
    
        % Define case
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        disp(caze)
    
        % Load data
        curvilinear = load(fullfile(curvilinear_path, strcat(caze, '_CURVILINEAR.mat')));
        curvilinear = curvilinear.output;
        
        % Wave Parameters
        wave_type  = wave;
        wavelength = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 2};
        amplitude  = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 3};
    
        % Save wavelength and amplitude
        integral.(caze).wavelength = wavelength * 1E-3;
        integral.(caze).amplitude = amplitude * 1E-3;
    
    
        %%% Boundary Layer Edge Detection
        for phase = 1:4
    
            fprintf("Phase %1.0f\n", phase)
        
            % Curvilinear coordinates
            vertical_lines = curvilinear.phase(phase).vertical_lines * 1E-3;
            horizontal_lines = curvilinear.phase(phase).horizontal_lines * 1E-3;
            
            % Waves
            % wave_profile = curvilinear.phase(phase).wave_profile * 1E-3;
            % max_wave_profile = curvilinear.phase(phase).max_wave_profile * 1E-3;
            
            % Find boundary layer edge
            u = curvilinear.phase(phase).u;
            u = juliaan_smooth(u, smooth_kernel);
        
            % Mask the data properly
            u(vertical_lines < left_mask | vertical_lines > right_mask) = nan;
            vertical_lines(vertical_lines < left_mask | vertical_lines > right_mask) = nan;
            horizontal_lines(vertical_lines < left_mask | vertical_lines > right_mask) = nan;
            
            % Find BL
            u(u > u_inf * boundary_layer_percent) = nan;
           
            %%% Boundary layer thickness
            thickness = nan(1, size(u,2));
            x_BL = nan(1, size(u,2));
        
            for i = 1:size(u,2)
                column = u(:,i);
                if all(isnan(column))
                    continue;
                end
        
                [~, idx] = max(column, [], 'omitnan');
                thickness(i) = horizontal_lines(idx, i);
                x_BL(i) = vertical_lines(idx, i);
            end
        
            % Kil bad edge values
            thickness_edge_mask = gradient(thickness) > threshold.thickness | gradient(thickness) < -threshold.thickness;
            thickness(thickness_edge_mask) = nan;
            x_BL(thickness_edge_mask) = nan;
        
        
    
            %%% Momentum and displacement thicknesses
            displacement = nan(1, size(u,2));
            momentum = nan(1, size(u,2));
            
            for i = 1:size(u,2)
            
                u_column = u(:,i);
                zeta_column = horizontal_lines(:,i);
            
                % Remove NaNs
                nan_mask = ~isnan(u_column) & ~isnan(zeta_column);
                u_column = u_column(nan_mask);
                zeta_column = zeta_column(nan_mask);
            
                % not enough points to integrate
                if length(u_column) < 2
                    continue; 
                end
            
                % Normalize velocity
                u_normalized = u_column / u_inf;
                    
                % Decreasing Î¶
                if zeta_column(end) < zeta_column(1) 
                    zeta_column = flipud(zeta_column);
                    u_normalized = flipud(u_normalized);
                end
            
                % Integrands
                displacement_integrand = (1 - u_normalized);
                momentum_integrand = u_normalized .* (1 - u_normalized);
        
                % Integrals
                displacement_cum = trapz(zeta_column, displacement_integrand);
                momentum_cum = trapz(zeta_column, momentum_integrand);
            
                % Final value = full integral over zeta
                displacement(i) = displacement_cum;
                momentum(i) = momentum_cum;
        
            end
        
            % Kill bad edge values
            displacement_edge_mask = gradient(displacement) > threshold.displacement | gradient(displacement) < -threshold.displacement;
            momentum_edge_mask = gradient(momentum) > threshold.momentum | gradient(momentum) < -threshold.momentum;
        
            displacement(displacement_edge_mask) = nan;
            momentum(momentum_edge_mask) = nan;
        
    
    
            %%% Save signals
            % Save reference wave profile
            integral.(caze)(phase).wave.wave_profile = curvilinear.phase(phase).wave_profile * 1E-3;
            integral.(caze)(phase).wave.x = curvilinear.X(1,:) * 1E-3;
             
            % Save raw signals
            integral.(caze)(phase).raw.thickness = thickness;
            integral.(caze)(phase).raw.displacement = displacement;
            integral.(caze)(phase).raw.momentum = momentum;
            integral.(caze)(phase).raw.shape = displacement ./ momentum;
            integral.(caze)(phase).x = x_BL;
        
            % Save filtered
            integral.(caze)(phase).filtered.thickness = FilterData(x_BL, thickness, 10, 10, sgolay_length);
            integral.(caze)(phase).filtered.displacement = FilterData(x_BL, displacement, gradient_tolerance, island_length, sgolay_length);
            integral.(caze)(phase).filtered.momentum = FilterData(x_BL, momentum, gradient_tolerance, island_length, sgolay_length);
            integral.(caze)(phase).filtered.shape = FilterData(x_BL, displacement ./ momentum, gradient_tolerance, island_length, sgolay_length);
        
            % Save polynomial fit
            integral.(caze)(phase).fitted.thickness = PolynomialFit(x_BL, integral.(caze)(phase).filtered.thickness, n);
            integral.(caze)(phase).fitted.displacement = PolynomialFit(x_BL, integral.(caze)(phase).filtered.displacement, n);
            integral.(caze)(phase).fitted.momentum = PolynomialFit(x_BL, integral.(caze)(phase).filtered.momentum, n);
            integral.(caze)(phase).fitted.shape = PolynomialFit(x_BL, integral.(caze)(phase).filtered.shape, n);
        
        end
    end
end

clear caze curvilinear horizontal_lines idx phase uv uv_profile 
clear vertical_lines w wave wave_profile z_profile wave_type s wind_speed 
clear curvilinear_path max_wave_profile means means_path project_path u_inf
clear wave_parameters wavelength amplitude cfs_tmp column i I M x_tmp
clear amplitude caze column curvilinear displacement displacement_cum displacement_edge_mask displacement_integrand
clear horizontal_lines i idx island_length left_mask max_wave_profile momentum momentum_cum momentum_edge_mask momentum_integrand
clear n nan_mask phase right_mask s sgolay_length smooth_kernerl thickness thickness_edge_mask threshold u u_column u_inf u_normalized
clear vertical_lines w wave wave_profile wave_parameters wave_type wave_length wind_speed x_BL zeta_column
clear gradient_tolerance project_path smooth_kernel wavelength curvilinear_path
clc;


%% ========================================================================
%  LEEWARD FACE SCALING PARAMETER ANALYSIS
%  ========================================================================
%  Compute and plot the scaling parameter: delta * (U_inf - c) / (lambda * U_inf)
%  to see if it collapses the phase 4 friction velocity data
%  ========================================================================

clc; close all;

% Phase to analyze (4 = lee face)
phase_to_analyze = 4;

% Storage for scatter plot data
scaling_data = struct();
scaling_data.parameter = [];
scaling_data.ustar = [];
scaling_data.ustar_normalized = [];
scaling_data.cf = [];
scaling_data.wave_idx = [];
scaling_data.wind_idx = [];
scaling_data.case_name = {};

% Also try other candidate parameters
scaling_data.wave_age = [];           % c / U_inf
scaling_data.delta_over_lambda = [];  % delta / lambda
scaling_data.ak = [];                 % wave steepness
scaling_data.Re_lambda = [];          % U_inf * lambda / nu
scaling_data.relative_velocity = [];  % (U_inf - c) / U_inf

% Kinematic viscosity of air
nu = 1.5e-5;  % m^2/s

% Loop through all cases
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    u_inf = freestreams.(wind_speed);
    
    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        
        % Get parameters
        lambda = wavelengths.(wave);
        c = phase_speeds.(wave);
        steepness = ak.(wave);
        
        % Get boundary layer thickness at phase 4
        % Use mean of the filtered thickness profile
        delta = mean(integral.(caze)(phase_to_analyze).filtered.thickness, 'omitnan');
        
        % Get friction velocity at phase 4
        ustar = friction_velocity.(caze)(phase_to_analyze);
        cf = skin_friction.(caze)(phase_to_analyze);
        
        % Compute the scaling parameter: delta * (U_inf - c) / (lambda * U_inf)
        scaling_parameter = delta * (u_inf - c) / (lambda * u_inf);
        
        % Store data
        scaling_data.parameter(end+1) = scaling_parameter;
        scaling_data.ustar(end+1) = ustar;
        scaling_data.ustar_normalized(end+1) = ustar / u_inf;
        scaling_data.cf(end+1) = cf;
        scaling_data.wave_idx(end+1) = w;
        scaling_data.wind_idx(end+1) = s;
        scaling_data.case_name{end+1} = caze;
        
        % Other candidate parameters
        scaling_data.wave_age(end+1) = c / u_inf;
        scaling_data.delta_over_lambda(end+1) = delta / lambda;
        scaling_data.ak(end+1) = steepness;
        scaling_data.Re_lambda(end+1) = u_inf * lambda / nu;
        scaling_data.relative_velocity(end+1) = (u_inf - c) / u_inf;
        
        % Print for debugging
        fprintf('%s: delta=%.4f m, lambda=%.4f m, c=%.3f m/s, U_inf=%.3f m/s\n', ...
                caze, delta, lambda, c, u_inf);
        fprintf('   Scaling param = %.4f, u* = %.4f m/s\n\n', scaling_parameter, ustar);
    end
end


%% Plot 1: u* vs the main scaling parameter (colored by wave type)

figure('color', 'white', 'Position', [100 100 600 500])
hold on

% Plot each point colored by wave type
for i = 1:length(scaling_data.parameter)
    w = scaling_data.wave_idx(i);
    scatter(scaling_data.parameter(i), scaling_data.ustar(i), 100, ...
            wave_colors_rgb(w,:), 'filled', 'MarkerEdgeColor', 'k')
end

% Add labels for wind speeds
for i = 1:length(scaling_data.parameter)
    text(scaling_data.parameter(i) + 0.002, scaling_data.ustar(i), ...
         wind_speeds{scaling_data.wind_idx(i)}, 'FontSize', 8)
end

xlabel('$\frac{\delta (U_\infty - c)}{\lambda U_\infty}$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$u^*_\zeta$ [m/s]', 'Interpreter', 'latex', 'FontSize', 14)
title(sprintf('Phase %d (Lee Face): Friction Velocity vs Scaling Parameter', phase_to_analyze))
grid on

% Add legend
legend_entries = {['\lambda_{', wavelength_names.A, '}, ak_{', steepnesses.A, '}'], ...
                  ['\lambda_{', wavelength_names.B, '}, ak_{', steepnesses.B, '}'], ...
                  ['\lambda_{', wavelength_names.C, '}, ak_{', steepnesses.C, '}'], ...
                  ['\lambda_{', wavelength_names.D, '}, ak_{', steepnesses.D, '}']};
% Create dummy plots for legend
for w = 1:4
    h(w) = scatter(nan, nan, 100, wave_colors_rgb(w,:), 'filled', 'MarkerEdgeColor', 'k');
end
legend(h, legend_entries, 'Location', 'best')

hold off


%% Plot 2: u*/U_inf vs scaling parameter (normalized friction velocity)

figure('color', 'white', 'Position', [100 100 600 500])
hold on

for i = 1:length(scaling_data.parameter)
    w = scaling_data.wave_idx(i);
    scatter(scaling_data.parameter(i), scaling_data.ustar_normalized(i), 100, ...
            wave_colors_rgb(w,:), 'filled', 'MarkerEdgeColor', 'k')
end

xlabel('$\frac{\delta (U_\infty - c)}{\lambda U_\infty}$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$u^*_\zeta / U_\infty$', 'Interpreter', 'latex', 'FontSize', 14)
title(sprintf('Phase %d (Lee Face): Normalized Friction Velocity', phase_to_analyze))
grid on

% Add legend
for w = 1:4
    h(w) = scatter(nan, nan, 100, wave_colors_rgb(w,:), 'filled', 'MarkerEdgeColor', 'k');
end
legend(h, legend_entries, 'Location', 'best')

hold off


%% Plot 3: Compare multiple candidate parameters (2x3 subplot)

figure('color', 'white', 'Position', [100 100 1200 700])
tiledlayout(2, 3, 'TileSpacing', 'compact')

% Parameter 1: Main scaling parameter
nexttile
hold on
for i = 1:length(scaling_data.parameter)
    w = scaling_data.wave_idx(i);
    scatter(scaling_data.parameter(i), scaling_data.ustar(i), 80, ...
            wave_colors_rgb(w,:), 'filled', 'MarkerEdgeColor', 'k')
end
xlabel('$\frac{\delta (U_\infty - c)}{\lambda U_\infty}$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$u^*_\zeta$ [m/s]', 'Interpreter', 'latex', 'FontSize', 12)
title('Main Scaling Parameter')
grid on
hold off

% Parameter 2: Wave age (c/U_inf)
nexttile
hold on
for i = 1:length(scaling_data.wave_age)
    w = scaling_data.wave_idx(i);
    scatter(scaling_data.wave_age(i), scaling_data.ustar(i), 80, ...
            wave_colors_rgb(w,:), 'filled', 'MarkerEdgeColor', 'k')
end
xlabel('$c / U_\infty$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$u^*_\zeta$ [m/s]', 'Interpreter', 'latex', 'FontSize', 12)
title('Wave Age')
grid on
hold off

% Parameter 3: delta/lambda
nexttile
hold on
for i = 1:length(scaling_data.delta_over_lambda)
    w = scaling_data.wave_idx(i);
    scatter(scaling_data.delta_over_lambda(i), scaling_data.ustar(i), 80, ...
            wave_colors_rgb(w,:), 'filled', 'MarkerEdgeColor', 'k')
end
xlabel('$\delta / \lambda$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$u^*_\zeta$ [m/s]', 'Interpreter', 'latex', 'FontSize', 12)
title('BL Thickness / Wavelength')
grid on
hold off

% Parameter 4: Wave steepness (ak)
nexttile
hold on
for i = 1:length(scaling_data.ak)
    w = scaling_data.wave_idx(i);
    scatter(scaling_data.ak(i), scaling_data.ustar(i), 80, ...
            wave_colors_rgb(w,:), 'filled', 'MarkerEdgeColor', 'k')
end
xlabel('$ak$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$u^*_\zeta$ [m/s]', 'Interpreter', 'latex', 'FontSize', 12)
title('Wave Steepness')
grid on
hold off

% Parameter 5: Relative velocity (U_inf - c)/U_inf
nexttile
hold on
for i = 1:length(scaling_data.relative_velocity)
    w = scaling_data.wave_idx(i);
    scatter(scaling_data.relative_velocity(i), scaling_data.ustar(i), 80, ...
            wave_colors_rgb(w,:), 'filled', 'MarkerEdgeColor', 'k')
end
xlabel('$(U_\infty - c) / U_\infty$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$u^*_\zeta$ [m/s]', 'Interpreter', 'latex', 'FontSize', 12)
title('Relative Velocity')
grid on
hold off

% Parameter 6: Re_lambda
nexttile
hold on
for i = 1:length(scaling_data.Re_lambda)
    w = scaling_data.wave_idx(i);
    scatter(scaling_data.Re_lambda(i), scaling_data.ustar(i), 80, ...
            wave_colors_rgb(w,:), 'filled', 'MarkerEdgeColor', 'k')
end
xlabel('$Re_\lambda = U_\infty \lambda / \nu$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$u^*_\zeta$ [m/s]', 'Interpreter', 'latex', 'FontSize', 12)
title('Wavelength Reynolds Number')
grid on
hold off

sgtitle(sprintf('Phase %d (Lee Face): Candidate Scaling Parameters', phase_to_analyze))


%% Plot 4: Colored by wind speed instead of wave type

figure('color', 'white', 'Position', [100 100 600 500])
hold on

wind_colors = [0.2, 0.4, 0.8;    % Blue for WT4
               0.8, 0.4, 0.2;    % Orange for WT6  
               0.4, 0.8, 0.2];   % Green for WT8

markers = {'o', 's', 'd', '^'};  % Different markers for each wave

for i = 1:length(scaling_data.parameter)
    w = scaling_data.wave_idx(i);
    s = scaling_data.wind_idx(i);
    scatter(scaling_data.parameter(i), scaling_data.ustar(i), 100, ...
            wind_colors(s,:), 'filled', markers{w}, 'MarkerEdgeColor', 'k')
end

xlabel('$\frac{\delta (U_\infty - c)}{\lambda U_\infty}$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$u^*_\zeta$ [m/s]', 'Interpreter', 'latex', 'FontSize', 14)
title(sprintf('Phase %d: Colored by Wind Speed, Marker by Wave', phase_to_analyze))
grid on

% Legend for wind speeds
h_wind = [];
for s = 1:3
    h_wind(s) = scatter(nan, nan, 100, wind_colors(s,:), 'filled', 'o', 'MarkerEdgeColor', 'k');
end
legend(h_wind, wind_speeds, 'Location', 'best')

hold off


%% Print summary table

fprintf('\n========== PHASE 4 SCALING ANALYSIS SUMMARY ==========\n\n')
fprintf('%-15s | %8s | %8s | %8s | %10s | %8s\n', ...
        'Case', 'u*', 'u*/U_inf', 'delta/lam', 'Param', 'c/U_inf')
fprintf('%s\n', repmat('-', 1, 75))

for i = 1:length(scaling_data.parameter)
    fprintf('%-15s | %8.4f | %8.4f | %8.4f | %10.4f | %8.3f\n', ...
            scaling_data.case_name{i}, ...
            scaling_data.ustar(i), ...
            scaling_data.ustar_normalized(i), ...
            scaling_data.delta_over_lambda(i), ...
            scaling_data.parameter(i), ...
            scaling_data.wave_age(i))
end


%% Compute correlation coefficients for each parameter

fprintf('\n========== CORRELATION ANALYSIS ==========\n\n')

params = {'parameter', 'wave_age', 'delta_over_lambda', 'ak', 'relative_velocity', 'Re_lambda'};
param_names = {'delta(U-c)/(lambda*U)', 'c/U_inf', 'delta/lambda', 'ak', '(U-c)/U', 'Re_lambda'};

for p = 1:length(params)
    x = scaling_data.(params{p});
    y = scaling_data.ustar;
    
    R = corrcoef(x, y);
    r = R(1,2);
    
    fprintf('%25s: R = %+.3f (R^2 = %.3f)\n', param_names{p}, r, r^2)
end



% Functions

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

function output = StatisticalGradientFilter(q, tol, max_island_width)
    % Take gradient
    grad = gradient(q);
    % Get standard deviation
    sd = std(grad, 0, 'omitnan');
    % Get average
    avg = mean(grad, 'all', 'omitnan');
    % Logical indexing
    q(grad > avg + tol * sd | grad < avg - tol * sd) = nan;
    % Remove Islands
    q = RemoveIslands(q,max_island_width);
    output = q;
end


% Chat GPT
function cleaned = RemoveIslands(signal, max_island_width)

    % Create logical mask of NaNs
    is_nan = isnan(signal);

    % Find islands of valid data (i.e., ~isnan) between NaNs
    cleaned = signal; 
    n = length(signal);
    i = 1;
    while i <= n
        % Skip NaNs
        if is_nan(i)
            i = i + 1;
            continue;
        end

        % Start of a non-NaN segment
        start_idx = i;
        while i <= n && ~is_nan(i)
            i = i + 1;
        end
        end_idx = i - 1;

        % Check if island is surrounded by NaNs
        is_surrounded = false;
        if start_idx > 1 && end_idx < n
            if is_nan(start_idx - 1) && is_nan(end_idx + 1)
                is_surrounded = true;
            end
        end

        % Remove if surrounded and short
        if is_surrounded && (end_idx - start_idx + 1) <= max_island_width
            cleaned(start_idx:end_idx) = NaN;
        end
    end
end


function interpolated = InterpolateInterior(x, y)
    % x: 1D array of independent variable (may contain NaNs)
    % y: 1D array of dependent variable (may contain NaNs)

    % Step 1: Find joint validity mask
    valid_mask = ~isnan(x) & ~isnan(y);

    % Step 2: Identify first and last *jointly valid* indices
    first_idx = find(valid_mask, 1, 'first');
    last_idx  = find(valid_mask, 1, 'last');

    if isempty(first_idx) || isempty(last_idx) || first_idx == last_idx
        % Nothing valid to interpolate
        interpolated = y;
        return
    end

    % Step 3: Trim x and y to internal valid range
    x_trim = x(first_idx:last_idx);
    y_trim = y(first_idx:last_idx);

    % Step 4: Mask internal NaNs in y, only where x is valid
    internal_valid = ~isnan(x_trim);
    nan_mask = isnan(y_trim) & internal_valid;

    % Step 5: Interpolate
    y_interp = y_trim;
    y_interp(nan_mask) = interp1(x_trim(~isnan(y_trim) & ~isnan(x_trim)), ...
                                 y_trim(~isnan(y_trim) & ~isnan(x_trim)), ...
                                 x_trim(nan_mask), 'spline');

    % Step 6: Rebuild full output with original NaN pads
    interpolated = y;
    interpolated(first_idx:last_idx) = y_interp;
end


function output = FilterData(x, q, tol, max_island_length, sgolay_length)
    output = smoothdata(InterpolateInterior(x, StatisticalGradientFilter(q, tol, max_island_length)), 'sgolay', sgolay_length);
end


function output = PolynomialFit(x, y, n)
    mask = ~isnan(x) & ~isnan(y);
    p = polyfit(x(mask), y(mask), n);
    output = polyval(p, x);
    output(isnan(y)) = nan;
end
