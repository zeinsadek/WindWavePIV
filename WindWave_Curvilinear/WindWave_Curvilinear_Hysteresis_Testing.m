% Plotting the hystersis observed between C_f and Re_{\theta} for the
% curvilinear velocity fields

addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
means_path = fullfile(project_path, 'means');
curvilinear_path = fullfile(project_path, 'curvilinear');
waves = {'A', 'B', 'C', 'D'};
wind_speeds = {'WT4', 'WT6', 'WT8'};

wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};

figure_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new';

% Approximate wavelengths in mm for labeling plots
wavelength_names.A = '410';
wavelength_names.B = '313';
wavelength_names.C = '189';
wavelength_names.D = '124';

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

    %%% Get Skin Friciton profiles
    for i = 1:size(uv,2)
        column = uv(:,i);
        if all(isnan(column))
            continue;
        end

        [M, I] = max(sqrt(-column), [], 'omitnan');
        cfs_tmp(i) = 2 * (M / u_inf).^2;
        % x_tmp(i) = vertical_lines(I, i);
    end


    % Save values
    % skin_friction_profiles.(caze)(phase).x = x_tmp;
    skin_friction_profiles.(caze).cf = cfs_tmp;
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
                    
                % Decreasing ζ
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

    % Cartesian coordinates
    X = means.X * 1E-3;
    Y = means.Y * 1E-3;
    
    % Find boundary layer edge
    u = means.ensemble.u;
    u = juliaan_smooth(u, smooth_kernel);

    % Mask the data properly
    u(X < left_mask | X > right_mask) = nan;
    X(X < left_mask | X > right_mask) = nan;
    Y(X < left_mask | X > right_mask) = nan;
    
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
        thickness(i) = Y(idx, i);
        x_BL(i) = X(idx, i);
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
        zeta_column = Y(:,i);
    
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
            
        % Decreasing ζ
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

    % Save filtered
    integral.(caze).filtered.thickness = FilterData(x_BL, thickness, 10, 10, sgolay_length);
    integral.(caze).filtered.displacement = FilterData(x_BL, displacement, gradient_tolerance, island_length, sgolay_length);
    integral.(caze).filtered.momentum = FilterData(x_BL, momentum, gradient_tolerance, island_length, sgolay_length);
    integral.(caze).filtered.shape = FilterData(x_BL, displacement ./ momentum, gradient_tolerance, island_length, sgolay_length);
    integral.(caze).filtered.x = X(1,:);

end

clc;

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


%% Check no wave values


parameter = 'momentum';
version = 'filtered';
linestyles = {'-.', '--', '-'};
lw = 2;


tickFontSize = 8;
labelFontSize = 16;
legendFontSize = 12;

close all; 
ax = figure('color', 'white');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)


hold on
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};

    if ismember(wind_speed(end), {'4'})
        u_inf = 2.4181;
    elseif ismember(wind_speed(end), {'6'})
        u_inf = 3.8709;
    elseif ismember(wind_speed(end), {'8'})
        u_inf = 5.4289;
    end

    caze = strcat(wind_speed, '_WV0_AGP');

    
    % label = sprintf('$\\lambda_{%s}, u_{\\infty} = %1.2f$', wavelengths.(wave), u_inf);
    
    data = integral.(caze).(version).(parameter);
    x = integral.(caze).(version).x;
    
    plot(x, data, ...
                'linewidth', lw, 'linestyle', linestyles{s}, ...
                'displayname', wind_speed);
end
hold off


%% Plot Re_{\theta} vs Cf for no-wave cases



parameter = 'momentum';
version = 'filtered';
linestyles = {'-.', '--', '-'};
lw = 2;
nu = 1.48E-5;

tickFontSize = 8;
labelFontSize = 16;
legendFontSize = 12;

smoothing_kernel = 16;

close all; 
ax = figure('color', 'white');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)


hold on
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};

    if ismember(wind_speed(end), {'4'})
        u_inf = 2.4181;
    elseif ismember(wind_speed(end), {'6'})
        u_inf = 3.8709;
    elseif ismember(wind_speed(end), {'8'})
        u_inf = 5.4289;
    end

    caze = strcat(wind_speed, '_WV0_AGP');

    % Get data
    buff = 30;
    x = integral.(caze).(version).x(buff:end-buff);
    momentum = smoothdata(integral.(caze).(version).(parameter)(buff:end-buff), 'movmean', smoothing_kernel);
    cf = smoothdata(skin_friction_profiles.(caze).cf(buff:end-buff), 'movmean', smoothing_kernel);
    % x = integral.(caze).(version).x;

    % Compute Re_{\theta}
    Re_theta = (u_inf * momentum) / nu;

    % Plot
    % plot(Re_theta, cf, ...
    %      'linewidth', lw, 'linestyle', linestyles{s}, ...
    %      'displayname', wind_speed);

    plot(x, cf, ...
         'linewidth', lw, 'linestyle', linestyles{s}, ...
         'displayname', wind_speed);
end
hold off






%% Plot Cf profiles to check

tickFontSize = 8;
labelFontSize = 16;
legendFontSize = 12;

wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
linestyles = {'-.', '--', '-'};
lw = 2;
buff = 5;
phase_4_Re_theta_scale = 1.25;

smoothing_kernel = 32;

x_crop = 1;

close all; 
figure('color', 'white');
tiledlayout(4,1)
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)

% Loop through phases
for phase = 1:4

    h(phase) = nexttile;
    hold on
    for s = 1:length(wind_speeds)
        wind_speed = wind_speeds{s};
    
        if ismember(wind_speed(end), {'4'})
            u_inf = 2.4181;
        elseif ismember(wind_speed(end), {'6'})
            u_inf = 3.8709;
        elseif ismember(wind_speed(end), {'8'})
            u_inf = 5.4289;
        end
    
        % for w = 1:length(waves)
        for w = 4
            wave = waves{w};
    
            caze = strcat(wind_speed, '_WV', wave, '_AG0');
            disp(caze)
    
            % Crop ends of data
            data = skin_friction_profiles.(caze)(phase).cf(buff:end-buff);
    
            % Smooth data
            data = smoothdata(data, 'movmean', smoothing_kernel);
    
            wavelength = wavelengths.(wave); 
            amplitude = amplitudes.(wave); 
            steepness = (2 * amplitude) / wavelength;
    
            x = skin_friction_profiles.(caze)(phase).x;
    
            label = sprintf('$\\lambda_{%s}, u_{\\infty} = %1.2f$', wavelength_names.(wave), u_inf);

            centered_normalized_x = (x(buff:end-buff) - mean(x, 'all', 'omitnan')) / wavelengths.(wave);

            % Crop data to one wavelength
            data(abs(centered_normalized_x) > x_crop) = nan;
    
            plot(centered_normalized_x, data, ...
                 'color', wave_colors{w}, 'linestyle', linestyles{s}, 'linewidth', lw, 'DisplayName', label)
    
            if wave == 'D' && s == 1
                tmp = plot((wave_profiles.(wave)(phase).x - mean(wave_profiles.(wave)(phase).x)) / wavelengths.(wave), phase_4_Re_theta_scale *  wave_profiles.(wave)(phase).profile + mean(data, 'all', 'omitnan'), ...
                            'color', 'black', 'linestyle', '--', 'linewidth', lw, 'handlevisibility', 'off');
                tmp.Color(4) = 0.5;
            end
    
        end
    end
    ylim([0, 0.02])
end

linkaxes(h, 'xy')

hold off
xlim([-1, 1])
ylim([0, 0.02])
yticks(0:0.005:0.02)
% legend('interpreter', 'latex', 'fontsize', 10, 'location',  'northeastoutside', 'box', 'off')
xlabel('$x / \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$C_f$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xticks(-1:0.25:1)
clc;






%% Plot to check momentum thickness curves


parameter = 'momentum';
version = 'filtered';
linestyles = {'-.', '--', '-'};
lw = 2;
alpha = 0.5;

x_crop = 1;

nu = 1.48E-5;

% Scale reference wave based on parameter
if strcmp(parameter, 'thickness')
    phase_4_Re_theta_scale = 8;
elseif strcmp(parameter, 'displacement')
    phase_4_Re_theta_scale = 3;
elseif strcmp(parameter, 'momentum')
    phase_4_Re_theta_scale = 1;
elseif strcmp(parameter, 'shape')
    phase_4_Re_theta_scale = 35;
end


tickFontSize = 8;
labelFontSize = 16;
legendFontSize = 12;

close all; 
ax = figure('color', 'white');
tiledlayout(2,1)
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)

for phase = 1:4
    h(phase) = nexttile;

    hold on
    for s = 1:length(wind_speeds)
        wind_speed = wind_speeds{s};
    
        if ismember(wind_speed(end), {'4'})
            u_inf = 2.4181;
        elseif ismember(wind_speed(end), {'6'})
            u_inf = 3.8709;
        elseif ismember(wind_speed(end), {'8'})
            u_inf = 5.4289;
        end
    
        % for w = 1:length(waves)
        for w = 4
            wave = waves{w};
            caze = strcat(wind_speed, '_WV', wave, '_AG0');
    
            wavelength = integral.(caze).wavelength; 
            amplitude = integral.(caze).amplitude;
            steepness = (2 * amplitude) / wavelength;
            
            label = sprintf('$\\lambda_{%s}, u_{\\infty} = %1.2f$', wavelengths.(wave), u_inf);
            
            data = integral.(caze)(phase).(version).(parameter);

            centered_normalized_x = (integral.(caze)(phase).x - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength;

            % Crop data
            data(abs(centered_normalized_x) > x_crop) = nan;

            % Compute Re_{\theta}
            Re_theta = (u_inf * data) / nu;

            H(w) = plot(centered_normalized_x, data, ...
                        'color', wave_colors{w}, 'linewidth', lw, 'linestyle', linestyles{s}, ...
                        'displayname', label);
    
            % Plot phase for reference
            if wave == 'D' && s == 1
                reference = plot((integral.(caze)(phase).wave.x  - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, (phase_4_Re_theta_scale * integral.(caze)(phase).wave.wave_profile) + mean(data, 'all', 'omitnan'), ...
                           'color', 'black', 'linewidth', lw, 'linestyle', '--', 'handlevisibility', 'off');
                reference.Color(4) = alpha;
            end
    
        end
    
    end
    hold off
    % legend('interpreter', 'latex', 'fontsize', 10, 'location',  'northeastoutside', 'box', 'off')
    
    % Axes labels
    if strcmp(parameter, 'thickness')
        vertLabel = '$\delta$ [m]';
    elseif strcmp(parameter, 'displacement')
        vertLabel = '$\delta^*$ [m]';
    elseif strcmp(parameter, 'momentum')
        vertLabel = '$\theta$ [m]';
    elseif strcmp(parameter, 'shape')
        vertLabel = '$H$';
        limits = [1, 1.7];
    end
end

linkaxes(h, 'xy')
xlabel('$x / \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(vertLabel, 'interpreter', 'latex', 'fontsize', labelFontSize)
xticks(-1:0.25:1)


%% Plot Hysteresis for wave D all phases

wind_speed = 'WT8';
wave = 'D';
caze = strcat(wind_speed, '_WV', wave, '_AG0');

if ismember(wind_speed(end), {'4'})
    u_inf = 2.4181;
elseif ismember(wind_speed(end), {'6'})
    u_inf = 3.8709;
elseif ismember(wind_speed(end), {'8'})
    u_inf = 5.4289;
end

smoothing_kernel = 8;

clc; close all
figure('color', 'white')
% tiledlayout(1,4)

hold on
for phase = 1:4
    % h(phase) = nexttile;
    % disp(phase)
    % title(num2str(phase))

    % Skin friction
    cf = skin_friction_profiles.(caze)(phase).cf;
    cf = smoothdata(cf, 'movmean', smoothing_kernel);
    
    % Re_{\theta}
    theta = integral.(caze)(phase).filtered.momentum;
    theta = smoothdata(theta, 'movmean', smoothing_kernel);
    Re_theta = (u_inf * theta) / nu;
    
    centered_normalized_x = (integral.(caze)(phase).x - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength;

    min_array_size = min([length(cf), length(Re_theta)]);
    centered_normalized_x = centered_normalized_x(1:min_array_size);
    cf = cf(1:min_array_size);
    Re_theta = Re_theta(1:min_array_size);

    % Crop data
    x_mask = centered_normalized_x;
    x_mask(abs(centered_normalized_x) > 1) = nan;
    x_mask(~isnan(x_mask)) = 1;

    cf = cf .* x_mask;
    Re_theta = Re_theta .* x_mask;
    masked_x = centered_normalized_x .* x_mask;


    % plot(Re_theta, cf, 'linewidth', 2, 'handlevisibility', 'off')
    scatter(Re_theta, cf, 50, 'filled', 'displayname', num2str(phase))

end
hold off
% legend()
% linkaxes(h, 'xy')








%% Plotting cropped phases 2 and 4 together + saving data

wind_speed = 'WT8';
wave = 'D';
caze = strcat(wind_speed, '_WV', wave, '_AG0');

if ismember(wind_speed(end), {'4'})
    u_inf = 2.4181;
elseif ismember(wind_speed(end), {'6'})
    u_inf = 3.8709;
elseif ismember(wind_speed(end), {'8'})
    u_inf = 5.4289;
end

smoothing_kernel = 2;

% clc; close all
% figure('color', 'white')
% 
% hold on
for phase = 1:4
    disp(phase)

    % Skin friction
    cf = skin_friction_profiles.(caze)(phase).cf;
    cf = smoothdata(cf, 'movmean', smoothing_kernel);
    
    % Re_{\theta}
    theta = integral.(caze)(phase).filtered.momentum;
    theta = smoothdata(theta, 'movmean', smoothing_kernel);
    Re_theta = (u_inf * theta) / nu;
    
    centered_normalized_x = (integral.(caze)(phase).x - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength;

    min_array_size = min([length(cf), length(Re_theta)]);
    centered_normalized_x = centered_normalized_x(1:min_array_size);
    cf = cf(1:min_array_size);
    Re_theta = Re_theta(1:min_array_size);

    % Crop data to 0.5 lambda
    x_mask = centered_normalized_x;
    x_mask(abs(centered_normalized_x) > 0.25) = nan;
    x_mask(~isnan(x_mask)) = 1;

    cf = cf .* x_mask;
    Re_theta = Re_theta .* x_mask;
    masked_x = centered_normalized_x .* x_mask;

    % Save data to handle stitching after
    hysteresis(phase).cf = cf(~isnan(cf));
    hysteresis(phase).Re_theta = Re_theta(~isnan(Re_theta));
    hysteresis(phase).x = masked_x(~isnan(masked_x));

    % Crop data to 0.25 lambda
    x_mask = centered_normalized_x;
    x_mask(abs(centered_normalized_x) > 0.125) = nan;
    x_mask(~isnan(x_mask)) = 1;

    cf = cf .* x_mask;
    Re_theta = Re_theta .* x_mask;
    % masked_x = centered_normalized_x .* x_mask;

    % Save data to handle stitching after
    hysteresis(phase).cf_small = cf(~isnan(cf));
    hysteresis(phase).Re_theta_small = Re_theta(~isnan(Re_theta));
    % hysteresis(phase).x = masked_x(~isnan(masked_x));

    % % Plot
    % if phase == 2
    %     color = 'red';
    % elseif phase == 4
    %     color = 'blue';
    % end
    % 
    % plot(Re_theta, cf, 'linewidth', 1, 'color', color)
    % scatter(Re_theta, cf, 40, 'filled', 'MarkerFaceColor', color)

end
hold off

% xlabel('Re_{\theta}')
% ylabel('Skin Friction')



%%% Stitching the two phases with offset and scaling
% Data per phase
phase_1_cf = hysteresis(1).cf;
phase_1_Re_theta = hysteresis(1).Re_theta;

phase_2_cf = hysteresis(2).cf;
phase_2_Re_theta = hysteresis(2).Re_theta;

phase_3_cf = hysteresis(3).cf;
phase_3_Re_theta = hysteresis(3).Re_theta;

phase_4_cf = hysteresis(4).cf;
phase_4_Re_theta = hysteresis(4).Re_theta;

% Compute shifts
phase_4_Re_theta_shift = phase_4_Re_theta(1) - phase_2_Re_theta(end);
phase_4_Re_theta_shifted = phase_4_Re_theta - phase_4_Re_theta_shift;

phase_4_cf_shift = phase_4_cf(1) - phase_2_cf(end);
phase_4_cf_shifted = phase_4_cf - phase_4_cf_shift;

% Compute scaling
phase_4_Re_theta_scale = phase_2_Re_theta(1) / phase_4_Re_theta_shifted(end);
phase_4_cf_scale = phase_2_cf(1) / phase_4_cf_shifted(end);

phase_4_Re_theta_fading = linspace(1, phase_4_Re_theta_scale, length(phase_4_Re_theta_shifted));
phase_4_cf_fading = linspace(1, phase_4_cf_scale, length(phase_4_cf_shifted));

% Shifted and scaled data
phase_4_Re_theta_shifted_scaled = phase_4_Re_theta_fading .* phase_4_Re_theta_shifted;
phase_4_cf_shifted_scaled = phase_4_cf_fading .* phase_4_cf_shifted;



%%% Plotting
tickFontSize = 8;
labelFontSize = 16;
legendFontSize = 12;
titleFontSize = 28;

Cf_LL = 2E-3;
Cf_UL = 20E-3;

Re_theta_LL = 2500;
Re_theta_UL = 5000;

Cf_ticks = Cf_LL:0.005:Cf_UL;
Re_theta_ticks = Re_theta_LL:500:Re_theta_UL;

clc; close all
lw = 3;
figure('color', 'white')
% tiledlayout(1,2)

% Plot hysteresis
hh(1) = nexttile;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
hold on
scatter(phase_1_Re_theta, phase_1_cf, 40, 'filled', 'MarkerFaceColor', 'red', 'displayname', 'phase 1')
scatter(phase_2_Re_theta, phase_2_cf, 40, 'filled', 'MarkerFaceColor', 'green', 'displayname', 'phase 2')
scatter(phase_3_Re_theta, phase_3_cf, 40, 'filled', 'MarkerFaceColor', 'cyan', 'displayname', 'phase 3')
scatter(phase_4_Re_theta, phase_4_cf, 40, 'filled', 'MarkerFaceColor', 'blue', 'displayname', 'phase 4')


% scatter(phase_2_Re_theta, phase_2_cf, 40, 'filled', 'MarkerFaceColor', 'black')
% scatter(phase_4_Re_theta_shifted_scaled, phase_4_cf_shifted_scaled, 40, 'filled', 'MarkerFaceColor', 'black')

hold off
legend('location', 'northwest')
title('Raw data')
xlabel('$Re_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$C_f$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xticks(Re_theta_ticks)
yticks(Cf_ticks)



hh(2) = nexttile;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
hold on
scatter(phase_2_Re_theta, phase_2_cf, 40, 'filled', 'MarkerFaceColor', 'red')
scatter(phase_4_Re_theta_shifted_scaled, phase_4_cf_shifted_scaled, 40, 'filled', 'MarkerFaceColor', 'green')

for i = 1:4
    scatter(hysteresis(i).Re_theta_small, hysteresis(i).cf_small, 30, 'filled')
end

hold off
title('Shifted and Scaled data')
xlabel('$Re_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlim([Re_theta_LL, Re_theta_UL])
ylim([Cf_LL, Cf_UL])
xticks(Re_theta_ticks)
yticks(Cf_ticks)


% Print shifts
fprintf('Re_theta shift and scaling\n\n')
fprintf('Shift applied: %.4f\n', phase_4_Re_theta_shift);
fprintf('Scale applied progressively up to: %.4f\n\n', phase_4_Re_theta_scale);

fprintf('Cf shift and scaling\n')
fprintf('Shift applied: %.4f\n', phase_4_cf_shift);
fprintf('Scale applied progressively up to: %.4f\n', phase_4_cf_scale);




% % Plot profiles
% hp(1) = nexttile([1,2]);
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% hold on
% yline(hysteresis(2).cf(1), 'HandleVisibility', 'off');
% yline(hysteresis(2).cf(end), 'HandleVisibility', 'off');
% plot((hysteresis(2).x) - 0.25, hysteresis(2).cf, 'linewidth', lw, 'color', 'red', 'displayname', 'Phase 2')
% plot((hysteresis(4).x + 0.5) - 0.25, hysteresis(4).cf, 'linewidth', lw, 'color', 'blue', 'linestyle', '-', 'Displayname', 'Phase 4 Raw')
% plot((hysteresis(4).x + 0.5) - 0.25, phase_4_cf_shifted_scaled, 'linewidth', lw, 'color', 'green', 'Displayname', 'Phase 4 Shifted + Scaled')
% 
% % Plot reference wave profile
% Cf_reference_scale = 0.5;
% Cf_reference_offset = mean([Cf_LL, Cf_UL]);
% reference = plot((integral.(caze)(1).wave.x  - mean(integral.(caze)(1).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, (Cf_reference_scale * integral.(caze)(1).wave.wave_profile) + Cf_reference_offset, ...
%                            'color', 'black', 'linewidth', lw, 'linestyle', '--', 'handlevisibility', 'off');
% reference.Color(4) = 0.5;
% 
% hold off
% ylabel('$C_f$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylim([Cf_LL, Cf_UL])
% yticks(Cf_ticks)
% xlim([-0.55, 0.55])
% xticks(-0.5:0.25:0.5)



% hp(2) = nexttile([1,2]);
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% hold on
% plot((hysteresis(2).x) - 0.25, hysteresis(2).Re_theta, 'linewidth', lw, 'color', 'red', 'displayname', 'Phase 2')
% plot((hysteresis(4).x + 0.5) - 0.25, hysteresis(4).Re_theta, 'linewidth', lw, 'color', 'blue', 'linestyle', '-', 'Displayname', 'Phase 4 Raw')
% plot((hysteresis(4).x + 0.5) - 0.25, phase_4_Re_theta_shifted_scaled, 'linewidth', lw, 'color', 'green', 'Displayname', 'Phase 4 Shifted + Scaled')
% 
% % Plot reference wave profile
% Re_theta_reference_scale = 1E5;
% Re_theta_reference_offset = mean([Re_theta_LL, Re_theta_UL]);
% reference = plot((integral.(caze)(1).wave.x  - mean(integral.(caze)(1).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, (Re_theta_reference_scale * integral.(caze)(1).wave.wave_profile) + Re_theta_reference_offset, ...
%                            'color', 'black', 'linewidth', lw, 'linestyle', '--', 'handlevisibility', 'off');
% reference.Color(4) = 0.5;
% 
% hold off
% yline(hysteresis(2).Re_theta(1))
% yline(hysteresis(2).Re_theta(end))
% ylabel('$Re_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% xlabel('$x / \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylim([Re_theta_LL, Re_theta_UL])
% yticks(Re_theta_ticks)
% xlim([-0.55, 0.55])
% xticks(-0.5:0.25:0.5)

% Align axes
% linkaxes(hp, 'x')
linkaxes(hh, 'xy')






%% Loop through and plot: Wind speed per subplot


wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
smoothing_kernel = 16;

wind_speeds = {'WT4', 'WT6', 'WT8'};
waves = {'A', 'B', 'C', 'D'};

figure('color', 'white')
tiledlayout(1, 3)

for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};

    if ismember(wind_speed(end), {'4'})
        u_inf = 2.4181;
    elseif ismember(wind_speed(end), {'6'})
        u_inf = 3.8709;
    elseif ismember(wind_speed(end), {'8'})
        u_inf = 5.4289;
    end

    h(s) = nexttile;
    hold on
    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        % Collect data
        for phase = 2:2:4
        
            % Skin friction
            cf = skin_friction_profiles.(caze)(phase).cf;
            cf = smoothdata(cf, 'movmean', smoothing_kernel);
            
            % Re_{\theta}
            theta = integral.(caze)(phase).filtered.momentum;
            theta = smoothdata(theta, 'movmean', smoothing_kernel);
            Re_theta = (u_inf * theta) / nu;
            
            centered_normalized_x = (integral.(caze)(phase).x - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength;
        
            min_array_size = min([length(cf), length(Re_theta)]);
            centered_normalized_x = centered_normalized_x(1:min_array_size);
            cf = cf(1:min_array_size);
            Re_theta = Re_theta(1:min_array_size);
        
            % Crop data
            x_mask = centered_normalized_x;
            x_mask(abs(centered_normalized_x) > 0.25) = nan;
            x_mask(~isnan(x_mask)) = 1;
        
            cf = cf .* x_mask;
            Re_theta = Re_theta .* x_mask;
            masked_x = centered_normalized_x .* x_mask;
        
            % Save data to handle stitching after
            hysteresis(phase).cf = cf(~isnan(cf));
            hysteresis(phase).Re_theta = Re_theta(~isnan(Re_theta));
            hysteresis(phase).x = masked_x(~isnan(masked_x));
        end

        % Stitch and scale phases
        phase_2_cf = hysteresis(2).cf;
        phase_2_Re_theta = hysteresis(2).Re_theta;
        
        phase_4_cf = hysteresis(4).cf;
        phase_4_Re_theta = hysteresis(4).Re_theta;
        
        % Compute shifts
        phase_4_Re_theta_shift = phase_4_Re_theta(1) - phase_2_Re_theta(end);
        phase_4_Re_theta_shifted = phase_4_Re_theta - phase_4_Re_theta_shift;
        
        phase_4_cf_shift = phase_4_cf(1) - phase_2_cf(end);
        phase_4_cf_shifted = phase_4_cf - phase_4_cf_shift;
        
        % Compute scaling
        phase_4_Re_theta_scale = phase_2_Re_theta(1) / phase_4_Re_theta_shifted(end);
        phase_4_cf_scale = phase_2_cf(1) / phase_4_cf_shifted(end);
        
        phase_4_Re_theta_fading = linspace(1, phase_4_Re_theta_scale, length(phase_4_Re_theta_shifted));
        phase_4_cf_fading = linspace(1, phase_4_cf_scale, length(phase_4_cf_shifted));
        
        % Shifted and scaled data
        phase_4_Re_theta_shifted_scaled = phase_4_Re_theta_fading .* phase_4_Re_theta_shifted;
        phase_4_cf_shifted_scaled = phase_4_cf_fading .* phase_4_cf_shifted;

        % Save data back to array to make seperate loop for plotting
        clean_data.(caze).cf_2 = phase_2_cf;
        clean_data.(caze).cf_4 = phase_4_cf_shifted_scaled;

        clean_data.(caze).Re_theta_2 = phase_2_Re_theta;
        clean_data.(caze).Re_theta_4 = phase_4_Re_theta_shifted_scaled;



        % Plot
        scatter(phase_2_Re_theta, phase_2_cf, 40, 'filled', 'MarkerFaceColor', wave_colors{w})
        scatter(phase_4_Re_theta_shifted_scaled, phase_4_cf_shifted_scaled, 40, 'MarkerEdgeColor', wave_colors{w})

        % % Try plotting a line between peak and trough
        % plot([phase_2_Re_theta(1), phase_2_Re_theta(end)], [phase_2_cf(1), phase_2_cf(end)], ...
        %      'linewidth', 2, 'color', wave_colors{w})

        clear hysteresis phase_4_Re_theta_shifted_scaled phase_4_cf_shifted_scaled 
    end
    hold off

end

linkaxes(h, 'y')




%% Loop through and plot: Wave speed per subplot


wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
smoothing_kernel = 18;

wind_speeds = {'WT4', 'WT6', 'WT8'};
waves = {'A', 'B', 'C', 'D'};

wind_speed_alphas = [0.1, 0.5, 1];

phase_markers = {'o', '^'};


clc; close all
figure('color', 'white')
tiledlayout(1, 4)

for w = 1:length(waves)
    wave = waves{w};

    h(w) = nexttile;
    hold on
    % for s = 1:length(wind_speeds)
    for s = 3
        wind_speed = wind_speeds{s};
    
        if ismember(wind_speed(end), {'4'})
            u_inf = 2.4181;
        elseif ismember(wind_speed(end), {'6'})
            u_inf = 3.8709;
        elseif ismember(wind_speed(end), {'8'})
            u_inf = 5.4289;
        end
    
        caze = strcat(wind_speed, '_WV', wave, '_AG0');


        % Plot saved data
        phase_2_Re_theta = clean_data.(caze).Re_theta_2;
        phase_2_cf = clean_data.(caze).cf_2;

        phase_4_Re_theta_shifted_scaled = clean_data.(caze).Re_theta_4;
        phase_4_cf_shifted_scaled = clean_data.(caze).cf_4;




        % Plot
        scatter(phase_2_Re_theta, phase_2_cf, 40, phase_markers{1}, 'filled', 'MarkerFaceColor', wave_colors{w}, ...
                'MarkerFaceAlpha', wind_speed_alphas(s), 'MarkerEdgeAlpha', wind_speed_alphas(s))

        scatter(phase_4_Re_theta_shifted_scaled, phase_4_cf_shifted_scaled, 40, phase_markers{2}, 'filled', 'MarkerFaceColor', wave_colors{w}, ...
                'MarkerFaceAlpha', wind_speed_alphas(s), 'MarkerEdgeAlpha', wind_speed_alphas(s))

        % Try plotting a line between peak and trough
        plot([phase_2_Re_theta(1), phase_2_Re_theta(end)], [phase_2_cf(1), phase_2_cf(end)], ...
             'linewidth', 2, 'color', wave_colors{w})

        clear hysteresis phase_4_Re_theta_shifted_scaled phase_4_cf_shifted_scaled 
    end
    hold off

end

linkaxes(h, 'xy')
xlim([2000, 8000])
ylim([1E-3 19E-3])


%% Compute \Delta C_f and \Delta Re_{\theta}

cf_ranges = nan(1,4);
Re_theta_ranges = nan(1,4);

wind_speeds = {'WT4', 'WT6', 'WT8'};
waves = {'A', 'B', 'C', 'D'};

for w = 1:length(waves)
    wave = waves{w};

    h(w) = nexttile;
    hold on
    for s = 3
        wind_speed = wind_speeds{s};
    
        if ismember(wind_speed(end), {'4'})
            u_inf = 2.4181;
        elseif ismember(wind_speed(end), {'6'})
            u_inf = 3.8709;
        elseif ismember(wind_speed(end), {'8'})
            u_inf = 5.4289;
        end
    
        caze = strcat(wind_speed, '_WV', wave, '_AG0');


        % Plot saved data
        phase_2_Re_theta = clean_data.(caze).Re_theta_2;
        phase_2_cf = clean_data.(caze).cf_2;

        phase_4_Re_theta_shifted_scaled = clean_data.(caze).Re_theta_4;
        phase_4_cf_shifted_scaled = clean_data.(caze).cf_4;


        cf_range = range([phase_2_cf, phase_4_cf_shifted_scaled]);
        Re_theta_range = range([phase_2_Re_theta, phase_4_Re_theta_shifted_scaled]);

        cf_ranges(w) = cf_range;
        Re_theta_ranges(w) = Re_theta_range;


        clear hysteresis phase_2_Re_theta phase_2_cf phase_4_Re_theta_shifted_scaled phase_4_cf_shifted_scaled 
    end
end

clc;close all
figure('color', 'white')
bar(cf_ranges)

figure('color', 'white')
bar(Re_theta_ranges)






















%% Functions

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
