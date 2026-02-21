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

figure_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/pdf_test4';

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



%% Plot Cf profiles to check

tickFontSize = 14;
labelFontSize = 16;
legendFontSize = 12;
titleFontSize = 18;

wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
linestyles = {'-.', '--', '-'};
lw = 2;
buff = 5;
phase_4_Re_theta_scale = 1.25;

smoothing_kernel = 32;

x_crop = 0.25;

close all; 
figure('color', 'white');
tiledlayout(2,1)
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)

% Loop through phases
for phase = 2:2:4

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
    
        for w = 1:length(waves)
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
xlabel('$\xi / \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$C_f$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xticks(-1:0.25:1)
clc;






%% Plot to check momentum thickness curves


parameter = 'momentum';
version = 'filtered';
linestyles = {'-.', '--', '-'};
lw = 2;
alpha = 0.5;

x_crop = 0.25;

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

for phase = 2:2:4
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
    
        for w = 1:length(waves)
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
xlabel('$\xi / \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(vertLabel, 'interpreter', 'latex', 'fontsize', labelFontSize)
xticks(-1:0.25:1)







%% Plotting cropped phases 2 and 4 together + saving data

wind_speed = 'WT6';
wave = 'D';
caze = strcat(wind_speed, '_WV', wave, '_AG0');

if ismember(wind_speed(end), {'4'})
    u_inf = 2.4181;
elseif ismember(wind_speed(end), {'6'})
    u_inf = 3.8709;
elseif ismember(wind_speed(end), {'8'})
    u_inf = 5.4289;
end

smoothing_kernel = 16;

clc; close all
figure('color', 'white')

hold on
for phase = 2:2:4
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

    % Plot
    if phase == 2
        color = 'red';
    elseif phase == 4
        color = 'blue';
    end

    plot(Re_theta, cf, 'linewidth', 1, 'color', color)
    scatter(Re_theta, cf, 40, 'filled', 'MarkerFaceColor', color)

end
hold off

xlabel('Re_{\theta}')
ylabel('Skin Friction')



%%% Stitching the two phases with offset and scaling
% Data per phase
phase_2_cf = hysteresis(2).cf;
phase_2_Re_theta = hysteresis(2).Re_theta;

phase_4_cf = hysteresis(4).cf;
phase_4_Re_theta = hysteresis(4).Re_theta;

% Compute shifts
phase_4_Re_theta_shift = phase_4_Re_theta(1) - phase_2_Re_theta(end);
phase_4_Re_theta_shifted = phase_4_Re_theta - phase_4_Re_theta_shift;

phase_4_cf_shift = phase_4_cf(1) - phase_2_cf(end);
phase_4_cf_shifted = phase_4_cf - phase_4_cf_shift;

% % Compute scaling
% phase_4_Re_theta_scale = phase_2_Re_theta(1) / phase_4_Re_theta_shifted(end);
% phase_4_cf_scale = phase_2_cf(1) / phase_4_cf_shifted(end);
% 
% phase_4_Re_theta_fading = linspace(1, phase_4_Re_theta_scale, length(phase_4_Re_theta_shifted));
% phase_4_cf_fading = linspace(1, phase_4_cf_scale, length(phase_4_cf_shifted));
% 
% % Shifted and scaled data
% phase_4_Re_theta_shifted_scaled = phase_4_Re_theta_fading .* phase_4_Re_theta_shifted;
% phase_4_cf_shifted_scaled = phase_4_cf_fading .* phase_4_cf_shifted;



%% Plot hysteresis

tickFontSize = 14;
labelFontSize = 16;
legendFontSize = 12;
titleFontSize = 18;

Cf_LL = 2E-3;
Cf_UL = 20E-3;

Re_theta_LL = 2800;
Re_theta_UL = 4600;

Cf_ticks = Cf_LL + 1E-3:0.005:Cf_UL;
Re_theta_ticks = Re_theta_LL:500:Re_theta_UL;

lw = 3;
sz = 20;

clc; close all
ax = figure('color', 'white', 'units', 'centimeters', 'position', [5,5,29,10]);
fig = tiledlayout(1,2, 'tilespacing', 'tight', 'padding', 'tight');

hh(1) = nexttile;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
hold on
scatter(phase_2_Re_theta, phase_2_cf, sz, 'filled', 'MarkerFaceColor', 'red')
plot(phase_2_Re_theta, phase_2_cf, 'Color', 'red', 'linewidth', 2)
scatter(phase_4_Re_theta, phase_4_cf, sz, 'filled', 'MarkerFaceColor', 'blue')
plot(phase_4_Re_theta, phase_4_cf, 'Color', 'blue', 'linewidth', 2)
hold off
% title('Raw Data')
% xlabel('$Re_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylabel('$C_f$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xticks(Re_theta_ticks)
yticks(Cf_ticks)

% Make y axis easier to read
ax = gca;
ax.YAxis.Exponent = -3;          % shows ×10^{-6} once
ax.YAxis.TickLabelFormat = '%2.0f';  % two digits on ticks



hh(2) = nexttile;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
hold on
scatter(phase_2_Re_theta, phase_2_cf, sz, 'filled', 'MarkerFaceColor', 'red')
plot(phase_2_Re_theta, phase_2_cf, 'Color', 'red', 'linewidth', 2)
scatter(phase_4_Re_theta_shifted, phase_4_cf_shifted, sz, 'filled', 'MarkerFaceColor', 'green')
plot(phase_4_Re_theta_shifted, phase_4_cf_shifted, 'Color', 'green', 'linewidth', 2)
hold off
% title('Enforced Boundary Condition')
% xlabel('$Re_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlim([Re_theta_LL, Re_theta_UL])
ylim([Cf_LL, Cf_UL])
xticks(Re_theta_ticks)
yticks(Cf_ticks)


linkaxes(hh, 'xy')
xlabel(fig, '$Re_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(fig, '$C_f$', 'interpreter', 'latex', 'fontsize', labelFontSize)


% Make y axis easier to read
ax = gca;
ax.YAxis.Exponent = -3;          % shows ×10^{-6} once
ax.YAxis.TickLabelFormat = '%2.0f';  % two digits on ticks

% Print shifts
fprintf('Re_theta shift and scaling\n\n')
fprintf('Shift applied: %.4f\n', phase_4_Re_theta_shift);

fprintf('Cf shift and scaling\n')
fprintf('Shift applied: %.4f\n', phase_4_cf_shift);


% Save figure
% pause(3)
% figure_name = 'HysteresisComparison.pdf';
% exportgraphics(fig, fullfile(figure_folder, 'Hysteresis', figure_name), 'Resolution', 600, 'ContentType', 'image');
% close all
% fprintf('Generated figure: %s\n\n', figure_name)




%% Plot Profiles

annotationFontSize = 12;

lw = 2;

ax = figure('color', 'white', 'units', 'centimeters', 'position', [5,5,29,12]);
fig = tiledlayout(2,1,"TileSpacing",'tight', 'padding', 'tight');

% C_f profiles
hp(1) = nexttile();
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
hold on
plot((hysteresis(2).x) + 0.25, hysteresis(2).cf, ...
     'linewidth', lw, 'color', 'red', 'displayname', '$\varphi = \lambda/4$')

plot((hysteresis(4).x + 0.5) + 0.25, hysteresis(4).cf, ...
     'linewidth', lw, 'color', 'blue', 'linestyle', '-', 'Displayname', '$\varphi = 3 \lambda/4$')

plot((hysteresis(4).x + 0.5) + 0.25, phase_4_cf_shifted, ...
     'linewidth', lw, 'color', 'green', 'Displayname', '$\varphi = 3 \lambda/4 + \Delta$')

% Plot reference wave profile
Cf_reference_scale = 1;
Cf_reference_offset = mean([Cf_LL, Cf_UL]);
reference = plot(0.5 + (integral.(caze)(1).wave.x  - mean(integral.(caze)(1).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, (Cf_reference_scale * integral.(caze)(1).wave.wave_profile) + Cf_reference_offset, ...
                           'color', 'black', 'linewidth', lw, 'linestyle', '--', 'handlevisibility', 'off');
reference.Color(4) = 0.5;

hold off
xline(0.5, 'color', 'black', 'linestyle', '-', 'linewidth', 2, 'HandleVisibility', 'off')
ylabel('$C_f$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylim([Cf_LL, Cf_UL])
yticks(Cf_ticks)
xlim([0, 1])
ax = gca;
ax.XAxis.Visible = 'off';

% Make y axis easier to read
ax.YAxis.Exponent = -3;          % shows ×10^{-6} once
ax.YAxis.TickLabelFormat = '%2.0f';  % two digits on ticks

% Show what the \Delta is
exponent = -3;
mantissa = phase_4_cf_shift / 10^exponent;
text(0.51, 0.5, sprintf('$\\Delta C_f = %0.1f \\times 10^{%d}$', mantissa, exponent), ...
    'Units', 'normalized', ...
    'Interpreter', 'latex', ...
    'FontSize', annotationFontSize)






% Re_{\theta} profiles
hp(2) = nexttile();
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
hold on
plot((hysteresis(2).x) + 0.25, hysteresis(2).Re_theta, ...
     'linewidth', lw, 'color', 'red', 'displayname', '$\varphi = \lambda/4$')

plot((hysteresis(4).x + 0.5) + 0.25, hysteresis(4).Re_theta, ...
     'linewidth', lw, 'color', 'blue', 'linestyle', '-', 'Displayname', '$\varphi = 3 \lambda/4$')

plot((hysteresis(4).x + 0.5) + 0.25, phase_4_Re_theta_shifted, ...
     'linewidth', lw, 'color', 'green', 'Displayname', '$\varphi = 3 \lambda/4 + \Delta$')

% Plot reference wave profile
Re_theta_reference_scale = 1.25E5;
Re_theta_reference_offset = mean([Re_theta_LL, Re_theta_UL]);
reference = plot(0.5 + (integral.(caze)(1).wave.x  - mean(integral.(caze)(1).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, (Re_theta_reference_scale * integral.(caze)(1).wave.wave_profile) + Re_theta_reference_offset, ...
                           'color', 'black', 'linewidth', lw, 'linestyle', '--', 'handlevisibility', 'on', 'Displayname', '$\eta$');
reference.Color(4) = 0.5;

hold off
xline(0.5, 'color', 'black', 'linestyle', '-', 'linewidth', 2, 'HandleVisibility', 'off')
ylabel('$Re_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlabel('$\xi / \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylim([Re_theta_LL, Re_theta_UL])
yticks(Re_theta_ticks)
xlim([0, 1])
xticks(0:0.25:1)

% Show what the \Delta is
text(0.51, 0.52, sprintf('$\\Delta Re_{\\theta} = %3.2f$', phase_4_Re_theta_shift), ...
    'Units', 'normalized', ...
    'Interpreter', 'latex', ...
    'FontSize', annotationFontSize)

leg = legend('Orientation', 'vertical', 'box', 'off', 'fontsize', legendFontSize, 'interpreter', 'latex');
leg.Layout.Tile = 'east';

% Align axes
linkaxes(hp, 'x')


% Save figure
% pause(3)
% figure_name = 'HysteresisStitchExamples.pdf';
% exportgraphics(fig, fullfile(figure_folder, 'Hysteresis', figure_name), 'Resolution', 600, 'ContentType', 'image');
% close all
% fprintf('Generated figure: %s\n\n', figure_name)





%% Loop through and plot: Wind speed per subplot NO SCALING


% wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
% smoothing_kernel = 16;
% 
% wind_speeds = {'WT4', 'WT6', 'WT8'};
% waves = {'A', 'B', 'C', 'D'};
% 
% figure('color', 'white')
% tiledlayout(1, 3)
% 
% for s = 1:length(wind_speeds)
%     wind_speed = wind_speeds{s};
% 
%     if ismember(wind_speed(end), {'4'})
%         u_inf = 2.4181;
%     elseif ismember(wind_speed(end), {'6'})
%         u_inf = 3.8709;
%     elseif ismember(wind_speed(end), {'8'})
%         u_inf = 5.4289;
%     end
% 
%     h(s) = nexttile;
%     hold on
%     for w = 1:length(waves)
%         wave = waves{w};
%         caze = strcat(wind_speed, '_WV', wave, '_AG0');
% 
%         % Collect data
%         for phase = 2:2:4
% 
%             % Skin friction
%             cf = skin_friction_profiles.(caze)(phase).cf;
%             cf = smoothdata(cf, 'movmean', smoothing_kernel);
% 
%             % Re_{\theta}
%             theta = integral.(caze)(phase).filtered.momentum;
%             theta = smoothdata(theta, 'movmean', smoothing_kernel);
%             Re_theta = (u_inf * theta) / nu;
% 
%             centered_normalized_x = (integral.(caze)(phase).x - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength;
% 
%             min_array_size = min([length(cf), length(Re_theta)]);
%             centered_normalized_x = centered_normalized_x(1:min_array_size);
%             cf = cf(1:min_array_size);
%             Re_theta = Re_theta(1:min_array_size);
% 
%             % Crop data
%             x_mask = centered_normalized_x;
%             x_mask(abs(centered_normalized_x) > 0.25) = nan;
%             x_mask(~isnan(x_mask)) = 1;
% 
%             cf = cf .* x_mask;
%             Re_theta = Re_theta .* x_mask;
%             masked_x = centered_normalized_x .* x_mask;
% 
%             % Save data to handle stitching after
%             hysteresis(phase).cf = cf(~isnan(cf));
%             hysteresis(phase).Re_theta = Re_theta(~isnan(Re_theta));
%             hysteresis(phase).x = masked_x(~isnan(masked_x));
%         end
% 
%         % Stitch and scale phases
%         phase_2_cf = hysteresis(2).cf;
%         phase_2_Re_theta = hysteresis(2).Re_theta;
% 
%         phase_4_cf = hysteresis(4).cf;
%         phase_4_Re_theta = hysteresis(4).Re_theta;
% 
%         % Compute shifts
%         phase_4_Re_theta_shift = phase_4_Re_theta(1) - phase_2_Re_theta(end);
%         phase_4_Re_theta_shifted = phase_4_Re_theta - phase_4_Re_theta_shift;
% 
%         phase_4_cf_shift = phase_4_cf(1) - phase_2_cf(end);
%         phase_4_cf_shifted = phase_4_cf - phase_4_cf_shift;
% 
%         % % Compute scaling
%         % phase_4_Re_theta_scale = phase_2_Re_theta(1) / phase_4_Re_theta_shifted(end);
%         % phase_4_cf_scale = phase_2_cf(1) / phase_4_cf_shifted(end);
%         % 
%         % phase_4_Re_theta_fading = linspace(1, phase_4_Re_theta_scale, length(phase_4_Re_theta_shifted));
%         % phase_4_cf_fading = linspace(1, phase_4_cf_scale, length(phase_4_cf_shifted));
%         % 
%         % % Shifted and scaled data
%         % phase_4_Re_theta_shifted_scaled = phase_4_Re_theta_fading .* phase_4_Re_theta_shifted;
%         % phase_4_cf_shifted_scaled = phase_4_cf_fading .* phase_4_cf_shifted;
% 
%         % Save data back to array to make seperate loop for plotting
%         clean_data.(caze).cf_2 = phase_2_cf;
%         clean_data.(caze).cf_4 = phase_4_cf_shifted;
% 
%         clean_data.(caze).Re_theta_2 = phase_2_Re_theta;
%         clean_data.(caze).Re_theta_4 = phase_4_Re_theta_shifted;
% 
% 
% 
%         % Plot
%         scatter(phase_2_Re_theta, phase_2_cf, 40, 'filled', 'MarkerFaceColor', wave_colors{w})
%         scatter(phase_4_Re_theta_shifted, phase_4_cf_shifted, 40, 'MarkerEdgeColor', wave_colors{w})
% 
%         % % Try plotting a line between peak and trough
%         % plot([phase_2_Re_theta(1), phase_2_Re_theta(end)], [phase_2_cf(1), phase_2_cf(end)], ...
%         %      'linewidth', 2, 'color', wave_colors{w})
% 
%         clear hysteresis phase_4_Re_theta_shifted phase_4_cf_shifted 
%     end
%     hold off
% 
% end
% 
% linkaxes(h, 'y')




%% Loop through and plot: Wave speed per subplot

tickFontSize = 14;
labelFontSize = 16;
legendFontSize = 14;
titleFontSize = 18;

wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
smoothing_kernel = 24;

wind_speeds = {'WT4', 'WT6', 'WT8'};
waves = {'A', 'B', 'C', 'D'};

wind_speed_alphas = [1, 1, 1];

phase_markers = {'o', '^'};

sz = 20;

clc; close all
ax = figure('color', 'white', 'units', 'centimeters', 'position', [3,5,34,10]);
t = tiledlayout(1, 4, 'Padding', 'tight');

for w = 1:length(waves)
    wave = waves{w};
    h(w) = nexttile;
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
    hold on

    for s = 2
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

        phase_4_Re_theta_shifted = clean_data.(caze).Re_theta_4;
        phase_4_cf_shifted = clean_data.(caze).cf_4;

        % Plot
        scatter(phase_2_Re_theta, phase_2_cf, sz, phase_markers{1}, 'filled', 'MarkerFaceColor', wave_colors{w}, ...
                'MarkerFaceAlpha', wind_speed_alphas(s), 'MarkerEdgeAlpha', wind_speed_alphas(s), 'HandleVisibility', 'off')
        plot(phase_2_Re_theta, phase_2_cf, 'color', wave_colors{w}, 'linewidth', 2, 'HandleVisibility', 'off')

        scatter(phase_4_Re_theta_shifted, phase_4_cf_shifted, sz, phase_markers{1}, 'filled', 'MarkerFaceColor', wave_colors{w}, ...
                'MarkerFaceAlpha', wind_speed_alphas(s), 'MarkerEdgeAlpha', wind_speed_alphas(s), 'HandleVisibility', 'off')
        plot(phase_4_Re_theta_shifted, phase_4_cf_shifted, 'color', wave_colors{w}, 'linewidth', 2, 'HandleVisibility', 'off')

        % Mark the starting point
        scatter(phase_2_Re_theta(1), phase_2_cf(1), 2*sz, '^', 'filled', 'MarkerFaceColor', 'white', ...
                'MarkerEdgeColor', 'black', 'LineWidth', 1.5', 'HandleVisibility', 'off')

        % Mark where peak is
        scatter(phase_4_Re_theta_shifted(1), phase_4_cf_shifted(1), 2*sz, 'o', 'filled', 'MarkerFaceColor', 'white', ...
                'MarkerEdgeColor', 'black', 'LineWidth', 1.5', 'HandleVisibility', 'off')

        % Mark where loop ends
        scatter(phase_4_Re_theta_shifted(end), phase_4_cf_shifted(end), 2*sz, 'v', 'filled', 'MarkerFaceColor', 'white', ...
                'MarkerEdgeColor', 'black', 'LineWidth', 1.5', 'HandleVisibility', 'off')

        clear hysteresis phase_4_Re_theta_shifted phase_4_cf_shifted 

        if w == 1
            ylabel('$C_f$', 'interpreter', 'latex', 'fontsize', labelFontSize)
        end

    end
    hold off
    
    % Make y axis easier to read
    ax = gca;
    ax.YAxis.Exponent = -3;          % shows ×10^{-6} once
    ax.YAxis.TickLabelFormat = '%2.0f';  % two digits on ticks

end

% Make a legend
hold on

% What colors mean
for w = 1:length(waves)
    disp(w)
    wave = waves{w};
    label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses.(wave));
    plot(nan, nan, 'color', wave_colors{w}, 'LineWidth', 3, 'DisplayName', label)
end

% Legend white space
plot(nan, nan, 'color', 'white', 'displayname', '')

% Legend for keypoint markers
% Mark the starting point
scatter(nan, nan, 2*sz, '^', 'filled', 'MarkerFaceColor', 'white', ...
        'MarkerEdgeColor', 'black', 'LineWidth', 1.5', 'displayname', '$\xi / \lambda = 0$')

% Mark the crest
scatter(nan, nan, 2*sz, 'o', 'filled', 'MarkerFaceColor', 'white', ...
        'MarkerEdgeColor', 'black', 'LineWidth', 1.5', 'displayname', '$\xi / \lambda = 0.5$')

% Mark the ending point
scatter(nan, nan, 2*sz, 'v', 'filled', 'MarkerFaceColor', 'white', ...
        'MarkerEdgeColor', 'black', 'LineWidth', 1.5', 'displayname', '$\xi / \lambda = 1$')
hold off
legend('interpreter', 'latex', 'box', 'off', 'location', 'eastoutside', 'fontsize', legendFontSize)


xlabel(t, '$Re_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)

linkaxes(h, 'xy')
xlim([2000, 6200])
ylim([1E-3 20E-3])


% Save figure
% pause(3)
% figure_name = strcat('Hysteresis_', wind_speeds{s}, '.pdf');
% exportgraphics(t, fullfile(figure_folder, 'Hysteresis', figure_name), 'resolution', 600, 'ContentType', 'image')
% close all




%% Compute Closure and ranges for \Delta C_f and \Delta Re_{\theta}


% Ranges
cf_ranges = nan(3,4);
Re_theta_ranges = nan(3,4);

% Closure
cf_closure = nan(3,4);
Re_theta_closure = nan(3,4);


wind_speeds = {'WT4', 'WT6', 'WT8'};
waves = {'A', 'B', 'C', 'D'};

for w = 1:length(waves)
    wave = waves{w};
    for s = 1:3
        wind_speed = wind_speeds{s};
    
        if ismember(wind_speed(end), {'4'})
            u_inf = 2.4181;
        elseif ismember(wind_speed(end), {'6'})
            u_inf = 3.8709;
        elseif ismember(wind_speed(end), {'8'})
            u_inf = 5.4289;
        end
    
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        % Get wavelengths
        wavelength = wavelengths.(wave);

        % Plot saved data
        phase_2_Re_theta = clean_data.(caze).Re_theta_2;
        phase_2_cf = clean_data.(caze).cf_2;

        phase_4_Re_theta_shifted = clean_data.(caze).Re_theta_4;
        phase_4_cf_shifted = clean_data.(caze).cf_4;

        % Compute + save range of data
        cf_range = range([phase_2_cf, phase_4_cf_shifted]);
        Re_theta_range = range([phase_2_Re_theta, phase_4_Re_theta_shifted]);

        cf_ranges(s,w) = cf_range;
        Re_theta_ranges(s,w) = Re_theta_range;

        % Compute + save closure
        cf_closure(s,w) = (phase_4_cf_shifted(end) - phase_2_cf(1));
        Re_theta_closure(s,w) = (phase_4_Re_theta_shifted(end) - phase_2_Re_theta(1));

        % cf_closure(s,w) = (phase_4_cf_shifted(end) - phase_2_cf(1)) / wavelength;
        % Re_theta_closure(s,w) = (phase_4_Re_theta_shifted(end) - phase_2_Re_theta(1)) / wavelength;

        clear hysteresis phase_2_Re_theta phase_2_cf phase_4_Re_theta_shifted phase_4_cf_shifted 
    end
end



%% Plot range of Cf against steepness and range of Re_{\theta} against wavelength

tickFontSize = 14;
labelFontSize = 16;
legendFontSize = 14;
titleFontSize = 18;

clc;
wind_speed_colors = {'#0075F2', '#FF8C42', '#D30C7B'};

% Range of Cf vs steepness
ax = figure('color', 'white', 'units', 'centimeters', 'position', [5,5,30,10]);
t = tiledlayout(1,2, 'padding', 'tight');

sz = 60;
lw = 3;

nexttile()
hold on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
for s = 1:3
    for w = 1:length(waves)
        wave = waves{w};
        wind_speed = wind_speeds{s};
    
        if ismember(wind_speed(end), {'4'})
            u_inf = 2.4181;
        elseif ismember(wind_speed(end), {'6'})
            u_inf = 3.8709;
        elseif ismember(wind_speed(end), {'8'})
            u_inf = 5.4289;
        end

        % Get wavelengths
        wavelength = wavelengths.(wave);
        amplitude = amplitudes.(wave);
        steepness = (2 * pi * amplitude) / wavelength;

        % Save steepnesses to plot a line
        steepnessez(w) = steepness;

        % Plot
        if w == 1 
            vis = 'on';
        else
            vis = 'off';
        end

        label = sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf);
        scatter(steepness, cf_ranges(s,w), sz, 'filled', 'markerfacecolor', wind_speed_colors{s}, ...
                'HandleVisibility', vis, 'displayname', label)
    end
    [sorted_steepnesses, sortIdx] = sort(steepnessez,'ascend');
    plot(sorted_steepnesses, cf_ranges(s, sortIdx), 'color', wind_speed_colors{s}, 'linewidth', lw, 'HandleVisibility', 'off')

end
hold off
leg = legend('Orientation', 'vertical', 'box', 'off', ...
             'fontsize', legendFontSize, 'interpreter', 'latex', 'location', 'northwest');

xlabel('$ak$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$\max \left( C_f \right) - \min \left( C_f \right)$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlim([0.15, 0.32])
ylim([0, 16E-3])
yticks(0:4E-3:16E-3)
% Make y axis easier to read
ax = gca;
ax.YAxis.Exponent = -3;          % shows ×10^{-6} once
ax.YAxis.TickLabelFormat = '%2.0f';  % two digits on ticks


% Range of Re_{\theta} vs wavelength
nexttile()
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
hold on
for s = 1:3
    for w = 1:length(waves)
        wave = waves{w};
        wind_speed = wind_speeds{s};
    
        if ismember(wind_speed(end), {'4'})
            u_inf = 2.4181;
        elseif ismember(wind_speed(end), {'6'})
            u_inf = 3.8709;
        elseif ismember(wind_speed(end), {'8'})
            u_inf = 5.4289;
        end

        % Get wavelengths
        wavelength = wavelengths.(wave);
        amplitude = amplitudes.(wave);
        steepness = (2 * pi * amplitude) / wavelength;

        % Save wavelengths to plot a line
        wvlenths(w) = wavelength;

        % Plot
        if w == 1 
            vis = 'on';
        else
            vis = 'off';
        end

        label = sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf);
        scatter(wavelength, Re_theta_ranges(s,w), sz, 'filled', 'markerfacecolor', wind_speed_colors{s}, ...
                'displayname', label, 'HandleVisibility', vis)
    end
    [sorted_wavelengths, sortIdx] = sort(wvlenths,'ascend');
    plot(sorted_wavelengths, Re_theta_ranges(s, sortIdx), 'color', wind_speed_colors{s}, 'linewidth', lw, 'HandleVisibility', 'off')

    % Try fitting a power law
    p = polyfit(log(sorted_wavelengths), log(Re_theta_ranges(s, sortIdx)), 1);
    m = p(1);
    b = exp(p(2));
    fprintf('b = %.4f,   m = %.4f\n', b, m);
    fitz = b * sorted_wavelengths.^m;
    % plot(sorted_wavelengths, fitz, 'linestyle', '--', 'linewidth', lw, 'Color', wind_speed_colors{s});
end
hold off

xlabel('$\lambda$ [m]', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$\max \left( {Re}_{\theta} \right) - \min \left( {Re}_{\theta} \right)$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylim([0, 6000])
yticks(0:2000:6000)


% Save figure
% pause(3)
% figure_name = 'Hysteresis_Ranges.pdf';
% exportgraphics(t, fullfile(figure_folder, 'Hysteresis', figure_name), 'resolution', 600, 'ContentType', 'image')
% close all



%% Plot closure of Cf against steepness and range of Re_{\theta} against wavelength

tickFontSize = 14;
labelFontSize = 16;
legendFontSize = 14;
titleFontSize = 18;

wind_speed_colors = {'#0075F2', '#FF8C42', '#D30C7B'};

% Range of Cf vs steepness
ax = figure('color', 'white', 'units', 'centimeters', 'position', [5,5,30,10]);
t = tiledlayout(1,2,'padding', 'tight');

sz = 60;
lw = 3;

h(1) = nexttile();
hold on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
for s = 1:3
    for w = 1:length(waves)
        wave = waves{w};
        wind_speed = wind_speeds{s};
    
        if ismember(wind_speed(end), {'4'})
            u_inf = 2.4181;
        elseif ismember(wind_speed(end), {'6'})
            u_inf = 3.8709;
        elseif ismember(wind_speed(end), {'8'})
            u_inf = 5.4289;
        end

        % Get wavelengths
        wavelength = wavelengths.(wave);
        amplitude = amplitudes.(wave);
        steepness = (2 * pi * amplitude) / wavelength;

        % Save steepnesses to plot a line
        steepnessez(w) = steepness;

        % Plot
        if w == 1 
            vis = 'on';
        else
            vis = 'off';
        end

        label = sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf);
        scatter(steepness, cf_closure(s,w), sz, 'filled', 'markerfacecolor', wind_speed_colors{s}, ...
                'HandleVisibility', vis, 'displayname', label)
    end
    [sorted_steepnesses, sortIdx] = sort(steepnessez,'ascend');
    plot(sorted_steepnesses, cf_closure(s, sortIdx), 'color', wind_speed_colors{s}, 'linewidth', lw, 'HandleVisibility', 'off')

end
hold off
leg = legend('Orientation', 'vertical', 'box', 'off', ...
             'fontsize', legendFontSize, 'interpreter', 'latex', 'location', 'northwest');

xlabel('$ak$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$\Delta C_f $', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlim([0.15, 0.32])
ylim([-10E-3, 10E-3])
% yticks(0:4E-3:16E-3)
% Make y axis easier to read
ax = gca;
ax.YAxis.Exponent = -3;          % shows ×10^{-6} once
ax.YAxis.TickLabelFormat = '%2.0f';  % two digits on ticks


% Range of Re_{\theta} vs wavelength
h(2) = nexttile();
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
hold on
for s = 1:3
    for w = 1:length(waves)
        wave = waves{w};
        wind_speed = wind_speeds{s};
    
        if ismember(wind_speed(end), {'4'})
            u_inf = 2.4181;
        elseif ismember(wind_speed(end), {'6'})
            u_inf = 3.8709;
        elseif ismember(wind_speed(end), {'8'})
            u_inf = 5.4289;
        end

        % Get wavelengths
        wavelength = wavelengths.(wave);
        amplitude = amplitudes.(wave);
        steepness = (2 * pi * amplitude) / wavelength;

        % Save wavelengths to plot a line
        wvlenths(w) = wavelength;

        % Plot
        if w == 1 
            vis = 'on';
        else
            vis = 'off';
        end

        label = sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf);
        scatter(wavelength, Re_theta_closure(s,w), sz, 'filled', 'markerfacecolor', wind_speed_colors{s}, ...
                'displayname', label, 'HandleVisibility', vis)
    end
    [sorted_wavelengths, sortIdx] = sort(wvlenths,'ascend');
    plot(sorted_wavelengths, Re_theta_closure(s, sortIdx), 'color', wind_speed_colors{s}, 'linewidth', lw, 'HandleVisibility', 'off')
end
hold off

xlabel('$\lambda$ [m]', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$\Delta {Re}_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylim([0, 2000])


% Save figure
% pause(3)
% figure_name = 'Hysteresis_Closures.pdf';
% exportgraphics(t, fullfile(figure_folder, 'Hysteresis',  figure_name), 'resolution', 600, 'ContentType', 'image')
% close all







%% Plot \Delta {Re}_{\theta} against u_{wave} / u_{\infty}


frequencies.A = 1.96;
frequencies.B = 2.27;
frequencies.C = 3;
frequencies.D = 3.93;


% Range of Re_{\theta} vs wavelength
figure('color', 'white')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
hold on
for s = 1:3
    for w = 1:length(waves)
        wave = waves{w};
        wind_speed = wind_speeds{s};
    
        if ismember(wind_speed(end), {'4'})
            u_inf = 2.4181;
        elseif ismember(wind_speed(end), {'6'})
            u_inf = 3.8709;
        elseif ismember(wind_speed(end), {'8'})
            u_inf = 5.4289;
        end

        % Get wavelengths
        wavelength = wavelengths.(wave);
        amplitude = amplitudes.(wave);
        steepness = (2 * amplitude) / wavelength;

        wave_speed = wavelength * frequencies.(wave);
        relative_wave_speed = wave_speed / u_inf;

        % Save wavelengths to plot a line
        tmp(w) = relative_wave_speed;

        % Plot
        if w == 1 
            vis = 'on';
        else
            vis = 'off';
        end

        label = sprintf('$u_{\\infty} = %1.2f m/s$', u_inf);
        scatter(relative_wave_speed, Re_theta_closure(s,w), sz, 'filled', 'markerfacecolor', wind_speed_colors{s}, ...
                'displayname', label, 'HandleVisibility', vis)
    end
    [sorted_tmp, sortIdx] = sort(tmp,'ascend');
    plot(sorted_tmp, Re_theta_closure(s, sortIdx), 'color', wind_speed_colors{s}, 'linewidth', lw, 'HandleVisibility', 'off')
end
hold off

xlabel('$\lambda f / u_{\infty}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$\Delta {Re}_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)








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
