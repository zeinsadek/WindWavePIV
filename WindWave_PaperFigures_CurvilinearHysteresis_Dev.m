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


%% Plotting cropped phases 2 and 4 together + saving data

nu = 1.48E-5;
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


%% Plot Profiles

wave_transparency = 0.25;


% Colors
phase_2_color = '#FF8552';
phase_4_color = '#78CAD2';
phase_4_shifted_color = '#094074';

% Fontsizes
tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;
annotationFontSize = 6;
linewidth = 1;
sz = 5;

% Plotting
clc; close all
ax = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,7]);
fig = tiledlayout(2,2,"TileSpacing",'tight', 'padding', 'tight');

% C_f profiles
% ax1 = nexttile([1 2]);
ax1 = nexttile(1);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
hold on
plot((hysteresis(2).x) + 0.25, hysteresis(2).cf, ...
     'linewidth', linewidth, 'color', phase_2_color)

plot((hysteresis(4).x + 0.5) + 0.25, hysteresis(4).cf, ...
     'linewidth', linewidth, 'color', phase_4_color, 'linestyle', '-')

plot((hysteresis(4).x + 0.5) + 0.25, phase_4_cf_shifted, ...
     'linewidth', linewidth, 'color', phase_4_shifted_color')

% Plot reference wave profile
Cf_reference_scale = 0.5;
Cf_reference_offset = 2E-3;
ref_x = 0.5 + (integral.(caze)(1).wave.x  - mean(integral.(caze)(1).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength;
ref_wave = (Cf_reference_scale * integral.(caze)(1).wave.wave_profile) + Cf_reference_offset;
plot(ref_x, ref_wave, 'color', 'black', 'linewidth', linewidth, 'linestyle', '-', 'handlevisibility', 'off');

% Shaded region below wave
hFill = patch( ...
[ref_x, fliplr(ref_x)], ...
[ref_wave, -10 * ones(size(ref_wave))], ...
'k', ...
'FaceAlpha', wave_transparency, ...      
'EdgeColor', 'none', ...    
'HandleVisibility', 'off'); 
uistack(hFill, 'bottom')

% Plot details
hold off
xline(0.5, 'color', 'black', 'linestyle', '-', 'linewidth', linewidth, 'HandleVisibility', 'off')
ylabel('$C_f$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylim([-1.5E-3, 0.02])
yticks(0:5E-3:25E-3)
xlim([0, 1])
ax = gca;
ax.XAxis.Visible = 'off';

% Make y axis easier to read
ax.YAxis.Exponent = -3;
ax.YAxis.TickLabelFormat = '%2.0f';

% Show what the \delta is
exponent = -3;
mantissa = phase_4_cf_shift / 10^exponent;
text(0.51, 0.55, sprintf('$\\delta C_f = %0.1f \\times 10^{%d}$', -mantissa, exponent), ...
    'Units', 'normalized', ...
    'Interpreter', 'latex', ...
    'FontSize', annotationFontSize)




%%% Re_{\theta} profiles
% ax2 = nexttile([1 2]);
ax2 = nexttile(3);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
hold on
plot((hysteresis(2).x) + 0.25, hysteresis(2).Re_theta, ...
     'linewidth', linewidth, 'color', phase_2_color, 'displayname', '$\varphi = \lambda/4$')

plot((hysteresis(4).x + 0.5) + 0.25, hysteresis(4).Re_theta, ...
     'linewidth', linewidth, 'color', phase_4_color, 'linestyle', '-', 'Displayname', '$\varphi = 3 \lambda/4$')

plot((hysteresis(4).x + 0.5) + 0.25, phase_4_Re_theta_shifted, ...
     'linewidth', linewidth, 'color', phase_4_shifted_color', 'Displayname', '$\varphi = 3 \lambda/4 + \delta$')

% Plot reference wave profile
Re_theta_reference_scale = 7E4;
Re_theta_reference_offset = 2500;
ref_x = 0.5 + (integral.(caze)(1).wave.x  - mean(integral.(caze)(1).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength;
ref_wave = (Re_theta_reference_scale * integral.(caze)(1).wave.wave_profile) + Re_theta_reference_offset;
plot(ref_x, ref_wave, 'color', 'black', 'linewidth', linewidth, 'linestyle', '-', 'handlevisibility', 'off', 'Displayname', '$\eta$');

% Shaded region below wave
hFill = patch( ...
[ref_x, fliplr(ref_x)], ...
[ref_wave, -10 * ones(size(ref_wave))], ...
'k', ...
'FaceAlpha', wave_transparency, ...        
'EdgeColor', 'none', ...     
'HandleVisibility', 'off'); 
uistack(hFill, 'bottom')

% Plot details
hold off
xline(0.5, 'color', 'black', 'linestyle', '-', 'linewidth', linewidth, 'HandleVisibility', 'off')
ylabel('$Re_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlabel('$\xi / \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylim([2000, 5000])
yticks(3000:1000:5000)
xlim([0, 1])
xticks(0:0.25:1)

% Show what the \Delta is
text(0.24, 0.39, sprintf('$\\delta Re_{\\theta} = %3.2f$', -phase_4_Re_theta_shift), ...
    'Units', 'normalized', ...
    'Interpreter', 'latex', ...
    'FontSize', annotationFontSize)

% Align axes
linkaxes([ax1, ax2], 'x')







%%% Add example hysteresis plot
Cf_LL = 2E-3;
Cf_UL = 20E-3;
Re_theta_LL = 2800;
Re_theta_UL = 4600;
Cf_ticks = Cf_LL + 1E-3:0.005:Cf_UL;
Re_theta_ticks = Re_theta_LL:500:Re_theta_UL;


% Set up tile
ax3 = nexttile([2, 1]);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
hold on

% Plot first phase
scatter(phase_2_Re_theta, phase_2_cf, sz, 'filled', ...
        'MarkerFaceColor', phase_2_color, 'HandleVisibility', 'off')
plot(phase_2_Re_theta, phase_2_cf, 'Color', phase_2_color, ...
     'linewidth', linewidth, 'displayname', '$\varphi = \lambda/4$')

% Plot raw phase
scatter(phase_4_Re_theta, phase_4_cf, sz, 'filled', ...
        'MarkerFaceColor', phase_4_color, 'HandleVisibility', 'off', ...
        'MarkerFaceAlpha', 0.25)
P = plot(phase_4_Re_theta, phase_4_cf, 'Color', phase_4_color, ...
         'linewidth', linewidth, 'displayname', '$\varphi = 3 \lambda/4$');
P.Color(4) = 0.25;

% Plot shifted phase
scatter(phase_4_Re_theta_shifted, phase_4_cf_shifted, sz, 'filled', ...
        'MarkerFaceColor', phase_4_shifted_color', 'HandleVisibility', 'off')
plot(phase_4_Re_theta_shifted, phase_4_cf_shifted, 'Color', phase_4_shifted_color', ...
     'linewidth', linewidth, 'displayname', '$\varphi = 3 \lambda/4 + \delta$')



% Legend for positions along wave
plot(nan, nan, 'color', 'white', 'displayname', '')

% Mark the starting point
scatter(phase_2_Re_theta(1), phase_2_cf(1), 2*sz, '^', 'filled', 'MarkerFaceColor', 'white', ...
        'MarkerEdgeColor', 'black', 'LineWidth', linewidth, 'HandleVisibility', 'off')
hLeg = plot(nan, nan, '^', ...
    'MarkerFaceColor', 'white', ...
    'MarkerEdgeColor', 'black', ...
    'MarkerSize', 4, ...
    'LineWidth', linewidth, ...
    'LineStyle', 'none', ...
    'DisplayName', '$\xi / \lambda = 0$');


% Mark the crest
scatter(phase_2_Re_theta(end), phase_2_cf(end), 2*sz, 'o', 'filled', 'MarkerFaceColor', 'white', ...
        'MarkerEdgeColor', 'black', 'LineWidth', linewidth, 'HandleVisibility', 'off')
hLeg = plot(nan, nan, 'o', ...
    'MarkerFaceColor', 'white', ...
    'MarkerEdgeColor', 'black', ...
    'MarkerSize', 4, ...
    'LineWidth', linewidth, ...
    'LineStyle', 'none', ...
    'DisplayName', '$\xi / \lambda = 0.5$');

% Mark the ending point
scatter(phase_4_Re_theta_shifted(end), phase_4_cf_shifted(end), 2*sz, 'v', 'filled', 'MarkerFaceColor', 'white', ...
        'MarkerEdgeColor', 'black', 'LineWidth', linewidth, 'HandleVisibility', 'off')
% scatter(nan, nan, sz, 'v', 'filled', 'MarkerFaceColor', 'white', ...
%         'MarkerEdgeColor', 'black', 'LineWidth', linewidth, 'displayname', '$\xi / \lambda = 1$')
hLeg = plot(nan, nan, 'v', ...
    'MarkerFaceColor', 'white', ...
    'MarkerEdgeColor', 'black', ...
    'MarkerSize', 4, ...
    'LineWidth', linewidth, ...
    'LineStyle', 'none', ...
    'DisplayName', '$\xi / \lambda = 1$');


% Plot lines to annotate \Delta between start and stop
x0 = phase_2_Re_theta(1);
y0 = phase_2_cf(1);
x1 = phase_4_Re_theta_shifted(end);
y1 = phase_4_cf_shifted(end);

% One polyline: vertical up to (x0,y1), then horizontal to (x1,y1)
P = plot([x0, x0, x1], [y0, y1, y1], ...
         'LineWidth', linewidth, ...
         'Color', 'black', ...
         'LineJoin', 'round', ...
         'HandleVisibility', 'off');
uistack(P,'bottom')
hold off

% Midpoints of each segment
midV = [x0, (y0+y1)/2];      % vertical segment midpoint
midH = [(x0+x1)/2, y1];      % horizontal segment midpoint

% Small offsets in DATA units (tune once)
dx = 0.01 * range(xlim);     % right/left nudge
dy = 0.01 * range(ylim);     % up/down nudge

% Put LaTeX text centered on the vertical segment (nudged left)
txtV = text(midV(1)+dx, midV(2), '$\Delta c_f$', ...
    'Interpreter','latex', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','middle', ...
    'FontSize', 6);

% Put LaTeX text centered on the horizontal segment (nudged up)
txtH = text(midH(1), midH(2)-dy, '$\Delta Re_{\theta}$', ...
    'Interpreter','latex', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','top', ...
    'FontSize', 6);


% Legend
leg = legend('box', 'off', 'interpreter', 'latex', 'location', 'northwest', 'fontsize', legendFontSize);
leg.IconColumnWidth = 10;




% Tick marks
xticks(Re_theta_ticks)
yticks(Cf_ticks)

% Make y axis easier to read
ax3.YAxis.Exponent = -3;
ax3.YAxis.TickLabelFormat = '%2.0f';

% Axes labels
xlabel('$Re_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$C_f$', 'interpreter', 'latex', 'fontsize', labelFontSize)

% Axes limits
xlim([Re_theta_LL, Re_theta_UL])
ylim([Cf_LL, Cf_UL])


% Print shifts
fprintf('Re_theta shift and scaling\n\n')
fprintf('Shift applied: %.4f\n', phase_4_Re_theta_shift);
fprintf('Cf shift and scaling\n')
fprintf('Shift applied: %.4f\n', phase_4_cf_shift);




%% See why WVA looks rough

nu = 1.48E-5;
wind_speed = 'WT6';
wave = 'A';
caze = strcat(wind_speed, '_WV', wave, '_AG0');

if ismember(wind_speed(end), {'4'})
    u_inf = 2.4181;
elseif ismember(wind_speed(end), {'6'})
    u_inf = 3.8709;
elseif ismember(wind_speed(end), {'8'})
    u_inf = 5.4289;
end

smoothing_kernel = 16;

previous_cf = 0;
previous_Re_theta = 0;

clc; close all
figure('color', 'white')

hold on
for phase = 2:2:4
    disp(phase)

    % Skin friction
    cf = skin_friction_profiles.(caze)(phase).cf;
    % cf = smoothdata(cf, 'movmean', smoothing_kernel);
    % cf = filloutliers(cf, 'spline');
    
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

    cf = cf(~isnan(cf));
    masked_x = masked_x(~isnan(masked_x));
    Re_theta = Re_theta(~isnan(Re_theta));
    % cf = smoothdata(cf, 'movmean', 12);
    cf = filloutliers(cf, 'previous', 'median');
    cf = sgolayfilt(cf, 3, 25);

    if strcmp(wind_speed, 'WT6') && strcmp(wave, 'A') && phase == 2
        % fprintf(':)\n')
        % Remove bad data around range
        badRegion = (masked_x >= 0.0) & (masked_x <= 0.24);
        cf_clean = cf;
        cf_clean(badRegion) = NaN;
        cf_filled = fillmissing(cf_clean, 'pchip');
        % cf_sgolay_filled = sgolayfilt(cf_filled, 3, 99);

        plot(Re_theta, cf_filled, 'color', 'red')

        previous_cf = cf_filled;
        previous_Re_theta = Re_theta;
    else

        delta_Re_theta = Re_theta(1) - previous_Re_theta(end);
        delta_cf = cf(1) - previous_cf(end);

        plot(Re_theta - delta_Re_theta, cf - delta_cf, 'color', 'blue')

    end

    

  

end
hold off



%% Loop through and plot: Wind speed per subplot NO SCALING


wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
% smoothing_kernel = 16;

smoothing_kernel = 16;
sg_order = 3;
sg_framelength = 25;

wind_speeds = {'WT4', 'WT6', 'WT8'};
waves = {'A', 'B', 'C', 'D'};

% figure('color', 'white')
% tiledlayout(1, 3)

clc;
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

            % Re_{\theta}
            theta = integral.(caze)(phase).filtered.momentum;
            Re_theta = (u_inf * theta) / nu;

            centered_normalized_x = (integral.(caze)(phase).x - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength;
            min_array_size = min([length(cf), length(Re_theta)]);
            centered_normalized_x = centered_normalized_x(1:min_array_size);
            cf = cf(1:min_array_size);
            Re_theta = Re_theta(1:min_array_size);

            if strcmp(wind_speed, 'WT6') && strcmp(wave, 'A') && phase == 2

                % Crop data
                x_mask = centered_normalized_x;
                x_mask(abs(centered_normalized_x) > 0.25) = nan;
                x_mask(~isnan(x_mask)) = 1;
    
                cf = cf .* x_mask;
                Re_theta = Re_theta .* x_mask;
                masked_x = centered_normalized_x .* x_mask;
    
                masked_x = masked_x(~isnan(masked_x));
                cf = cf(~isnan(cf));
                Re_theta = Re_theta(~isnan(Re_theta));

                cf = smoothdata(cf, 'movmean', smoothing_kernel);
                cf = filloutliers(cf, 'previous', 'median');
                cf = sgolayfilt(cf, sg_order, sg_framelength);

                % Remove bad data around range
                badRegion = (masked_x >= 0.0) & (masked_x <= 0.24);
                cf_clean = cf;
                cf_clean(badRegion) = NaN;
                cf = fillmissing(cf_clean, 'pchip');
            else

                % Crop data
                x_mask = centered_normalized_x;
                x_mask(abs(centered_normalized_x) > 0.25) = nan;
                x_mask(~isnan(x_mask)) = 1;
    
                cf = cf .* x_mask;
                Re_theta = Re_theta .* x_mask;
                masked_x = centered_normalized_x .* x_mask;

                masked_x = masked_x(~isnan(masked_x));
                cf = cf(~isnan(cf));
                Re_theta = Re_theta(~isnan(Re_theta));

                cf = smoothdata(cf, 'movmean', 12);
                % cf = filloutliers(cf, 'previous', 'median');
                % cf = sgolayfilt(cf, sg_order, sg_framelength);
    
            end


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

        % Save data back to array to make seperate loop for plotting
        clean_data.(caze).cf_2 = phase_2_cf;
        clean_data.(caze).cf_4 = phase_4_cf_shifted;

        clean_data.(caze).Re_theta_2 = phase_2_Re_theta;
        clean_data.(caze).Re_theta_4 = phase_4_Re_theta_shifted;

        clear hysteresis phase_4_Re_theta_shifted phase_4_cf_shifted 
    end
    % hold off

end







% Loop through and plot: Wave speed per subplot
tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;

linewidth = 1;

wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};

wind_speeds = {'WT4', 'WT6', 'WT8'};
waves = {'A', 'B', 'C', 'D'};

wind_speed_alphas = [1, 1, 1];

% phase_markers = {'o', '^'};

sz = 5;

close all
ax = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,5]);
t = tiledlayout(1, 4, 'Padding', 'tight', 'TileSpacing', 'compact');

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
        phase_2_cf = clean_data.(caze).cf_2 * 1E3;

        phase_4_Re_theta_shifted = clean_data.(caze).Re_theta_4;
        phase_4_cf_shifted = clean_data.(caze).cf_4 * 1E3;

        % Plot
        % scatter(phase_2_Re_theta, phase_2_cf, sz, phase_markers{1}, 'filled', 'MarkerFaceColor', wave_colors{w}, ...
        %         'MarkerFaceAlpha', wind_speed_alphas(s), 'MarkerEdgeAlpha', wind_speed_alphas(s), 'HandleVisibility', 'off')
        plot(phase_2_Re_theta, phase_2_cf, 'color', wave_colors{w}, 'linewidth', 2, 'HandleVisibility', 'off')

        % scatter(phase_4_Re_theta_shifted, phase_4_cf_shifted, sz, phase_markers{1}, 'filled', 'MarkerFaceColor', wave_colors{w}, ...
        %         'MarkerFaceAlpha', wind_speed_alphas(s), 'MarkerEdgeAlpha', wind_speed_alphas(s), 'HandleVisibility', 'off')
        plot(phase_4_Re_theta_shifted, phase_4_cf_shifted, 'color', wave_colors{w}, 'linewidth', 2, 'HandleVisibility', 'off')

        % Mark the starting point
        scatter(phase_2_Re_theta(1), phase_2_cf(1), 2*sz, 'square', 'filled', 'MarkerFaceColor', 'white', ...
                'MarkerEdgeColor', 'black', 'LineWidth', 1', 'HandleVisibility', 'off')

        % Mark where peak is
        scatter(phase_4_Re_theta_shifted(1), phase_4_cf_shifted(1), 2*sz, 'o', 'filled', 'MarkerFaceColor', 'white', ...
                'MarkerEdgeColor', 'black', 'LineWidth', 1', 'HandleVisibility', 'off')

        % Mark where loop ends
        scatter(phase_4_Re_theta_shifted(end), phase_4_cf_shifted(end), 2*sz, 'diamond', 'filled', 'MarkerFaceColor', 'white', ...
                'MarkerEdgeColor', 'black', 'LineWidth', 1', 'HandleVisibility', 'off')

        clear hysteresis phase_4_Re_theta_shifted phase_4_cf_shifted 

        if w == 1
            ylabel('$C_f \times 10^{3}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
        end

        if w ~= 1
            yticklabels([]);
        end

    end
    hold off
    xticks(3000:2000:6000)
    label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses.(wave));
    title(label, 'interpreter', 'latex', 'fontsize', 10)

    
    % Make y axis easier to read
    % ax = gca;
    % ax.YAxis.Exponent = -3;          % shows ×10^{-6} once
    % ax.YAxis.TickLabelFormat = '%2.0f';  % two digits on ticks

end

% Make a legend
hold on

% What colors mean
% for w = 1:length(waves)
%     disp(w)
%     wave = waves{w};
%     label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses.(wave));
%     plot(nan, nan, 'color', wave_colors{w}, 'LineWidth', 1, 'DisplayName', label)
% end
% 
% % Legend white space
% plot(nan, nan, 'color', 'white', 'displayname', '')

% Legend for keypoint markers
% Mark the starting point
hLeg = plot(nan, nan, 'square', ...
    'MarkerFaceColor', 'white', ...
    'MarkerEdgeColor', 'black', ...
    'MarkerSize', 4, ...
    'LineWidth', linewidth, ...
    'LineStyle', 'none', ...
    'DisplayName', '$\xi / \lambda = 0$');

% Mark the crest
hLeg = plot(nan, nan, 'o', ...
    'MarkerFaceColor', 'white', ...
    'MarkerEdgeColor', 'black', ...
    'MarkerSize', 4, ...
    'LineWidth', linewidth, ...
    'LineStyle', 'none', ...
    'DisplayName', '$\xi / \lambda = 0.5$');

% Mark the ending point
hLeg = plot(nan, nan, 'diamond', ...
    'MarkerFaceColor', 'white', ...
    'MarkerEdgeColor', 'black', ...
    'MarkerSize', 4, ...
    'LineWidth', linewidth, ...
    'LineStyle', 'none', ...
    'DisplayName', '$\xi / \lambda = 1$');

hold off
leg = legend('interpreter', 'latex', 'box', 'off', 'orientation', 'horizontal', 'fontsize', legendFontSize);
% leg.IconColumnWidth = 19;
leg.Layout.Tile = 'north';
leg.ItemTokenSize = [10 8];


xlabel(t, '$Re_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)

linkaxes(h, 'xy')
% xlim([2000, 6200])
xlim([2000, 6500])
% ylim([0 20E-3])
ylim([0 20])


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

% tickFontSize = 14;
% labelFontSize = 16;
% legendFontSize = 14;
% titleFontSize = 18;
% 
% clc;
% wind_speed_colors = {'#0075F2', '#FF8C42', '#D30C7B'};
% 
% % Range of Cf vs steepness
% ax = figure('color', 'white', 'units', 'centimeters', 'position', [5,5,30,10]);
% t = tiledlayout(1,2, 'padding', 'tight');
% 
% sz = 60;
% linewidth = 3;
% 
% nexttile()
% hold on
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% for s = 1:3
%     for w = 1:length(waves)
%         wave = waves{w};
%         wind_speed = wind_speeds{s};
% 
%         if ismember(wind_speed(end), {'4'})
%             u_inf = 2.4181;
%         elseif ismember(wind_speed(end), {'6'})
%             u_inf = 3.8709;
%         elseif ismember(wind_speed(end), {'8'})
%             u_inf = 5.4289;
%         end
% 
%         % Get wavelengths
%         wavelength = wavelengths.(wave);
%         amplitude = amplitudes.(wave);
%         steepness = (2 * pi * amplitude) / wavelength;
% 
%         % Save steepnesses to plot a line
%         steepnessez(w) = steepness;
% 
%         % Plot
%         if w == 1 
%             vis = 'on';
%         else
%             vis = 'off';
%         end
% 
%         label = sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf);
%         scatter(steepness, cf_ranges(s,w), sz, 'filled', 'markerfacecolor', wind_speed_colors{s}, ...
%                 'HandleVisibility', vis, 'displayname', label)
%     end
%     [sorted_steepnesses, sortIdx] = sort(steepnessez,'ascend');
%     plot(sorted_steepnesses, cf_ranges(s, sortIdx), 'color', wind_speed_colors{s}, 'linewidth', linewidth, 'HandleVisibility', 'off')
% 
% end
% hold off
% leg = legend('Orientation', 'vertical', 'box', 'off', ...
%              'fontsize', legendFontSize, 'interpreter', 'latex', 'location', 'northwest');
% 
% xlabel('$ak$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylabel('$\max \left( C_f \right) - \min \left( C_f \right)$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% xlim([0.15, 0.32])
% ylim([0, 16E-3])
% yticks(0:4E-3:16E-3)
% % Make y axis easier to read
% ax = gca;
% ax.YAxis.Exponent = -3;          % shows ×10^{-6} once
% ax.YAxis.TickLabelFormat = '%2.0f';  % two digits on ticks
% 
% 
% % Range of Re_{\theta} vs wavelength
% nexttile()
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% hold on
% for s = 1:3
%     for w = 1:length(waves)
%         wave = waves{w};
%         wind_speed = wind_speeds{s};
% 
%         if ismember(wind_speed(end), {'4'})
%             u_inf = 2.4181;
%         elseif ismember(wind_speed(end), {'6'})
%             u_inf = 3.8709;
%         elseif ismember(wind_speed(end), {'8'})
%             u_inf = 5.4289;
%         end
% 
%         % Get wavelengths
%         wavelength = wavelengths.(wave);
%         amplitude = amplitudes.(wave);
%         steepness = (2 * pi * amplitude) / wavelength;
% 
%         % Save wavelengths to plot a line
%         wvlenths(w) = wavelength;
% 
%         % Plot
%         if w == 1 
%             vis = 'on';
%         else
%             vis = 'off';
%         end
% 
%         label = sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf);
%         scatter(wavelength, Re_theta_ranges(s,w), sz, 'filled', 'markerfacecolor', wind_speed_colors{s}, ...
%                 'displayname', label, 'HandleVisibility', vis)
%     end
%     [sorted_wavelengths, sortIdx] = sort(wvlenths,'ascend');
%     plot(sorted_wavelengths, Re_theta_ranges(s, sortIdx), 'color', wind_speed_colors{s}, 'linewidth', linewidth, 'HandleVisibility', 'off')
% 
%     % Try fitting a power law
%     p = polyfit(log(sorted_wavelengths), log(Re_theta_ranges(s, sortIdx)), 1);
%     m = p(1);
%     b = exp(p(2));
%     fprintf('b = %.4f,   m = %.4f\n', b, m);
%     fitz = b * sorted_wavelengths.^m;
%     % plot(sorted_wavelengths, fitz, 'linestyle', '--', 'linewidth', lw, 'Color', wind_speed_colors{s});
% end
% hold off
% 
% xlabel('$\lambda$ [m]', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylabel('$\max \left( {Re}_{\theta} \right) - \min \left( {Re}_{\theta} \right)$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylim([0, 6000])
% yticks(0:2000:6000)


% Save figure
% pause(3)
% figure_name = 'Hysteresis_Ranges.pdf';
% exportgraphics(t, fullfile(figure_folder, 'Hysteresis', figure_name), 'resolution', 600, 'ContentType', 'image')
% close all



%% Plot closure of Cf against steepness and range of Re_{\theta} against wavelength

% tickFontSize = 14;
% labelFontSize = 16;
% legendFontSize = 14;
% titleFontSize = 18;
% 
% wind_speed_colors = {'#0075F2', '#FF8C42', '#D30C7B'};
% 
% % Range of Cf vs steepness
% ax = figure('color', 'white', 'units', 'centimeters', 'position', [5,5,30,10]);
% t = tiledlayout(1,2,'padding', 'tight');
% 
% sz = 60;
% linewidth = 3;
% 
% h(1) = nexttile();
% hold on
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% for s = 1:3
%     for w = 1:length(waves)
%         wave = waves{w};
%         wind_speed = wind_speeds{s};
% 
%         if ismember(wind_speed(end), {'4'})
%             u_inf = 2.4181;
%         elseif ismember(wind_speed(end), {'6'})
%             u_inf = 3.8709;
%         elseif ismember(wind_speed(end), {'8'})
%             u_inf = 5.4289;
%         end
% 
%         % Get wavelengths
%         wavelength = wavelengths.(wave);
%         amplitude = amplitudes.(wave);
%         steepness = (2 * pi * amplitude) / wavelength;
% 
%         % Save steepnesses to plot a line
%         steepnessez(w) = steepness;
% 
%         % Plot
%         if w == 1 
%             vis = 'on';
%         else
%             vis = 'off';
%         end
% 
%         label = sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf);
%         scatter(steepness, cf_closure(s,w), sz, 'filled', 'markerfacecolor', wind_speed_colors{s}, ...
%                 'HandleVisibility', vis, 'displayname', label)
%     end
%     [sorted_steepnesses, sortIdx] = sort(steepnessez,'ascend');
%     plot(sorted_steepnesses, cf_closure(s, sortIdx), 'color', wind_speed_colors{s}, 'linewidth', linewidth, 'HandleVisibility', 'off')
% 
% end
% hold off
% leg = legend('Orientation', 'vertical', 'box', 'off', ...
%              'fontsize', legendFontSize, 'interpreter', 'latex', 'location', 'northwest');
% 
% xlabel('$ak$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylabel('$\Delta C_f $', 'interpreter', 'latex', 'fontsize', labelFontSize)
% xlim([0.15, 0.32])
% ylim([-10E-3, 10E-3])
% % yticks(0:4E-3:16E-3)
% % Make y axis easier to read
% ax = gca;
% ax.YAxis.Exponent = -3;          % shows ×10^{-6} once
% ax.YAxis.TickLabelFormat = '%2.0f';  % two digits on ticks
% 
% 
% % Range of Re_{\theta} vs wavelength
% h(2) = nexttile();
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% hold on
% for s = 1:3
%     for w = 1:length(waves)
%         wave = waves{w};
%         wind_speed = wind_speeds{s};
% 
%         if ismember(wind_speed(end), {'4'})
%             u_inf = 2.4181;
%         elseif ismember(wind_speed(end), {'6'})
%             u_inf = 3.8709;
%         elseif ismember(wind_speed(end), {'8'})
%             u_inf = 5.4289;
%         end
% 
%         % Get wavelengths
%         wavelength = wavelengths.(wave);
%         amplitude = amplitudes.(wave);
%         steepness = (2 * pi * amplitude) / wavelength;
% 
%         % Save wavelengths to plot a line
%         wvlenths(w) = wavelength;
% 
%         % Plot
%         if w == 1 
%             vis = 'on';
%         else
%             vis = 'off';
%         end
% 
%         label = sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf);
%         scatter(wavelength, Re_theta_closure(s,w), sz, 'filled', 'markerfacecolor', wind_speed_colors{s}, ...
%                 'displayname', label, 'HandleVisibility', vis)
%     end
%     [sorted_wavelengths, sortIdx] = sort(wvlenths,'ascend');
%     plot(sorted_wavelengths, Re_theta_closure(s, sortIdx), 'color', wind_speed_colors{s}, 'linewidth', linewidth, 'HandleVisibility', 'off')
% end
% hold off
% 
% xlabel('$\lambda$ [m]', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylabel('$\Delta {Re}_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylim([0, 2000])


% Save figure
% pause(3)
% figure_name = 'Hysteresis_Closures.pdf';
% exportgraphics(t, fullfile(figure_folder, 'Hysteresis',  figure_name), 'resolution', 600, 'ContentType', 'image')
% close all





%% Create 1x3 plot of ranges and Re_theta closure

tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;

clc; close all
wind_speed_colors = {'#0075F2', '#FF8C42', '#D30C7B'};

% Range of Cf vs steepness
ax = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,5]);
t = tiledlayout(1,3, 'padding', 'compact', 'TileSpacing', 'loose');

sz = 10;
linewidth = 1.2;

h(1) = nexttile;
hold on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)

all_steepness = [];
all_cf        = [];

for s = 1:3
    tmpX = nan(1, length(waves));
    tmpY = nan(1, length(waves));

    steepness_missing = NaN;

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

        wavelength = wavelengths.(wave);
        amplitude  = amplitudes.(wave);
        steepness  = (2*pi*amplitude)/wavelength;

        % ---- special missing point ----
        if (s == 3) && (w == 4)
            tmpX(w) = steepness;
            tmpY(w) = NaN;
            steepness_missing = steepness;
            continue
        end

        % collect for global fit
        all_steepness(end+1) = steepness;
        all_cf(end+1)        = 1E3 * cf_ranges(s,w);

        % scatter (original data)
        vis = 'on';
        if w > 1, vis = 'off'; end

        label = sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf);
        scatter(steepness, 1E3 * cf_ranges(s,w), sz, 'filled', ...
            'MarkerFaceColor', wind_speed_colors{s}, ...
            'HandleVisibility', vis, ...
            'DisplayName', label)

        tmpX(w) = steepness;
        tmpY(w) = 1E3* cf_ranges(s,w);
    end

    % ---- line plot with interpolation ----
    [sortedX, sortIdx] = sort(tmpX, 'ascend');
    sortedY = tmpY(sortIdx);

    valid = ~isnan(sortedX);
    sortedX = sortedX(valid);
    sortedY = sortedY(valid);

    yLine = fillmissing(sortedY, 'spline', ...
                        'SamplePoints', sortedX);

    plot(sortedX, yLine, ...
        'Color', wind_speed_colors{s}, ...
        'LineWidth', linewidth, ...
        'HandleVisibility', 'off')

    % ---- plot interpolated scatter point ----
    if s == 3 && ~isnan(steepness_missing)
        y_interp = interp1(sortedX, yLine, steepness_missing);

        scatter(steepness_missing, y_interp, sz, ...
            'MarkerFaceColor', wind_speed_colors{s}, ...
            'MarkerEdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end
end


% ---- global fit ----
p    = polyfit(all_steepness, all_cf, 1);
xfit = linspace(0.1, 0.4, 10);
yfit = polyval(p, xfit);
P = plot(xfit, yfit, 'color', 'black', 'LineStyle', '--', 'LineWidth', 0.5, ...
         'DisplayName', 'Linear fit (all data)');
uistack(P, 'bottom');


% Show equation
m = p(1);
b = p(2);

exp_m = floor(log10(abs(m)));
mant_m  = m / 10^exp_m;
exp_b = floor(log10(abs(b)));
mant_b  = b / 10^exp_b;

R = corrcoef(all_steepness, all_cf);
R2 = R(1,2)^2;
eqn = sprintf('$C_f = %.3f\\, ak + %.3f$\\quad \n($R^2 = %.2f$)', m, b, R2);

% text(0.05, 0.95, eqn, ...
%      'Units', 'normalized', ...
%      'Interpreter', 'latex', ...
%      'FontSize', annotationFontSize, ...
%      'VerticalAlignment', 'top')

hold off

axis square
xlabel('$ak$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylabel('$\max \left( C_f \right) - \min \left( C_f \right)$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$\mathrm{range} \left( C_f \right) \times 10^{3}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlim([0.17, 0.32])
ylim([0, 16])
yticks(0:4:16)
ax = gca;
% ax.YAxis.Exponent = -3;          % shows ×10^{-6} once
% ax.YAxis.TickLabelFormat = '%2.0f';  % two digits on ticks




% Range of Re_{\theta} vs wavelength
h(2) = nexttile;
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
        if w == 100 
            vis = 'on';
        else
            vis = 'off';
        end

        label = sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf);
        scatter(wavelength, Re_theta_ranges(s,w), sz, 'filled', 'markerfacecolor', wind_speed_colors{s}, ...
                'displayname', label, 'HandleVisibility', vis)
    end
    [sorted_wavelengths, sortIdx] = sort(wvlenths,'ascend');
    plot(sorted_wavelengths, Re_theta_ranges(s, sortIdx), 'color', wind_speed_colors{s}, 'linewidth', linewidth, 'HandleVisibility', 'off')

    x = sorted_wavelengths(:);
    y = Re_theta_ranges(s, sortIdx).';   % make column
    
    % Fit y = b*x (through origin)
    % b = (x' * y) / (x' * x);
    p = polyfit(x,y,1);
    m = p(1);
    b = p(2);    
    
    % Smooth line for plotting
    xf = linspace(min(x), max(x), 10).';
    yf = m * xf + b;
    
    P = plot(xf, yf, '--', 'LineWidth', 0.5, 'Color', 'black', 'HandleVisibility','off');
    uistack(P, 'bottom')
    

    % After computing b
    exp10 = floor(log10(abs(m)));
    % exp10 = 3;
    mant  = m / 10^exp10;
    
    label = sprintf('$\\sim %.2f \\times 10^{%d}\\,\\lambda$', mant, exp10);
    
    % Place near mid-range of line
    % text(max(xf), max(yf), label, ...
    %      'Interpreter','latex', ...
    %      'FontSize', annotationFontSize, ...
    %      'Color', 'black', ...
    %      'BackgroundColor','white', ...
    %      'Margin', 2, ...
    %      'HorizontalAlignment','left', ...
    %      'VerticalAlignment','middle');


end
hold off

axis square
xlabel('$\lambda$ [m]', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$\mathrm{range} (Re_\theta)$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylabel('$\max \left( {Re}_{\theta} \right) - \min \left( {Re}_{\theta} \right)$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylim([0, 6000])
yticks(0:2000:6000)
xlim([0.1, 0.45])





% Range of Re_{\theta} vs wavelength
h(3) = nexttile();
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
        if w == 100
            vis = 'on';
        else
            vis = 'off';
        end

        label = sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf);
        scatter(wavelength, Re_theta_closure(s,w), sz, 'filled', 'markerfacecolor', wind_speed_colors{s}, ...
                'displayname', label, 'HandleVisibility', vis)
    end
    [sorted_wavelengths, sortIdx] = sort(wvlenths,'ascend');
    plot(sorted_wavelengths, Re_theta_closure(s, sortIdx), 'color', wind_speed_colors{s}, 'linewidth', linewidth, 'HandleVisibility', 'off')


    x = sorted_wavelengths(:);
    y = Re_theta_closure(s, sortIdx).';   % make column
    
    % Fit y = b*x (through origin)
    % b = (x' * y) / (x' * x);
    p = polyfit(x,y,1);
    m = p(1);
    b = p(2);    
    
    % Smooth line for plotting
    % xf = linspace(min(x), max(x), 200).';
    xf = linspace(min(x), max(x), 10).';
    yf = m * xf + b;
    
    P = plot(xf, yf, '--', 'LineWidth', 0.5, 'Color', 'black', 'HandleVisibility','off');
    uistack(P, 'bottom')
    
    % fprintf('through-origin slope b = %g\n', b);

    % After computing b
    exp10 = floor(log10(abs(m)));
    % exp10 = 3;
    mant  = m / 10^exp10;
    
    label = sprintf('$\\sim %.2f \\times 10^{%d}\\,\\lambda$', mant, exp10);
    
    % Place near mid-range of line
    % text(max(xf), max(yf), label, ...
    %      'Interpreter','latex', ...
    %      'FontSize', annotationFontSize, ...
    %      'Color', 'black', ...
    %      'BackgroundColor','white', ...
    %      'Margin', 2, ...
    %      'HorizontalAlignment','left', ...
    %      'VerticalAlignment','middle');

end

% Make legend
for s = 1:3
    wind_speed = wind_speeds{s};
    if ismember(wind_speed(end), {'4'})
        u_inf = 2.4181;
    elseif ismember(wind_speed(end), {'6'})
        u_inf = 3.8709;
    elseif ismember(wind_speed(end), {'8'})
        u_inf = 5.4289;
    end

    label = sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf);
    plot(nan, nan, 'linewidth', linewidth, 'color', wind_speed_colors{s}, ...
            'displayname', label, 'HandleVisibility', 'on')
end

hold off

axis square
xlabel('$\lambda$ [m]', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$\Delta {Re}_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylim([0, 2000])
xlim([0.08, 0.4])
leg = legend('Orientation', 'horizontal', 'box', 'off', ...
             'fontsize', legendFontSize, 'interpreter', 'latex');
leg.Layout.Tile = 'north';
leg.IconColumnWidth = 19;
linkaxes(h(2:3), 'x')

% Add panel labels
addPanelLabels(h, {'a', 'b', 'c'}, 'FontSize', 10, 'Offset', [0.05, 1])

% addPanelLabels(ax, labels) adds (a),(b),... just OUTSIDE top-left of each axes.
% ax     : array of axes handles (e.g., from tiledlayout / findall)
% labels : cellstr like {'a','b','c'} or string array ["a" "b" "c"]
%
% Optional name-value:
% 'Offset'   : [dx dy] in normalized axes units (default [-0.10 1.02])
% 'FontSize' : default 12
% 'FontName' : default 'Times New Roman'


% Save figure
% pause(3)
% figure_name = 'Hysteresis_RangesClosures_Combined.pdf';
% exportgraphics(t, fullfile(figure_folder, 'Hysteresis',  figure_name), 'resolution', 600, 'ContentType', 'image')
% close all





%% Create 1x3 plot of ranges and Re_theta closure AGAINST STEEPNESS

tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;

clc; close all
wind_speed_colors = {'#0075F2', '#FF8C42', '#D30C7B'};

% Range of Cf vs steepness
ax = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,5]);
t = tiledlayout(1,3, 'padding', 'compact', 'TileSpacing', 'loose');

sz = 10;
linewidth = 1.2;

h(1) = nexttile;
hold on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)

all_steepness = [];
all_cf        = [];

for s = 1:3
    tmpX = nan(1, length(waves));
    tmpY = nan(1, length(waves));

    steepness_missing = NaN;

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

        wavelength = wavelengths.(wave);
        amplitude  = amplitudes.(wave);
        steepness  = (2*pi*amplitude)/wavelength;

        % ---- special missing point ----
        if (s == 3) && (w == 4)
            tmpX(w) = steepness;
            tmpY(w) = NaN;
            steepness_missing = steepness;
            continue
        end

        % collect for global fit
        all_steepness(end+1) = steepness;
        all_cf(end+1)        = 1E3 * cf_ranges(s,w);

        % scatter (original data)
        vis = 'on';
        if w > 1, vis = 'off'; end

        label = sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf);
        scatter(steepness, 1E3 * cf_ranges(s,w), sz, 'filled', ...
            'MarkerFaceColor', wind_speed_colors{s}, ...
            'HandleVisibility', vis, ...
            'DisplayName', label)

        tmpX(w) = steepness;
        tmpY(w) = 1E3* cf_ranges(s,w);
    end

    % ---- line plot with interpolation ----
    [sortedX, sortIdx] = sort(tmpX, 'ascend');
    sortedY = tmpY(sortIdx);

    valid = ~isnan(sortedX);
    sortedX = sortedX(valid);
    sortedY = sortedY(valid);

    yLine = fillmissing(sortedY, 'spline', ...
                        'SamplePoints', sortedX);

    plot(sortedX, yLine, ...
        'Color', wind_speed_colors{s}, ...
        'LineWidth', linewidth, ...
        'HandleVisibility', 'off')

    % ---- plot interpolated scatter point ----
    if s == 3 && ~isnan(steepness_missing)
        y_interp = interp1(sortedX, yLine, steepness_missing);

        scatter(steepness_missing, y_interp, sz, ...
            'MarkerFaceColor', wind_speed_colors{s}, ...
            'MarkerEdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end
end


% ---- global fit ----
p    = polyfit(all_steepness, all_cf, 1);
xfit = linspace(0.1, 0.4, 10);
yfit = polyval(p, xfit);
P = plot(xfit, yfit, 'color', 'black', 'LineStyle', '--', 'LineWidth', 0.5, ...
         'DisplayName', 'Linear fit (all data)');
uistack(P, 'bottom');


% Show equation
m = p(1);
b = p(2);

exp_m = floor(log10(abs(m)));
mant_m  = m / 10^exp_m;
exp_b = floor(log10(abs(b)));
mant_b  = b / 10^exp_b;

R = corrcoef(all_steepness, all_cf);
R2 = R(1,2)^2;
eqn = sprintf('$C_f = %.3f\\, ak + %.3f$\\quad \n($R^2 = %.2f$)', m, b, R2);

% text(0.05, 0.95, eqn, ...
%      'Units', 'normalized', ...
%      'Interpreter', 'latex', ...
%      'FontSize', annotationFontSize, ...
%      'VerticalAlignment', 'top')

hold off

axis square
xlabel('$ak$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$\mathrm{range} \left( C_f \right) \times 10^{3}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlim([0.17, 0.32])
ylim([0, 16])
yticks(0:4:16)
ax = gca;




% Range of Re_{\theta} vs wavelength
h(2) = nexttile;
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
        wvlenths(w) = steepness;

        % Plot
        if w == 100 
            vis = 'on';
        else
            vis = 'off';
        end

        label = sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf);
        scatter(steepness, Re_theta_ranges(s,w), sz, 'filled', 'markerfacecolor', wind_speed_colors{s}, ...
                'displayname', label, 'HandleVisibility', vis)
    end
    [sorted_wavelengths, sortIdx] = sort(wvlenths,'ascend');
    plot(sorted_wavelengths, Re_theta_ranges(s, sortIdx), 'color', wind_speed_colors{s}, 'linewidth', linewidth, 'HandleVisibility', 'off')

    x = sorted_wavelengths(:);
    y = Re_theta_ranges(s, sortIdx).';   % make column
    
    % Fit y = b*x (through origin)
    % b = (x' * y) / (x' * x);
    p = polyfit(x,y,1);
    m = p(1);
    b = p(2);    
    
    % Smooth line for plotting
    xf = linspace(min(x), max(x), 10).';
    yf = m * xf + b;
    
    P = plot(xf, yf, '--', 'LineWidth', 0.5, 'Color', 'black', 'HandleVisibility','off');
    uistack(P, 'bottom')
    

    % After computing b
    exp10 = floor(log10(abs(m)));
    % exp10 = 3;
    mant  = m / 10^exp10;
    
    label = sprintf('$\\sim %.2f \\times 10^{%d}\\,\\lambda$', mant, exp10);
    
    % Place near mid-range of line
    % text(max(xf), max(yf), label, ...
    %      'Interpreter','latex', ...
    %      'FontSize', annotationFontSize, ...
    %      'Color', 'black', ...
    %      'BackgroundColor','white', ...
    %      'Margin', 2, ...
    %      'HorizontalAlignment','left', ...
    %      'VerticalAlignment','middle');


end
hold off

axis square
xlabel('$ak$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$\mathrm{range} (Re_\theta)$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylim([0, 6000])
yticks(0:2000:6000)
% xlim([0.1, 0.45])





% Range of Re_{\theta} vs wavelength
h(3) = nexttile();
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
        wvlenths(w) = steepness;

        % Plot
        if w == 100
            vis = 'on';
        else
            vis = 'off';
        end

        label = sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf);
        scatter(steepness, Re_theta_closure(s,w), sz, 'filled', 'markerfacecolor', wind_speed_colors{s}, ...
                'displayname', label, 'HandleVisibility', vis)
    end
    [sorted_wavelengths, sortIdx] = sort(wvlenths,'ascend');
    plot(sorted_wavelengths, Re_theta_closure(s, sortIdx), 'color', wind_speed_colors{s}, 'linewidth', linewidth, 'HandleVisibility', 'off')


    x = sorted_wavelengths(:);
    y = Re_theta_closure(s, sortIdx).';   % make column
    
    % Fit y = b*x (through origin)
    % b = (x' * y) / (x' * x);
    p = polyfit(x,y,1);
    m = p(1);
    b = p(2);    
    
    % Smooth line for plotting
    % xf = linspace(min(x), max(x), 200).';
    xf = linspace(min(x), max(x), 10).';
    yf = m * xf + b;
    
    P = plot(xf, yf, '--', 'LineWidth', 0.5, 'Color', 'black', 'HandleVisibility','off');
    uistack(P, 'bottom')
    
    % fprintf('through-origin slope b = %g\n', b);

    % After computing b
    exp10 = floor(log10(abs(m)));
    % exp10 = 3;
    mant  = m / 10^exp10;
    
    label = sprintf('$\\sim %.2f \\times 10^{%d}\\,\\lambda$', mant, exp10);
    
    % Place near mid-range of line
    % text(max(xf), max(yf), label, ...
    %      'Interpreter','latex', ...
    %      'FontSize', annotationFontSize, ...
    %      'Color', 'black', ...
    %      'BackgroundColor','white', ...
    %      'Margin', 2, ...
    %      'HorizontalAlignment','left', ...
    %      'VerticalAlignment','middle');

end

% Make legend
for s = 1:3
    wind_speed = wind_speeds{s};
    if ismember(wind_speed(end), {'4'})
        u_inf = 2.4181;
    elseif ismember(wind_speed(end), {'6'})
        u_inf = 3.8709;
    elseif ismember(wind_speed(end), {'8'})
        u_inf = 5.4289;
    end

    label = sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf);
    plot(nan, nan, 'linewidth', linewidth, 'color', wind_speed_colors{s}, ...
            'displayname', label, 'HandleVisibility', 'on')
end

hold off

axis square
xlabel('$ak$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$\Delta {Re}_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylim([0, 2000])
% xlim([0.08, 0.4])
leg = legend('Orientation', 'horizontal', 'box', 'off', ...
             'fontsize', legendFontSize, 'interpreter', 'latex');
leg.Layout.Tile = 'north';
leg.IconColumnWidth = 19;
linkaxes(h(2:3), 'x')

% Add panel labels
addPanelLabels(h, {'a', 'b', 'c'}, 'FontSize', 10, 'Offset', [0.05, 1])

% addPanelLabels(ax, labels) adds (a),(b),... just OUTSIDE top-left of each axes.
% ax     : array of axes handles (e.g., from tiledlayout / findall)
% labels : cellstr like {'a','b','c'} or string array ["a" "b" "c"]
%
% Optional name-value:
% 'Offset'   : [dx dy] in normalized axes units (default [-0.10 1.02])
% 'FontSize' : default 12
% 'FontName' : default 'Times New Roman'


% Save figure
% pause(3)
% figure_name = 'Hysteresis_RangesClosures_Combined.pdf';
% exportgraphics(t, fullfile(figure_folder, 'Hysteresis',  figure_name), 'resolution', 600, 'ContentType', 'image')
% close all











%% Plot Re_thet vs Cf for wave D to show the spiral behavior

% tickFontSize = 14;
% labelFontSize = 16;
% legendFontSize = 14;
% titleFontSize = 18;
% 
% 
% phase_names = {'$\varphi = 0$', ...
%                '$\varphi = \lambda / 4$', ...
%                '$\varphi = \lambda / 2$', ...
%                '$\varphi = 3 \lambda / 4$'};
% 
% wind_speed = 'WT4';
% wave = 'D';
% caze = strcat(wind_speed, '_WV', wave, '_AG0');
% 
% if ismember(wind_speed(end), {'4'})
%     u_inf = 2.4181;
% elseif ismember(wind_speed(end), {'6'})
%     u_inf = 3.8709;
% elseif ismember(wind_speed(end), {'8'})
%     u_inf = 5.4289;
% end
% 
% smoothing_kernel = 20;
% 
% colors = slanCM(35, 4);
% 
% clc; close all
% figure('color', 'white', 'units', 'centimeters', 'position', [5,5,30,15]);
% fig = tiledlayout(2,2,"TileSpacing",'tight', 'padding', 'tight');
% 
% ax1 = nexttile(1);
% ax2 = nexttile(3);
% ax3 = nexttile([2,1]);
% 
% 
% ax1.XAxis.Visible = 'off';
% 
% % Make y axis easier to read
% ax1.YAxis.Exponent = -3;
% ax1.YAxis.TickLabelFormat = '%2.0f';
% ax3.YAxis.Exponent = -3;
% ax3.YAxis.TickLabelFormat = '%2.0f';
% 
% 
% set([ax1, ax2, ax3], 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% 
% X_all = [];
% Re_all = [];
% 
% hold(ax1,'on'); hold(ax2,'on'); hold(ax3,'on'); 
% for phase = 1:4
%     disp(phase)
% 
%     % Skin friction
%     cf = skin_friction_profiles.(caze)(phase).cf;
%     cf = smoothdata(cf, 'movmean', smoothing_kernel);
% 
%     % Re_{\theta}
%     theta = integral.(caze)(phase).filtered.momentum;
%     theta = smoothdata(theta, 'movmean', smoothing_kernel);
%     Re_theta = (u_inf * theta) / nu;
% 
%     wave_profile = integral.(caze)(phase).wave.wave_profile;
%     centered_normalized_x = (integral.(caze)(phase).x - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength;
%     min_array_size = min([length(cf), length(Re_theta)]);
%     centered_normalized_x = centered_normalized_x(1:min_array_size);
%     cf = cf(1:min_array_size);
%     Re_theta = Re_theta(1:min_array_size);
%     wave_profile = wave_profile(1:min_array_size);
% 
%     % Crop data
%     x_mask = centered_normalized_x;
%     % Need to crop left and right differently for this to work
%     x_mask(centered_normalized_x < -0.5) = nan;
%     x_mask(centered_normalized_x > 0.25) = nan;
%     x_mask(~isnan(x_mask)) = 1;
% 
%     % Mask curves
%     cf = cf .* x_mask;
%     Re_theta = Re_theta .* x_mask;
%     wave_profile = wave_profile .* x_mask;
%     masked_x = centered_normalized_x .* x_mask;
%     cf = cf(~isnan(cf));
%     Re_theta = Re_theta(~isnan(Re_theta));
%     wave_profile = wave_profile(~isnan(wave_profile));
%     masked_x = masked_x(~isnan(masked_x));
% 
%     % Shift profiles appropriately
%     if phase ~= 1 
%         Re_theta_shift = Re_theta_old(end) - Re_theta(1);
%         cf_shift = cf_old(end) - cf(1);
%     else
%         Re_theta_shift = 0;
%         cf_shift = 0;
%     end
% 
% 
%     % Set color
%     color = colors(phase,:);
% 
%     % Plot stitched curves over 4x phases
%     plot(ax2, masked_x + (0.75 * (phase - 1)) + 0.5, Re_theta + Re_theta_shift, ...
%          'color', color, 'linewidth', 3, 'HandleVisibility', 'off')
%     plot(ax1, masked_x + (0.75 * (phase - 1)) + 0.5, cf + cf_shift, ...
%          'color', color, 'linewidth', 3, 'HandleVisibility', 'off')
% 
%     % Plot an X where there is a seam
%     if phase ~= 1
%         scatter(ax2, masked_x(1) + (0.75 * (phase - 1)) + 0.5,  Re_theta(1) + Re_theta_shift, ...
%                 100, 'x', 'markeredgecolor', 'black', 'linewidth', 2)
%         scatter(ax1, masked_x(1) + (0.75 * (phase - 1)) + 0.5,  cf(1) + cf_shift, ...
%                 100, 'x', 'markeredgecolor', 'black', 'linewidth', 2)
%     end
% 
%     % Mark seams
%     xline(ax1, (0.75 * (phase - 1)), 'color', 'black', 'linewidth', 2, 'alpha', 0.1, 'HandleVisibility', 'off')
%     xline(ax2, (0.75 * (phase - 1)), 'color', 'black', 'linewidth', 2, 'alpha', 0.1, 'HandleVisibility', 'off')
% 
% 
%     % Plot hysteresis loops over 4x phases
%     plot(ax3, Re_theta + Re_theta_shift, cf + cf_shift, 'color', color, ...
%          'linewidth', 2, 'DisplayName', phase_names{phase})
%     scatter(ax3, Re_theta + Re_theta_shift, cf + cf_shift, 40, 'filled', ...
%             'markerfacecolor', color, 'HandleVisibility', 'off')
% 
%     % Plot an X where there is a seam
%     if phase ~= 1
%         scatter(ax3, Re_theta(1) + Re_theta_shift, cf(1) + cf_shift, ...
%                 100, 'x', 'markeredgecolor', 'black', 'linewidth', 2, ...
%                 'HandleVisibility', 'off')
%     end
% 
%     % Mark start/stop
%     if phase == 1
%         scatter(Re_theta(1) + Re_theta_shift, cf(1) + cf_shift, 2*sz, '^', 'filled', 'MarkerFaceColor', 'white', ...
%         'MarkerEdgeColor', 'black', 'LineWidth', 2, 'HandleVisibility', 'off')
%     end
% 
%     % Mark the ending point
%     if phase == 4
%         scatter(Re_theta(end) + Re_theta_shift, cf(end) + cf_shift, 2*sz, 'v', 'filled', 'MarkerFaceColor', 'white', ...
%                 'MarkerEdgeColor', 'black', 'LineWidth', 2, 'HandleVisibility', 'off')
% 
% 
%     end
% 
%     % Save Re_theta values for fitting
%     x_plot  = masked_x + (0.75 * (phase - 1)) + 0.5;
%     Re_plot = Re_theta + Re_theta_shift;
%     % collect
%     X_all  = [X_all;  x_plot(:)];
%     Re_all = [Re_all; Re_plot(:)];
% 
% 
% 
%     % Update previous curve
%     Re_theta_old = Re_theta + Re_theta_shift;
%     cf_old = cf + cf_shift;
% 
% end
% 
% % Finish legend
% plot(ax3, nan, nan, 'color', 'white', 'DisplayName', ' ')
% scatter(ax3, nan, nan, 2*sz, '^', 'filled', 'MarkerFaceColor', 'white', ...
%         'MarkerEdgeColor', 'black', 'LineWidth', 1.5', 'displayname', '$\xi / \lambda = 0$')
% scatter(ax3, nan, nan, 2*sz, 'v', 'filled', 'MarkerFaceColor', 'white', ...
%                 'MarkerEdgeColor', 'black', 'LineWidth', 1.5', 'displayname', '$\xi / \lambda = 3$')
% 
% legend(ax3, 'interpreter', 'latex', 'box', 'off', 'fontsize', legendFontSize, 'location', 'northwest')
% 
% % clean + sort
% good = isfinite(X_all) & isfinite(Re_all);
% X = X_all(good);
% Y = Re_all(good);
% [X, I] = sort(X);
% Y = Y(I);
% p = polyfit(X, Y, 1);
% m = p(1);
% b = p(2);
% xfit = linspace(min(X), max(X), 200);
% yfit = polyval(p, xfit);
% 
% P = plot(ax2, xfit, yfit, 'k--', 'LineWidth', 2, 'DisplayName', 'Linear fit');
% uistack(P, 'bottom')
% str = sprintf('$Re_{\\theta} \\approx %.3g\\, x/\\lambda + %4.0f$', m,b);
% 
% text(ax2, 0.05, 0.93, str, ...
%     'Units','normalized', ...           
%     'Interpreter','latex', ...
%     'FontSize', 14, ...
%     'VerticalAlignment','top', ...
%     'BackgroundColor','white', ...
%     'Margin', 4, ...
%     'EdgeColor', [1 1 1], ...
%     'Clipping','on');
% 
% 
% 
% 
% % Plot reference profile
% x_ref = 0:0.001:3;
% ref_wave = 1 * cos((2 * pi) * (x_ref - 0.5));
% cf_wave = 8E-4 * ref_wave + 0.8E-3;
% plot(ax1, x_ref, cf_wave, 'linewidth', 2, 'color', 'black')
% 
% % Shaded region below wave
% hFill = patch( ax1, ...
% [x_ref, fliplr(x_ref)], ...
% [cf_wave, -10 * ones(size(cf_wave))], ...
% 'k', ...
% 'FaceAlpha', 0.5, ...      
% 'EdgeColor', 'none', ...    
% 'HandleVisibility', 'off'); 
% uistack(hFill, 'bottom')
% 
% 
% % Plot reference profile
% x_ref = 0:0.001:3;
% ref_wave = 1 * cos((2 * pi) * (x_ref - 0.5));
% Re_theta_wave = 50 * ref_wave + 1700;
% plot(ax2, x_ref, Re_theta_wave, 'linewidth', 2, 'color', 'black')
% 
% % Shaded region below wave
% hFill = patch( ax2, ...
% [x_ref, fliplr(x_ref)], ...
% [Re_theta_wave, -10 * ones(size(Re_theta_wave))], ...
% 'k', ...
% 'FaceAlpha', 0.5, ...      
% 'EdgeColor', 'none', ...    
% 'HandleVisibility', 'off'); 
% uistack(hFill, 'bottom')
% 
% 
% hold(ax2,'off'); hold(ax2,'off'); hold(ax3,'off'); 
% linkaxes([ax1, ax2], 'x')
% xticks([ax1, ax2], 0:0.75:4)
% ylabel(ax2, '$Re_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylabel(ax1, '$C_f$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% xlabel(ax2, '$\xi / \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% xlabel(ax3, '$Re_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylabel(ax3, '$C_f$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% 
% ylim(ax1, [-1E-3, 16E-3])
% ylim(ax2, [1600, 2600])
% ylim(ax3, [0, 18E-3])
% 
% 
% % Save figure
% % pause(3)
% % figure_name = sprintf('Hysteresis_%s_WV%s_Spiral.pdf', wind_speed, wave);
% % exportgraphics(fig, fullfile(figure_folder, 'Hysteresis', figure_name), 'resolution', 600, 'ContentType', 'image')
% % close all








%% Functions

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
