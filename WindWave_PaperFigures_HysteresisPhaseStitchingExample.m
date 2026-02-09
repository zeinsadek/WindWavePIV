% Plotting the hystersis observed between C_f and Re_{\theta} for the
% curvilinear velocity fields

addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps/slanCM')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
means_path = fullfile(project_path, 'means');
curvilinear_path = fullfile(project_path, 'curvilinear');
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
fullfigure = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,7]);
fig = tiledlayout(2,2,"TileSpacing",'compact', 'padding', 'tight');

% C_f profiles
% ax1 = nexttile([1 2]);
ax1 = nexttile(1);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
hold on
plot((hysteresis(2).x) + 0.25, hysteresis(2).cf * 1E3, ...
     'linewidth', linewidth, 'color', phase_2_color)

plot((hysteresis(4).x + 0.5) + 0.25, hysteresis(4).cf * 1E3, ...
     'linewidth', linewidth, 'color', phase_4_color, 'linestyle', '-')

plot((hysteresis(4).x + 0.5) + 0.25, phase_4_cf_shifted * 1E3, ...
     'linewidth', linewidth, 'color', phase_4_shifted_color')

% Starting Point
scatter(0,  hysteresis(2).cf(1) * 1E3, 2*sz, 'square', 'filled', 'MarkerFaceColor', 'white', ...
        'MarkerEdgeColor', 'black', 'LineWidth', linewidth, 'HandleVisibility', 'off')

% Middle
mp = scatter(0.5,  hysteresis(2).cf(end) * 1E3, 2*sz, 'o', 'filled', 'MarkerFaceColor', 'white', ...
        'MarkerEdgeColor', 'black', 'LineWidth', linewidth, 'HandleVisibility', 'off');

% End point
scatter(1,  phase_4_cf_shifted(end) * 1E3, 2*sz, 'diamond', 'filled', 'MarkerFaceColor', 'white', ...
        'MarkerEdgeColor', 'black', 'LineWidth', linewidth, 'HandleVisibility', 'off')





% Plot reference wave profile
Cf_reference_scale = 0.5 * 1E3;
Cf_reference_offset = 2E-3 * 1E3;
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

% Plot vetical line
plot([0.5, 0.5], [-10, 1E6], 'color', 'black', 'linewidth', linewidth)

% Plot details
hold off
uistack(mp, 'top')
ylabel('$C_f \times 10^{3}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylim([-2E-3 * 1E3, 0.025 * 1E3])
yticks(0:10:30)
xlim([0, 1])
ax = gca;
ax.XAxis.Visible = 'off';

% Make y axis easier to read

% Show what the \delta is
exponent = -3;
mantissa = phase_4_cf_shift / 10^exponent;
text(0.51, 0.49, sprintf('$\\delta C_f = %0.1f \\times 10^{%d}$', -mantissa, exponent), ...
    'Units', 'normalized', ...
    'Interpreter', 'latex', ...
    'FontSize', annotationFontSize)






%%% Re_{\theta} profiles
% ax2 = nexttile([1 2]);
ax2 = nexttile(3);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
hold on
% xline(0.5, 'color', 'black', 'linestyle', '-', 'linewidth', linewidth, 'HandleVisibility', 'off');

plot((hysteresis(2).x) + 0.25, hysteresis(2).Re_theta, ...
     'linewidth', linewidth, 'color', phase_2_color, 'displayname', '$\varphi = \lambda/4$')

plot((hysteresis(4).x + 0.5) + 0.25, hysteresis(4).Re_theta, ...
     'linewidth', linewidth, 'color', phase_4_color, 'linestyle', '-', 'Displayname', '$\varphi = 3 \lambda/4$')

plot((hysteresis(4).x + 0.5) + 0.25, phase_4_Re_theta_shifted, ...
     'linewidth', linewidth, 'color', phase_4_shifted_color', 'Displayname', '$\varphi = 3 \lambda/4 + \delta$')

% Starting Point
scatter(0,  hysteresis(2).Re_theta(1), 2*sz, 'square', 'filled', 'MarkerFaceColor', 'white', ...
        'MarkerEdgeColor', 'black', 'LineWidth', linewidth, 'HandleVisibility', 'off')

% Middle
mp = scatter(0.5,  hysteresis(2).Re_theta(end), 2*sz, 'o', 'filled', 'MarkerFaceColor', 'white', ...
        'MarkerEdgeColor', 'black', 'LineWidth', linewidth, 'HandleVisibility', 'off');

% End point
scatter(1,  phase_4_Re_theta_shifted(end), 2*sz, 'diamond', 'filled', 'MarkerFaceColor', 'white', ...
        'MarkerEdgeColor', 'black', 'LineWidth', linewidth, 'HandleVisibility', 'off')




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

% Plot vetical line
plot([0.5, 0.5], [0, 1E6], 'color', 'black', 'linewidth', linewidth)

% Plot details
hold off
uistack(mp, 'top')
ylabel('$Re_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlabel('$\xi \mathbin{/} \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylim([2000, 5000])
yticks(3000:1000:5000)
xlim([0, 1])
xticks(0:0.25:1)

% Show what the \Delta is
text(0.24, 0.38, sprintf('$\\delta Re_{\\theta} = %3.2f$', -phase_4_Re_theta_shift), ...
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
scatter(phase_2_Re_theta, phase_2_cf * 1E3, sz, 'filled', ...
        'MarkerFaceColor', phase_2_color, 'HandleVisibility', 'off')
plot(phase_2_Re_theta, phase_2_cf * 1E3, 'Color', phase_2_color, ...
     'linewidth', linewidth, 'displayname', '$\varphi = \lambda/4$')

% Plot raw phase
scatter(phase_4_Re_theta, phase_4_cf * 1E3, sz, 'filled', ...
        'MarkerFaceColor', phase_4_color, 'HandleVisibility', 'off')
P = plot(phase_4_Re_theta, phase_4_cf * 1E3, 'Color', phase_4_color, ...
         'linewidth', linewidth, 'displayname', '$\varphi = 3 \lambda/4$');


% Plot shifted phase
scatter(phase_4_Re_theta_shifted, phase_4_cf_shifted * 1E3, sz, 'filled', ...
        'MarkerFaceColor', phase_4_shifted_color', 'HandleVisibility', 'off')
plot(phase_4_Re_theta_shifted, phase_4_cf_shifted * 1E3, 'Color', phase_4_shifted_color', ...
     'linewidth', linewidth, 'displayname', '$\varphi = 3 \lambda/4 + \delta$')


% Mark the starting point
scatter(phase_2_Re_theta(1), phase_2_cf(1) * 1E3, 2*sz, 'square', 'filled', 'MarkerFaceColor', 'white', ...
        'MarkerEdgeColor', 'black', 'LineWidth', linewidth, 'HandleVisibility', 'off')

% Mark the crest
scatter(phase_2_Re_theta(end), phase_2_cf(end) * 1E3, 2*sz, 'o', 'filled', 'MarkerFaceColor', 'white', ...
        'MarkerEdgeColor', 'black', 'LineWidth', linewidth, 'HandleVisibility', 'off')

% Mark the ending point
scatter(phase_4_Re_theta_shifted(end), phase_4_cf_shifted(end) * 1E3, 2*sz, 'diamond', 'filled', 'MarkerFaceColor', 'white', ...
        'MarkerEdgeColor', 'black', 'LineWidth', linewidth, 'HandleVisibility', 'off')



% Plot lines to annotate \Delta between start and stop
x0 = phase_2_Re_theta(1);
y0 = phase_2_cf(1) * 1E3;
x1 = phase_4_Re_theta_shifted(end);
y1 = phase_4_cf_shifted(end) * 1E3;

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
leg = legend('box', 'off', 'interpreter', 'latex', 'orientation', 'horizontal', 'fontsize', legendFontSize);
leg.IconColumnWidth = 10;
leg.Layout.Tile = 'north';




% Tick marks
xticks(Re_theta_ticks)
yticks(Cf_ticks * 1E3)

% Axes labels
xlabel('$Re_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$C_f \times 10^{3}$', 'interpreter', 'latex', 'fontsize', labelFontSize)

% Axes limits
xlim([3000, Re_theta_UL])
ylim([3, 18])


% Print shifts
fprintf('Re_theta shift and scaling\n\n')
fprintf('Shift applied: %.4f\n', phase_4_Re_theta_shift);
fprintf('Cf shift and scaling\n')
fprintf('Shift applied: %.4f\n', phase_4_cf_shift);

% Add a,b,c labels
addPanelLabelsFixed(fullfigure, [ax1, ax2, ax3], {'a', 'b' ,'c'}, 'FontSize', 10, 'OffsetPts', [-2,-18])


% Save figure
pause(3)
figure_name = 'HysteresisStitchExamples.pdf';
exportgraphics(fig, fullfile(figure_folder, 'Hysteresis', figure_name), 'Resolution', 600, 'ContentType', 'image');
close all
fprintf('Generated figure: %s\n\n', figure_name)




%%

% %% Plot Re_thet vs Cf for wave D to show the spiral behavior
% 
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

function hAnn = addPanelLabelsFixed(fig, ax, labels, varargin)
% Places (a),(b),... just outside top-left with FIXED physical offsets.
%
% fig    : figure handle (e.g., totalFigure)
% ax     : array of axes handles
% labels : {'a','b',...} or ["a","b",...]
%
% Name-value:
% 'OffsetPts' : [dx dy] in points (default [-10 6])  (left, up)
% 'FontSize'  : default 12
% 'FontName'  : default 'Times New Roman'

p = inputParser;
addParameter(p,'OffsetPts',[-10 6]);
addParameter(p,'FontSize',12);
addParameter(p,'FontName','Times New Roman');
parse(p,varargin{:});

labels = string(labels);
offPts = p.Results.OffsetPts;

ppi = get(0,'ScreenPixelsPerInch');
offPx = offPts/72 * ppi;   % points -> pixels

hAnn = gobjects(numel(ax),1);

for k = 1:numel(ax)
    if ~isgraphics(ax(k),'axes'), continue; end

    % Axes outer position in figure pixels
    op = getpixelposition(ax(k), true);  % [x y w h] in pixels relative to figure

    % Anchor point: outside top-left of axes
    x = op(1) + offPx(1);
    y = op(2) + op(4) + offPx(2);

    % TeX gives roman parentheses + italic letter
    str = sprintf('(\\it%s\\rm)', labels(k));

    hAnn(k) = annotation(fig,'textbox', ...
        'Units','pixels', ...
        'Position',[x y 40 20], ...   % small box; big enough for "(a)"
        'String',str, ...
        'Interpreter','tex', ...
        'FontName',p.Results.FontName, ...
        'FontSize',p.Results.FontSize, ...
        'LineStyle','none', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom');
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
