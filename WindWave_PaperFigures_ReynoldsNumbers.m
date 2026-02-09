% Compute different Reynolds numbers for large table in paper
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

% Wind speeds
freestreams.('WT4') = 2.4181;
freestreams.('WT6') = 3.8709;
freestreams.('WT8') = 5.4289;

% Wave parameters
wavelengths.A = 410.8;
wavelengths.B = 313.3;
wavelengths.C = 189.6;
wavelengths.D = 124.3;

amplitudes.A = 11.78;
amplitudes.B = 10.53;
amplitudes.C = 9.21;
amplitudes.D = 5.29;

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


%% Compute wave age and wave Reynolds number per wind+wave combination

% Kinematic viscosity
nu = 1.46E-5;

clc;
for w = 1:length(waves)
    wave = waves{w};
    wavelength = wavelengths.(wave);
    amplitude = amplitudes.(wave);
    frequency = frequencies.(wave);
    wave_speed = wavelength * frequency;

    for s = 1:length(wind_speeds)
        wind_speed = wind_speeds{s};
        freestream = freestreams.(wind_speed);
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        % Wave age
        wave_age = wave_speed / freestream;

        % Wave Reynolds number
        Re_wave = (freestream * wavelength) / nu;

        % Loop through phases to find max values of friction and roughness
        % Reynolds numbers
        Re_frictions = nan(4,171);
        Re_roughnesses = nan(4,171);
        Re_momentums = nan(4,171);
        for phase = 1:4
            u_star = friction_velocity.(caze)(phase);
            delta = integral.(caze)(phase).filtered.thickness;
            momentum = integral.(caze)(phase).filtered.momentum;

            Re_frictions(phase, :) = (u_star * delta) / nu;
            Re_roughnesses(phase, :) = (u_star * (2 * amplitude)) / nu;
            Re_momentums(phase, :) = (freestream * momentum) / nu;
        end

        % Determine min/max ranges of Re_friction (~1000)
        min_Re_friction = min(Re_frictions, [], 'all', 'omitnan');
        max_Re_friction = max(Re_frictions, [], 'all', 'omitnan');

        % Determine min/max ranges of Re_roughness (~100)
        min_Re_roughness = min(Re_roughnesses, [], 'all', 'omitnan');
        max_Re_roughness = max(Re_roughnesses, [], 'all', 'omitnan');

         % Determine min/max ranges of Re_momentum (~1000)
        min_Re_momentum = min(Re_momentums, [], 'all', 'omitnan');
        max_Re_momentum = max(Re_momentums, [], 'all', 'omitnan');

        
        % Print results
        disp(caze)
        fprintf('u_inf = %1.2f\n', freestream)
        fprintf('Wave age (c/u_inf) = %2.3f\n', wave_age)
        fprintf('Re_wave = %3.2f X 10^3\n', Re_wave * 1E-3)
        fprintf('Range of Re_friction = %2.3f - %2.3f X 10^3\n', min_Re_friction * 1E-3, max_Re_friction * 1E-3)
        fprintf('Range of Re_roughness = %3.2f - %3.2f\n', min_Re_roughness, max_Re_roughness)
        fprintf('Range of Re_momentum = %3.2f - %3.2f\n', min_Re_momentum * 1E-3, max_Re_momentum * 1E-3)
        fprintf('\n')

    end
    fprintf('\n\n\n\n\n')
end






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
