%%% Curvilinear Integral Boundary Layer Parameters

addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
curvilinear_path = fullfile(project_path, 'curvilinear');
means_path = fullfile(project_path, 'means');
wave_parameters = readcell("Offshore_Waves.xlsx");

% Boundary layer detection parameters
boundary_layer_percent = 0.96;
left_mask = 4E-3;
right_mask = 234E-3;

% Edge killing values
% OG was 0.008
threshold.thickness = 0.005;
threshold.displacement = 0.005;
threshold.momentum = 0.005;

% Data filtering and smoothing
smooth_kernel = 9;




% Loop through all waves and wind speeds
wind_speeds = {'WT4', 'WT6', 'WT8'};
waves = {'A', 'B', 'C', 'D'};
wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
line_styles = {'-.', '--', '-'};
lw = 2;
alpha = 0.25;
fs = 24;


figure('color', 'white')
hold on

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
            wave_profile = curvilinear.phase(phase).wave_profile * 1E-3;
            
            
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

            %%% Get u profile
            % Get profiles
            u = curvilinear.phase(phase).u;
            idx = round(length(u)/2);
            y_profile = horizontal_lines(:, idx) - wave_profile(idx); 
            u_profile = u(:, idx);


            %%% CLEAN UP U PROFILE
            u_profile = StatisticalGradientFilter(u_profile, 2.4, 10);
            u_profile(find(~isnan(u_profile), 1, 'first')) = nan;



            %%% PLOT
            % Normalizing by displacement over thickness
            difference_u = u_inf - u_profile;
            normalized_difference_u = (u_inf - u_profile) /u_inf;
            corrected_normalized_u = (u_inf - u_profile) / (u_inf * (displacement(idx) / thickness(idx)));
            normalized_y = y_profile / thickness(idx);

            P = plot(corrected_normalized_u, normalized_y, ...
                     'color', wave_colors{w}, 'linestyle', line_styles{s}, 'linewidth', lw, 'displayname', caze);
            P.Color(4) = alpha;
        

        end
    end
end



%%% No Wave Profiles
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


    caze = strcat(wind_speed, '_WV0_AGP');
    disp(caze)

    % Load data
    data = load(fullfile(means_path, strcat(caze, '_MEANS.mat')));
    data = data.output;
 
      
    % Curvilinear coordinates
    x = data.X(1,:) * 1E-3;
    y = data.Y(:,1) * 1E-3;

    % Find boundary layer edge
    u = data.ensemble.u;
    u = juliaan_smooth(u, smooth_kernel);

    % Mask the data properly
    u(x < left_mask | x > right_mask) = nan;
    x(x < left_mask | x > right_mask) = nan;
    y(x < left_mask | x > right_mask) = nan;
            
    % Find BL
    u(u > u_inf * boundary_layer_percent) = nan;
   
    %%% Boundary layer thickness
    thickness = nan(1, size(u,2));
    % x_BL = nan(1, size(u,2));
        
    for i = 1:size(u,2)
        column = u(:,i);
        if all(isnan(column))
            continue;
        end

        [~, idx] = max(column, [], 'omitnan');
        thickness(i) = horizontal_lines(idx, i);
        % x_BL(i) = x(idx, i);
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

    %%% Get u profile
    % Get profiles
    u = data.ensemble.u;
    idx = round(length(u)/2);
    y_profile = horizontal_lines(:, idx) - wave_profile(idx); 
    u_profile = u(:, idx);


    %%% CLEAN UP U PROFILE
    u_profile = StatisticalGradientFilter(u_profile, 2.4, 10);
    u_profile(find(~isnan(u_profile), 1, 'first')) = nan;


    % Normalizing by displacement over thickness
    difference_u = u_inf - u_profile;
    normalized_difference_u = (u_inf - u_profile) /u_inf;
    corrected_normalized_u = (u_inf - u_profile) / (u_inf * (displacement(idx) / thickness(idx)));
    normalized_y = y_profile / thickness(idx);


    %%% PLOT
    P = plot(corrected_normalized_u, normalized_y, ...
             'color', 'black', 'linestyle', line_styles{s}, 'linewidth', lw, 'displayname', caze);
    P.Color(4) = 1.0;
        
    ylabel('$y / \delta$', 'interpreter', 'latex', 'fontsize', 24)

end

hold off
% legend('fontsize', 2)
% ylim([0, 2])
% xlim([-0.5, 3])
% xlabel("$\frac{u_{\infty} - u}{\left( u_{\infty} \frac{\delta^*}{\delta} \right)}$", 'Interpreter', 'latex', 'fontsize', fs)
% ylabel("$\frac{y}{\delta}$", 'Interpreter', 'latex', 'fontsize', fs)



clear amplitude caze column curvilinear displacement displacement_cum displacement_edge_mask displacement_integrand
clear horizontal_lines i idx island_length left_mask max_wave_profile momentum momentum_cum momentum_edge_mask momentum_integrand
clear n nan_mask phase right_mask s sgolay_length smooth_kernerl thickness thickness_edge_mask threshold u u_column u_inf u_normalized
clear x w wave wave_profile wave_parameters wave_type wave_length wind_speed x_BL zeta_column
clear gradient_tolerance project_path smooth_kernel wavelength curvilinear_path y_profile u_profile corrected_normalized_u

clc;









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

