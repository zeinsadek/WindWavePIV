%%% Curvilinear Integral Boundary Layer Parameters

addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
means_path = fullfile(project_path, 'means');
curvilinear_path = fullfile(project_path, 'curvilinear');

% Case
caze = 'WT8_WVA_AG0';

% Free stream
tunnel_freq = caze(strfind(caze, 'WT') + 2);

if ismember(tunnel_freq, {'4'})
    u_inf = 2.4181;
elseif ismember(tunnel_freq, {'6'})
    u_inf = 3.8709;
elseif ismember(tunnel_freq, {'8'})
    u_inf = 5.4289;
end

% Load data
curvilinear = load(fullfile(curvilinear_path, strcat(caze, '_CURVILINEAR.mat')));
curvilinear = curvilinear.output;

% Wave Parameters
wave_type       = caze(strfind(caze, 'WV') + 2);
wave_parameters = readcell("Offshore_Waves.xlsx");
wavelength      = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 2};
amplitude       = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 3};

clear means_path curvilinear_path project_path wave_parameters wave_type

%% Boundary Layer Edge Detection

percent = 0.96;
left_mask = 4E-3;
right_mask = 234E-3;

threshold.thickness = 0.008;
threshold.displacement = 0.005;
threshold.momentum = 0.005;

smooth_kernel = 5;
gradient_tolerance = 1;
island_length = 2;
sgolay_length = 35;
n = 9;


for phase = 1:4

    % X (xi) in mm
    % X = curvilinear.X * 1E-3;
    % x = X(1,:);
    % dx = mean(diff(X(1,:)));
    vertical_lines = curvilinear.phase(phase).vertical_lines * 1E-3;
    
    % Y (zeta) in mm
    % Y = curvilinear.Y * 1E-3;
    % y = Y(:,1);
    % dy = mean(diff(Y(:,1)));
    horizontal_lines = curvilinear.phase(phase).horizontal_lines * 1E-3;
    
    % Waves
    wave_profile = curvilinear.phase(phase).wave_profile * 1E-3;
    max_wave_profile = curvilinear.phase(phase).max_wave_profile * 1E-3;
    
    % Find boundary layer edge
    u = curvilinear.phase(phase).u;
    u = juliaan_smooth(u, smooth_kernel);

    %%% GPT
    % Mask the data properly
    u(vertical_lines < left_mask | vertical_lines > right_mask) = nan;
    vertical_lines(vertical_lines < left_mask | vertical_lines > right_mask) = nan;
    horizontal_lines(vertical_lines < left_mask | vertical_lines > right_mask) = nan;
    
    % Find BL
    u(u > u_inf * percent) = nan;
   
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
    % Initialize
    displacement = nan(1, size(u,2));   % displacement thickness
    momentum = nan(1, size(u,2));       % momentum thickness
    
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
            fprintf('Im Flipping Out!!\n')
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

    % Save raw signals
    integral(phase).raw.thickness = thickness;
    integral(phase).raw.displacement = displacement;
    integral(phase).raw.momentum = momentum;
    integral(phase).raw.shape = displacement ./ momentum;
    integral(phase).x = x_BL;

    % Save filtered
    integral(phase).filtered.thickness = FilterData(x_BL, thickness, 10, 10, sgolay_length);
    integral(phase).filtered.displacement = FilterData(x_BL, displacement, gradient_tolerance, island_length, sgolay_length);
    integral(phase).filtered.momentum = FilterData(x_BL, momentum, gradient_tolerance, island_length, sgolay_length);
    integral(phase).filtered.shape = FilterData(x_BL, displacement ./ momentum, gradient_tolerance, island_length, sgolay_length);

    % Save polynomial fit
    integral(phase).fitted.thickness = PolynomialFit(x_BL, integral(phase).filtered.thickness, n);
    integral(phase).fitted.displacement = PolynomialFit(x_BL, integral(phase).filtered.displacement, n);
    integral(phase).fitted.momentum = PolynomialFit(x_BL, integral(phase).filtered.momentum, n);
    integral(phase).fitted.shape = PolynomialFit(x_BL, integral(phase).filtered.shape, n);

    % Plot
    % h(phase) = nexttile;
    % hold on
    % contourf(vertical_lines, horizontal_lines, u, 20)
    % plot(x, wave_profile, 'color', 'black', 'linewidth', 3)
    % plot(x, max_wave_profile, 'color', 'red', 'linewidth', 3)
    % plot(x_BL, thickness, 'color', 'cyan', 'linewidth', 2)
    % hold off
    % axis equal
    % title(sprintf("Phase %1.0f", phase))
end

clc;
% linkaxes(h, 'xy')


%%% Plot
clc;
lw = 2;
left_bound = 3E-3;
right_bound = 230E-3;

close all
figure()
tiledlayout(1,4)
h(1) = nexttile;
hold on
for phase = 1:4
    plot(integral(phase).x, integral(phase).raw.thickness, 'linewidth', lw, 'color', 'black')
    plot(integral(phase).x, integral(phase).filtered.thickness, 'linewidth', lw, 'color', 'red')
    plot(integral(phase).x, integral(phase).fitted.thickness, 'linewidth', lw, 'color', 'blue')
end
hold off
title('Curvilinear Bounary Layer Thickness')

h(2) = nexttile;
hold on
for phase = 1:4
    plot(integral(phase).x, integral(phase).raw.displacement, 'linewidth', lw, 'color', 'black')
    plot(integral(phase).x, integral(phase).filtered.displacement, 'linewidth', lw, 'color', 'red')
    plot(integral(phase).x, integral(phase).fitted.displacement, 'linewidth', lw, 'color', 'blue')
end
hold off
title('Curvilinear Displacement Thickness')

h(3) = nexttile;
hold on
for phase = 1:4
    plot(integral(phase).x, integral(phase).raw.momentum, 'linewidth', lw, 'color', 'black')
    plot(integral(phase).x, integral(phase).filtered.momentum, 'linewidth', lw, 'color', 'red')
    plot(integral(phase).x, integral(phase).fitted.momentum, 'linewidth', lw, 'color', 'blue')
end
hold off
title('Curvilinear Momentum Thickness')

h(4) = nexttile;
hold on
for phase = 1:4
    plot(integral(phase).x, integral(phase).raw.shape, 'linewidth', lw, 'color', 'black')
    plot(integral(phase).x, integral(phase).filtered.shape, 'linewidth', lw, 'color', 'red')
    plot(integral(phase).x, integral(phase).fitted.shape, 'linewidth', lw, 'color', 'blue')
end
hold off
title('Curvilinear Shape Factor')

% linkaxes(h, 'xy')

%% Pseudo boundary layer stitching

WL = wavelength * 1E-3;
phase_shifts = [0, WL/4, WL/2, 3*WL/4];

figure()
hold on
for i = 1:4
    plot((integral(i).x - phase_shifts(i))/WL, integral(i).filtered.momentum, 'linewidth', lw)
end
hold off


%% Testing filtering and smoothing data

gradient_tolerance = 1;
island_length = 2;
sgolay_length = 25;

close all; clc;
figure()
tiledlayout(2,2)
for phase = 1:4

    fprintf('Phase %1.0f\n\n', phase);
    x = integral(phase).x;
    q = integral(phase).raw.momentum;
    % Remove spikes and data inbetween them
    q_filt = StatisticalGradientFilter(q, 1, 10);
    % Fill-in removed spiked
    q_interp = InterpolateInterior(x,q_filt);
    % Smooth out
    q_smooth = smoothdata(q_interp, 'sgolay', sgolay_length);

    % Test all-in-one function
    test = FilterData(x, q, gradient_tolerance, island_length, sgolay_length);

    %%% Polynomial fit?
    % valid = ~isnan(x) & ~isnan(q_interp);
    % p = polyfit(x(valid), q_interp(valid), 9);
    % q_fitted = polyval(p, x);
    % q_fitted(isnan(q_interp)) = nan;
    test_poly_fit = PolynomialFit(x, q_interp, 9);
     
    h(phase) = nexttile;
    hold on
    plot(x, q, 'linewidth', lw, 'color', 'black')
    plot(x, q_interp, 'color', 'blue', 'linewidth', 2);
    % plot(x, q_fitted, 'color', 'red', 'linewidth', 2);
    plot(x, test_poly_fit, 'color', 'green', 'linewidth', 2);
    hold off
end

linkaxes(h, 'xy')









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





