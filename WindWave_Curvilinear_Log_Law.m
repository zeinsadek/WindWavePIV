%%% Curvilinear Log Law

addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
% means_path = fullfile(project_path, 'means');
curvilinear_path = fullfile(project_path, 'curvilinear');

% Case 
caze = 'WT8_WVD_AG0';

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


%% UV profiles across different phases

lw = 2;

clc; close all;
figure()
hold on
for phase = 1:4

    % X (xi) in mm
    vertical_lines = curvilinear.phase(phase).vertical_lines * 1E-3;
    
    % Y (zeta) in mm
    horizontal_lines = curvilinear.phase(phase).horizontal_lines * 1E-3;
    
    % Waves
    wave_profile = curvilinear.phase(phase).wave_profile * 1E-3;
    max_wave_profile = curvilinear.phase(phase).max_wave_profile * 1E-3;
    
    % Get components
    uv = curvilinear.phase(phase).uv;
    
    % Get profiles
    idx = round(length(uv)/2);
    z_profile = horizontal_lines(:, idx) - wave_profile(idx); 
    uv_profile = uv(:, idx); 

    % Polynomial fit
    valid = ~isnan(z_profile) & ~isnan(uv_profile);
    n = 9;
    p = polyfit(z_profile(valid), uv_profile(valid), n);
    uv_fit = polyval(p, z_profile);
    uv_fit(isnan(uv_profile)) = nan;


    %%% Get Friction Velocity
    u_star_raw(phase) = max(sqrt(-uv_profile), [], 'all', 'omitnan');
    u_star_fit(phase) = max(sqrt(-uv_fit), [], 'all', 'omitnan');
    fprintf("Phase %1.0f:\nraw u* = %2.4f m/s\n", phase, u_star_raw(phase))
    fprintf("fit u* = %2.4f m/s\n\n", u_star_fit(phase))


    % Plot
    plot(uv_profile, z_profile, 'linewidth', lw, 'handlevisibility', 'off', 'color', 'black')
    plot(uv_fit, z_profile, 'linewidth', lw, 'displayname', sprintf("Phase %1.0f", phase))
    

end
hold off
legend('location', 'northwest')
title(caze, 'interpreter', 'none')
ylim([0, 0.25])
xlim([-0.25, 0.05])


figure()
bar([1,2,3,4], u_star_raw)
yline(mean(u_star_raw))
fprintf("AVG u* = %1.3f m/s\n\n", mean(u_star_raw))
ylim([0,1])
title(caze, 'interpreter', 'none')



%% Skin Friction Profiles

lw = 2;

clc; close all;
figure()
hold on
for phase = 1:4

    % X (xi) in mm
    vertical_lines = curvilinear.phase(phase).vertical_lines * 1E-3;
    
    % Y (zeta) in mm
    horizontal_lines = curvilinear.phase(phase).horizontal_lines * 1E-3;
    
    % Waves
    wave_profile = curvilinear.phase(phase).wave_profile * 1E-3;
    max_wave_profile = curvilinear.phase(phase).max_wave_profile * 1E-3;
    
    % Get components
    uv = curvilinear.phase(phase).uv;
    
    % Get profiles
    for i = 1:size(uv,2)
        column = uv(:,i);
        if all(isnan(column))
            continue;
        end

        [M, I] = max(sqrt(-column), [], 'omitnan');
        u_stars_tmp(i) = M;
        cfs_tmp(i) = 2 * (M / u_inf).^2;
        x_tmp(i) = vertical_lines(I, i);
    end


    % Save values
    friction(phase).u_star = u_stars_tmp;
    friction(phase).x = x_tmp;
    friction(phase).cf = cfs_tmp;
    
    % Plot
    plot(x_tmp, cfs_tmp, 'linewidth', lw', 'DisplayName', sprintf("Phase %1.0f", phase))

end
hold off
legend('location', 'northwest')
title(caze, 'interpreter', 'none')
% ylim([0, 0.25])
% xlim([-0.25, 0.05])

%% Velocity profiles across different phases

lw = 2;
nu = 1.46E-5;
sz = 20;
colors = {'#AF1B3F', '#FC814A', '#3EC300', '#084B83'};

clc; close all;
figure()
hold on
for phase = 1:4

    % X (xi) in mm
    vertical_lines = curvilinear.phase(phase).vertical_lines * 1E-3;
    
    % Y (zeta) in mm
    horizontal_lines = curvilinear.phase(phase).horizontal_lines * 1E-3;
    
    % Waves
    wave_profile = curvilinear.phase(phase).wave_profile * 1E-3;
    max_wave_profile = curvilinear.phase(phase).max_wave_profile * 1E-3;
    
    % Get components
    u = curvilinear.phase(phase).u;
    uv = curvilinear.phase(phase).uv;
    
    % Get profiles
    idx = round(length(uv)/2);
    z_profile = horizontal_lines(:, idx) - wave_profile(idx); 
    uv_profile = uv(:, idx); 
    u_profile = u(:, idx);

    
    % Test cleaning up u profile
    test = StatisticalGradientFilter(u_profile, 2.4, 10);

    
    % Plot
    % plot(u_profile, z_profile)
    plot(u_profile, z_profile, 'linewidth', lw, 'displayname', sprintf("Phase %1.0f", phase), 'color', colors{phase})
    scatter(u_profile, z_profile, sz, 'filled', 'HandleVisibility', 'off', 'MarkerFaceColor', colors{phase})
    

end
hold off
legend('location', 'northwest')
title(caze, 'interpreter', 'none')

ylim([0, 0.25])
xlim([1, 6])


%% Log law profiles across different phases

lw = 2;
nu = 1.46E-5;
sz = 20;
colors = {'#AF1B3F', '#FC814A', '#3EC300', '#084B83'};

clc; close all;
figure()
hold on
for phase = 1:4

    % X (xi) in mm
    vertical_lines = curvilinear.phase(phase).vertical_lines * 1E-3;
    
    % Y (zeta) in mm
    horizontal_lines = curvilinear.phase(phase).horizontal_lines * 1E-3;
    
    % Waves
    wave_profile = curvilinear.phase(phase).wave_profile * 1E-3;
    max_wave_profile = curvilinear.phase(phase).max_wave_profile * 1E-3;
    
    % Get components
    u = curvilinear.phase(phase).u;
    uv = curvilinear.phase(phase).uv;
    
    % Get profiles
    idx = round(length(uv)/2);
    z_profile = horizontal_lines(:, idx) - wave_profile(idx); 
    uv_profile = uv(:, idx); 
    u_profile = u(:, idx);

    
    % Test cleaning up u profile
    test = StatisticalGradientFilter(u_profile, 2.4, 10);

    
    %%% Get Friction Velocity
    u_star_raw(phase) = max(sqrt(-uv_profile), [], 'all', 'omitnan');
    u_star_fit(phase) = max(sqrt(-uv_fit), [], 'all', 'omitnan');
    fprintf("Phase %1.0f:\nraw u* = %2.4f m/s\n", phase, u_star_raw(phase))
    fprintf("fit u* = %2.4f m/s\n\n", u_star_fit(phase))

    % Law of the wall
    % u_plus = test / u_star_raw(phase);
    u_plus = test / mean(u_star_raw);
    % y_plus = (z_profile * u_star_raw(phase)) / nu;
    y_plus = (z_profile * mean(u_star_raw)) / nu;

    % Plot
    plot(y_plus, u_plus, 'linewidth', lw, 'displayname', sprintf("Phase %1.0f", phase), 'color', colors{phase})
    scatter(y_plus, u_plus, sz, 'filled', 'HandleVisibility', 'off', 'MarkerFaceColor', colors{phase})
    

end
hold off
legend('location', 'northwest')
title(caze, 'interpreter', 'none')
xscale('log')
ylim([6, 25])
xlim([1E2, 1E4])
grid on









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
