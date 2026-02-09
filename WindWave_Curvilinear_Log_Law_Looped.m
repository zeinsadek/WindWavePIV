%%% Curvilinear Log Law

addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
means_path = fullfile(project_path, 'means');
curvilinear_path = fullfile(project_path, 'curvilinear');


clc; close all;
% figure()
% hold on

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

        for phase = 1:4
        
            % X (xi) in mm
            vertical_lines = curvilinear.phase(phase).vertical_lines * 1E-3;
            
            % Y (zeta) in mm
            horizontal_lines = curvilinear.phase(phase).horizontal_lines * 1E-3;
            
            % Waves
            wave_profile = curvilinear.phase(phase).wave_profile * 1E-3;
            % max_wave_profile = curvilinear.phase(phase).max_wave_profile * 1E-3;
            
            % Get components
            u = curvilinear.phase(phase).u;
            uv = curvilinear.phase(phase).uv;
            
            % Get profiles
            idx = round(length(uv)/2);
            y_profile = horizontal_lines(:, idx) - wave_profile(idx); 
            uv_profile = uv(:, idx); 
            u_profile = u(:, idx);


            %%% CLEAN UP U PROFILE
            u_profile = StatisticalGradientFilter(u_profile, 2.4, 10);
            u_profile(find(~isnan(u_profile), 1, 'first')) = nan;
        
        
            %%% Get Friction Velocity
            u_star.(caze)(phase).raw = max(sqrt(-uv_profile), [], 'all', 'omitnan');
            % u_star.(caze)(phase).fit = max(sqrt(-uv_fit), [], 'all', 'omitnan');
            u_star.(caze)(phase).u_profile = u_profile;
            u_star.(caze)(phase).u_profile_normalized = u_profile / (u_inf);
            u_star.(caze)(phase).uv_profile = uv_profile;
            u_star.(caze)(phase).uv_profile_normalized = uv_profile / (u_inf^2);
            u_star.(caze)(phase).y = y_profile;
        
            % Plot
            % plot(uv_profile / (u_inf^2), y_profile, 'linewidth', lw, 'handlevisibility', 'off', 'color', 'black')
            % plot(uv_fit / (u_inf^2), y_profile, 'linewidth', lw, 'displayname', sprintf("Phase %1.0f", phase))
            
        
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
    u = means.ensemble.u;
    uv = means.ensemble.uv;

    % Get profiles
    idx = round(length(uv)/2);
    y_profile = means.Y(:, 1) * 1E-3; 
    uv_profile = uv(:, idx); 
    u_profile = u(:, idx);

    % CLEAN UP U PROFILE
    u_profile = StatisticalGradientFilter(u_profile, 2.4, 10);
    u_profile(find(~isnan(u_profile), 1, 'first')) = nan;

    % Polynomial fit
    % valid = ~isnan(y_profile) & ~isnan(uv_profile);
    % n = 9;
    % p = polyfit(y_profile(valid), uv_profile(valid), n);
    % uv_fit = polyval(p, y_profile);
    % uv_fit(isnan(uv_profile)) = nan;


    %%% Get Friction Velocity
    u_star.(caze).raw = max(sqrt(-uv_profile), [], 'all', 'omitnan');
    % u_star.(caze).fit = max(sqrt(-uv_fit), [], 'all', 'omitnan');
    u_star.(caze).u_profile = u_profile;
    u_star.(caze).u_profile_normalized = u_profile / (u_inf);
    u_star.(caze).uv_profile = uv_profile;
    u_star.(caze).uv_profile_normalized = uv_profile / (u_inf^2);
    u_star.(caze).y = y_profile;
end

% hold off
% title(caze, 'interpreter', 'none')
% ylim([0, 0.25])


clear caze curvilinear horizontal_lines idx lw max_wave_profile n p phase s u_inf 
clear uv uv_fit uv_profile valid vertical_lines w wave wave_profile y_profile wind_speed
clear cfs_tmp column i I M u u_profile u_stars_tmp x_tmp means curvilinear_path means_path


%% Velocity profiles across different phases: non-normalized

lw = 2;
nu = 1.46E-5;
sz = 20;
wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};

clc; close all;
figure('color', 'white')
hold on

for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};

    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        for phase = 1:4
        
            % Get components
            u = u_star.(caze)(phase).u_profile;
            y = u_star.(caze)(phase).y;
            
            % Plot
            plot(u, y, 'linewidth', lw, 'displayname', sprintf("Phase %1.0f", phase), 'color', wave_colors{phase})
            scatter(u, y, sz, 'filled', 'HandleVisibility', 'off', 'MarkerFaceColor', wave_colors{phase})
        
        end
    end
    
    % No wave
    no_wave_caze = strcat(wind_speed, '_WV0_AGP');

    % Get components
    u = u_star.(no_wave_caze).u_profile;
    y = u_star.(no_wave_caze).y;
    
    % Plot
    plot(u, y, 'linewidth', lw, 'displayname', 'No Waves', 'color', 'black')
    scatter(u, y, sz, 'filled', 'HandleVisibility', 'off', 'MarkerFaceColor', 'black')


end

hold off
title(caze, 'interpreter', 'none')
ylim([0, 0.25])
xlim([1, 6])

clear caze lw no_wave_caze phase s sz u w wave wind_speed y



%% Velocity profiles across different phases: normalized

lw = 2;
nu = 1.46E-5;
sz = 20;
wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};

clc; close all;
figure('color', 'white')
hold on

for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};

    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        for phase = 1:4
        
            % Get components
            u = u_star.(caze)(phase).u_profile_normalized;
            y = u_star.(caze)(phase).y;
            
            % Plot
            plot(u, y, 'linewidth', lw, 'displayname', sprintf("Phase %1.0f", phase), 'color', wave_colors{phase})
            scatter(u, y, sz, 'filled', 'HandleVisibility', 'off', 'MarkerFaceColor', wave_colors{phase})
        
        end
    end
    
    % No wave
    no_wave_caze = strcat(wind_speed, '_WV0_AGP');

    % Get components
    u = u_star.(no_wave_caze).u_profile_normalized;
    y = u_star.(no_wave_caze).y;
    
    % Plot
    plot(u, y, 'linewidth', lw, 'displayname', 'No Waves', 'color', 'black')
    scatter(u, y, sz, 'filled', 'HandleVisibility', 'off', 'MarkerFaceColor', 'black')


end

hold off
ylim([0, 0.25])
xlim([0.4, 1.1])
xlabel('$u / u_{\infty}$', 'interpreter', 'latex', 'FontSize', 16)
ylabel('$y$ [m]', 'interpreter', 'latex', 'FontSize', 16)

clear caze lw no_wave_caze phase s sz u w wave wind_speed y



%% Log law profiles across different phases

lw = 2;
nu = 1.46E-5;
sz = 20;

wind_speed = 'WT8';

clc; close all;
ax = figure('color', 'white');
tiledlayout(1, length(waves))
sgtitle(wind_speed)

for w = 1:length(waves)
    wave = waves{w};
    caze = strcat(wind_speed, '_WV', wave, '_AG0');

    h(w) = nexttile();
    title(wave)
    hold on
    for phase = 1:4
    
        y = u_star.(caze)(phase).y;
        u_profile = u_star.(caze)(phase).u_profile;
        friction_velocity = u_star.(caze)(phase).raw;
       
        % Law of the wall
        u_plus = u_profile / friction_velocity;
        y_plus = (y * friction_velocity) / nu;
    
        % Plot
        plot(y_plus, u_plus, 'linewidth', lw, 'displayname', sprintf("Phase %1.0f", phase), 'color', wave_colors{phase})
        scatter(y_plus, u_plus, sz, 'filled', 'HandleVisibility', 'off', 'MarkerFaceColor', wave_colors{phase})
        
    
    end

    % No-Wave
    no_wave_caze = strcat(wind_speed, '_WV0_AGP');

    y = u_star.(no_wave_caze).y;
    u_profile = u_star.(no_wave_caze).u_profile;
    friction_velocity = u_star.(no_wave_caze).raw;
   
    % Law of the wall
    u_plus = u_profile / friction_velocity;
    y_plus = (y * friction_velocity) / nu;

    % Plot
    plot(y_plus, u_plus, 'linewidth', lw, 'displayname', 'No Wave', 'color', 'black')
    scatter(y_plus, u_plus, sz, 'filled', 'HandleVisibility', 'off', 'MarkerFaceColor', 'black')

    hold off
    xscale('log')
    grid on
    xlabel('$y^+$', 'interpreter', 'latex', 'fontsize', 16)

    if w == 1
        ylabel('$u^+$', 'interpreter', 'latex', 'fontsize', 16)
    end

end

linkaxes(h, 'xy')
xlim([0.5E2, 1E4])
ylim([0, 30])
leg = legend('Orientation', 'Horizontal');
leg.Layout.Tile = 'north';

clear ax caze friction_velocity h leg lw phase sz u_plus y_plus u_profile w wave wind_speed y


%% Log law profiles across different phases: using 'average-of-phases' u*

lw = 2;
nu = 1.46E-5;
sz = 20;

wind_speed = 'WT8';

clc; close all;
ax = figure('color', 'white');
tiledlayout(1, length(waves))
sgtitle(wind_speed)

for w = 1:length(waves)
    wave = waves{w};
    caze = strcat(wind_speed, '_WV', wave, '_AG0');

    h(w) = nexttile();
    title(wave)
    hold on


    %%% Average-of-phases u*
    avg_friction_velocity = nan(1,4);
    for i = 1:4
        avg_friction_velocity(1,i) = u_star.(caze)(i).raw;
    end

    for phase = 1:4
    
        y = u_star.(caze)(phase).y;
        u_profile = u_star.(caze)(phase).u_profile;
        friction_velocity = mean(avg_friction_velocity, 'all', 'omitnan');
       
        % Law of the wall
        u_plus = u_profile / friction_velocity;
        y_plus = (y * friction_velocity) / nu;
    
        % Plot
        plot(y_plus, u_plus, 'linewidth', lw, 'displayname', sprintf("Phase %1.0f", phase), 'color', wave_colors{phase})
        scatter(y_plus, u_plus, sz, 'filled', 'HandleVisibility', 'off', 'MarkerFaceColor', wave_colors{phase})
        
    
    end

    % No-Wave
    no_wave_caze = strcat(wind_speed, '_WV0_AGP');

    y = u_star.(no_wave_caze).y;
    u_profile = u_star.(no_wave_caze).u_profile;
    friction_velocity = u_star.(no_wave_caze).raw;
   
    % Law of the wall
    u_plus = u_profile / friction_velocity;
    y_plus = (y * friction_velocity) / nu;

    % Plot
    plot(y_plus, u_plus, 'linewidth', lw, 'displayname', 'No Wave', 'color', 'black')
    scatter(y_plus, u_plus, sz, 'filled', 'HandleVisibility', 'off', 'MarkerFaceColor', 'black')

    hold off
    xscale('log')
    grid on
    xlabel('$y^+$', 'interpreter', 'latex', 'fontsize', 16)

    if w == 1
        ylabel('$u^+$', 'interpreter', 'latex', 'fontsize', 16)
    end

end

linkaxes(h, 'xy')
xlim([0.5E2, 1E4])
ylim([0, 30])
leg = legend('Orientation', 'Horizontal');
leg.Layout.Tile = 'north';

clear ax caze friction_velocity h leg lw phase sz u_plus y_plus u_profile w wave wind_speed y








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
