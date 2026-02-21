%%% Curvilinear Log Law

addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
means_path = fullfile(project_path, 'means');
curvilinear_path = fullfile(project_path, 'curvilinear');
waves = {'A', 'B', 'C', 'D'};
wind_speeds = {'WT4', 'WT6', 'WT8'};

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
            normalized_uv_profiles.(wind_speed).(wave)(phase).data = uv_profile / (u_inf^2);
            normalized_uv_profiles.(wind_speed).(wave)(phase).y = horizontal_lines(:, idx) - wave_profile(idx);

        
            %%% Get Friction Velocity
            friction_velocity.(wind_speed).(wave)(phase) = max(sqrt(-uv_profile), [], 'all', 'omitnan');
            skin_friction.(wind_speed).(wave)(phase) = 2 * (friction_velocity.(wind_speed).(wave)(phase) / u_inf)^2;

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
            skin_friction_profiles.(wind_speed).(wave)(phase).x = x_tmp;
            skin_friction_profiles.(wind_speed).(wave)(phase).cf = cfs_tmp;
        
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
    normalized_uv_profiles.(wind_speed).No_Wave.data = uv_profile / (u_inf^2);
    normalized_uv_profiles.(wind_speed).No_Wave.y = means.Y(:,1) * 1E-3;

    %%% Get Friction Velocity
    friction_velocity.(wind_speed).No_Wave = max(sqrt(-uv_profile), [], 'all', 'omitnan');
    skin_friction.(wind_speed).No_Wave = 2 * (friction_velocity.(wind_speed).(wave)(phase) / u_inf)^2;
end

clc;
clear caze curvilinear horizontal_lines idx phase uv uv_profile 
clear vertical_lines w wave wave_profile z_profile wave_type s wind_speed 
clear curvilinear_path max_wave_profile means means_path project_path u_inf
clear wave_parameters wavelength amplitude cfs_tmp column i I M x_tmp


%% Plot shear stress profiles used to get u*
% rows ~ fixed wind speeds
% cols ~ fixed phases

lw = 2; 

ax = figure();
tiledlayout(length(wind_speeds), 4)

c = 1;
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    for phase = 1:4
        h(c) = nexttile;
        hold on
        title(sprintf("%s Phase %1.0f", wind_speed, phase))
        for w = 1:length(waves)
            wave = waves{w};
            % Phase
            plot(normalized_uv_profiles.(wind_speed).(wave)(phase).data, normalized_uv_profiles.(wind_speed).(wave)(phase).y, ...
                 'linewidth', lw, 'DisplayName', wave)
        end
        % No-Wave
        plot(normalized_uv_profiles.(wind_speed).No_Wave.data, normalized_uv_profiles.(wind_speed).No_Wave.y, ...
             'color', 'black', 'linewidth', lw, 'DisplayName', 'No Waves')

        c = c + 1;
        hold off
    end
end

linkaxes(h, 'xy')
leg = legend('Orientation', 'Horizontal');
leg.Layout.Tile = 'north';

clear c s phase w wave wind_speed h ax

%% Plot of phase average uv profiles, averagr of them, and no-wave profile

ax = figure();
tiledlayout(length(wind_speeds), length(waves))

c = 1;
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    for w = 1:length(waves)
        wave = waves{w};
        h(c) = nexttile;
        hold on
        title(sprintf("%s %s", wind_speed, wave))
        for phase = 1:4
            % Phase
            plot(normalized_uv_profiles.(wind_speed).(wave)(phase).data, normalized_uv_profiles.(wind_speed).(wave)(phase).y, ...
                 'linewidth', lw, 'DisplayName', sprintf("%1.0f", phase))
        end

        % Average of phases
        temp_uv = nan(171,4);
        temp_y = nan(171,4);
        for i = 1:4
            temp_uv(:,i) = normalized_uv_profiles.(wind_speed).(wave)(i).data;
            temp_y(:,i) = normalized_uv_profiles.(wind_speed).(wave)(i).y;
        end

        plot(mean(temp_uv, 2, 'omitnan'), mean(temp_y, 2, 'omitnan'), ...
                 'linewidth', lw, 'DisplayName', 'Average of Phases', 'color', 'cyan')

        % No-Wave
        plot(normalized_uv_profiles.(wind_speed).No_Wave.data, normalized_uv_profiles.(wind_speed).No_Wave.y, ...
             'color', 'black', 'linewidth', lw, 'DisplayName', 'No Waves')

        c = c + 1;
        hold off
    end
end 

linkaxes(h, 'xy')
leg = legend('Orientation', 'Horizontal');
leg.Layout.Tile = 'north';

clear c s phase w wave wind_speed h ax



%% Plot u* at same phase for different waves

wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
u_star_data = nan(length(waves), length(waves));
cf_data = nan(length(waves), length(waves));

wind_speed = 'WT8';

for phase = 1:4
    for w = 1:length(waves)
        wave = waves{w};
        u_star_data(phase, w) = friction_velocity.(wind_speed).(wave)(phase);
        cf_data(phase, w) = skin_friction.(wind_speed).(wave)(phase);
    end
end


close all; figure()
b = bar({'Phase 1', 'Phase 2', 'Phase 3', 'Phase 4'}, cf_data);
legend(b, waves)
ylim([0,0.03])
ylabel('Skin Friction Factor $C_f$', 'interpreter', 'latex')
title(wind_speed)

for i = 1:4
    b(i).FaceColor = wave_colors{i};
end

clear u_star_data cf_data wind_speed

%% Plot u* and Cf averaged across all 4 phases for each wave and each wind speed

u_star_data = nan(length(waves), length(wind_speeds));
cf_data = nan(length(waves), length(wind_speeds));

for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    for w = 1:length(waves)
        wave = waves{w};
        % Avg quantity
        % u_star_data(w, s) = mean(friction_velocity.(wind_speed).(wave)(:));
        % cf_data(w, s) = mean(skin_friction.(wind_speed).(wave)(:));

        % Avg quantity - no wave
        u_star_data(w, s) = mean(friction_velocity.(wind_speed).(wave)(:)) - friction_velocity.(wind_speed).No_Wave;
        cf_data(w, s) = mean(skin_friction.(wind_speed).(wave)(:)) - skin_friction.(wind_speed).No_Wave;
    end
end



figure()
b = bar(waves, cf_data);
% ylim([0,0.03])
ylim([-0.02,0.02])
ylabel('Skin Friction Factor $C_f$', 'interpreter', 'latex')


clear u_star_data cf_data


%% Plot Cf profiles

phase = 4;
wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
linestyles = {'-.', '--', '-'};
lw = 2;
buff = 5;
scale = 1.25;

close all; figure();
hold on
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};

    for w = 1:length(waves)
        wave = waves{w};

        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        disp(caze)

        data = skin_friction_profiles.(wind_speed).(wave)(phase).cf(buff:end-buff);

        x = skin_friction_profiles.(wind_speed).(wave)(phase).x;
        plot((x(buff:end-buff) - mean(x, 'all', 'omitnan')) / wavelengths.(wave), data, ...
             'color', wave_colors{w}, 'linestyle', linestyles{s}, 'linewidth', lw, 'DisplayName', strcat(wind_speed, '_WV', wave))

        if wave == 'D' && s == 1
            tmp = plot((wave_profiles.(wave)(phase).x - mean(wave_profiles.(wave)(phase).x)) / wavelengths.(wave), scale *  wave_profiles.(wave)(phase).profile + mean(data, 'all', 'omitnan'), ...
                        'color', 'black', 'linestyle', '--', 'linewidth', lw, 'handlevisibility', 'off');
            tmp.Color(4) = 0.5;
        end

    end
end
hold off
xlim([-1, 1])
ylim([0, 0.03])
legend('interpreter', 'none', 'fontsize', 6)
xlabel('$x / \lambda$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$C_f$', 'interpreter', 'latex', 'fontsize', 16)
title(sprintf("Phase %1.0f", phase))
xticks(-1:0.25:1)





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
