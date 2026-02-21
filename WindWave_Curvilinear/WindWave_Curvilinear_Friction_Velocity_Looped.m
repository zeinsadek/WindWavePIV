%%% Curvilinear Friction Velocity

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
            caze = strcat(wind_speed, '_WV', wave, '_AG0');

            % Phase
            plot(normalized_uv_profiles.(caze)(phase).uv, normalized_uv_profiles.(caze)(phase).y, ...
                 'linewidth', lw, 'DisplayName', wave)
        end

        % No-Wave
        no_wave_caze = strcat(wind_speed, '_WV0_AGP');
        plot(normalized_uv_profiles.(no_wave_caze).uv, normalized_uv_profiles.(no_wave_caze).y, ...
             'color', 'black', 'linewidth', lw, 'DisplayName', 'No Waves')

        c = c + 1;
        hold off
    end
end

linkaxes(h, 'xy')
leg = legend('Orientation', 'Horizontal');
leg.Layout.Tile = 'north';

clear c s phase w wave wind_speed h ax caze no_wave_caze leg lw 

%% Plot of phase average uv profiles, averagr of them, and no-wave profile

lw = 2;
ax = figure();
tiledlayout(length(wind_speeds), length(waves))

c = 1;
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        h(c) = nexttile;
        hold on
        title(sprintf("%s %s", wind_speed, wave))
        for phase = 1:4
            % Phase
            plot(normalized_uv_profiles.(caze)(phase).uv, normalized_uv_profiles.(caze)(phase).y, ...
                 'linewidth', lw, 'DisplayName', sprintf("%1.0f", phase))
        end

        % Average of phases
        temp_uv = nan(171,4);
        temp_y = nan(171,4);
        for i = 1:4
            temp_uv(:,i) = normalized_uv_profiles.(caze)(i).uv;
            temp_y(:,i) = normalized_uv_profiles.(caze)(i).y;
        end

        plot(mean(temp_uv, 2, 'omitnan'), mean(temp_y, 2, 'omitnan'), ...
                 'linewidth', lw, 'DisplayName', 'Average of Phases', 'color', 'cyan')

        % No-Wave
        no_wave_caze = strcat(wind_speed, '_WV0_AGP');
        plot(normalized_uv_profiles.(no_wave_caze).uv, normalized_uv_profiles.(no_wave_caze).y, ...
             'color', 'black', 'linewidth', lw, 'DisplayName', 'No Waves')

        c = c + 1;
        hold off
    end
end 

linkaxes(h, 'xy')
leg = legend('Orientation', 'Horizontal');
leg.Layout.Tile = 'north';

clear c s phase w wave wind_speed h ax caze i leg lw no_wave_caze


%% Compare average of phases for each wave, to the no wave case

lw = 3;
clc; close all;
ax = figure();
tiledlayout(1, length(wind_speeds))

c = 1;
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    h(c) = nexttile;

    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        
        hold on
        title(sprintf("%s %s", wind_speed, wave))
        for phase = 1:4
            % Phase
            p = plot(normalized_uv_profiles.(caze)(phase).uv, normalized_uv_profiles.(caze)(phase).y, ...
                     'linewidth', 5, 'handlevisibility', 'off', 'color', wave_colors{w});
            p.Color(4) = 0.1;
        end

        if s == 1
            vis = 'on';
            ylabel('$y$ [m]', 'interpreter', 'latex', 'fontsize', 16)
        else 
            vis = 'off';
        end

        % Average of phases
        temp_uv = nan(171,4);
        temp_y = nan(171,4);
        for i = 1:4
            temp_uv(:,i) = normalized_uv_profiles.(caze)(i).uv;
            temp_y(:,i) = normalized_uv_profiles.(caze)(i).y;
        end

        plot(mean(temp_uv, 2, 'omitnan'), mean(temp_y, 2, 'omitnan'), ...
                 'linewidth', lw, 'DisplayName', wave, 'color', wave_colors{w})

        % No-Wave
        no_wave_caze = strcat(wind_speed, '_WV0_AGP');
        plot(normalized_uv_profiles.(no_wave_caze).uv, normalized_uv_profiles.(no_wave_caze).y, ...
             'color', 'black', 'linewidth', lw, 'DisplayName', 'No Waves', 'HandleVisibility', vis)

        c = c + 1;
        hold off
        xlabel('$\overline{uv} / {u_{\infty}}^2$', 'interpreter', 'latex', 'fontsize', 16)
    end
end 

linkaxes(h, 'xy')
leg = legend('Orientation', 'Horizontal');
leg.Layout.Tile = 'north';

clear ax c caze h i leg lw no_wave_caze p phase s temp_uv temp_y vis w wind_speed wave


%% Plot u* at same phase for different waves


u_star_data = nan(length(waves), length(waves));
cf_data = nan(length(waves), length(waves));

wind_speed = 'WT8';

for phase = 1:4
    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        u_star_data(phase, w) = friction_velocity.(caze)(phase);
        cf_data(phase, w) = skin_friction.(caze)(phase);
    end
end


close all; figure('color', 'white')
b = bar({'Phase 1', 'Phase 2', 'Phase 3', 'Phase 4'}, cf_data);
legend(b, waves)
ylim([0,0.03])
ylabel('Skin Friction Factor $C_f$', 'interpreter', 'latex')
title(wind_speed)
yline(skin_friction.(strcat(wind_speed, '_WV0_AGP')), 'displayname', 'No Waves', 'linewidth', 2)

for i = 1:4
    b(i).FaceColor = wave_colors{i};
end

clear u_star_data cf_data wind_speed b caze i phase temp_uv temp_y w wave

%% Plot u* and Cf averaged across all 4 phases for each wave and each wind speed

u_star_data = nan(length(waves), length(wind_speeds));
cf_data = nan(length(waves), length(wind_speeds));

for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        no_wave_caze = strcat(wind_speed, '_WV0_AGP');
        % Avg quantity
        % u_star_data(w, s) = mean(friction_velocity.(caze)(:));
        % cf_data(w, s) = mean(skin_friction.(caze)(:));

        % Avg quantity - no wave
        u_star_data(w, s) = mean(friction_velocity.(caze)(:)) - friction_velocity.(no_wave_caze);
        cf_data(w, s) = (mean(skin_friction.(caze)(:)) - skin_friction.(no_wave_caze)) / skin_friction.(no_wave_caze);
    end
end



figure()
b = bar(waves, cf_data);
legend(b, {'WT4', 'WT6', 'WT8'})
% ylim([0,0.02])
% ylim([-0.02,0.02])
ylabel('Skin Friction Factor $C_f$', 'interpreter', 'latex')


clear u_star_data cf_data ans s b caze no_wave_caze w wave wind_speed


%% Plot Cf profiles

phase = 1;
wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
linestyles = {'-.', '--', '-'};
lw = 3;
buff = 5;
scale = 1.25;

close all; figure('color', 'white');
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

        data = skin_friction_profiles.(caze)(phase).cf(buff:end-buff);

        wavelength = wavelengths.(wave); 
        amplitude = amplitudes.(wave); 
        steepness = (2 * amplitude) / wavelength;
        
        label = sprintf('$u_{\\infty} = %1.2f, H / \\lambda = %1.3f$', u_inf, steepness);

        x = skin_friction_profiles.(caze)(phase).x;
        plot((x(buff:end-buff) - mean(x, 'all', 'omitnan')) / wavelengths.(wave), data, ...
             'color', wave_colors{w}, 'linestyle', linestyles{s}, 'linewidth', lw, 'DisplayName', label)

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
legend('interpreter', 'latex', 'fontsize', 10, 'location' ,'northeastoutside')
xlabel('$x / \lambda$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$C_f$', 'interpreter', 'latex', 'fontsize', 16)
title(sprintf("Phase %1.0f", phase))
xticks(-1:0.25:1)



%% Plot average of phases Cf against wave steepness: single wind speed

wind_speed = 'WT4';
close all; figure('color', 'white')
hold on
for w = 1:length(waves)

    wave = waves{w};
    wavelength = wavelengths.(wave); 
    amplitude = amplitudes.(wave); 
    steepness = (2 * amplitude) / wavelength;
    disp(steepness)

    for phase = 1:4
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        u_star_data(phase) = friction_velocity.(caze)(phase);
        cf_data(phase) = skin_friction.(caze)(phase);

        scatter(steepness, skin_friction.(caze)(phase), 50, 'filled', 'markerfacecolor', wave_colors{w})
    end

    scatter(steepness, mean(cf_data, 'all', 'omitnan'), 100, 'filled', 'markerfacecolor', 'black')

    clear u_star_data cf_data
end
hold off

clear u_star_data cf_data wind_speed b caze i phase temp_uv temp_y w wave



%% Plot average of phases Cf against wave steepness: all wind speeds

wind_speed_shapes = {'o', '^', 'square'};

close all; figure('color', 'white')
hold on

for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    
    for w = 1:length(waves)
    
        wave = waves{w};
        wavelength = wavelengths.(wave); 
        amplitude = amplitudes.(wave); 
        steepness = (2 * amplitude) / wavelength;
        disp(steepness)
    
        for phase = 1:4
            caze = strcat(wind_speed, '_WV', wave, '_AG0');
            u_star_data(phase) = friction_velocity.(caze)(phase);
            cf_data(phase) = skin_friction.(caze)(phase);
    
            scatter(steepness, skin_friction.(caze)(phase), 100, 'filled', 'markerfacecolor', wave_colors{w})
        end
    
        scatter(steepness, mean(cf_data, 'all', 'omitnan'), 200, 'filled', wind_speed_shapes{s}, 'markerfacecolor', 'black')
    
        clear u_star_data cf_data
    end
end
hold off




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
