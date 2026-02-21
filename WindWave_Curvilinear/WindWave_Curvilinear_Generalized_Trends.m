%%% Curvilinear Integral Boundary Layer Parameters and Friction Velocity

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
        frequency  = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 4};
    
        % Save wavelength and amplitude
        integral.(caze).wavelength = wavelength * 1E-3;
        integral.(caze).amplitude = amplitude * 1E-3;
        integral.(caze).frequncy = frequency;
        integral.(caze).phase_velocity = (wavelength * 1E-3) * frequency;
    
    
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
            % integral.(caze)(phase).fitted.thickness = PolynomialFit(x_BL, integral.(caze)(phase).filtered.thickness, n);
            % integral.(caze)(phase).fitted.displacement = PolynomialFit(x_BL, integral.(caze)(phase).filtered.displacement, n);
            % integral.(caze)(phase).fitted.momentum = PolynomialFit(x_BL, integral.(caze)(phase).filtered.momentum, n);
            % integral.(caze)(phase).fitted.shape = PolynomialFit(x_BL, integral.(caze)(phase).filtered.shape, n);



            %%% Get friction velocity
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
            u_star.(caze)(phase).u_profile = u_profile;
            u_star.(caze)(phase).u_profile_normalized = u_profile / (u_inf);
            u_star.(caze)(phase).uv_profile = uv_profile;
            u_star.(caze)(phase).uv_profile_normalized = uv_profile / (u_inf^2);
            u_star.(caze)(phase).y = y_profile;
        
        end
    end
end



% No-Wave Cases
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

    x = data.X(1,:) * 1E-3;
    y = data.Y(:,1) * 1E-3;
    u = data.ensemble.u;
    u = juliaan_smooth(u, smooth_kernel);
    uv = data.ensemble.uv;


    % Mask the data properly
    u(x < left_mask | x > right_mask) = nan;
    x(x < left_mask | x > right_mask) = nan;
    % y(x < left_mask | x > right_mask) = nan;
    
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
    % x_BL(thickness_edge_mask) = nan;



    %%% Momentum and displacement thicknesses
    displacement = nan(1, size(u,2));
    momentum = nan(1, size(u,2));
    
    for i = 1:size(u,2)
    
        u_column = u(:,i);
        y = data.Y(:,1) * 1E-3;
        % zeta_column = horizontal_lines(:,i);
    
        % Remove NaNs
        nan_mask = ~isnan(u_column);% & ~isnan(zeta_column);
        u_column = u_column(nan_mask);
        y = y(nan_mask);
    
        % not enough points to integrate
        if length(u_column) < 2
            continue; 
        end
    
        % Normalize velocity
        u_normalized = u_column / u_inf;
            
        % Decreasing ζ
        if y(end) < y(1) 
            y = flipud(y);
            u_normalized = flipud(u_normalized);
        end
    
        % Integrands
        displacement_integrand = (1 - u_normalized);
        momentum_integrand = u_normalized .* (1 - u_normalized);

        % Integrals
        displacement_cum = trapz(y, displacement_integrand);
        momentum_cum = trapz(y, momentum_integrand);
    
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
    % integral.(caze).wave.wave_profile = curvilinear.phase(phase).wave_profile * 1E-3;
    % integral.(caze).wave.x = curvilinear.X(1,:) * 1E-3;
     
    % Save raw signals
    integral.(caze).raw.thickness = thickness;
    integral.(caze).raw.displacement = displacement;
    integral.(caze).raw.momentum = momentum;
    integral.(caze).raw.shape = displacement ./ momentum;
    integral.(caze).x = x;

    % Save filtered
    integral.(caze).filtered.thickness = FilterData(x, thickness, 10, 10, sgolay_length);
    integral.(caze).filtered.displacement = FilterData(x, displacement, gradient_tolerance, island_length, sgolay_length);
    integral.(caze).filtered.momentum = FilterData(x, momentum, gradient_tolerance, island_length, sgolay_length);
    integral.(caze).filtered.shape = FilterData(x, displacement ./ momentum, gradient_tolerance, island_length, sgolay_length);


    %%%% Friction Velocity
    % Get profiles
    idx = round(length(uv)/2);
    y_profile = data.Y(:, 1) * 1E-3; 
    uv_profile = uv(:, idx); 
    u_profile = u(:, idx);

    % CLEAN UP U PROFILE
    u_profile = StatisticalGradientFilter(u_profile, 2.4, 10);
    u_profile(find(~isnan(u_profile), 1, 'first')) = nan;


    %%% Get Friction Velocity
    u_star.(caze).raw = max(sqrt(-uv_profile), [], 'all', 'omitnan');
    u_star.(caze).u_profile = u_profile;
    u_star.(caze).u_profile_normalized = u_profile / (u_inf);
    u_star.(caze).uv_profile = uv_profile;
    u_star.(caze).uv_profile_normalized = uv_profile / (u_inf^2);
    u_star.(caze).y = y_profile;


end



clear amplitude caze column curvilinear displacement displacement_cum displacement_edge_mask displacement_integrand
clear horizontal_lines i idx island_length left_mask max_wave_profile momentum momentum_cum momentum_edge_mask momentum_integrand
clear n nan_mask phase right_mask s sgolay_length smooth_kernerl thickness thickness_edge_mask threshold u u_column u_inf u_normalized
clear vertical_lines w wave wave_profile wave_parameters wave_type wave_length wind_speed x_BL zeta_column frequency
clear gradient_tolerance project_path smooth_kernel wavelength curvilinear_path u_profile uv uv_profile y_profile data
clc;


% %% Check no-wave parameters
% 
% figure
% hold on
% for s = 1:length(wind_speeds)
%     wind_speed = wind_speeds{s};
% 
%     caze = strcat(wind_speed, '_WV0_AGP');
% 
%     thickness = integral.(caze).filtered.thickness;
%     displacement = integral.(caze).filtered.displacement;
%     momentum = integral.(caze).filtered.momentum;
%     x = integral.(caze).x;
% 
%     plot(x, momentum)
% end
% hold off

%% Plot wave-informed free stream against x / lambda

% Quantities to be looped
parameters = {'thickness', 'displacement', 'momentum', 'shape'};
wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
linestyles = {'-.', '--', '-'};

lw = 3;
alpha = 0.5;

profile_smoothing_kernel = 50;


clc; close all
figure('color', 'white')
tiledlayout(1, 4);
% sgtitle(wind_speed)

% Interate through waves
for p = 1:4
    
    h(p) = nexttile;
    hold on

    % Iterate through different wind speeds
    for s = 1:length(wind_speeds)
        wind_speed = wind_speeds{s};
        
        if ismember(wind_speed(end), {'4'})
            u_inf = 2.4181;
        elseif ismember(wind_speed(end), {'6'})
            u_inf = 3.8709;
        elseif ismember(wind_speed(end), {'8'})
            u_inf = 5.4289;
        end

        % Phase angle
        phase_angle = (pi / 2) * (p - 1);

        % Reference wave profile
        tst = linspace(-1, 1, 50);
        plot(tst, -0.07 * cos(2 * pi * (tst - (p - 1)/4)) + 0.17, 'color', 'black', 'linewidth', lw)

        % Iterate through waves
        for w = 1:length(waves)
            wave = waves{w};
            caze = strcat(wind_speed, '_WV', wave, '_AG0');
            disp(caze)
    
            % Plot data
            wavelength = integral.(caze)(1).wavelength;
            waveheight = 2 * integral.(caze)(1).amplitude;
            geometric_steepness = waveheight / wavelength;

            displacement = smoothdata(integral.(caze)(p).filtered.displacement, 'gaussian', profile_smoothing_kernel);
            thickness = smoothdata(integral.(caze)(p).filtered.thickness, 'gaussian', profile_smoothing_kernel);
            wave_informed_freestream = (displacement ./ thickness);
    
            H(w) = plot((integral.(caze)(p).x - mean(integral.(caze)(p).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, wave_informed_freestream, ...
                        'color', wave_colors{w}, 'linewidth', lw, 'linestyle', linestyles{s}, ...
                        'displayname', sprintf("Wave %s", wave));
    
           title(sprintf("Phase %1.0f", p))

           ylim([0, 0.35])
        end
    end
    hold off
end

% leg = legend(H, 'Orientation', 'Horizontal');
% leg.Layout.Tile = 'north';
% set(leg,'Box','off')

% Share axes for same parameter across different plots
for i = 1:4
    linkaxes(h, 'xy')
end
clc;

clear alpha caze count data h H i leg lw p v reference scales c version_colors verrions w wave wind_speed ylim_max
clear displacement geometric_steepness linestyles s thickness u_inf versions wave_informed_freestream waveheight wavelength



%% Plotting displacement/thickness against friction velocity / wave velocity

fs = 18;
sz = 75;
markers = {"^", "square", "o"};
wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};



phase_sizes = [50, 75, 100, 125];

% shape = wind speed
% color = wave
% size = phase

figure('color', 'white')
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

        wavelength = integral.(caze)(1).wavelength;
        amplitude = integral.(caze)(1).amplitude;
        phase_velocity = integral.(caze)(1).phase_velocity;

        wave_steepness = (2 * amplitude) / wavelength;


        % Phase 4 is throwing things off. Cf values may be off in this
        % region of the wave
        for p = 1:3
            friction_velocity = u_star.(caze)(p).raw;
            skin_friction = 2 * (friction_velocity / u_inf)^2;

            thickness = integral.(caze)(p).filtered.thickness;
            displacement = integral.(caze)(p).filtered.displacement;

            idx = round(length(thickness) / 2);

            scatter(skin_friction, displacement(idx) / thickness(idx), phase_sizes(p), markers{s}, 'filled', ...
                    'MarkerEdgeColor', 'none', 'MarkerFaceColor', wave_colors{w})
        end
    end

    title('Phases 1 - 3')
    clear friction_velocity skin_friction thickness displacement

    %%% No Wave Case
    no_wave_caze = strcat(wind_speed, '_WV0_AGP');

    friction_velocity = u_star.(no_wave_caze).raw;
    skin_friction = 2 * (friction_velocity / u_inf)^2;

    thickness = integral.(no_wave_caze).filtered.thickness;
    displacement = integral.(no_wave_caze).filtered.displacement;

    idx = round(length(thickness) / 2);

    scatter(skin_friction, displacement(idx) / thickness(idx), 100, markers{s}, 'filled', ...
            'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'black')

    % scatter(skin_friction, mean(displacement, 'all', 'omitnan') / mean(thickness, 'all', 'omitnan'), 100, markers{s}, 'filled', ...
    %         'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'black')
    

end
hold off

xlabel('$C_f$', 'Interpreter', 'latex', 'fontsize', fs)
ylabel('$\frac{\delta^*}{\delta}$', 'Interpreter', 'latex', 'fontsize', fs)
xlim([0.002, 0.018])
ylim([0.06, 0.26])

clear amplitude caze displacement friction_velocity idx markers p phase_velocity s sz thickness u_inf w 
clear wave wave_steepness wavelength wind_speed fs no_wave_caze skin_friction y
clc;







%% Plotting displacement/thickness against wave steepness
% All the phase will be plotted against the same wave steepness

fs = 18;
sz = 75;
markers = {"^", "square", "o"};
wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};

profile_smoothing_kernel = 10;

phase_sizes = [50, 75, 100, 125];

% shape = wind speed
% color = wave
% size = phase

figure('color', 'white')
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

        % Wave properties
        wavelength = integral.(caze)(1).wavelength;
        amplitude = integral.(caze)(1).amplitude;
        phase_velocity = integral.(caze)(1).phase_velocity;
        wave_steepness = (2 * amplitude) / wavelength;


        % Also average the displacment / thickness across the 4 phases
        tmp = [];

        % Loop through phases
        for p = 1:4
            friction_velocity = u_star.(caze)(p).raw;
            skin_friction = 2 * (friction_velocity / u_inf)^2;

            thickness = smoothdata(integral.(caze)(p).filtered.thickness, 'gaussian', profile_smoothing_kernel);
            displacement = smoothdata(integral.(caze)(p).filtered.displacement, 'gaussian', profile_smoothing_kernel);
            idx = round(length(thickness) / 2);

            % Important ratio
            displacement_over_thickness = displacement(idx) / thickness(idx);
            
            % Phase angle
            phase_angle = (pi*(p - 1) / 2);

            % Scaling term
            scaling_term = 1*cos(phase_angle) + 1;

            % Phase scaled parameter
            scaled_displacement_over_thickness = displacement_over_thickness;% * scaling_term;

            % Save to array to average
            tmp = [scaled_displacement_over_thickness, tmp];

            scatter(wave_steepness, scaled_displacement_over_thickness, phase_sizes(p), markers{s}, 'filled', ...
                    'MarkerEdgeColor', 'none', 'MarkerFaceColor', wave_colors{w})
        end

        % Plot average of phases
        scatter(wave_steepness, mean(tmp), 200, markers{s}, 'filled', ...
                    'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'black')
    end

    clear friction_velocity skin_friction thickness displacement

    %%% No Wave Case
    % no_wave_caze = strcat(wind_speed, '_WV0_AGP');
    % 
    % friction_velocity = u_star.(no_wave_caze).raw;
    % skin_friction = 2 * (friction_velocity / u_inf)^2;
    % 
    % thickness = integral.(no_wave_caze).filtered.thickness;
    % displacement = integral.(no_wave_caze).filtered.displacement;
    % 
    % idx = round(length(thickness) / 2);
    % 
    % scatter(skin_friction, displacement(idx) / thickness(idx), 100, markers{s}, 'filled', ...
    %         'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'black')

end
hold off

xlabel('$\frac{H}{\lambda}$', 'Interpreter', 'latex', 'fontsize', fs)
ylabel('$\frac{\delta^*}{\delta}$', 'Interpreter', 'latex', 'fontsize', fs)
xlim([0.05, 0.105])
ylim([0.05, 0.26])


clear amplitude caze displacement friction_velocity idx markers p phase_velocity s sz thickness u_inf w 
clear wave wave_steepness wavelength wind_speed fs no_wave_caze skin_friction y
clc;



%% Plotting displacement/thickness against Reynolds number from lambda

fs = 18;
sz = 75;
markers = {"^", "square", "o"};
wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};


phase_sizes = [25, 100, 200, 400];

nu = 1.6E-5;

% shape = wind speed
% color = wave
% size = phase

close all;
figure('color', 'white')
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

        % Wave properties
        wavelength = integral.(caze)(1).wavelength;
        amplitude = integral.(caze)(1).amplitude;
        phase_velocity = integral.(caze)(1).phase_velocity;
        wave_steepness = (2 * amplitude) / wavelength;

        % Relative velocities
        % relative_velocity = u_inf - phase_velocity;
        normalized_relative_velocity = (u_inf / phase_velocity);

        % Wavelength reynolds number
        Re_lambda =  (phase_velocity * wavelength) / nu;


        % Also average the displacment / thickness across the 4 phases
        tmp = [];

        % Loop through phases
        for p = 1:4
            friction_velocity = u_star.(caze)(p).raw;
            skin_friction = 2 * (friction_velocity / u_inf)^2;
            thickness = integral.(caze)(p).filtered.thickness;
            displacement = integral.(caze)(p).filtered.displacement;
            idx = round(length(thickness) / 2);

            % Important ratio
            displacement_over_thickness = displacement(idx) / thickness(idx);

            % Phase scaled parameter
            scaled_displacement_over_thickness = displacement_over_thickness; % * scaling_term;

            % Save to array to average
            tmp = [scaled_displacement_over_thickness, tmp];

            scatter(Re_lambda, scaled_displacement_over_thickness, phase_sizes(p), markers{s}, 'filled', ...
                    'MarkerEdgeColor', 'none', 'MarkerFaceColor', wave_colors{w})
        end

        % Plot average of phases
        scatter(Re_lambda, mean(tmp), 200, markers{s}, 'filled', ...
                    'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'black')

    end

    clear friction_velocity skin_friction thickness displacement

end
hold off

% xlabel('$\frac{H}{\lambda}$', 'Interpreter', 'latex', 'fontsize', fs)
% ylabel('$\frac{\delta^*}{\delta}$', 'Interpreter', 'latex', 'fontsize', fs)
% xlim([0.05, 0.105])
% ylim([0.05, 0.26])

xscale('log')

clear amplitude caze displacement friction_velocity idx markers p phase_velocity s sz thickness u_inf w 
clear wave wave_steepness wavelength wind_speed fs no_wave_caze skin_friction y
clc;




%% Considering phase as a cyclic variable

fs = 18;
sz = 75;
markers = {"^", "square", "o"};
wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};


phase_sizes = [50, 75, 100, 125];

% shape = wind speed
% color = wave
% size = phase

figure('color', 'white')
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

        wavelength = integral.(caze)(1).wavelength;
        amplitude = integral.(caze)(1).amplitude;
        phase_velocity = integral.(caze)(1).phase_velocity;

        wave_steepness = (2 * amplitude) / wavelength;


        % Also average the displacment / thickness across the 4 phases
        tmp = [];
        for p = 1:4
            friction_velocity = u_star.(caze)(p).raw;
            skin_friction = 2 * (friction_velocity / u_inf)^2;

            thickness = integral.(caze)(p).filtered.thickness;
            displacement = integral.(caze)(p).filtered.displacement;

            idx = round(length(thickness) / 2);

            tmp = [displacement(idx) / thickness(idx), tmp];

            % Convert phase into a phase angle
            phase_angle = (pi*(p - 1) / 2);

            % Displacement / Thickness
            displacement_over_thickness = displacement(idx) / thickness(idx);

            % wave_steepness * sin(phase_angle)

            scatter(sin(phase_angle), displacement_over_thickness, phase_sizes(p), markers{s}, 'filled', ...
                    'MarkerEdgeColor', 'none', 'MarkerFaceColor', wave_colors{w})
        end

        % Plot average of phases
        scatter(0, mean(tmp), phase_sizes(p), markers{s}, 'filled', ...
                    'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'black')
    end

    clear friction_velocity skin_friction thickness displacement

end
hold off

xlabel('Phase Index')
% xlabel('$\frac{H}{\lambda}$', 'Interpreter', 'latex', 'fontsize', fs)
ylabel('$\left(\frac{\delta^*}{\delta}\right)$', 'Interpreter', 'latex', 'fontsize', fs)

% xlim([0.5, 4.5])
% xticks(1:4)

clear amplitude caze displacement friction_velocity idx markers p phase_velocity s sz thickness u_inf w 
clear wave wave_steepness wavelength wind_speed fs no_wave_caze skin_friction y
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
