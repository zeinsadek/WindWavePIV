%%% Curvilinear Velocity Profile Collapse

addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
curvilinear_path = fullfile(project_path, 'curvilinear');
means_path = fullfile(project_path, 'means');
wave_parameters = readcell("Offshore_Waves.xlsx");

figure_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/pdf_test4';

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
alpha = 1.0;
fs = 24;

wavelengths.A = '410';
wavelengths.B = '313';
wavelengths.C = '189';
wavelengths.D = '124';


%% Non-normalized profiles

tickFontSize = 14;
labelFontSize = 16;
legendFontSize = 12;

phase_line_styles = {'-', '--', '-.', ':'};

ax = figure('color', 'white', 'units', 'centimeters', 'position', [2,4,30,10]);

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)

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
        steepness = (2 * pi * amplitude) / wavelength;
    
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


            % Legend shit; for 1 wind speed and first phase
            if s == 1 && phase == 1
                vis = 'on';
            else
                vis = 'off';
            end


            %%% PLOT
            label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%1.3f}$', wavelengths.(wave), steepness);
            % Normalizing by displacement over thickness
            difference_u = u_inf - u_profile;
            normalized_difference_u = (u_inf - u_profile) /u_inf;
            corrected_normalized_u = (u_inf - u_profile) / (u_inf * (displacement(idx) / thickness(idx)));
            normalized_y = y_profile / thickness(idx);

            P = plot(u_profile, normalized_y, ...
                     'color', wave_colors{w}, 'linestyle', phase_line_styles{phase}, 'linewidth', lw, ...
                     'displayname', label, 'HandleVisibility', vis);
            P.Color(4) = alpha;
        

        end
    end

    % Annotate wind speeds
    x_label_positions = [0.1, 0.33, 0.57];
    y_label_positions = [0.84, 0.82, 0.77];
    annotation(ax, 'textbox', [x_label_positions(s), y_label_positions(s), 0.5, 0.04], ...
            'String', sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf), ...
            'Interpreter', 'latex', ...
            'FontSize', 14, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'EdgeColor', 'none');
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

    % Legend visability
    if s == 1
        vis = 'on';
    else
        vis = 'off';
    end

    %%% PLOT
    P = plot(u_profile, normalized_y, ...
             'color', 'black', 'linestyle', '-', 'linewidth', lw, 'HandleVisibility', vis, 'DisplayName', 'No Wave');
    P.Color(4) = 1.0;
        
end

% Dummy line to add a break
plot(nan, nan, 'color', 'white', 'displayname', '')
% plot(nan, nan, 'color', 'white', 'displayname', '')

% Dummy lines to get legend
% phase_names = {'Peak', 'Ascending', 'Trough', 'Descending'};
phase_names = {'$\varphi = 0$', '$\varphi = \lambda/4$', '$\varphi = \lambda/2$', '$\varphi = 3 \lambda / 4$'};
for phase = 1:4
    P = plot(nan, nan, 'linestyle', phase_line_styles{phase}, ...
         'linewidth', lw, 'color', 'black', 'DisplayName', phase_names{phase});
    % P.Color(4) = 0.5;
end

hold off



ylabel('$\zeta / \delta$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlabel('$\overline{u}, \hspace{1mm} \langle u_{\xi} \rangle$ [m/s]', 'interpreter', 'latex', 'fontsize', labelFontSize)
xticks(1:6)
yticks(0:0.5:2.5)
xlim([1, 6])
ylim([0, 2.5])
legend('location', 'northwest', 'interpreter', 'latex', 'fontsize', legendFontSize, 'box', 'off')





% clear amplitude caze column curvilinear displacement displacement_cum displacement_edge_mask displacement_integrand
% clear horizontal_lines i idx island_length left_mask max_wave_profile momentum momentum_cum momentum_edge_mask momentum_integrand
% clear n nan_mask phase right_mask s sgolay_length smooth_kernerl thickness thickness_edge_mask threshold u u_column u_inf u_normalized
% clear x w wave wave_profile wave_parameters wave_type wave_length wind_speed x_BL zeta_column
% clear gradient_tolerance project_path smooth_kernel wavelength curvilinear_path y_profile u_profile corrected_normalized_u

clc;

% Save figure
pause(3)
figure_name = 'CollapsedProfiles_CurvilinearNonNormalized.pdf';
exportgraphics(ax, fullfile(figure_folder, 'Collapse', figure_name), 'Resolution', 600, 'ContentType', 'image');
close all
fprintf('Generated figure: %s\n\n', figure_name)





%% Traditionally normalized profiles vs Collapsed profiles

tickFontSize = 14;
labelFontSize = 16;
legendFontSize = 12;
titleFontSize = 18;

annotation_dim = [0.5, 0.8, 0.4, 0.4];

clc; close all
ax = figure('color', 'white', 'units', 'centimeters', 'position', [2,2,30,7]);
tile = tiledlayout(1, 2, 'TileSpacing', 'compact');

annotation_x = 0.6;
annotation_y = 0.7;

%%% Traditionally Normalized Profiles
h(1) = nexttile;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
hold on
%%% Wave Profiles
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

            P = plot(normalized_difference_u, normalized_y, ...
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
    P = plot(normalized_difference_u, normalized_y, ...
             'color', 'black', 'linestyle', line_styles{s}, 'linewidth', lw, 'displayname', caze);
    P.Color(4) = 1.0;
        
end
hold off
ylim([-0.1, 2.4])
xlim([-0.05, 0.55])
xticks(0:0.2:0.5)
yticks(0:0.5:2)
ylabel('$\zeta / \delta$', 'interpreter', 'latex', 'fontsize', labelFontSize)

text(annotation_x, annotation_y, '$\frac{u_{\infty} - u}{u_{\infty}}$', 'Interpreter', 'latex', ...
    'Units', 'normalized', 'HorizontalAlignment', 'center', 'fontsize', titleFontSize);


%%% Collapsed Profiles
h(2) = nexttile;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
hold on
%%% Wave Profiles
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
        
end
hold off
ylim([-0.1, 2.2])
xlim([-0.4, 2.7])
xticks(0:0.5:2.7)
yticks(0:0.5:2)

text(annotation_x, annotation_y, '$\frac{u_{\infty} - u}{u_{\infty} \frac{\delta^*}{\delta}}$', 'Interpreter', 'latex', ...
    'Units', 'normalized', 'HorizontalAlignment', 'center', 'fontsize', titleFontSize);

linkaxes(h, 'y')
ylim([0, 2.3])

% Name figure
figure_name = 'CollapsedProfiles_CurvilinearNormalizedDeficit_vs_Collapsed.pdf';

% Save figure
pause(3)
exportgraphics(tile, fullfile(figure_folder, 'Collapse', figure_name), 'Resolution', 600, 'ContentType', 'image'); 
close all
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

