%%% Curvilinear Log Law

addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
means_path = fullfile(project_path, 'means');
curvilinear_path = fullfile(project_path, 'curvilinear_new');
wave_parameters = readcell("Offshore_Waves.xlsx");

figure_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/pdf_test7';

% Approximate wavelengths in mm for labeling plots
wavelength_names.A = '410';
wavelength_names.B = '313';
wavelength_names.C = '189';
wavelength_names.D = '124';

steepnesses_names.A = '0.180';
steepnesses_names.B = '0.211';
steepnesses_names.C = '0.305';
steepnesses_names.D = '0.267';

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


            if phase == 1
                stat_grad_tol = 2.9;
            else
                stat_grad_tol = 5;
            end

            %%% CLEAN UP U PROFILE
            u_profile = StatisticalGradientFilter(u_profile, stat_grad_tol, 9);
            u_profile(find(~isnan(u_profile), 1, 'first')) = nan;
        
        
            %%% Get Friction Velocity
            u_star.(caze)(phase).raw = max(sqrt(-uv_profile), [], 'all', 'omitnan');
            % u_star.(caze)(phase).fit = max(sqrt(-uv_fit), [], 'all', 'omitnan');
            u_star.(caze)(phase).u_profile = u_profile;
            u_star.(caze)(phase).u_profile_normalized = u_profile / (u_inf);
            u_star.(caze)(phase).uv_profile = uv_profile;
            u_star.(caze)(phase).uv_profile_normalized = uv_profile / (u_inf^2);
            u_star.(caze)(phase).y = y_profile;
        
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
    %%% NEW UPDATE; DONT CLEAN UP PROFILE FOR NO WAVE CASES
    % u_profile = StatisticalGradientFilter(u_profile, 2.4, 10);
    u_profile(find(~isnan(u_profile), 1, 'first')) = nan;

    %%% Get Friction Velocity
    u_star.(caze).raw = max(sqrt(-uv_profile), [], 'all', 'omitnan');
    % u_star.(caze).fit = max(sqrt(-uv_fit), [], 'all', 'omitnan');
    u_star.(caze).u_profile = u_profile;
    u_star.(caze).u_profile_normalized = u_profile / (u_inf);
    u_star.(caze).uv_profile = uv_profile;
    u_star.(caze).uv_profile_normalized = uv_profile / (u_inf^2);
    u_star.(caze).y = y_profile;
end

% clear amplitude caze column curvilinear displacement displacement_cum displacement_edge_mask displacement_integrand
% clear horizontal_lines i idx island_length left_mask max_wave_profile momentum momentum_cum momentum_edge_mask momentum_integrand
% clear n nan_mask phase right_mask s sgolay_length smooth_kernerl thickness thickness_edge_mask threshold u u_column u_inf u_normalized
% clear vertical_lines w wave wave_profile wave_parameters wave_type wave_length wind_speed x_BL zeta_column
% clear gradient_tolerance project_path smooth_kernel wavelength curvilinear_path
% clc;


% clear caze curvilinear horizontal_lines idx lw max_wave_profile n p phase s u_inf 
% clear uv uv_fit uv_profile valid vertical_lines w wave wave_profile y_profile wind_speed
% clear cfs_tmp column i I M u u_profile u_stars_tmp x_tmp means means_path
% clc;



%% Check u* value for WVB phase 4, WT6

caze = 'WT6_WVD_AG0';
phase = 4;


% Load data
curvilinear = load(fullfile(curvilinear_path, strcat(caze, '_CURVILINEAR.mat')));
curvilinear = curvilinear.output;

 % X (xi) in mm
vertical_lines = curvilinear.phase(phase).vertical_lines * 1E-3;

% Y (zeta) in mm
horizontal_lines = curvilinear.phase(phase).horizontal_lines * 1E-3;
wave_profile = curvilinear.phase(phase).wave_profile * 1E-3;

% Get components
uv = curvilinear.phase(phase).uv;

% Get profiles
idx = round(length(uv)/2);
y_profile = horizontal_lines(:, idx) - wave_profile(idx); 
uv_profile = uv(:, idx); 

figure('color', 'white')
plot(-uv_profile, y_profile)
xlim([0, 0.35])













%% Example of constant slope region

wind_speed_colors = {'#0075F2', '#FF8C42', '#D30C7B'};

% Replace the fixed bounds section with:
karman_constant = 0.41;
c_plus = 5;
nu = 1.46E-5;
smooth_y_plus = 0:10:1E4;
smooth_u_plus = (1/karman_constant) * log(smooth_y_plus) + c_plus;

% Wind speeds
freestreams.('WT4') = 2.4181;
freestreams.('WT6') = 3.8709;
freestreams.('WT8') = 5.4289;

% Plot
clc; close all
figure('color', 'white')
t = tiledlayout(1,3);

% Loop through wind speeds
for s = 1:length(wind_speeds)

    % Get wind speed
    wind_speed = wind_speeds{s};

    % No-Wave
    no_wave_caze = strcat(wind_speed, '_WV0_AGP');
    y = u_star.(no_wave_caze).y;
    u_profile = u_star.(no_wave_caze).u_profile;
    friction_velocity = u_star.(no_wave_caze).raw;

    % Law of the wall
    u_plus = u_profile / friction_velocity;
    y_plus = (y * friction_velocity) / nu;

    % Find constant slope region (log region)
    [bounds, slope_info] = find_constant_slope_region(y_plus, u_plus, ...
                                                      'plateau_tol', 0.2, ...
                                                      'smooth_window', 7);

    % Compute \Delta U^+
    in_region = (y_plus >= bounds(2)) & (y_plus <= bounds(1));
    y_fit = y_plus(in_region);
    u_fit = u_plus(in_region);
    
    % Smooth wall reference
    u_smooth = (1/karman_constant) * log(y_fit) + c_plus;
    
    % Delta U+
    deltaU = mean(u_smooth - u_fit, 'omitnan');
    
    % Printfs
    disp(wind_speed)
    fprintf('Log Region: %3.1f - %3.1f\n', bounds(2), bounds(1))
    fprintf('Delta U = %1.3f\n\n', deltaU)

    h(s) = nexttile;
    title(sprintf('$u_{\\infty} = %1.2f$ m/s', freestreams.(wind_speed)), 'interpreter', 'latex')
    hold on
    P = plot(y_plus, u_plus, 'color', wind_speed_colors{s}, 'linewidth', 2);
    P.Color(4) = 0.5;
    plot(y_plus, u_plus + deltaU, 'LineWidth', 2, 'color', wind_speed_colors{s})
    plot(smooth_y_plus, smooth_u_plus, 'LineStyle', '--', 'color', 'black')
    hold off
    xline(bounds)
    xscale('log')
    xlim([5, 5000])
end

xlabel(t, '$y^+$', 'Interpreter', 'latex')
ylabel(t, '$u^+$', 'Interpreter', 'latex')
close all

%% No-waves Log Law

% lw = 2;
% nu = 1.46E-5;
% sz = 30;
% 
% tickFontSize = 14;
% labelFontSize = 16;
% legendFontSize = 12;
% titleFontSize = 18;
% 
% wind_speed_markers = {'o', 'o', 'o'};
% wind_speed_colors = {'#0075F2', '#FF8C42', '#D30C7B'};
% 
% clc; close all;
% ax = figure('color', 'white', 'units', 'centimeters', 'position', [2,2,20,10]);
% tiledlayout(1,1,'Padding', 'tight')
% ax = nexttile;
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% 
% % Loop through wind speeds
% hold on
% 
% % Smooth wall equation
% karman_constant = 0.41;
% c_plus = 5;
% 
% % Smooth wall reference
% smooth_y_plus = 0:10:1E4;
% smooth_u_plus = (1/karman_constant) * log(smooth_y_plus) + c_plus;
% plot(smooth_y_plus, smooth_u_plus, 'linestyle', '--', 'HandleVisibility', 'off', 'color', 'black', 'linewidth', 2)
% 
% % Make array to hold onto log-region bounds
% no_wave_log_regions = nan(3,2);
% 
% % Loop through wind speeds
% for s = 1:length(wind_speeds)
% 
%     % Get wind speed
%     wind_speed = wind_speeds{s};
% 
%     if ismember(wind_speed(end), {'4'})
%         u_inf = 2.4181;
%     elseif ismember(wind_speed(end), {'6'})
%         u_inf = 3.8709;
%     elseif ismember(wind_speed(end), {'8'})
%         u_inf = 5.4289;
%     end
% 
%     % No-Wave
%     no_wave_caze = strcat(wind_speed, '_WV0_AGP');
%     y = u_star.(no_wave_caze).y;
%     u_profile = u_star.(no_wave_caze).u_profile;
%     friction_velocity = u_star.(no_wave_caze).raw;
% 
%     % Law of the wall
%     u_plus = u_profile / friction_velocity;
%     y_plus = (y * friction_velocity) / nu;
% 
%     %%% NEW METHOD
%     % Find constant slope region (log region)
%     [bounds, slope_info] = find_constant_slope_region(y_plus, u_plus, ...
%                                                       'plateau_tol', 0.2, ...
%                                                       'smooth_window', 7);
% 
%     % Save bounds
%     no_wave_log_regions(s,:) = bounds;
% 
%     % Compute \Delta U^+
%     in_region = (y_plus >= bounds(2)) & (y_plus <= bounds(1));
%     y_fit = y_plus(in_region);
%     u_fit = u_plus(in_region);
% 
%     % Smooth wall reference
%     u_smooth = (1/karman_constant) * log(y_fit) + c_plus;
% 
%     % Delta U+
%     deltaU = mean(u_smooth - u_fit, 'omitnan');
% 
%     % Prints
%     disp(wind_speed)
%     fprintf('Log Region: %3.1f - %3.1f\n', bounds(2), bounds(1))
%     fprintf('Delta U = %1.3f\n\n', deltaU)
% 
%     % Plot
%     spacing = 1;
%     label = sprintf('$u_{\\infty} = %1.2f$ m/s, $\\Delta U^+ = %1.2f$',  u_inf, deltaU);
%     % Plot lines
%     plot(y_plus(1:spacing:end), u_plus(1:spacing:end) + deltaU, ...
%          'linewidth', lw, 'displayname', 'No Wave', ...
%          'color', wind_speed_colors{s}, 'HandleVisibility', 'off')
%     % Plot points
%     scatter(y_plus(1:spacing:end), u_plus(1:spacing:end) + deltaU, sz, wind_speed_markers{s}, 'filled', ...
%         'MarkerFaceColor', wind_speed_colors{s}, ...
%         'displayname', label)
%     % Plot log region per wind speed
%     % xline(bounds, 'color', wind_speed_colors{s}, 'alpha', 0.5, 'HandleVisibility', 'off')
% 
% end
% l = xline(mean(no_wave_log_regions, 1, 'omitnan'), 'color', 'black', ...
%           'linewidth', 1, 'HandleVisibility', 'off');
% uistack(l, 'bottom')
% hold off
% 
% legend('location', 'northwest', 'interpreter', 'latex', 'fontsize', legendFontSize, 'box', 'off')
% xscale('log')
% grid on
% xlabel('$y^+$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylabel('$u^+ + \Delta U^+$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% xlim([5, 5E3])


% clear ax caze friction_velocity h leg lw phase sz u_plus y_plus u_profile w wave wind_speed y
% Save figure
% pause(3)
% figure_name = 'LogLaw_NoWave_Ensemble.pdf';
% exportgraphics(ax, fullfile(figure_folder, 'LogLaw', figure_name), 'resolution', 600, 'ContentType', 'image')
% close all


%% Curvilinear log law sorted per phase: SHIFTED

% Font sizes
tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;
titleFontSize = 10;

% What to plot
wind_speed = 'WT6';
phase_names = {'$\varphi = 0$', '$\varphi = \lambda / 4$', '$\varphi = \lambda / 2$', '$\varphi = 3 \lambda / 4$'};
% wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
wave_colors = {'#FE6202', '#DC2680', '#775EEF', '#648FFF'};
lw = 1;
nu = 1.46E-5;
sz = 2;

grid_trans = 0.05;

% Smooth wall reference
smooth_y_plus = 0:10:1E4;
smooth_u_plus = (1/karman_constant) * log(smooth_y_plus) + c_plus;

% Plotting
clc; close all;
ax = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,9]); %#ok<*NASGU>
t = tiledlayout(2, 4, 'padding' ,'tight', 'TileSpacing', 'compact');

% Make array to hold all Log-region bounds
c = 1;
tmp_wave_log_regions = nan(16, 2);

% Loop through phases
for phase = 1:4

    % Start tile
    h(phase) = nexttile();
    hold on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
    title(phase_names{phase}, 'interpreter', 'latex', 'fontsize', titleFontSize)

    % Smooth wall reference
    plot(smooth_y_plus, smooth_u_plus, 'linestyle', '--', 'HandleVisibility', 'off', 'color', 'black', 'linewidth', 1)
    
    % Loop through waves
    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        % Get data
        y = u_star.(caze)(phase).y;
        u_profile = u_star.(caze)(phase).u_profile;
        friction_velocity = u_star.(caze)(phase).raw;
       
        % Law of the wall
        u_plus = u_profile / friction_velocity;
        y_plus = (y * friction_velocity) / nu;


        %%% NEW METHOD
        % Find constant slope region (log region)
        [bounds, ~] = find_constant_slope_region(y_plus, u_plus, ...
                                                          'plateau_tol', 0.2, ...
                                                          'smooth_window', 7);
    
        % Save bounds
        tmp_wave_log_regions(c,:) = bounds;
        wave_log_regions(phase).(caze) = bounds;
    
        % Compute \Delta U^+
        in_region = (y_plus >= bounds(2)) & (y_plus <= bounds(1));
        y_fit = y_plus(in_region);
        u_fit = u_plus(in_region);
        
        % Smooth wall reference
        u_smooth = (1/karman_constant) * log(y_fit) + c_plus;
        
        % Delta U+
        deltaU = mean(u_smooth - u_fit, 'omitnan');

        % Save Delta U+
        wave_DeltaUs(phase).(caze) = deltaU;

        % Plot
        label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses_names.(wave));
        plot(y_plus, u_plus + deltaU, 'linewidth', lw, 'color', wave_colors{w}, 'HandleVisibility', 'off')
        scatter(y_plus, u_plus + deltaU, sz, 'filled', 'HandleVisibility', 'on', 'MarkerFaceColor', wave_colors{w}, 'DisplayName', label)

        if phase > 1
            set(gca, 'YTickLabel', [])
        end
        axis square
        
        % Increment case counter
        c = c + 1;
    end

    % No-Wave
    no_wave_caze = strcat(wind_speed, '_WV0_AGP');
    y = u_star.(no_wave_caze).y;
    u_profile = u_star.(no_wave_caze).u_profile;
    friction_velocity = u_star.(no_wave_caze).raw;

    % Law of the wall
    u_plus = u_profile / friction_velocity;
    y_plus = (y * friction_velocity) / nu;

    %%% NEW METHOD
    % Find constant slope region (log region)
    [bounds, slope_info] = find_constant_slope_region(y_plus, u_plus, ...
                                                      'plateau_tol', 0.2, ...
                                                      'smooth_window', 7);

    % Compute \Delta U^+
    in_region = (y_plus >= bounds(2)) & (y_plus <= bounds(1));
    y_fit = y_plus(in_region);
    u_fit = u_plus(in_region);
    
    % Smooth wall reference
    u_smooth = (1/karman_constant) * log(y_fit) + c_plus;
    
    % Delta U+
    deltaU = mean(u_smooth - u_fit, 'omitnan');

    if phase == 1
        ylabel('$\langle u_{\xi} \rangle^+ + \Delta U_{\xi}^+$', 'interpreter', 'latex', 'fontsize', labelFontSize)
    end

    hold off
    xscale('log')
    grid on
    tmp = gca;
    tmp.GridAlpha = grid_trans;
    tmp.MinorGridAlpha = grid_trans;
end


linkaxes(h, 'xy')
xlim([1E1, 1E4])
ylim([10, 27])
leg = legend('Orientation', 'horizontal', 'box', 'off', 'interpreter', 'latex', 'fontsize', legendFontSize);
leg.Layout.Tile = 'north';

% Add a centered xlabel across the top row
annotation(ax, 'textbox', [0.3, 0.4, 0.44, 0.04], ...
           'String', '$\zeta^+$', ...
           'Interpreter', 'latex', ...
           'FontSize', labelFontSize, ...
           'HorizontalAlignment', 'center', ...
           'VerticalAlignment', 'bottom', ...
           'EdgeColor', 'none');


% Print useful metrics
clc;
fprintf('Avg Log Region Bounds = %3.2f - %4.2f\n', fliplr(mean(tmp_wave_log_regions, 1, 'omitnan')))
fprintf('Min Bound = %3.2f\nMax bound = %4.2f\n', min(tmp_wave_log_regions, [], 'all'), max(tmp_wave_log_regions, [], 'all'))



%%% Curvilinear local u* bar chart
friction_velocity_bar_chart = nan(4,length(waves));

for phase = 1:4
    hold on
    title(phase_names{phase})
    
    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        % Get data and save to array
        friction_velocity_bar_chart(phase, w) = u_star.(caze)(phase).raw;
    end
end

% No-Wave
no_wave_caze = strcat(wind_speed, '_WV0_AGP');
no_wave_friction_velocity = u_star.(no_wave_caze).raw;
ax2 = nexttile([1,4]);
b = bar({'$\varphi = 0$', '$\varphi = \lambda / 4$', '$\varphi = \lambda / 2$', '$\varphi = 3 \lambda / 4$'}, friction_velocity_bar_chart);
yline(no_wave_friction_velocity, '-', 'No Wave', 'Interpreter', 'latex', 'linewidth', 1, 'fontsize', 6)
yticks(0:0.05:0.5)

% Set bar colors
for k = 1:4
    b(k).FaceColor = wave_colors{k};
    b(k).EdgeColor = 'none';
end

tmp = gca;
tmp.TickLabelInterpreter = 'latex';
tmp.XAxis.FontSize = labelFontSize;
tmp.YAxis.FontSize = tickFontSize;
ylabel('$u_{\xi}^*$ [m/s]', 'interpreter', 'latex', 'fontsize', labelFontSize)
box off
ylim([0, 0.45])
yticks(0:0.1:0.4)


% Add a,b labels
addPanelLabelsFixed(ax, [h(1), ax2], {'a', 'b'}, 'FontSize', 10, 'OffsetPts', [-30,0])

% Save figure
% pause(3)
% figure_name = ['LogLaw_Curvilinear_', wind_speed, '_LogLawCombinedPlot_Shifted.pdf'];
% exportgraphics(ax, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 600, 'ContentType', 'image'); 
% close all




%% Curvilinear log law sorted per phase: UNSHIFTED

% Font sizes
tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;
titleFontSize = 10;

% What to plot
wind_speed = 'WT6';
phase_names = {'$\varphi = 0$', '$\varphi = \lambda / 4$', '$\varphi = \lambda / 2$', '$\varphi = 3 \lambda / 4$'};
% wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
wave_colors = {'#FE6202', '#DC2680', '#775EEF', '#648FFF'};
lw = 1;
nu = 1.46E-5;
sz = 2;

grid_trans = 0.05;

% Smooth wall reference
smooth_y_plus = 0:10:1E4;
smooth_u_plus = (1/karman_constant) * log(smooth_y_plus) + c_plus;

% Plotting
clc; close all;
ax = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,7.5]); %#ok<*NASGU>
t = tiledlayout(2, 4, 'padding' ,'tight', 'TileSpacing', 'compact');

% Make array to hold all Log-region bounds
c = 1;
tmp_wave_log_regions = nan(16, 2);

% Loop through phases
for phase = 1:4

    % Start tile
    h(phase) = nexttile();
    hold on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
    title(phase_names{phase}, 'interpreter', 'latex', 'fontsize', titleFontSize)

    % Smooth wall reference
    plot(smooth_y_plus, smooth_u_plus, 'linestyle', '--', 'HandleVisibility', 'off', 'color', 'black', 'linewidth', 1)
    
    % Loop through waves
    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        % Get data
        y = u_star.(caze)(phase).y;
        u_profile = u_star.(caze)(phase).u_profile;
        friction_velocity = u_star.(caze)(phase).raw;
       
        % Law of the wall
        u_plus = u_profile / friction_velocity;
        y_plus = (y * friction_velocity) / nu;


        %%% NEW METHOD
        % Find constant slope region (log region)
        [bounds, ~] = find_constant_slope_region(y_plus, u_plus, ...
                                                          'plateau_tol', 0.2, ...
                                                          'smooth_window', 7);
    
        % Save bounds
        tmp_wave_log_regions(c,:) = bounds;
        wave_log_regions(phase).(caze) = bounds;
    
        % Compute \Delta U^+
        in_region = (y_plus >= bounds(2)) & (y_plus <= bounds(1));
        y_fit = y_plus(in_region);
        u_fit = u_plus(in_region);
        
        % Smooth wall reference
        u_smooth = (1/karman_constant) * log(y_fit) + c_plus;
        
        % Delta U+
        deltaU = mean(u_smooth - u_fit, 'omitnan');
        if strcmp(wave, 'B') && phase == 4
            disp(deltaU)
        end

        % Save Delta U+
        wave_DeltaUs(phase).(caze) = deltaU;

        % Plot
        label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses_names.(wave));
        plot(y_plus, u_plus, 'linewidth', lw, 'color', wave_colors{w}, 'HandleVisibility', 'off')
        scatter(y_plus, u_plus, sz, 'filled', 'HandleVisibility', 'off', 'MarkerFaceColor', wave_colors{w}, 'DisplayName', label)
        
        if phase > 1
            set(gca, 'YTickLabel', [])
        end
        axis square

        % Increment case counter
        c = c + 1;
    end

    % No-Wave
    no_wave_caze = strcat(wind_speed, '_WV0_AGP');
    y = u_star.(no_wave_caze).y;
    u_profile = u_star.(no_wave_caze).u_profile;
    friction_velocity = u_star.(no_wave_caze).raw;

    % Law of the wall
    u_plus = u_profile / friction_velocity;
    y_plus = (y * friction_velocity) / nu;

    %%% NEW METHOD
    % Find constant slope region (log region)
    [bounds, slope_info] = find_constant_slope_region(y_plus, u_plus, ...
                                                      'plateau_tol', 0.2, ...
                                                      'smooth_window', 7);

    % Compute \Delta U^+
    in_region = (y_plus >= bounds(2)) & (y_plus <= bounds(1));
    y_fit = y_plus(in_region);
    u_fit = u_plus(in_region);
    
    % Smooth wall reference
    u_smooth = (1/karman_constant) * log(y_fit) + c_plus;
    
    % Delta U+
    deltaU = mean(u_smooth - u_fit, 'omitnan');

    if phase == 1
        ylabel('$\langle u_{\xi} \rangle^+$', 'interpreter', 'latex', 'fontsize', labelFontSize)
    end

    hold off
    xscale('log')
    grid on
    tmp = gca;
    tmp.GridAlpha = grid_trans;
    tmp.MinorGridAlpha = grid_trans;
end

% Legend
hold on
for w = 1:length(waves)
    wave = waves{w};
    label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses_names.(wave));

    hLeg = plot(nan, nan, 'o', ...
    'MarkerFaceColor',  wave_colors{w}, ...
    'MarkerEdgeColor','none', ...
    'MarkerSize', 4, ...        % <-- CONTROL LEGEND SIZE HERE
    'LineWidth', 1, ...
    'LineStyle','none', ...
    'DisplayName',label);
end
hold off


linkaxes(h, 'xy')
xlim([1E1, 1E4])
ylim([0, 30])
leg = legend('Orientation', 'horizontal', 'box', 'off', 'interpreter', 'latex', 'fontsize', legendFontSize);
leg.Layout.Tile = 'north';
leg.ItemTokenSize(1) = 10;

% Add a centered xlabel across the top row
annotation(ax, 'textbox', [0.32, 0.35, 0.44, 0.04], ...
           'String', '$\zeta^+$', ...
           'Interpreter', 'latex', ...
           'FontSize', labelFontSize, ...
           'HorizontalAlignment', 'center', ...
           'VerticalAlignment', 'bottom', ...
           'EdgeColor', 'none');


% Print useful metrics
clc;
fprintf('Avg Log Region Bounds = %3.2f - %4.2f\n', fliplr(mean(tmp_wave_log_regions, 1, 'omitnan')))
fprintf('Min Bound = %3.2f\nMax bound = %4.2f\n', min(tmp_wave_log_regions, [], 'all'), max(tmp_wave_log_regions, [], 'all'))



%%% Curvilinear local u* bar chart
friction_velocity_bar_chart = nan(4,length(waves));

for phase = 1:4
    hold on
    title(phase_names{phase})
    
    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        % Get data and save to array
        friction_velocity_bar_chart(phase, w) = u_star.(caze)(phase).raw;
    end
end

% No-Wave
no_wave_caze = strcat(wind_speed, '_WV0_AGP');
no_wave_friction_velocity = u_star.(no_wave_caze).raw;
ax2 = nexttile([1,4]);
b = bar({'$\varphi = 0$', '$\varphi = \lambda / 4$', '$\varphi = \lambda / 2$', '$\varphi = 3 \lambda / 4$'}, friction_velocity_bar_chart);
yline(no_wave_friction_velocity, '-', 'No Wave', 'Interpreter', 'latex', 'linewidth', 1, 'fontsize', 6)
yticks(0:0.05:0.5)

% Set bar colors
for k = 1:4
    b(k).FaceColor = wave_colors{k};
    b(k).EdgeColor = 'none';
end

tmp = gca;
tmp.TickLabelInterpreter = 'latex';
tmp.XAxis.FontSize = labelFontSize;
tmp.YAxis.FontSize = tickFontSize;
ylabel('$u_{\xi}^*$ [m/s]', 'interpreter', 'latex', 'fontsize', labelFontSize)
box off
ylim([0, 0.45])
yticks(0:0.1:0.4)

% Add a,b labels
addPanelLabelsFixed(ax, [h(1), ax2], {'a', 'b'}, 'FontSize', 10, 'OffsetPts', [-30,0])


% Save figure
% pause(3)
% figure_name = ['LogLaw_Curvilinear_', wind_speed, '_LogLawCombinedPlot_Unshifted.pdf'];
% exportgraphics(ax, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 600, 'ContentType', 'image'); 
% close all







%% Curvilinear log law sorted per phase: UNSHIFTED + Rearranged phases

% Font sizes
tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;
titleFontSize = 10;

% What to plot
wind_speed = 'WT6';
phase_names = {'$\varphi = 0$', '$\varphi = \lambda / 4$', '$\varphi = \lambda / 2$', '$\varphi = 3 \lambda / 4$'};
% wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
wave_colors = {'#FE6202', '#DC2680', '#775EEF', '#648FFF'};
lw = 1;
nu = 1.46E-5;
sz = 2;

grid_trans = 0.05;

% Smooth wall reference
smooth_y_plus = 0:10:1E4;
smooth_u_plus = (1/karman_constant) * log(smooth_y_plus) + c_plus;


% Move phases around
reordered_phases = [1,4,3,2];

% Plotting
clc; close all;
ax = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,7.5]); %#ok<*NASGU>
t = tiledlayout(2, 4, 'padding' ,'tight', 'TileSpacing', 'compact');

% Make array to hold all Log-region bounds
c = 1;
tmp_wave_log_regions = nan(16, 2);

% Loop through phases
for p = 1:4

    phase = reordered_phases(p);
    % Start tile
    h(p) = nexttile();
    hold on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
    title(phase_names{p}, 'interpreter', 'latex', 'fontsize', titleFontSize)

    % Smooth wall reference
    plot(smooth_y_plus, smooth_u_plus, 'linestyle', '--', 'HandleVisibility', 'off', 'color', 'black', 'linewidth', 1)
    
    % Loop through waves
    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        % Get data
        y = u_star.(caze)(phase).y;
        u_profile = u_star.(caze)(phase).u_profile;
        friction_velocity = u_star.(caze)(phase).raw;
       
        % Law of the wall
        u_plus = u_profile / friction_velocity;
        y_plus = (y * friction_velocity) / nu;


        %%% NEW METHOD
        % Find constant slope region (log region)
        [bounds, ~] = find_constant_slope_region(y_plus, u_plus, ...
                                                          'plateau_tol', 0.2, ...
                                                          'smooth_window', 7);
    
        % Save bounds
        tmp_wave_log_regions(c,:) = bounds;
        wave_log_regions(phase).(caze) = bounds;
    
        % Compute \Delta U^+
        in_region = (y_plus >= bounds(2)) & (y_plus <= bounds(1));
        y_fit = y_plus(in_region);
        u_fit = u_plus(in_region);
        
        % Smooth wall reference
        u_smooth = (1/karman_constant) * log(y_fit) + c_plus;
        
        % Delta U+
        deltaU = mean(u_smooth - u_fit, 'omitnan');
        if strcmp(wave, 'B') && phase == 4
            disp(deltaU)
        end

        % Save Delta U+
        wave_DeltaUs(phase).(caze) = deltaU;

        % Plot
        label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses_names.(wave));
        plot(y_plus, u_plus, 'linewidth', lw, 'color', wave_colors{w}, 'HandleVisibility', 'off')
        scatter(y_plus, u_plus, sz, 'filled', 'HandleVisibility', 'off', 'MarkerFaceColor', wave_colors{w}, 'DisplayName', label)
        
        if phase > 1
            set(gca, 'YTickLabel', [])
        end
        axis square

        % Increment case counter
        c = c + 1;
    end

    % No-Wave
    no_wave_caze = strcat(wind_speed, '_WV0_AGP');
    y = u_star.(no_wave_caze).y;
    u_profile = u_star.(no_wave_caze).u_profile;
    friction_velocity = u_star.(no_wave_caze).raw;

    % Law of the wall
    u_plus = u_profile / friction_velocity;
    y_plus = (y * friction_velocity) / nu;

    %%% NEW METHOD
    % Find constant slope region (log region)
    [bounds, slope_info] = find_constant_slope_region(y_plus, u_plus, ...
                                                      'plateau_tol', 0.2, ...
                                                      'smooth_window', 7);

    % Compute \Delta U^+
    in_region = (y_plus >= bounds(2)) & (y_plus <= bounds(1));
    y_fit = y_plus(in_region);
    u_fit = u_plus(in_region);
    
    % Smooth wall reference
    u_smooth = (1/karman_constant) * log(y_fit) + c_plus;
    
    % Delta U+
    deltaU = mean(u_smooth - u_fit, 'omitnan');

    if p == 1
        ylabel('$\langle u_{\xi} \rangle^+$', 'interpreter', 'latex', 'fontsize', labelFontSize)
    end

    hold off
    xscale('log')
    grid on
    tmp = gca;
    tmp.GridAlpha = grid_trans;
    tmp.MinorGridAlpha = grid_trans;
end

% Legend
hold on
for w = 1:length(waves)
    wave = waves{w};
    label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses_names.(wave));

    hLeg = plot(nan, nan, 'o', ...
    'MarkerFaceColor',  wave_colors{w}, ...
    'MarkerEdgeColor','none', ...
    'MarkerSize', 4, ...        % <-- CONTROL LEGEND SIZE HERE
    'LineWidth', 1, ...
    'LineStyle','none', ...
    'DisplayName',label);
end
hold off


linkaxes(h, 'xy')
xlim([1E1, 1E4])
ylim([0, 30])
leg = legend('Orientation', 'horizontal', 'box', 'off', 'interpreter', 'latex', 'fontsize', legendFontSize);
leg.Layout.Tile = 'north';
leg.ItemTokenSize(1) = 10;

% Add a centered xlabel across the top row
annotation(ax, 'textbox', [0.32, 0.35, 0.44, 0.04], ...
           'String', '$\zeta^+$', ...
           'Interpreter', 'latex', ...
           'FontSize', labelFontSize, ...
           'HorizontalAlignment', 'center', ...
           'VerticalAlignment', 'bottom', ...
           'EdgeColor', 'none');


% Print useful metrics
clc;
fprintf('Avg Log Region Bounds = %3.2f - %4.2f\n', fliplr(mean(tmp_wave_log_regions, 1, 'omitnan')))
fprintf('Min Bound = %3.2f\nMax bound = %4.2f\n', min(tmp_wave_log_regions, [], 'all'), max(tmp_wave_log_regions, [], 'all'))



%%% Curvilinear local u* bar chart
friction_velocity_bar_chart = nan(4,length(waves));

for p = 1:4
    phase = reordered_phases(p);
    hold on
    title(phase_names{phase})
    
    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        % Get data and save to array
        friction_velocity_bar_chart(p, w) = u_star.(caze)(phase).raw;
    end
end

% No-Wave
no_wave_caze = strcat(wind_speed, '_WV0_AGP');
no_wave_friction_velocity = u_star.(no_wave_caze).raw;
ax2 = nexttile([1,4]);
hold on
b = bar(phase_names, friction_velocity_bar_chart);

yline(no_wave_friction_velocity, '-', 'No Wave', 'Interpreter', 'latex', 'linewidth', 1, 'fontsize', 6)
yticks(0:0.05:0.5)

% Set bar colors
for k = 1:4
    b(k).FaceColor = wave_colors{k};
    b(k).EdgeColor = 'none';
end

% Plot a wave line?
% wave_x = 0.5:0.1:5;
% P = plot(wave_x, 0.04 * cos((2*pi/4) * (wave_x - 1)) + 0.073, 'color', 'black', 'linewidth', 2);
% P.Color(4) = 0.5;
hold off

tmp = gca;

cats = tmp.XAxis.Categories;          % the category labels on the x-axis
tmp.XLim = categorical(cats([1 4]));  % show categories 1 through 4

tmp.XRuler.TickLength = [0 0];   % removes x-axis tick marks only
tmp.TickLabelInterpreter = 'latex';
tmp.XAxis.FontSize = labelFontSize;
tmp.YAxis.FontSize = tickFontSize;
tmp.XAxis.TickLabelInterpreter = 'latex';
ylabel('$u_{\xi}^*$ [m/s]', 'interpreter', 'latex', 'fontsize', labelFontSize)
box off
ylim([0, 0.40])
yticks(0:0.1:0.4)

% Add a,b labels
addPanelLabelsFixed(ax, [h(1), ax2], {'a', 'b'}, 'FontSize', 10, 'OffsetPts', [-30,0])


% Save figure
% pause(3)
% figure_name = ['LogLaw_Curvilinear_', wind_speed, '_LogLawCombinedPlot_Unshifted_Rearranged.pdf'];
% exportgraphics(ax, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 600, 'ContentType', 'image'); 
% close all






%% Try plotting \Delta U^+ against wave phase (1-4)

% Get \Delta U^+ for each wind speed
for phase = 1:4

    for s = 1:length(wind_speeds)
        wind_speed = wind_speeds{s};
        % Loop through waves
        for w = 1:length(waves)
            wave = waves{w};
            caze = strcat(wind_speed, '_WV', wave, '_AG0');
    
            % Get data
            y = u_star.(caze)(phase).y;
            u_profile = u_star.(caze)(phase).u_profile;
            friction_velocity = u_star.(caze)(phase).raw;
           
            % Law of the wall
            u_plus = u_profile / friction_velocity;
            y_plus = (y * friction_velocity) / nu;
    
    
            %%% NEW METHOD
            % Find constant slope region (log region)
            [bounds, ~] = find_constant_slope_region(y_plus, u_plus, ...
                                                              'plateau_tol', 0.2, ...
                                                              'smooth_window', 7);
        
            % Compute \Delta U^+
            in_region = (y_plus >= bounds(2)) & (y_plus <= bounds(1));
            y_fit = y_plus(in_region);
            u_fit = u_plus(in_region);
            
            % Smooth wall reference
            u_smooth = (1/karman_constant) * log(y_fit) + c_plus;
            
            % Delta U+
            deltaU = mean(u_smooth - u_fit, 'omitnan');
    
            % Save Delta U+
            wave_DeltaUs(phase).(caze) = deltaU;
    
        end
    end
end


wind_speed_markers = {'^', 'square', 'o'};
% Wind speeds
freestreams.('WT4') = 2.4181;
freestreams.('WT6') = 3.8709;
freestreams.('WT8') = 5.4289;

steepnesses.A = 0.180;
steepnesses.B = 0.211;
steepnesses.C = 0.305;
steepnesses.D = 0.267;

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


%% Same plot as above but reordered to be shown as (peak, descending, trough, ascending)

wave_transparency = 0.25;

% Font sizes
tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;
titleFontSize = 10;

% Move phases around
reordered_phases = [1,4,3,2];

% Plots sizes
sz = 20;
lw = 1.5;

phase_axes_names = {'$0$', '$\lambda / 4$', '$\lambda / 2$', '$3 \lambda / 4$'};




% Phase \Delta U^+ vs wave phase
clc; close all
figure('color', 'white', 'units', 'centimeters', 'position', [10,10,12,5.5]); %#ok<*NASGU>
tiledlayout(1,1,'padding', 'loose')
ax = nexttile;
hold on
% Plot at each phase averaged over all three wind speeds
for p = 1:4
    phase = reordered_phases(p);
    disp(phase)
    for w = 1:length(waves)
        wave = waves{w};

        tmp = nan(1,3);
        for s = 1:length(wind_speeds)
            wind_speed = wind_speeds{s};
            caze = strcat(wind_speed, '_WV', wave, '_AG0');
            tmp(s) = wave_DeltaUs(phase).(caze);
        end

        
        scatter(p, mean(tmp, 'all', 'omitnan'), sz, 'filled', ...
                'markerfacecolor', wave_colors{w}, ...
                'HandleVisibility', 'off')

        % Start and end plot at wave peak (repeat peak values)
        if p == 1
             scatter(5, mean(tmp, 'all', 'omitnan'), sz, 'filled', ...
                'markerfacecolor', wave_colors{w}, ...
                'HandleVisibility', 'off')
        end

        % Save values to plot as a line
        mean_DeltaUs(phase).(wave) =  mean(tmp, 'all', 'omitnan');
        max_DeltaUs(phase).(wave) = max(tmp, [], 'all', 'omitnan');
        min_DeltaUs(phase).(wave) = min(tmp, [], 'all', 'omitnan');
    end
end


for w = 1:length(waves)
    wave = waves{w};
    for p = 1:4
        phase = reordered_phases(p);
        tmp(p) = mean_DeltaUs(phase).(wave);
        tmp_max(p) = max_DeltaUs(phase).(wave);
        tmp_min(p) = min_DeltaUs(phase).(wave);
    end

    % Repeat peak values
    tmp(5) = tmp(1);
    tmp_max(5) = tmp_max(1); 
    tmp_min(5) = tmp_min(1);

    plot(1:5, tmp, 'linewidth', lw, 'color', wave_colors{w}, 'HandleVisibility', 'off')

    % % Shaded region
    % hFill = patch( ...
    % [1:5, fliplr(1:5)], ...
    % [tmp_max, fliplr(tmp_min)], ...
    % hex2rgb(wave_colors{w}), ...
    % 'FaceAlpha', 0.15, ...        % transparency
    % 'EdgeColor', 'none', ...      % no outline
    % 'HandleVisibility', 'off');   % keep legend clean
    % uistack(hFill, 'bottom')
end

% Mark were smooth wall is
yline(0, 'linestyle', '--', 'Label', 'Smooth Wall', 'interpreter', 'latex', ...
      'fontsize', 6, 'HandleVisibility', 'off', 'linewidth', 1)

% Plot representtive wave profile
fake_x = 0:0.01:6;
amplitude = 2;
offset = -5;
reference_wave = amplitude * cos(((2*pi) / 4) * (fake_x - 1)) + offset;
plot(fake_x, reference_wave, 'color', 'black', 'linewidth', 1.5, 'HandleVisibility', 'off')
scatter(1:5, [offset + amplitude, offset, offset - amplitude, offset, offset + amplitude], sz, 'filled', ...
        'MarkerFaceColor', 'black', 'HandleVisibility', 'off')

% Shade below
X = [fake_x(:); flipud(fake_x(:))];
Y = [reference_wave(:); (-10)*ones(size(reference_wave(:)))];

patch(X, Y, 'k', ...
      'FaceAlpha', wave_transparency, ...
      'EdgeColor', 'none', ...
      'handlevisibility', 'off');


% Add legend showing what colors mean
for w = 1:length(waves)
    wave = waves{w};
    label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses_names.(wave));
    % plot(nan, nan, 'linewidth', lw, 'color', wave_colors{w}, ...
    %      'DisplayName', label)
    hLeg = plot(nan, nan, 'o', ...
    'MarkerFaceColor',  wave_colors{w}, ...
    'MarkerEdgeColor','none', ...
    'MarkerSize', 4, ...        % <-- CONTROL LEGEND SIZE HERE
    'LineWidth', 1, ...
    'LineStyle','none', ...
    'DisplayName',label);
end

leg = legend('interpreter', 'latex', 'box', 'off', 'orientation', 'horizontal', 'fontsize', legendFontSize);
% leg.IconColumnWidth = 19;
leg.Layout.Tile = 'north';
leg.ItemTokenSize(1) = 10;

hold off
xticks(1:4)
xlim([0.5, 6])
ylim([-10 15])
yticks(0:5:15)
ylabel('$\overline{\left(\Delta U_{\xi}^+ \right)}_{u_{\infty}}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlabel('$\varphi$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xticks([1 2 3 4 5])
xticklabels([phase_axes_names, '$\lambda$'])
tmp = gca;
tmp.FontSize = tickFontSize;
set(tmp, 'TickLabelInterpreter', 'latex', 'fontsize', tickFontSize)



% Save figure
% pause(3)
% figure_name = 'LogLaw_Curvilinear_RoughnessFunction_WavePhase.pdf';
% exportgraphics(ax, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 600, 'ContentType', 'image'); 
% close all





%% Average phase \Delta U^+ vs wave age (c/u_{\infty})

% Font sizes
tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;
titleFontSize = 10;

% Plot line sizes
sz = 20;
lw = 1;

% Plot
clc; close all
figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,5]); %#ok<*NASGU>
tiledlayout(1, 1, 'padding', 'tight')
ax = nexttile;
hold on
for w = 1:length(waves)
    wave = waves{w};
    steepness = steepnesses_names.(wave);
    wavelength = wavelengths.(wave);
    frequency = frequencies.(wave);
    wave_speed = (wavelength* 1E-3 * frequency);

    % Hold data to plot line
    xplot = nan(1,3);
    yplot = nan(1,3);

    % Hold data to shade variability
    maxs = nan(1,3);
    mins = nan(1,3);

    % Loop through wind speeds
    for s = 1:length(wind_speeds)
        wind_speed = wind_speeds{s};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        freestream = freestreams.(wind_speed);
        wave_age = wave_speed / freestream;
    
        % Get data across all 4 phases
        tmpY = nan(1,4);
        for phase = 1:4
            tmpY(phase) = wave_DeltaUs(phase).(caze);
        end

        % Plot scatter
        scatter(wave_age, mean(tmpY, 'all', 'omitnan'), sz, 'filled', wind_speed_markers{s}, ...
                'markerfacecolor', wave_colors{w}, 'HandleVisibility', 'off')

        % Save to plot line
        xplot(s) = wave_age;
        yplot(s) = mean(tmpY, 'all', 'omitnan');

        % Save max and min values
        maxs(s) = max(tmpY, [], 'all', 'omitnan');
        mins(s) = min(tmpY, [], 'all', 'omitnan');
    end
    p = plot(xplot, yplot, 'linewidth', lw, 'color', wave_colors{w}, 'HandleVisibility', 'off');

    % Shaded region
    hFill = patch( ...
    [xplot, fliplr(xplot)], ...
    [maxs, fliplr(mins)], ...
    hex2rgb(wave_colors{w}), ...
    'FaceAlpha', 0.15, ...        % transparency
    'EdgeColor', 'none', ...      % no outline
    'HandleVisibility', 'off');   % keep legend clean
    uistack(hFill, 'bottom')
end

% Make legend for line color
for w = 1:length(waves)
    wave = waves{w};
    label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses_names.(wave));
    plot(nan, nan, 'linewidth', 2, 'color', wave_colors{w}, ...
         'DisplayName', label)
end

plot(nan, nan, 'color', 'white', 'displayname', ' ')

% Make a legend for marker type
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    freestream = freestreams.(wind_speed);
    label = sprintf('$u_{\\infty} = %1.2f$ m/s', freestream);
    hLeg = plot(nan, nan, wind_speed_markers{s}, ...
    'MarkerFaceColor','black', ...
    'MarkerEdgeColor','black', ...
    'MarkerSize', 4, ...        % <-- CONTROL LEGEND SIZE HERE
    'LineWidth', 1, ...
    'LineStyle','none', ...
    'DisplayName', label);
end


hold off
leg = legend('interpreter', 'latex', 'box', 'off', 'location', 'eastoutside', 'fontsize', legendFontSize);
leg.IconColumnWidth = 19;

% Mark were smooth wall is
yline(0, 'linestyle', '--', 'Label', 'Smooth Wall', 'interpreter', 'latex', ...
      'fontsize', 6, 'HandleVisibility', 'off', 'linewidth', 1, 'labelhorizontalAlignment', 'left')
axis square
xlim([0.05, 0.35])
ylim([-4, 17])
yticks(-4:4:16)
ylabel('$\overline{\left(\Delta U_{\xi}^+\right)}_{\varphi}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlabel('$c \mathbin{/} u_{\infty}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
tmp = gca;
tmp.FontSize = labelFontSize;
set(tmp, 'TickLabelInterpreter', 'latex')


% Save figure
% pause(3)
% figure_name = 'LogLaw_Curvilinear_RoughnessFunction_WaveAge.pdf';
% exportgraphics(ax, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 600, 'ContentType', 'image'); 
% close all









%% Average phase \Delta U^+ vs friction wave age (c/u*)


% Font sizes
tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;
titleFontSize = 10;

% Plot line sizes
sz = 20;
lw = 1;

% Plot
clc; close all
figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,7]); %#ok<*NASGU>
tiledlayout(1, 1, 'padding', 'tight')
ax = nexttile;
hold on
for w = 1:length(waves)
    wave = waves{w};
    steepness = steepnesses_names.(wave);
    wavelength = wavelengths.(wave);
    frequency = frequencies.(wave);
    wave_speed = (wavelength* 1E-3 * frequency);

    % Hold data to plot line
    xplot = nan(1,3);
    yplot = nan(1,3);

    % Hold data to shade variability
    maxs = nan(1,3);
    mins = nan(1,3);

    % Loop through wind speeds
    for s = 1:length(wind_speeds)
        wind_speed = wind_speeds{s};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        freestream = freestreams.(wind_speed);


        % avg_friction_velocity = mean(LL_Waves.(caze).phase(:).u_star);

        
    
        % Get data across all 4 phases
        friction_velocities = nan(1,4);
        tmpY = nan(1,4);
        for phase = 1:4
            tmpY(phase) = wave_DeltaUs(phase).(caze);
            friction_velocities(phase) = u_star.(caze)(phase).raw;
        end

        avg_friction_velocity = mean(friction_velocities, 'all', 'omitnan');
        wave_age = wave_speed / avg_friction_velocity;

        % Plot scatter
        scatter(wave_age, mean(tmpY, 'all', 'omitnan'), sz, 'filled', wind_speed_markers{s}, ...
                'markerfacecolor', wave_colors{w}, 'HandleVisibility', 'off')

        % Save to plot line
        xplot(s) = wave_age;
        yplot(s) = mean(tmpY, 'all', 'omitnan');

        % Save max and min values
        maxs(s) = max(tmpY, [], 'all', 'omitnan');
        mins(s) = min(tmpY, [], 'all', 'omitnan');
    end
    p = plot(xplot, yplot, 'linewidth', lw, 'color', wave_colors{w}, 'HandleVisibility', 'off');

    % Shaded region
    hFill = patch( ...
    [xplot, fliplr(xplot)], ...
    [maxs, fliplr(mins)], ...
    hex2rgb(wave_colors{w}), ...
    'FaceAlpha', 0.15, ...        % transparency
    'EdgeColor', 'none', ...      % no outline
    'HandleVisibility', 'off');   % keep legend clean
    uistack(hFill, 'bottom')
end

% Make legend for line color
for w = 1:length(waves)
    wave = waves{w};
    label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses_names.(wave));
    plot(nan, nan, 'linewidth', 2, 'color', wave_colors{w}, ...
         'DisplayName', label)
end

plot(nan, nan, 'color', 'white', 'displayname', ' ')

% Make a legend for marker type
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    freestream = freestreams.(wind_speed);
    label = sprintf('$u_{\\infty} = %1.2f$ m/s', freestream);
    hLeg = plot(nan, nan, wind_speed_markers{s}, ...
    'MarkerFaceColor','black', ...
    'MarkerEdgeColor','black', ...
    'MarkerSize', 4, ...        % <-- CONTROL LEGEND SIZE HERE
    'LineWidth', 1, ...
    'LineStyle','none', ...
    'DisplayName', label);
end


hold off
leg = legend('interpreter', 'latex', 'box', 'off', 'location', 'eastoutside', 'fontsize', legendFontSize);
leg.IconColumnWidth = 19;

% Mark were smooth wall is
yline(0, 'linestyle', '--', 'Label', 'Smooth Wall', 'interpreter', 'latex', ...
      'fontsize', 6, 'HandleVisibility', 'off', 'linewidth', 1, 'labelhorizontalAlignment', 'left')
axis square
xlim([1, 7])
ylim([-4, 17])
yticks(-4:4:16)
ylabel('$\overline{\left(\Delta U_{\xi}^+\right)}_{\varphi}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlabel('$c \mathbin{/} \overline{ \left( u_{\xi}^* \right)}_{\varphi}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
tmp = gca;
tmp.FontSize = labelFontSize;
set(tmp, 'TickLabelInterpreter', 'latex')


% Save figure
% pause(3)
% figure_name = 'LogLaw_Curvilinear_RoughnessFunction_FrictionWaveAge.pdf';
% exportgraphics(ax, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 600, 'ContentType', 'image'); 
% close all



%% Plot with power law fit over-top


% Font sizes
tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;
titleFontSize = 10;
% Plot line sizes
sz = 40;
lw = 1;

% Collect all data for power law fit
all_wave_ages = [];
all_deltaU = [];

% Plot
clc; close all
figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,6]); %#ok<*NASGU>
tiledlayout(1, 1, 'padding', 'tight')
ax = nexttile;
hold on
for w = 1:length(waves)
    wave = waves{w};
    steepness = steepnesses_names.(wave);
    wavelength = wavelengths.(wave);
    frequency = frequencies.(wave);
    wave_speed = (wavelength* 1E-3 * frequency);
    % Hold data to plot line
    xplot = nan(1,3);
    yplot = nan(1,3);
    % Hold data to shade variability
    maxs = nan(1,3);
    mins = nan(1,3);
    % Loop through wind speeds
    for s = 1:length(wind_speeds)
        wind_speed = wind_speeds{s};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        freestream = freestreams.(wind_speed);
        % Get data across all 4 phases
        friction_velocities = nan(1,4);
        tmpY = nan(1,4);
        for phase = 1:4
            tmpY(phase) = wave_DeltaUs(phase).(caze);
            friction_velocities(phase) = u_star.(caze)(phase).raw;
        end
        avg_friction_velocity = mean(friction_velocities, 'all', 'omitnan');
        wave_age = wave_speed / avg_friction_velocity;
        % Plot scatter
        scatter(wave_age, mean(tmpY, 'all', 'omitnan'), sz, 'filled', wind_speed_markers{s}, ...
            'markerfacecolor', wave_colors{w}, 'HandleVisibility', 'off')
        % Save to plot line
        xplot(s) = wave_age;
        yplot(s) = mean(tmpY, 'all', 'omitnan');
        % Save max and min values
        maxs(s) = max(tmpY, [], 'all', 'omitnan');
        mins(s) = min(tmpY, [], 'all', 'omitnan');
        
        % Collect for power law fit
        all_wave_ages = [all_wave_ages, wave_age];
        all_deltaU = [all_deltaU, mean(tmpY, 'all', 'omitnan')];
    end
    % p = plot(xplot, yplot, 'linewidth', lw, 'color', wave_colors{w}, 'HandleVisibility', 'off');
    % % Shaded region
    % hFill = patch( ...
    %     [xplot, fliplr(xplot)], ...
    %     [maxs, fliplr(mins)], ...
    %     hex2rgb(wave_colors{w}), ...
    %     'FaceAlpha', 0.15, ...
    %     'EdgeColor', 'none', ...
    %     'HandleVisibility', 'off');
    % uistack(hFill, 'bottom')
end

% Fit power law: deltaU = a * (c/u*)^b
% Take log of both sides: log(deltaU) = log(a) + b*log(c/u*)
% Linear regression in log-space
log_x = log(all_wave_ages);
log_y = log(all_deltaU);

% Remove any NaN or Inf values
valid_idx = isfinite(log_x) & isfinite(log_y);
log_x = log_x(valid_idx);
log_y = log_y(valid_idx);

% Linear fit
p_fit = polyfit(log_x, log_y, 1);
b = p_fit(1);  % exponent
a = exp(p_fit(2));  % coefficient

% Calculate R^2
y_pred = polyval(p_fit, log_x);
SS_res = sum((log_y - y_pred).^2);
SS_tot = sum((log_y - mean(log_y)).^2);
R2 = 1 - SS_res/SS_tot;

% Plot fitted curve
x_fit = linspace(min(all_wave_ages)*0.9, max(all_wave_ages)*1.1, 100);
y_fit = a * x_fit.^b;
P = plot(x_fit, y_fit, 'k-', 'linewidth', 1, 'HandleVisibility', 'off');
uistack(P, 'bottom')

% Print fit results to console
fprintf('Power law fit: ΔU⁺ = %.2f × (c/u*)^(%.2f)\n', a, b);
fprintf('R² = %.3f\n', R2);

% Make legend for line color
for w = 1:length(waves)
    wave = waves{w};
    label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses_names.(wave));
    % plot(nan, nan, 'linewidth', 2, 'color', wave_colors{w}, ...
    %     'DisplayName', label)
    hLeg = plot(nan, nan, 'o', ...
        'MarkerFaceColor', wave_colors{w}, ...
        'MarkerEdgeColor','none', ...
        'MarkerSize', 6, ...
        'LineWidth', 1, ...
        'LineStyle','none', ...
        'DisplayName', label);
end

plot(nan, nan, 'color', 'white', 'displayname', ' ')

% Make a legend for marker type
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    freestream = freestreams.(wind_speed);
    label = sprintf('$u_{\\infty} = %1.2f$ m/s', freestream);
    hLeg = plot(nan, nan, wind_speed_markers{s}, ...
        'MarkerFaceColor','black', ...
        'MarkerEdgeColor','black', ...
        'MarkerSize', 4, ...
        'LineWidth', 1, ...
        'LineStyle','none', ...
        'DisplayName', label);
end

% Add fit equation to legend
% fit_label = sprintf('$%.1f \\cdot (c/u_*^\\xi)^{%.2f}$', a, b);
% plot(nan, nan, 'k-', 'linewidth', 1, 'DisplayName', fit_label)

% Add fit equation and R^2 to plot
text(2.7, 14, sprintf('$\\overline{\\left(\\Delta U_{\\xi}^+\\right)}_{\\varphi} \\approx %.1f \\left( c \\mathbin{/} \\overline{ \\left( u_{\\xi}^* \\right)}_{\\varphi} \\right)^{%.2f}$', a, b), ...
     'interpreter', 'latex', 'fontsize', 8)
text(3.8, 12, sprintf('$R^2 = %.2f$', R2), ...
     'interpreter', 'latex', 'fontsize', 8)

hold off
leg = legend('interpreter', 'latex', 'box', 'off', 'location', 'eastoutside', 'fontsize', legendFontSize);
leg.IconColumnWidth = 19;
% Mark were smooth wall is
% yline(0, 'linestyle', '--', 'Label', 'Smooth Wall', 'interpreter', 'latex', ...
%     'fontsize', 6, 'HandleVisibility', 'off', 'linewidth', 1, 'labelhorizontalAlignment', 'left')
axis square
xlim([1, 7])
% ylim([-4, 17])
ylim([0, 17])
yticks(-4:4:16)
ylabel('$\overline{\left(\Delta U_{\xi}^+\right)}_{\varphi}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlabel('$c \mathbin{/} \overline{ \left( u_{\xi}^* \right)}_{\varphi}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
tmp = gca;
tmp.FontSize = labelFontSize;
set(tmp, 'TickLabelInterpreter', 'latex')


% Save figure
pause(3)
figure_name = 'LogLaw_Curvilinear_RoughnessFunction_FrictionWaveAge_wPowerLaw.pdf';
exportgraphics(ax, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 600, 'ContentType', 'image'); 
close all





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


function [log_mask, diagnostics] = find_log_region(y_plus, u_plus, karman_constant, varargin)
    % Parse optional parameters
    p = inputParser;
    addParameter(p, 'min_yplus', 50, @isnumeric);      % Buffer layer ends ~30-50
    addParameter(p, 'max_yplus_frac', 0.3, @isnumeric); % Fraction of BL thickness
    addParameter(p, 'kappa_tol', 0.15, @isnumeric);     % Tolerance from 1/kappa
    addParameter(p, 'smooth_window', 5, @isnumeric);    % Smoothing for gradient
    parse(p, varargin{:});
    opts = p.Results;
    
    % Compute diagnostic function: Xi = y+ * du+/dy+
    % In log region, Xi = 1/kappa
    log_yplus = log(y_plus);
    
    % Smooth gradient computation (less sensitive to noise)
    if opts.smooth_window > 1
        u_plus_smooth = movmean(u_plus, opts.smooth_window);
    else
        u_plus_smooth = u_plus;
    end
    
    % du+/d(ln y+) = y+ * du+/dy+ in log region should equal 1/kappa
    du_dlog = gradient(u_plus_smooth, log_yplus);
    Xi = du_dlog;  % This is the diagnostic function
    
    target = 1/karman_constant;
    
    % Find where Xi is within tolerance of 1/kappa
    within_kappa = abs(Xi - target) < opts.kappa_tol * target;
    
    % Apply y+ bounds
    yplus_valid = (y_plus > opts.min_yplus) & (y_plus < opts.max_yplus_frac * max(y_plus));
    
    % Combined mask
    log_mask = within_kappa & yplus_valid;
    
    % Find largest contiguous region (avoid scattered points)
    log_mask = find_largest_contiguous(log_mask);
    
    % Diagnostics for inspection
    diagnostics.Xi = Xi;
    diagnostics.target = target;
    diagnostics.within_kappa = within_kappa;
    diagnostics.yplus_valid = yplus_valid;
end

function mask_out = find_largest_contiguous(mask_in)
    % Find the largest contiguous true region
    d = diff([0; mask_in(:); 0]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;
    
    if isempty(starts)
        mask_out = mask_in;
        return
    end
    
    lengths = ends - starts + 1;
    [~, idx] = max(lengths);
    
    mask_out = false(size(mask_in));
    mask_out(starts(idx):ends(idx)) = true;
end


function [log_bounds, slope_info] = find_constant_slope_region(y_plus, u_plus, varargin)
    % Find region where d(u+)/d(ln y+) is approximately constant
    %
    % Returns bounds [ymin, ymax] where the profile has constant logarithmic slope
    
    p = inputParser;
    addParameter(p, 'min_yplus', 30, @isnumeric);       % Exclude buffer layer
    addParameter(p, 'max_yplus_frac', 0.5, @isnumeric); % Exclude outer wake
    addParameter(p, 'smooth_window', 7, @isnumeric);    % Smoothing for gradient
    addParameter(p, 'plateau_tol', 0.10, @isnumeric);   % 10% deviation from median slope
    parse(p, varargin{:});
    opts = p.Results;
    
    % Apply initial bounds
    valid = (y_plus > opts.min_yplus) & (y_plus < opts.max_yplus_frac * max(y_plus));
    
    if sum(valid) < 10
        log_bounds = [opts.min_yplus, opts.max_yplus_frac * max(y_plus)];
        slope_info = struct('kappa_eff', NaN, 'slope_std', NaN);
        return
    end
    
    % Compute slope: du+/d(ln y+)
    log_yplus = log(y_plus);
    
    % Smooth first to reduce noise
    u_smooth = movmean(u_plus, opts.smooth_window);
    
    % Local slope
    slope = gradient(u_smooth, log_yplus);
    
    % Only consider valid region
    slope_valid = slope(valid);
    yplus_valid = y_plus(valid);
    
    % Find where slope is close to its median (i.e., plateau)
    median_slope = median(slope_valid, 'omitnan');
    near_median = abs(slope_valid - median_slope) < opts.plateau_tol * median_slope;
    
    % Find largest contiguous region near median slope
    d = diff([0; near_median(:); 0]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;
    
    if isempty(starts)
        % Fallback: use middle 50% of valid region
        log_bounds = [prctile(yplus_valid, 25), prctile(yplus_valid, 75)];
        slope_info.kappa_eff = 1 / median_slope;
        slope_info.slope_std = std(slope_valid, 'omitnan');
        slope_info.method = 'fallback';
        return
    end
    
    lengths = ends - starts + 1;
    [~, idx] = max(lengths);
    
    region_idx = starts(idx):ends(idx);
    log_bounds = [yplus_valid(region_idx(1)), yplus_valid(region_idx(end))];
    
    % Effective kappa from the constant-slope region
    slope_in_region = slope_valid(region_idx);
    slope_info.kappa_eff = 1 / mean(slope_in_region, 'omitnan');
    slope_info.slope_std = std(slope_in_region, 'omitnan');
    slope_info.median_slope = median_slope;
    slope_info.method = 'constant_slope';
end


function [deltaU, info] = compute_deltaU_robust(y_plus, u_plus, karman_constant, c_plus)
    % Find where slope is constant
    [bounds, slope_info] = find_constant_slope_region(y_plus, u_plus);
    
    % Select that region
    in_region = (y_plus >= bounds(1)) & (y_plus <= bounds(2));
    y_fit = y_plus(in_region);
    u_fit = u_plus(in_region);
    
    % Smooth wall reference
    u_smooth = (1/karman_constant) * log(y_fit) + c_plus;
    
    % Delta U+
    deltaU = mean(u_smooth - u_fit, 'omitnan');
    
    % Package diagnostics
    info.bounds = bounds;
    info.n_points = sum(in_region);
    info.kappa_eff = slope_info.kappa_eff;
    info.kappa_deviation = abs(slope_info.kappa_eff - karman_constant) / karman_constant;
    info.residual_std = std(u_smooth - u_fit, 'omitnan');
end

