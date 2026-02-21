%%% Curvilinear Log Law

addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
means_path = fullfile(project_path, 'means');
curvilinear_path = fullfile(project_path, 'curvilinear');
wave_parameters = readcell("Offshore_Waves.xlsx");

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

% hold off
% title(caze, 'interpreter', 'none')
% ylim([0, 0.25])

clear amplitude caze column curvilinear displacement displacement_cum displacement_edge_mask displacement_integrand
clear horizontal_lines i idx island_length left_mask max_wave_profile momentum momentum_cum momentum_edge_mask momentum_integrand
clear n nan_mask phase right_mask s sgolay_length smooth_kernerl thickness thickness_edge_mask threshold u u_column u_inf u_normalized
clear vertical_lines w wave wave_profile wave_parameters wave_type wave_length wind_speed x_BL zeta_column
clear gradient_tolerance project_path smooth_kernel wavelength curvilinear_path
clc;


clear caze curvilinear horizontal_lines idx lw max_wave_profile n p phase s u_inf 
clear uv uv_fit uv_profile valid vertical_lines w wave wave_profile y_profile wind_speed
clear cfs_tmp column i I M u u_profile u_stars_tmp x_tmp means means_path
clc;


%% No-waves Log Law

lw = 2;
nu = 1.46E-5;
sz = 30;

tickFontSize = 14;
labelFontSize = 16;
legendFontSize = 12;
titleFontSize = 18;

wind_speed_markers = {'o', 'o', 'o'};
wind_speed_colors = {'#0075F2', '#FF8C42', '#D30C7B'};

clc; close all;
ax = figure('color', 'white', 'units', 'centimeters', 'position', [2,2,20,10]);
tiledlayout(1,1,'Padding', 'tight')
ax = nexttile;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)

% Loop through wind speeds
hold on

% Smooth wall equation
karman_constant = 0.41;
c_plus = 5;

smooth_y_plus = 0:10:1E4;
smooth_u_plus = (1/karman_constant) * log(smooth_y_plus) + c_plus;
plot(smooth_y_plus, smooth_u_plus, 'linestyle', '--', 'HandleVisibility', 'off', 'color', 'black', 'linewidth', 2)

% Loop through wind speeds
for s = 1:length(wind_speeds)

    % Get wind speed
    wind_speed = wind_speeds{s};

    if ismember(wind_speed(end), {'4'})
        u_inf = 2.4181;
    elseif ismember(wind_speed(end), {'6'})
        u_inf = 3.8709;
    elseif ismember(wind_speed(end), {'8'})
        u_inf = 5.4289;
    end

    % No-Wave
    no_wave_caze = strcat(wind_speed, '_WV0_AGP');
    y = u_star.(no_wave_caze).y;
    u_profile = u_star.(no_wave_caze).u_profile;
    friction_velocity = u_star.(no_wave_caze).raw;

   
    % Law of the wall
    u_plus = u_profile / friction_velocity;
    y_plus = (y * friction_velocity) / nu;


    % Calculate vertical shift from smooth wall
    ymin = 600;
    ymax = min(0.8 * max(y_plus), 5000); % up to 30% of max, capped at 5000
    log_region = (y_plus > ymin & y_plus < ymax);

    % Difference (measured - smooth)
    same_grid_smooth_u_plus = (1/karman_constant) * log(y_plus) + c_plus;
    deltaU = mean(same_grid_smooth_u_plus(log_region) - u_plus(log_region), 'all', 'omitnan');
    disp(deltaU)
    % deltaU = 0;


    % Plot
    spacing = 1;
    label = sprintf('$u_{\\infty} = %1.2f$ m/s, $\\Delta U^+ = %1.2f$',  u_inf, deltaU);
    plot(y_plus(1:spacing:end), u_plus(1:spacing:end) + deltaU, 'linewidth', lw, 'displayname', 'No Wave', 'color', wind_speed_colors{s}, 'HandleVisibility', 'off')
    scatter(y_plus(1:spacing:end), u_plus(1:spacing:end) + deltaU, sz, wind_speed_markers{s}, 'filled', 'MarkerFaceColor', wind_speed_colors{s}, 'displayname', label)

end
hold off

legend('location', 'northwest', 'interpreter', 'latex', 'fontsize', legendFontSize, 'box', 'off')
xscale('log')
grid on
xlabel('$y^+$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$u^+ + \Delta U^+$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylabel('$u^+$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlim([5, 5E3])

% clear ax caze friction_velocity h leg lw phase sz u_plus y_plus u_profile w wave wind_speed y


% Save figure
% pause(3)
% figure_name = 'LogLaw_NoWave_Ensemble.pdf';
% exportgraphics(ax, fullfile(figure_folder, 'LogLaw', figure_name), 'resolution', 600, 'ContentType', 'image')
% close all


%% Curvilinear log law sorted per phase


tickFontSize = 14;
labelFontSize = 16;
legendFontSize = 12;
titleFontSize = 18;

wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};

lw = 1;
nu = 1.46E-5;
sz = 10;

wind_speed = 'WT6';
phase_names = {'$\varphi = 0$', '$\varphi = \lambda / 4$', '$\varphi = \lambda / 2$', '$\varphi = 3 \lambda / 4$'};

clc; close all;
ax = figure('color', 'white', 'units', 'centimeters', 'position', [2,2,30,8]);
t = tiledlayout(1, 4);

for phase = 1:4

    h(phase) = nexttile();
    hold on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
    title(phase_names{phase}, 'interpreter', 'latex', 'fontsize', 16)

    smooth_y_plus = 0:10:1E4;
    smooth_u_plus = (1/karman_constant) * log(smooth_y_plus) + c_plus;
    plot(smooth_y_plus, smooth_u_plus, 'linestyle', '--', 'HandleVisibility', 'off', 'color', 'black')
    
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

        % Calculate vertical shift from smooth wall
        ymin = 900;
        ymax = min(0.9 * max(y_plus), 5000); % up to 30% of max, capped at 5000
        log_region = (y_plus > ymin & y_plus < ymax);

        % Difference (measured - smooth)
        same_grid_smooth_u_plus = (1/karman_constant) * log(y_plus) + c_plus;
        deltaU = mean(u_plus(log_region) - same_grid_smooth_u_plus(log_region), 'all', 'omitnan');
        disp(deltaU)
        deltaU = 0;
    
        % Plot
        label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses.(wave));
        plot(y_plus, u_plus - deltaU, 'linewidth', lw, 'color', wave_colors{w}, 'HandleVisibility', 'off')
        scatter(y_plus, u_plus - deltaU, sz, 'filled', 'HandleVisibility', 'on', 'MarkerFaceColor', wave_colors{w}, 'DisplayName', label)
        
    
    end

    % No-Wave
    no_wave_caze = strcat(wind_speed, '_WV0_AGP');

    y = u_star.(no_wave_caze).y;
    u_profile = u_star.(no_wave_caze).u_profile;
    friction_velocity = u_star.(no_wave_caze).raw;

    % Law of the wall
    u_plus = u_profile / friction_velocity;
    y_plus = (y * friction_velocity) / nu;

    % Calculate vertical shift from smooth wall
    ymin = 1200;
    ymax = min(1 * max(y_plus), 5000);
    log_region = (y_plus > ymin & y_plus < ymax);

    % Difference (measured - smooth)
    same_grid_smooth_u_plus = (1/karman_constant) * log(y_plus) + c_plus;
    deltaU = mean(u_plus(log_region) - same_grid_smooth_u_plus(log_region), 'all', 'omitnan');
    disp(deltaU)
    deltaU = 0;


    P1 = plot(y_plus, u_plus - deltaU, 'linewidth', lw, 'displayname', 'No Wave', 'color', 'black', 'handlevisibility', 'off');
    scatter(y_plus, u_plus - deltaU, sz, 'filled', 'HandleVisibility', 'on', 'MarkerFaceColor', 'black', 'displayname', 'No Wave');

    hold off
    xscale('log')
    grid on

    if phase == 1
        ylabel('$\langle u_{\xi} \rangle^+ + \Delta U^+$', 'interpreter', 'latex', 'fontsize', 16)
    end
end

xlabel(t, '$\zeta^+$', 'interpreter', 'latex', 'fontsize', 16)

linkaxes(h, 'xy')
% ylim([10, 26])
xlim([1E1, 1E4])
leg = legend('Orientation', 'Vertical', 'box', 'off', 'interpreter', 'latex', 'fontsize', 18);
leg.Layout.Tile = 'East';


% Name figure
% figure_name = ['LogLaw_Curvilinear_', wind_speed, '_Profiles.pdf'];
% pause(3)
% exportgraphics(ax, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 600, 'ContentType', 'image'); 
% close all




%% Curvilinear local u* bar chart


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

clc; close all;
ax = figure('color', 'white', 'units', 'centimeters', 'position', [2,2,30,10]);
b = bar({'$\varphi = 0$', '$\varphi = \lambda / 4$', '$\varphi = \lambda / 2$', '$\varphi = 3 \lambda / 4$'}, friction_velocity_bar_chart);
yline(no_wave_friction_velocity, '--', 'No Wave', 'Interpreter', 'latex')
yticks(0:0.05:0.5)

% Set bar colors
for k = 1:4
    b(k).FaceColor = wave_colors{k};
    b(k).EdgeColor = 'none';
end

ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = tickFontSize;
ylabel('$u^*$ [m/s]', 'interpreter', 'latex', 'fontsize', labelFontSize)
box off
ylim([0, 0.45])
yticks(0:0.1:0.4)



% Name figure
% figure_name = ['LogLaw_Curvilinear_', wind_speed, '_FrictionVelocity.pdf'];
% pause(3)
% exportgraphics(ax, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 600, 'ContentType', 'image'); 
% close all















%%
% %% Scaling friction velocity by (\delta^* / \delta)
% 
% 
% wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
% 
% lw = 2;
% nu = 1.46E-5;
% sz = 20;
% 
% idx = 86;
% 
% wind_speed = 'WT4';
% 
% clc; close all;
% ax = figure('color', 'white');
% tiledlayout(1, 4)
% sgtitle(wind_speed)
% 
% for phase = 1:4
% 
%     h(phase) = nexttile();
%     hold on
% 
%     smooth_y_plus = 0:10:1E4;
%     smooth_u_plus = (1/karman_constant) * log(smooth_y_plus) + c_plus;
%     plot(smooth_y_plus, smooth_u_plus, 'linestyle', '--', 'HandleVisibility', 'off', 'color', 'black')
% 
%     for w = 1:length(waves)
%         wave = waves{w};
%         caze = strcat(wind_speed, '_WV', wave, '_AG0');
% 
%         % Get data
%         y = u_star.(caze)(phase).y;
%         u_profile = u_star.(caze)(phase).u_profile;
%         friction_velocity = u_star.(caze)(phase).raw;
% 
%         % Get boundary layer parameters
%         thickness = integral.(caze)(phase).filtered.thickness(idx);
%         displacement = integral.(caze)(phase).filtered.displacement(idx);
% 
%         % Try scaling friction velocity
%         friction_velocity = friction_velocity * (displacement / thickness);
% 
%         % Law of the wall
%         u_plus = u_profile / friction_velocity;
%         y_plus = (y * friction_velocity) / nu;
% 
%         % Calculate vertical shift from smooth wall
%         ymin = 100;
%         ymax = min(0.9 * max(y_plus), 5000); % up to 30% of max, capped at 5000
%         log_region = (y_plus > ymin & y_plus < ymax);
% 
%         % Difference (measured - smooth)
%         same_grid_smooth_u_plus = (1/karman_constant) * log(y_plus) + c_plus;
%         deltaU = mean(u_plus(log_region) - same_grid_smooth_u_plus(log_region), 'all', 'omitnan');
%         disp(deltaU)
%         % deltaU = 0;
% 
%         % Plot
%         plot(y_plus, u_plus - deltaU, 'linewidth', lw, 'displayname', sprintf("Phase %1.0f", phase), 'color', wave_colors{w})
%         scatter(y_plus, u_plus - deltaU, sz, 'filled', 'HandleVisibility', 'off', 'MarkerFaceColor', wave_colors{w})
% 
% 
%     end
% 
%     % % No-Wave
%     % no_wave_caze = strcat(wind_speed, '_WV0_AGP');
%     % 
%     % y = u_star.(no_wave_caze).y;
%     % u_profile = u_star.(no_wave_caze).u_profile;
%     % friction_velocity = u_star.(no_wave_caze).raw;
%     % 
%     % % Law of the wall
%     % u_plus = u_profile / friction_velocity;
%     % y_plus = (y * friction_velocity) / nu;
%     % 
%     % % Calculate vertical shift from smooth wall
%     % ymin = 600;
%     % ymax = min(0.8 * max(y_plus), 5000); % up to 30% of max, capped at 5000
%     % log_region = (y_plus > ymin & y_plus < ymax);
%     % 
%     % % Difference (measured - smooth)
%     % same_grid_smooth_u_plus = (1/karman_constant) * log(y_plus) + c_plus;
%     % deltaU = mean(u_plus(log_region) - same_grid_smooth_u_plus(log_region), 'all', 'omitnan');
%     % disp(deltaU)
%     % 
%     % % Plot
%     % plot(y_plus, u_plus - deltaU, 'linewidth', lw, 'displayname', 'No Wave', 'color', 'black')
%     % scatter(y_plus, u_plus - deltaU, sz, 'filled', 'HandleVisibility', 'off', 'MarkerFaceColor', 'black')
% 
%     hold off
%     xscale('log')
%     grid on
%     % xlabel('$y^+$', 'interpreter', 'latex', 'fontsize', 16)
%     % 
%     % if w == 1
%     %     ylabel('$u^+$', 'interpreter', 'latex', 'fontsize', 16)
%     % end
% 
% end
% % 
% linkaxes(h, 'xy')
% % xlim([0.5E2, 1E4])
% % ylim([0, 30])
% leg = legend('Orientation', 'Horizontal');
% leg.Layout.Tile = 'north';





















%% Log law profiles across different phases

% wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
% 
% lw = 2;
% nu = 1.46E-5;
% sz = 20;
% 
% wind_speed = 'WT8';
% 
% clc; close all;
% ax = figure('color', 'white');
% tiledlayout(1, length(waves))
% sgtitle(wind_speed)
% 
% for w = 1:length(waves)
%     wave = waves{w};
%     caze = strcat(wind_speed, '_WV', wave, '_AG0');
% 
%     h(w) = nexttile();
%     title(wave)
%     hold on
%     for phase = 1:4
% 
%         y = u_star.(caze)(phase).y;
%         u_profile = u_star.(caze)(phase).u_profile;
%         friction_velocity = u_star.(caze)(phase).raw;
% 
%         % Law of the wall
%         u_plus = u_profile / friction_velocity;
%         y_plus = (y * friction_velocity) / nu;
% 
%         % Plot
%         plot(y_plus, u_plus, 'linewidth', lw, 'displayname', sprintf("Phase %1.0f", phase), 'color', wave_colors{phase})
%         scatter(y_plus, u_plus, sz, 'filled', 'HandleVisibility', 'off', 'MarkerFaceColor', wave_colors{phase})
% 
% 
%     end
% 
%     % No-Wave
%     no_wave_caze = strcat(wind_speed, '_WV0_AGP');
% 
%     y = u_star.(no_wave_caze).y;
%     u_profile = u_star.(no_wave_caze).u_profile;
%     friction_velocity = u_star.(no_wave_caze).raw;
% 
%     % Law of the wall
%     u_plus = u_profile / friction_velocity;
%     y_plus = (y * friction_velocity) / nu;
% 
%     % Plot
%     plot(y_plus, u_plus, 'linewidth', lw, 'displayname', 'No Wave', 'color', 'black')
%     scatter(y_plus, u_plus, sz, 'filled', 'HandleVisibility', 'off', 'MarkerFaceColor', 'black')
% 
%     hold off
%     xscale('log')
%     grid on
%     xlabel('$y^+$', 'interpreter', 'latex', 'fontsize', 16)
% 
%     if w == 1
%         ylabel('$u^+$', 'interpreter', 'latex', 'fontsize', 16)
%     end
% 
% end
% 
% linkaxes(h, 'xy')
% xlim([0.5E2, 1E4])
% ylim([0, 30])
% leg = legend('Orientation', 'Horizontal');
% leg.Layout.Tile = 'north';
% 
% clear ax caze friction_velocity h leg lw phase sz u_plus y_plus u_profile w wave wind_speed y
% 
% 
% %% Log law profiles across different phases: using 'average-of-phases' u*
% 
% lw = 2;
% nu = 1.46E-5;
% sz = 20;
% 
% wind_speed = 'WT8';
% 
% clc; close all;
% ax = figure('color', 'white');
% tiledlayout(1, length(waves))
% sgtitle(wind_speed)
% 
% for w = 1:length(waves)
%     wave = waves{w};
%     caze = strcat(wind_speed, '_WV', wave, '_AG0');
% 
%     h(w) = nexttile();
%     title(wave)
%     hold on
% 
% 
%     %%% Average-of-phases u*
%     avg_friction_velocity = nan(1,4);
%     for i = 1:4
%         avg_friction_velocity(1,i) = u_star.(caze)(i).raw;
%     end
% 
%     for phase = 1:4
% 
%         y = u_star.(caze)(phase).y;
%         u_profile = u_star.(caze)(phase).u_profile;
%         friction_velocity = mean(avg_friction_velocity, 'all', 'omitnan');
% 
%         % Law of the wall
%         u_plus = u_profile / friction_velocity;
%         y_plus = (y * friction_velocity) / nu;
% 
%         % Plot
%         plot(y_plus, u_plus, 'linewidth', lw, 'displayname', sprintf("Phase %1.0f", phase), 'color', wave_colors{phase})
%         scatter(y_plus, u_plus, sz, 'filled', 'HandleVisibility', 'off', 'MarkerFaceColor', wave_colors{phase})
% 
% 
%     end
% 
%     % No-Wave
%     no_wave_caze = strcat(wind_speed, '_WV0_AGP');
% 
%     y = u_star.(no_wave_caze).y;
%     u_profile = u_star.(no_wave_caze).u_profile;
%     friction_velocity = u_star.(no_wave_caze).raw;
% 
%     % Law of the wall
%     u_plus = u_profile / friction_velocity;
%     y_plus = (y * friction_velocity) / nu;
% 
%     % Plot
%     plot(y_plus, u_plus, 'linewidth', lw, 'displayname', 'No Wave', 'color', 'black')
%     scatter(y_plus, u_plus, sz, 'filled', 'HandleVisibility', 'off', 'MarkerFaceColor', 'black')
% 
%     hold off
%     xscale('log')
%     grid on
%     xlabel('$y^+$', 'interpreter', 'latex', 'fontsize', 16)
% 
%     if w == 1
%         ylabel('$u^+$', 'interpreter', 'latex', 'fontsize', 16)
%     end
% 
% end
% 
% linkaxes(h, 'xy')
% xlim([0.5E2, 1E4])
% ylim([0, 30])
% leg = legend('Orientation', 'Horizontal');
% leg.Layout.Tile = 'north';
% 
% clear ax caze friction_velocity h leg lw phase sz u_plus y_plus u_profile w wave wind_speed y








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


