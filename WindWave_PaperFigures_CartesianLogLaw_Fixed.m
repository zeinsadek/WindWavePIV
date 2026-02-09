%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions/');
experiment_log  = readcell('Offshore_Inflow_Log.xlsx');
wave_parameters = readcell('Offshore_Waves.xlsx');

s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";

top_bound_value   = 205;       % relative to Y centered at still water
left_bound_value  = -121;      % relative to X centered at DaVis default
right_bound_value = 115;       % relative to X centered at DaVis default
range = abs(left_bound_value) + abs(right_bound_value);

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

%% Import Data

% Cases with waves
WTs = {'4', '6', '8'};
WVs = {'A', 'B', 'C', 'D'};

BL_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/boundary_layer';
LL_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/log_law_fixed/';
MN_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means/';
figure_path = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures';
% Load cases with waves
for s = 1:length(WTs)
    speed = WTs{s};
    for i = 1:length(WVs)
        % Case
        caze = strcat('WT', num2str(speed), '_WV', WVs{i}, '_AG0');
        fprintf('Loading Case: %s...\n', caze)
    
        % Paths
        BL_path = fullfile(BL_folder, strcat(caze, '_BL.mat'));
        LL_path = fullfile(LL_folder, strcat(caze, '_LL_Fixed.mat'));
        MN_path = fullfile(MN_folder, strcat(caze, '_MEANS.mat'));
        
        % Save LL to Structure
        LL_temp = load(LL_path);
        LL_temp = LL_temp.output;
        LL_Waves.(caze) = LL_temp;

        % Save BL to Structure
        BL_temp = load(BL_path);
        BL_temp = BL_temp.output;
        BL_Waves.(caze) = BL_temp;
       
        % Save MN to Structure
        means_temp = load(MN_path);
        means_temp = means_temp.output;
        MN_Waves.(caze) = means_temp;
    end
    fprintf('\n')
end

% Load cases without waves
for s = 1:length(WTs)
    % Case
    no_wave_case = strcat('WT', num2str(WTs{s}), '_WV0_AGP');
    fprintf('Loading Case: %s...\n', no_wave_case)
    
    % Paths
    BL_path = fullfile(BL_folder, strcat(no_wave_case, '_BL.mat'));
    LL_path = fullfile(LL_folder, strcat(no_wave_case, '_LL_Fixed.mat'));
    MN_path = fullfile(MN_folder, strcat(no_wave_case, '_MEANS.mat'));
    
    % Save LL to Structure
    LL_temp = load(LL_path);
    LL_temp = LL_temp.output;
    LL_No_Waves.(no_wave_case) = LL_temp;

    % Save BL to Structure
    BL_temp = load(BL_path);
    BL_temp = BL_temp.output;
    BL_No_Waves.(no_wave_case) = BL_temp;
    
    % Save MN to Structure
    means_temp = load(MN_path);
    means_temp = means_temp.output;
    MN_No_Waves.(no_wave_case) = means_temp;
end

% Get Case Names
wave_cases    = fields(MN_Waves);
no_wave_cases = fields(MN_No_Waves);


%% All Log Laws

% WT_colors = turbo(length(WTs));
% WV_colors = cool(length(WVs));
% markers    = {'o', 'square', 'diamond', '^'};
% alpha = 1;
% size  = 5;
% lw    = 1;
% 
% clc;
% figure('Name', 'Ensemble Average Log Law')
% hold on
% for i = 1:length(wave_cases)
% 
%     % marker = markers{i};
%     % color = WV_colors(i,:);   
%     caze  = wave_cases{i};
%     disp(caze)
% 
%     % Color based on wave
%     if contains(caze, 'WVA') == 1
%         color = WV_colors(1,:);
%     elseif contains(caze, 'WVB') == 1
%         color = WV_colors(2,:);
%     elseif contains(caze, 'WVC') == 1
%         color = WV_colors(3,:);
%     elseif contains(caze, 'WVD') == 1
%         color = WV_colors(4,:);
%     end
% 
%      % Linestyle Based on wind speed
%     if contains(caze, 'WT4') == 1
%         ls = '-';
%     elseif contains(caze, 'WT6') == 1
%         ls = '--';
%     elseif contains(caze, 'WT8') == 1
%         ls = ':';
%     end
% 
%     % Ensemble
%     ensemble_u_plus = LL_Waves.(caze).ensemble.u_plus;
%     ensemble_y_plus = LL_Waves.(caze).ensemble.y_plus;
% 
%     scatter(ensemble_y_plus, ensemble_u_plus, size, 'o', 'filled', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha)
%     plot(ensemble_y_plus, ensemble_u_plus, 'color', color, 'linewidth', lw, 'LineStyle', ls)
% 
%     % Phase
%     for j = 1:4
%         % Phase Averge Profile with Ensemble u*
%         ens_phase_u_plus = LL_Waves.(caze).phase(j).ens_u_plus;
%         ens_phase_y_plus = LL_Waves.(caze).phase(j).ens_y_plus;
%         scatter(ens_phase_y_plus, ens_phase_u_plus, size, 'o', 'filled', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha)
%         plot(ens_phase_y_plus, ens_phase_u_plus, 'color', color, 'linewidth', lw, 'LineStyle', ls)
% 
%         % Phase Averge Profile with Phase u*
%         phase_u_plus = LL_Waves.(caze).phase(j).u_plus;
%         phase_y_plus = LL_Waves.(caze).phase(j).y_plus;
%         scatter(phase_y_plus, phase_u_plus, size, 'o', 'filled', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha)
%         plot(phase_y_plus, phase_u_plus, 'color', color, 'linewidth', lw, 'LineStyle', ls)
%     end
% end
% 
% for i = 1:length(no_wave_cases)
%     % marker = markers{i};
%     % color  = WV_colors(i,:);   
%     caze   = no_wave_cases{i};
%     disp(caze)
% 
%      % Linestyle Based on wind speed
%     if contains(caze, 'WT4') == 1
%         ls = '-';
%     elseif contains(caze, 'WT6') == 1
%         ls = '--';
%     elseif contains(caze, 'WT8') == 1
%         ls = ':';
%     end
% 
%     % Ensemble
%     ensemble_u_plus = LL_No_Waves.(caze).ensemble.u_plus;
%     ensemble_y_plus = LL_No_Waves.(caze).ensemble.y_plus;
%     scatter(ensemble_y_plus, ensemble_u_plus, size, 'o', 'filled', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha)
%     plot(ensemble_y_plus, ensemble_u_plus, 'color', 'black', 'linewidth', lw, 'LineStyle', ls)
% end
% hold off
% 
% title('All Log Laws')
% ylabel('$u^+$', 'Interpreter', 'Latex')
% xlabel('$y^+$', 'Interpreter', 'Latex')
% set(gca,'xscale','log')
% xlim([5, 1E4])
% ylim([5,30])
% grid on



%% Cartesian phase average log laws binned by phase

tickFontSize = 14;
labelFontSize = 16;
legendFontSize = 12;
titleFontSize = 18;

wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};

lw = 1;
sz = 10;

clc; close all
ax = figure('color', 'white', 'units', 'centimeters', 'position', [2,2,30,8]);
t = tiledlayout(1, 4);

wind = '6';
% phase_names = {'Peak', 'Ascending', 'Trough', 'Descending'};
phase_names = {'$\varphi = 0$', '$\varphi = \lambda / 4$', '$\varphi = \lambda / 2$', '$\varphi = 3 \lambda / 4$'};

for p = 1:4
    h(p) = nexttile;

    hold on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
    title(phase_names{p}, 'interpreter', 'latex', 'fontsize', 16)

    % Plot smooth wall
    karman_constant = 0.41;
    c_plus = 5;
    
    smooth_y_plus = 0:10:1E4;
    smooth_u_plus = (1/karman_constant) * log(smooth_y_plus) + c_plus;
    plot(smooth_y_plus, smooth_u_plus, 'linestyle', '--', 'HandleVisibility', 'off', 'color', 'black')

    for WV = 1:length(WVs)

        wave = WVs{WV};
        caze = strcat('WT', wind, '_WV', wave, '_AG0');
        disp(caze)

        u_plus = LL_Waves.(caze).phase(p).u_plus;
        y_plus = LL_Waves.(caze).phase(p).y_plus;

        % Calculate vertical shift from smooth wall
        ymin = 1000;
        ymax = min(1 * max(y_plus), 5000); % up to 30% of max, capped at 5000
        log_region = (y_plus > ymin & y_plus < ymax);
    
        % Difference (measured - smooth)
        same_grid_smooth_u_plus = (1/karman_constant) * log(y_plus) + c_plus;
        deltaU = mean(u_plus(log_region) - same_grid_smooth_u_plus(log_region), 'all', 'omitnan');
        disp(deltaU)
        deltaU = 0;

        label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses.(wave));
        plot(y_plus, u_plus - deltaU, 'linewidth', lw, 'color', wave_colors{WV}, 'HandleVisibility', 'off')
        scatter(y_plus, u_plus - deltaU, sz, 'filled', 'HandleVisibility', 'on', 'MarkerFaceColor', wave_colors{WV}, 'DisplayName', label)

    end

    no_wave_caze = strcat('WT', wind, '_WV0_AGP');
    ens_u_plus = LL_No_Waves.(no_wave_caze).ensemble.u_plus;
    ens_y_plus = LL_No_Waves.(no_wave_caze).ensemble.y_plus;

    % Calculate vertical shift from smooth wall
    ymin = 1000;
    ymax = min(1 * max(ens_y_plus), 5000); % up to 30% of max, capped at 5000
    log_region = (ens_y_plus > ymin & ens_y_plus < ymax);

    % Difference (measured - smooth)
    same_grid_smooth_u_plus = (1/karman_constant) * log(ens_y_plus) + c_plus;
    deltaU = mean(ens_u_plus(log_region) - same_grid_smooth_u_plus(log_region), 'all', 'omitnan');
    disp(deltaU)
    deltaU = 0

    % label = sprintf('$\\lambda_{%s}$', wavelength_names.(wave));
    P = plot(ens_y_plus, ens_u_plus - deltaU, 'color', 'black', 'linewidth', lw, 'handlevisibility', 'off');
    % P.Color(4) = 0.5;
    scatter(ens_y_plus, ens_u_plus - deltaU, sz, 'filled', 'HandleVisibility', 'on', 'MarkerFaceColor', 'black', 'displayname', 'No Waves');

    hold off
    grid on
    set(gca,'xscale','log')
    ylim([10, 26])
    xlim([1E1, 1E4])

    if p == 1
        ylabel('$\langle u \rangle^+ + \Delta U^+$', 'interpreter', 'latex', 'fontsize', labelFontSize)
    end
end

xlabel(t, '$y^+$', 'interpreter', 'latex', 'fontsize', labelFontSize)
linkaxes(h, 'xy');
leg = legend('Orientation', 'Vertical', 'box', 'off', 'interpreter', 'latex', 'fontsize', 18);
leg.Layout.Tile = 'East';


% Save figure
% figure_name = ['LogLaw_Cartesian_WT', wind, '_Profiles_Fixed.pdf'];
% pause(3)
% exportgraphics(ax, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 600, 'ContentType', 'image'); 
% close all


%% Bar chart showing different friction velocities

friction_velocity_bar_chart = nan(4,length(WVs));

for phase = 1:4

    hold on
    title(phase_names{phase})
    
    for w = 1:length(WVs)
        wave = WVs{w};
        caze = strcat('WT', wind, '_WV', wave, '_AG0');

        % Get data and save to array
        friction_velocity_bar_chart(phase, w) = LL_Waves.(caze).phase(phase).u_star;

    end
end

% No-Wave
no_wave_caze = strcat('WT', wind, '_WV0_AGP');
no_wave_friction_velocity = LL_No_Waves.(no_wave_caze).ensemble.u_star.u_star;

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


% Save figure
% pause(3)
% figure_name = ['LogLaw_Cartesian_WT', wind, '_FrictionVelocity.pdf'];
% exportgraphics(ax, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 600, 'ContentType', 'image'); 
% close all




%% Cf vs x / \lambda


% WT = 6;
% gauss_window = 4;
% phase = 2;
% 
% tickFontSize = 8;
% labelFontSize = 16;
% legendFontSize = 12;
% 
% wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
% linestyles = {'-.', '--', '-'};
% lw = 2;
% buff = 5;
% scale = 1.5;
% 
% 
% close all; 
% ax = figure('color', 'white', 'units', 'centimeters', 'position', [20,20,29,10]);
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% 
% hold on
% 
% for wt = 1:length(WTs)
% 
%     WT = WTs{wt};
%     for w = 1:length(WVs)
%         caze = strcat('WT', WT, '_WV', WVs{w}, '_AG0');
% 
%         if contains(caze, 'WVA') == 1
%             color = WV_colors(1,:);
%             wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
%             wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
%             % label = 'WVA';
%         elseif contains(caze, 'WVB') == 1
%             color = WV_colors(2,:);
%             wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
%             wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
%             % label = 'WVB';
%         elseif contains(caze, 'WVC') == 1
%             color = WV_colors(3,:);
%             wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
%             wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
%             % label = 'WVC';
%         elseif contains(caze, 'WVD') == 1
%             color = WV_colors(4,:);
%             wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
%             wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
%             % label = 'WVD';
%         end
% 
%         % Data
%         u_inf = BL_Waves.(caze).u_inf;
%         uv    = MN_Waves.(caze).phase(phase).uv;
%         x     = BL_Waves.(caze).x;
% 
%         % Compute u*
%         u_star_profile = sqrt(max(-1 * uv, [], 1, 'omitnan'));
%         u_star_profile = smoothdata(u_star_profile, 'gaussian', gauss_window);
% 
%         % Compute Cf
%         cf = 2 * (u_star_profile / u_inf).^2;
% 
%         % Waves
%         centered_normalized_x = (x - (range/2)) / wave_length;
%         reference_wave = MN_Waves.(caze).phase(phase).reference_wave / wave_amplitude;
% 
%         % Scale Wave
%         wave_scale = 0.007;
%         wave_offset = 0.012;
%         reference_wave = reference_wave * wave_scale;
%         reference_wave = reference_wave + wave_offset;
% 
%         % Plot
%         label = sprintf('$\\lambda_{%s}, u_{\\infty} = %1.2f$', wavelength_names.(WVs{w}), u_inf);
%         plot(centered_normalized_x, cf, 'color', wave_colors{w}, 'linewidth', lw, 'displayname', label, 'linestyle', linestyles{wt});
% 
%     end
% end
% 
% % Reference wave
% a = plot(centered_normalized_x, reference_wave, 'color', 'black', 'linewidth', lw, 'HandleVisibility', 'off', 'linestyle', '--');
% a.Color(4) = 0.5;
% uistack(a, 'bottom')
% hold off
% 
% ylim([0, 0.025])
% yticks(0:0.005:0.025)
% xticks(-1:0.25:1)
% 
% legend('interpreter', 'latex', 'fontsize', 10, 'location',  'northeastoutside', 'box', 'off')
% xlabel('$x / \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylabel('$C_f$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% 
% % % Name figure
% % figure_name = strcat('SkinFriction_Cartesian_Phase', num2str(phase), '.png');
% % 
% % % Save figure
% % exportgraphics(ax, fullfile(figure_folder, figure_name), 'Resolution', 300); 
% % close all









%% Log Law Different Profiles [Ensemble Phases]

% phase_colors = cool(4);
% lw = 1.5;
% 
% 
% clc;
% ax = figure('name', 'LL Profiles');
% t = tiledlayout(length(WTs), length(WVs));
% 
% for WT = 1:length(WTs)
%     for WV = 1:length(WVs)
% 
%         wind = WTs{WT};
%         wave = WVs{WV};
%         caze = strcat('WT', wind, '_WV', wave, '_AG0');
%         disp(caze)
% 
%         ens_u_plus = LL_Waves.(caze).ensemble.u_plus;
%         ens_y_plus = LL_Waves.(caze).ensemble.y_plus;
% 
%         ax(WV) = nexttile();
%         title(caze, 'interpreter', 'none')
%         hold on
%         for j = 1:4
%             color = phase_colors(j,:);
%             ens_phase_u_plus = LL_Waves.(caze).phase(j).ens_u_plus;
%             ens_phase_y_plus = LL_Waves.(caze).phase(j).ens_y_plus;
%             plot(ens_phase_y_plus, ens_phase_u_plus, 'linewidth', lw, 'color', color)
%         end
%         plot(ens_y_plus, ens_u_plus, 'color', 'black', 'linewidth', lw)
%         hold off
%         set(gca,'xscale','log')
%         xlim([5, 1E4])
%         % ylim([5,30])
%         grid on
%     end
%     linkaxes(ax,'y')
% end
% 
% %% Get AVG Freestream from all cases
% 
% clc;
% for s = 1:length(WTs)
% 
%     c = 1;
%     all_u_infs = zeros(1,5);
%     WT = WTs{s};
%     % Waves
%     for w = 1:length(WVs)
%         WV = WVs{w};
%         caze = strcat('WT', WT, '_WV', WV, '_AG0');
%         u_inf = BL_Waves.(caze).u_inf;
%         all_u_infs(c) = u_inf;
%         c = c + 1;
%     end
% 
%     % No Waves
%     caze = strcat('WT', WT, '_WV0_AGP');
%     u_inf = BL_No_Waves.(caze).u_inf;
%     all_u_infs(5) = u_inf;
%     avg_u_inf = mean(all_u_infs, 'all');
%     fprintf("WT %s = %2.4f m/s\n", WT, avg_u_inf);
%     u_infs.(strcat("WT", WT)) = avg_u_inf;
% end
% 
% 
% %% Log Law Ensembles: Waves and No Waves
% 
% wave_colors = cool(4);
% lw = 1.5;
% 
% clc;
% ax = figure('name', 'LL Profiles', 'Position', [200,200,1000,300]);
% t = tiledlayout(1, length(WTs), 'TileSpacing', 'compact', 'Padding', 'tight');
% 
% for WT = 1:length(WTs)
%     h(WT) = nexttile();
%     hold on
% 
%     % No Waves
%     caze = strcat('WT', wind, '_WV0_AGP');
%     disp(caze)
%     ens_u_plus = LL_No_Waves.(caze).ensemble.u_plus;
%     ens_y_plus = LL_No_Waves.(caze).ensemble.y_plus;
%     d(1) = plot(ens_y_plus, ens_u_plus, 'color', 'black', 'linewidth', lw);
% 
%     % With Waves
%     for WV = 1:length(WVs)
% 
%         wind = WTs{WT};
%         wave = WVs{WV};
%         caze = strcat('WT', wind, '_WV', wave, '_AG0');
%         disp(caze)
% 
%         if wave == "A"
%             color = wave_colors(1,:);
%         elseif wave == "B"
%             color = wave_colors(2,:);
%         elseif wave == "C"
%             color = wave_colors(3,:);
%         elseif wave == "D"
%             color = wave_colors(4,:);
%         end
% 
%         ens_u_plus = LL_Waves.(caze).ensemble.u_plus;
%         ens_y_plus = LL_Waves.(caze).ensemble.y_plus;
% 
%         d(WV+1) = plot(ens_y_plus, ens_u_plus, 'color', color, 'linewidth', lw);
%     end
% 
%     hold off
%     grid on
%     set(gca,'xscale','log')
%     xlim([7, 0.5E4])
%     tit = sprintf('$u_{\\infty}$ =  %2.2f m/s', u_infs.(strcat("WT", WTs{WT})));
%     title(tit, 'interpreter', 'latex', 'fontsize', 14)
%     ylabel('$u^+$', 'interpreter', 'latex', 'fontsize', 18)
%     xlabel('$y^+$', 'interpreter', 'latex', 'fontsize', 18)
% end
% linkaxes(h, 'y');
% leg = legend(d, 'WV0', 'WVA', 'WVB', 'WVC', 'WVD', 'orientation', 'horizontal', 'box', 'off');
% leg.Layout.Tile = 'north';
% 
% 
% % exportgraphics(ax, fullfile(figure_path, strcat('Ensemble_Mean_LL_Profiles.png')), 'resolution', 300);
% 
% 
% 
% 
% %% Log Law Phase: Waves and No Waves
% 
% wave_colors = cool(4);
% lw = 1.5;
% phase_styles = {"-","--","-.",":"};
% 
% clc;
% clear d;
% ax = figure('name', 'LL Profiles', 'Position', [200,200,1000,300]);
% t = tiledlayout(1, length(WTs), 'TileSpacing', 'compact', 'Padding', 'tight');
% 
% for WT = 1:length(WTs)
%     h(WT) = nexttile();
%     hold on
% 
%     % With Waves
%     % for WV = 1:length(WVs)
%     for WV = 3
% 
%         wind = WTs{WT};
%         wave = WVs{WV};
%         caze = strcat('WT', wind, '_WV', wave, '_AG0');
%         disp(caze)
% 
%         color = wave_colors(WV,:);
% 
%         for p = 1:4
% 
%             style = phase_styles{p};
%             ens_u_plus = LL_Waves.(caze).phase(p).ens_u_plus;
%             ens_y_plus = LL_Waves.(caze).phase(p).ens_y_plus;
% 
%             if WT == 1
%                 d(WV) = plot(ens_y_plus, ens_u_plus, 'color', color, 'linewidth', lw, 'linestyle', style);
%             else
%                 plot(ens_y_plus, ens_u_plus, 'color', color, 'linewidth', lw, 'linestyle', style);
%             end
%         end
%         ens_u_plus = LL_Waves.(caze).ensemble.u_plus;
%         ens_y_plus = LL_Waves.(caze).ensemble.y_plus;
%         plot(ens_y_plus, ens_u_plus, 'color', 'black', 'linewidth', lw);
% 
%     end
% 
%     hold off
%     grid on
%     set(gca,'xscale','log')
%     tit = sprintf('$u_{\\infty}$ =  %2.2f m/s', u_infs.(strcat("WT", WTs{WT})));
%     title(tit, 'interpreter', 'latex', 'fontsize', 14)
%     ylabel('$u^+$', 'interpreter', 'latex', 'fontsize', 18)
%     xlabel('$y^+$', 'interpreter', 'latex', 'fontsize', 18)
% end
% linkaxes(h, 'xy');
% leg = legend(d, 'WVA', 'WVB', 'WVC', 'WVD', 'orientation', 'horizontal', 'box', 'off');
% leg.Layout.Tile = 'north';
% 
% 
% % exportgraphics(ax, fullfile(figure_path, strcat('Phase_Mean_LL_Profiles.png')), 'resolution', 300);



%% Cf Phase Subplots

WT = 4;
gauss_window = 5;
WV_colors = cool(length(WVs));

ax = figure('position', [200,200,800,400]);
t = tiledlayout(2,2,'padding', 'tight', 'TileSpacing','tight');
% sgtitle(strcat('WT', num2str(WT)))
for i = 1:4
    h(i) = nexttile();
    for w = 1:length(WVs)
        caze = strcat('WT', num2str(WT), '_WV', WVs{w}, '_AG0');

        if contains(caze, 'WVA') == 1
            color = WV_colors(1,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
            label = 'WVA';
        elseif contains(caze, 'WVB') == 1
            color = WV_colors(2,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
            label = 'WVB';
        elseif contains(caze, 'WVC') == 1
            color = WV_colors(3,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
            label = 'WVC';
        elseif contains(caze, 'WVD') == 1
            color = WV_colors(4,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
            label = 'WVD';
        end
    
        % Data
        u_inf = BL_Waves.(caze).u_inf;
        uv    = MN_Waves.(caze).phase(i).uv;
        x     = BL_Waves.(caze).x;
    
        % Compute u*
        u_star_profile = sqrt(max(-1 * uv, [], 1, 'omitnan'));
        u_star_profile = smoothdata(u_star_profile, 'gaussian', gauss_window);

        % Compute Cf
        cf = 2 * (u_star_profile / u_inf).^2;

        % Waves
        centered_normalized_x = (x - (range/2)) / wave_length;
        reference_wave = MN_Waves.(caze).phase(i).reference_wave / wave_amplitude;
        
        % Scale Wave
        wave_scale = 0.005;
        wave_offset = 0.008;
        reference_wave = reference_wave * wave_scale;
        reference_wave = reference_wave + wave_offset;
         
        % Plot
        hold on
        a = plot(centered_normalized_x, reference_wave, 'color', 'black', 'linewidth', 2, 'HandleVisibility', 'off');
        b(w) = plot(centered_normalized_x, cf, 'color', color, 'linewidth', 2);
        a.Color(4) = 1;
        uistack(a, 'bottom')
        % Plot Vertical Line At Max Cf
        % [max_cf, max_cf_idx] = max(cf);
        % xline(centered_normalized_x(max_cf_idx), 'linestyle', '--', 'color', color)
        hold off
        tit = strcat('Phase', {' '}, num2str(i));
        title(tit{1})
        xlabel('$x / \lambda$', 'Interpreter', 'latex')
        ylabel('$C_f$', 'Interpreter', 'latex')
        
    end
end

leg = legend(b, 'WVA', 'WVB', 'WVC', 'WVD', 'orientation', 'horizontal');
leg.Layout.Tile = 'north';
linkaxes(h, 'xy')

% exportgraphics(ax, fullfile(figure_path, strcat("WT_", num2str(WT), "_Cf_Phases.png")), resolution=300)


%% Stitched Cf Across All Waves

WT = 4;
gauss_window = 1;
WV_colors = cool(length(WVs));

ax = figure();

hold on
for i = 1:4
    for w = 1:length(WVs)
        caze = strcat('WT', num2str(WT), '_WV', WVs{w}, '_AG0');

        if contains(caze, 'WVA') == 1
            color = WV_colors(1,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
            label = 'WVA';
        elseif contains(caze, 'WVB') == 1
            color = WV_colors(2,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
            label = 'WVB';
        elseif contains(caze, 'WVC') == 1
            color = WV_colors(3,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
            label = 'WVC';
        elseif contains(caze, 'WVD') == 1
            color = WV_colors(4,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
            label = 'WVD';
        end
    
        % Data
        u_inf = BL_Waves.(caze).u_inf;
        uv    = MN_Waves.(caze).phase(i).uv;
        x     = BL_Waves.(caze).x;
    
        % Compute u*
        u_star_profile = sqrt(max(-1 * uv, [], 1, 'omitnan'));
        u_star_profile = smoothdata(u_star_profile, 'gaussian', gauss_window);

        % Compute Cf
        cf = 2 * (u_star_profile / u_inf).^2;

        % Waves
        centered_normalized_x = ((x - (range/2)) / wave_length) - ((i - 1)/4);
        reference_wave = MN_Waves.(caze).phase(i).reference_wave / wave_amplitude;
        
        % Scale Wave
        wave_scale = 0.005;
        wave_offset = 0.008;
        reference_wave = reference_wave * wave_scale;
        reference_wave = reference_wave + wave_offset;
        % 
        % Plot
        a = plot(centered_normalized_x, reference_wave, 'color', 'black', 'linewidth', 2, 'HandleVisibility', 'off');
        b(w) = plot(centered_normalized_x, cf, 'color', color, 'linewidth', 2);
        a.Color(4) = 1;
        uistack(a, 'bottom')

        % Plot Vertical Line At Max Cf
        [max_cf, max_cf_idx] = max(cf);
        xline(centered_normalized_x(max_cf_idx), 'linestyle', '--', 'color', color)
        
    end
end
hold off
legend(b, 'WVA', 'WVB', 'WVC', 'WVD');
xlabel('$x / \lambda$', 'Interpreter', 'latex')
ylabel('$C_f$', 'Interpreter', 'latex')
xlim([-1.85 1.15])
ylim([0.0020 0.0220])
title(strcat('WT', num2str(WT), ' All Waves: Stitched Phases'));

%% Re Phase Subplots

WT = 4;
gauss_window = 1;
WV_colors = cool(length(WVs));
kinematic_viscosity = 1.48E-5;

ax = figure();
t = tiledlayout(2,2);
sgtitle(strcat('WT', num2str(WT)))
for i = 1:4
    disp(i)
    h(i) = nexttile();
    for w = 1:length(WVs)
        caze = strcat('WT', num2str(WT), '_WV', WVs{w}, '_AG0');

        if contains(caze, 'WVA') == 1
            color = WV_colors(1,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
            label = 'WVA';
        elseif contains(caze, 'WVB') == 1
            color = WV_colors(2,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
            label = 'WVB';
        elseif contains(caze, 'WVC') == 1
            color = WV_colors(3,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
            label = 'WVC';
        elseif contains(caze, 'WVD') == 1
            color = WV_colors(4,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
            label = 'WVD';
        end
    
        % Data
        u_inf = BL_Waves.(caze).u_inf;
        uv    = MN_Waves.(caze).phase(i).uv;
        theta = BL_Waves.(caze).phase(i).momentum * 1E-3;
        x     = BL_Waves.(caze).x;
    
        % Compute Re
        Re = (u_inf * theta) / kinematic_viscosity;
        Re = smoothdata(Re, 'gaussian', gauss_window);

        % Waves
        centered_normalized_x = (x - (range/2)) / wave_length;
        reference_wave = MN_Waves.(caze).phase(i).reference_wave / wave_amplitude;
        
        % Scale Wave
        wave_scale = 800;
        wave_offset = 1900;
        reference_wave = reference_wave * wave_scale;
        reference_wave = reference_wave + wave_offset;
         
        % Plot
        hold on
        a = plot(centered_normalized_x, reference_wave, 'color', 'black', 'linewidth', 2, 'HandleVisibility', 'off');
        b(w) = plot(centered_normalized_x, Re, 'color', color, 'linewidth', 2);
        a.Color(4) = 1;
        uistack(a, 'bottom')
        % Plot Vertical Line At Max Re
        [max_Re, max_Re_idx] = max(Re);
        xline(centered_normalized_x(max_Re_idx), 'linestyle', '--', 'color', color)
        hold off
        tit = strcat('Phase', {' '}, num2str(i));
        title(tit{1})
        xlabel('$x / \lambda$', 'Interpreter', 'latex')
        ylabel('$Re_{\theta}$', 'Interpreter', 'latex')
        
    end
end

leg = legend(b, 'WVA', 'WVB', 'WVC', 'WVD', 'orientation', 'horizontal');
leg.Layout.Tile = 'north';
linkaxes(h, 'xy')


%% Stitched Re_Theta Across All Waves

WT = 4;
gauss_window = 1;
WV_colors = cool(length(WVs));
kinematic_viscosity = 1.48E-5;

ax = figure();

hold on
for i = 1:4
    for w = 1:length(WVs)
        caze = strcat('WT', num2str(WT), '_WV', WVs{w}, '_AG0');

        if contains(caze, 'WVA') == 1
            color = WV_colors(1,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
            label = 'WVA';
        elseif contains(caze, 'WVB') == 1
            color = WV_colors(2,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
            label = 'WVB';
        elseif contains(caze, 'WVC') == 1
            color = WV_colors(3,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
            label = 'WVC';
        elseif contains(caze, 'WVD') == 1
            color = WV_colors(4,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
            label = 'WVD';
        end
    
        % Data
        u_inf = BL_Waves.(caze).u_inf;
        uv    = MN_Waves.(caze).phase(i).uv;
        theta = BL_Waves.(caze).phase(i).momentum * 1E-3;
        x     = BL_Waves.(caze).x;
    
        % Compute u*
        u_star_profile = sqrt(max(-1 * uv, [], 1, 'omitnan'));
        u_star_profile = smoothdata(u_star_profile, 'gaussian', gauss_window);

        % Compute Re
        Re = (u_inf * theta) / kinematic_viscosity;
        Re = smoothdata(Re, 'gaussian', gauss_window);

        % Waves
        centered_normalized_x = ((x - (range/2)) / wave_length) - ((i - 1)/4);
        reference_wave = MN_Waves.(caze).phase(i).reference_wave / wave_amplitude;
        
        % Scale Wave
        wave_scale = 800;
        wave_offset = 1900;
        reference_wave = reference_wave * wave_scale;
        reference_wave = reference_wave + wave_offset;
       
        % Plot
        a = plot(centered_normalized_x, reference_wave, 'color', 'black', 'linewidth', 2, 'HandleVisibility', 'off');
        b(w) = plot(centered_normalized_x, Re, 'color', color, 'linewidth', 2);
        a.Color(4) = 1;
        uistack(a, 'bottom')
        % Plot Vertical Line At Max Cf
        [max_Re, max_Re_idx] = max(Re);
        xline(centered_normalized_x(max_Re_idx), 'linestyle', '--', 'color', color)
        
    end
end
hold off
legend(b, 'WVA', 'WVB', 'WVC', 'WVD');
xlabel('$x / \lambda$', 'Interpreter', 'latex')
ylabel('$Re_{\theta}$', 'Interpreter', 'latex')
xlim([-1.85 1.15])
% ylim([0.0020 0.0220])
title(strcat('WT', num2str(WT), ' All Waves: Stitched Phases'));



%% Cf vs Re_theta Phases 2 + 4

WT = 4;
gauss_window = 20;
WV_colors = cool(length(WVs));
size = 20;
kinematic_viscosity = 1.48E-5;

ax = figure();
tit = strcat('WT', num2str(WT), ': Phases 2 and 4');
sgtitle(strcat('WT', num2str(WT), ': Phases 2 and 4'))

h = 0;
hold on
for phase = 2:2:4
    disp(phase)
    for w = 1:4
        caze = strcat('WT', num2str(WT), '_WV', WVs{w}, '_AG0');
    
        if contains(caze, 'WVA') == 1
            color = WV_colors(1,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
            label = 'WVA';
        elseif contains(caze, 'WVB') == 1
            color = WV_colors(2,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
            label = 'WVB';
        elseif contains(caze, 'WVC') == 1
            color = WV_colors(3,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
            label = 'WVC';
        elseif contains(caze, 'WVD') == 1
            color = WV_colors(4,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
            label = 'WVD';
        end
    
        if phase == 2
            marker = '^';
        end

        if phase == 4
            marker = 'o';
        end

        % Data
        u_inf = BL_Waves.(caze).u_inf;
        uv    = MN_Waves.(caze).phase(phase).uv;
        theta = BL_Waves.(caze).phase(phase).momentum * 1E-3;
        x     = BL_Waves.(caze).x;
    
        % Compute u*
        u_star_profile = sqrt(max(-1 * uv, [], 1, 'omitnan'));
        % u_star_profile = smoothdata(u_star_profile, 'gaussian', gauss_window);
    
        % Compute Cf
        cf = 2 * (u_star_profile / u_inf).^2;
    
        % Compute Re
        Re_theta = (u_inf * theta) / kinematic_viscosity;
    
        % Smooth
        cf = smoothdata(cf, 'gaussian', gauss_window);
        Re_theta = smoothdata(Re_theta, 'gaussian', gauss_window);
    
        % Crop
        bound = 0.25;
        centered_normalized_x = ((x - (range/2)) / wave_length);
        cf(centered_normalized_x < -bound) = nan;
        cf(centered_normalized_x > bound) = nan;
        Re_theta(centered_normalized_x < -bound) = nan;
        Re_theta(centered_normalized_x > bound) = nan;
    
        % Plot
        % hold on
        h(w) = plot(Re_theta, cf, 'color', color, 'linewidth', 1);
        scatter(Re_theta, cf, size, marker, 'filled', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none')
        p = plot(Re_theta, cf, 'color', color, 'linewidth', 10);
        % hold off
        p.Color(4) = 0.2;
        
    end
end
hold off
legend(h, 'WVA', 'WVB', 'WVC', 'WVD', 'location', 'southeast')
xlabel('$Re_{\theta}$', 'Interpreter', 'latex')
ylabel('$C_f$', 'Interpreter', 'latex')


%% Cf vs Re_theta Phases 2 + 4

WT = 4;
gauss_window = 1;
WV_colors = cool(length(WVs));
size = 20;
kinematic_viscosity = 1.48E-5;

ax = figure();
t = tiledlayout(1,4);
sgtitle(strcat('WT', num2str(WT), ': Phases 2 and 4'))

h = 0;
for w = 1:4
    h(w) = nexttile();
    for phase = 2:2:4
        disp(phase)
        caze = strcat('WT', num2str(WT), '_WV', WVs{w}, '_AG0');
    
        if contains(caze, 'WVA') == 1
            color = WV_colors(1,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
            label = 'WVA';
        elseif contains(caze, 'WVB') == 1
            color = WV_colors(2,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
            label = 'WVB';
        elseif contains(caze, 'WVC') == 1
            color = WV_colors(3,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
            label = 'WVC';
        elseif contains(caze, 'WVD') == 1
            color = WV_colors(4,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
            label = 'WVD';
        end
    
        if phase == 2
            marker = '^';
        end

        if phase == 4
            marker = 'o';
        end

        % Data
        u_inf = BL_Waves.(caze).u_inf;
        uv    = MN_Waves.(caze).phase(phase).uv;
        theta = BL_Waves.(caze).phase(phase).momentum * 1E-3;
        x     = BL_Waves.(caze).x;
    
        % Compute u*
        u_star_profile = sqrt(max(-1 * uv, [], 1, 'omitnan'));
    
        % Compute Cf
        cf = 2 * (u_star_profile / u_inf).^2;
    
        % Compute Re
        Re_theta = (u_inf * theta) / kinematic_viscosity;
    
        % Smooth
        cf = smoothdata(cf, 'gaussian', gauss_window);
        Re_theta = smoothdata(Re_theta, 'gaussian', gauss_window);
    
        % Crop
        bound = 0.25;
        centered_normalized_x = ((x - (range/2)) / wave_length);
        cf(centered_normalized_x < -bound) = nan;
        cf(centered_normalized_x > bound) = nan;
        Re_theta(centered_normalized_x < -bound) = nan;
        Re_theta(centered_normalized_x > bound) = nan;
    
        % Path Colors
        num_values = sum(~isnan(cf));
        temp_color = turbo(2 * num_values);
        disp(sum(~isnan(cf)))

        mask      = ~isnan(cf);
        cf_no_nan = cf(mask);
        Re_no_nan = Re_theta(mask);

        % Plot
        if phase == 2
            hold on
            plot(Re_theta, cf, 'color', 'black', 'linewidth', 1);
            for i = 1:length(cf_no_nan)
                scatter(Re_no_nan(i), cf_no_nan(i), size, marker, 'filled', 'MarkerFaceColor', temp_color(i,:), 'MarkerEdgeColor', 'none')
            end
            hold off
        end

        if phase == 4
            hold on
            plot(Re_theta, cf, 'color', 'black', 'linewidth', 1);
            for i = 1:length(cf_no_nan)
                scatter(Re_no_nan(i), cf_no_nan(i), size, marker, 'filled', 'MarkerFaceColor', temp_color(i + num_values,:), 'MarkerEdgeColor', 'none')
            end
            hold off
        end

        title(strcat('WV', WVs{w}))
        xlabel('$Re_{\theta}$', 'Interpreter', 'latex')
        ylabel('$C_f$', 'Interpreter', 'latex')
        
    end
end
linkaxes(h, 'xy')



%% UV Profiles at Max u* Location

gauss_window = 1;
WT = 4;
WV_colors = cool(length(WVs));
ref_colors = parula(length(WVs));

ax = figure();
t = tiledlayout(1,4);
for i = 1:4
    h(i) = nexttile();
    hold on
    for w = 1:length(WVs)
    
        caze = strcat('WT', num2str(WT), '_WV', WVs{w}, '_AG0');

        if contains(caze, 'WVA') == 1
            color = WV_colors(1,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
        elseif contains(caze, 'WVB') == 1
            color = WV_colors(2,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
        elseif contains(caze, 'WVC') == 1
            color = WV_colors(3,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
        elseif contains(caze, 'WVD') == 1
            color = WV_colors(4,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
        end
    
        % Data
        u_inf = BL_Waves.(caze).u_inf;
        uv    = MN_Waves.(caze).phase(i).uv;
        y     = flipud(unique(MN_Waves.(caze).Y));
    
        % Compute u*
        u_star_profile = sqrt(max(-1 * uv, [], 1, 'omitnan'));
        u_star_profile = smoothdata(u_star_profile, 'gaussian', gauss_window);
        [max_u_star, max_u_star_index] = max(u_star_profile);

        % Get Profile
        uv_profile = -1 * uv(:,max_u_star_index);
        uv_profile = smoothdata(uv_profile, 'gaussian', 5);
        norm_uv_profile = uv_profile / (u_inf^2);

        % Plot
        plot(norm_uv_profile, y, 'color', color, 'linewidth', 2)
        xlabel('$-\bar{uv} / u_{\infty}$', 'Interpreter', 'latex')
        ylabel('$y [mm]$', 'Interpreter', 'latex')
        ylim([min(y), 200])

        tit = strcat('Phase', {' '}, num2str(i));
        title(tit{1})

    end
    hold off
end
linkaxes(h,'xy')


%% Profiles + Contours

caze = "WT4_WVD_AG0";
step = 10;
scale = 600;
gauss_smooth = 5;
phase = 3;
lw = 1.5;


X = MN_Waves.(caze).X;
Y = MN_Waves.(caze).Y;
y = flipud(unique(Y));


ax = figure();
sgtitle(caze, 'interpreter', 'none')
hold on
contourf(X,Y,MN_Waves.(caze).phase(phase).uv, 500, 'linestyle', 'none');
for i = 1:step:length(x)
    minus_uv_profile = -1 * MN_Waves.(caze).phase(phase).uv(:,i);
    minus_uv_profile = smoothdata(minus_uv_profile, 'gaussian', gauss_smooth);

    plot(scale * minus_uv_profile + x(i), y, 'color', 'black', 'linewidth', lw)
    xl = xline(x(i), 'color', 'black', 'linestyle', '--');
    xl.Color(4) = 0.5;
end
plot(x, MN_Waves.(caze).phase(phase).max_wave_profile, 'color', 'red', 'linewidth', 2)
plot(x, MN_Waves.(caze).phase(phase).reference_wave, 'color', 'green', 'linewidth', 2)
hold off
axis equal
xlim([0, max(x)])
ylim([-20, 199])


%% Test
 
caze = "WT4_WVC_AG0";
phase = 4;

figure()
contourf(MN_Waves.(caze).X, MN_Waves.(caze).Y, -1 * MN_Waves.(caze).phase(phase).uv, 100, 'linestyle', 'none')
axis equal
colorbar()












