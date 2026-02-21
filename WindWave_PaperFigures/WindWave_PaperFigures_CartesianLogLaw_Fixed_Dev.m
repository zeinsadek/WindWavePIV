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

figure_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/pdf_test6';

% Approximate wavelengths in mm for labeling plots
wavelength_names.A = '410';
wavelength_names.B = '313';
wavelength_names.C = '189';
wavelength_names.D = '124';

steepnesses_names.A = '0.180';
steepnesses_names.B = '0.211';
steepnesses_names.C = '0.305';
steepnesses_names.D = '0.267';


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


%% Cartesian phase average log laws binned by phase: SHIFTED

% Font sizes
tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;
titleFontSize = 10;

% Plot details
wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
wind = '6';
phase_names = {'$\varphi = 0$', '$\varphi = \lambda / 4$', '$\varphi = \lambda / 2$', '$\varphi = 3 \lambda / 4$'};
lw = 1;
sz = 2;

grid_trans = 0.05;

% Smooth wall reference
karman_constant = 0.41;
c_plus = 5;
smooth_y_plus = 0:10:1E4;
smooth_u_plus = (1/karman_constant) * log(smooth_y_plus) + c_plus;

% Plot
clc; close all
ax = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,9]); %#ok<*NASGU>
t = tiledlayout(2, 4, 'padding', 'tight', 'TileSpacing', 'compact');

% Make array to hold all Log-region bounds
c = 1;
tmp_wave_log_regions = nan(16, 2); 

% Loop through phases
for phase = 1:4
    % Start tile
    h(phase) = nexttile;
    hold on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
    title(phase_names{phase}, 'interpreter', 'latex', 'fontsize', titleFontSize)

    % Plot smooth wall
    plot(smooth_y_plus, smooth_u_plus, 'linestyle', '--', 'HandleVisibility', 'off', 'color', 'black', 'linewidth', 1)

    % Loop through waves
    for WV = 1:length(WVs)

        % Get numbers
        wave = WVs{WV};
        caze = strcat('WT', wind, '_WV', wave, '_AG0');
        disp(caze)
        u_plus = LL_Waves.(caze).phase(phase).u_plus;
        y_plus = LL_Waves.(caze).phase(phase).y_plus;

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

        label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses_names.(wave));
        plot(y_plus, u_plus + deltaU, 'linewidth', lw, 'color', wave_colors{WV}, 'HandleVisibility', 'off')
        scatter(y_plus, u_plus + deltaU, sz, 'filled', 'HandleVisibility', 'on', 'MarkerFaceColor', wave_colors{WV}, 'DisplayName', label)

        if phase > 1
            set(gca, 'YTickLabel', [])
        end
        axis square

        % Increment case counter
        c = c + 1;
    end

    no_wave_caze = strcat('WT', wind, '_WV0_AGP');
    ens_u_plus = LL_No_Waves.(no_wave_caze).ensemble.u_plus;
    ens_y_plus = LL_No_Waves.(no_wave_caze).ensemble.y_plus;

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
    deltaU = mean(u_smooth + u_fit, 'omitnan');


    hold off
    hold off
    xscale('log')
    grid on
    tmp = gca;
    tmp.GridAlpha = grid_trans;
    tmp.MinorGridAlpha = grid_trans;

    xlim([1E1, 1E4])
    ylim([10, 27])

    if phase == 1
        ylabel('$\langle u \rangle^+ + \Delta U^+$', 'interpreter', 'latex', 'fontsize', labelFontSize)
    end
end


linkaxes(h, 'xy');
leg = legend('Orientation', 'Horizontal', 'box', 'off', 'interpreter', 'latex', 'fontsize', legendFontSize);
leg.Layout.Tile = 'North';

% Add a centered xlabel across the top row
annotation(ax, 'textbox', [0.3, 0.4, 0.44, 0.04], ...
           'String', '$y^+$', ...
           'Interpreter', 'latex', ...
           'FontSize', labelFontSize, ...
           'HorizontalAlignment', 'center', ...
           'VerticalAlignment', 'bottom', ...
           'EdgeColor', 'none');

% Print useful metrics
clc;
fprintf('Avg Log Region Bounds = %3.2f - %4.2f\n', fliplr(mean(tmp_wave_log_regions, 1, 'omitnan')))
fprintf('Min Bound = %3.2f\nMax bound = %4.2f\n', min(tmp_wave_log_regions, [], 'all'), max(tmp_wave_log_regions, [], 'all'))




%%% Bar chart showing different friction velocities
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

% clc; close all;
% ax = figure('color', 'white', 'units', 'centimeters', 'position', [2,2,30,10]);
ax2 = nexttile([1,4]);
b = bar({'$\varphi = 0$', '$\varphi = \lambda / 4$', '$\varphi = \lambda / 2$', '$\varphi = 3 \lambda / 4$'}, friction_velocity_bar_chart);
yline(no_wave_friction_velocity, '-', 'No Wave', 'Interpreter', 'latex', 'HandleVisibility', 'off', 'linewidth', 1, 'FontSize', 6)
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
ylabel('$u^*$ [m/s]', 'interpreter', 'latex', 'fontsize', labelFontSize)
box off
ylim([0, 0.45])
yticks(0:0.1:0.4)


% Add a,b labels
addPanelLabelsFixed(ax, [h(1), ax2], {'a', 'b'}, 'FontSize', 10, 'OffsetPts', [-30,0])


% Save figure
% pause(3)
% figure_name = ['LogLaw_Cartesian_', ['WT', wind], '_LogLawCombinedPlot_Shifted.pdf'];
% exportgraphics(ax, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 600, 'ContentType', 'image'); 
% close all




%% Cartesian phase average log laws binned by phase: UNSHIFTED

% Font sizes
tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;
titleFontSize = 10;

% Plot details
wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
wind = '6';
phase_names = {'$\varphi = 0$', '$\varphi = \lambda / 4$', '$\varphi = \lambda / 2$', '$\varphi = 3 \lambda / 4$'};
lw = 1;
sz = 2;

grid_trans = 0.05;

% Smooth wall reference
karman_constant = 0.41;
c_plus = 5;
smooth_y_plus = 0:10:1E4;
smooth_u_plus = (1/karman_constant) * log(smooth_y_plus) + c_plus;

% Plot
clc; close all
ax = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,7.5]); %#ok<*NASGU>
t = tiledlayout(2, 4, 'padding', 'tight', 'TileSpacing', 'compact');

% Make array to hold all Log-region bounds
c = 1;
tmp_wave_log_regions = nan(16, 2); 

% Loop through phases
for phase = 1:4
    % Start tile
    h(phase) = nexttile;
    hold on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
    title(phase_names{phase}, 'interpreter', 'latex', 'fontsize', titleFontSize)

    % Plot smooth wall
    plot(smooth_y_plus, smooth_u_plus, 'linestyle', '--', 'HandleVisibility', 'off', 'color', 'black', 'linewidth', 1)

    % Loop through waves
    for WV = 1:length(WVs)

        % Get numbers
        wave = WVs{WV};
        caze = strcat('WT', wind, '_WV', wave, '_AG0');
        disp(caze)
        u_plus = LL_Waves.(caze).phase(phase).u_plus;
        y_plus = LL_Waves.(caze).phase(phase).y_plus;

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

        label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses_names.(wave));
        plot(y_plus, u_plus, 'linewidth', lw, 'color', wave_colors{WV}, 'HandleVisibility', 'off')
        scatter(y_plus, u_plus, sz, 'filled', 'HandleVisibility', 'on', 'MarkerFaceColor', wave_colors{WV}, 'DisplayName', label)

        if phase > 1
            set(gca, 'YTickLabel', [])
        end
        axis square

        % Increment case counter
        c = c + 1;
    end

    no_wave_caze = strcat('WT', wind, '_WV0_AGP');
    ens_u_plus = LL_No_Waves.(no_wave_caze).ensemble.u_plus;
    ens_y_plus = LL_No_Waves.(no_wave_caze).ensemble.y_plus;

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
    deltaU = mean(u_smooth + u_fit, 'omitnan');


    hold off
    hold off
    xscale('log')
    grid on
    tmp = gca;
    tmp.GridAlpha = grid_trans;
    tmp.MinorGridAlpha = grid_trans;

    xlim([1E1, 1E4])
    ylim([0, 30])

    if phase == 1
        ylabel('$\langle u \rangle^+$', 'interpreter', 'latex', 'fontsize', labelFontSize)
    end
end


linkaxes(h, 'xy');
leg = legend('Orientation', 'Horizontal', 'box', 'off', 'interpreter', 'latex', 'fontsize', legendFontSize);
leg.Layout.Tile = 'North';

% Add a centered xlabel across the top row
annotation(ax, 'textbox', [0.32 0.35, 0.44, 0.04], ...
           'String', '$y^+$', ...
           'Interpreter', 'latex', ...
           'FontSize', labelFontSize, ...
           'HorizontalAlignment', 'center', ...
           'VerticalAlignment', 'bottom', ...
           'EdgeColor', 'none');

% Print useful metrics
clc;
fprintf('Avg Log Region Bounds = %3.2f - %4.2f\n', fliplr(mean(tmp_wave_log_regions, 1, 'omitnan')))
fprintf('Min Bound = %3.2f\nMax bound = %4.2f\n', min(tmp_wave_log_regions, [], 'all'), max(tmp_wave_log_regions, [], 'all'))




%%% Bar chart showing different friction velocities
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


ax2 = nexttile([1,4]);
b = bar({'$\varphi = 0$', '$\varphi = \lambda / 4$', '$\varphi = \lambda / 2$', '$\varphi = 3 \lambda / 4$'}, friction_velocity_bar_chart);
yline(no_wave_friction_velocity, '-', 'No Wave', 'Interpreter', 'latex', 'HandleVisibility', 'off', 'linewidth', 1, 'FontSize', 6)
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
ylabel('$u^*$ [m/s]', 'interpreter', 'latex', 'fontsize', labelFontSize)
box off
ylim([0, 0.45])
yticks(0:0.1:0.4)

% Add a,b labels
addPanelLabelsFixed(ax, [h(1), ax2], {'a', 'b'}, 'FontSize', 10, 'OffsetPts', [-30,0])


% Save figure
pause(3)
figure_name = ['LogLaw_Cartesian_', ['WT', wind], '_LogLawCombinedPlot_Unshifted.pdf'];
exportgraphics(ax, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 600, 'ContentType', 'image'); 
close all




%% Try plotting \Delta U^+ against wave phase (1-4)

wind_speeds = {'WT4', 'WT6', 'WT8'};
waves = {'A', 'B', 'C', 'D'};

% Get \Delta U^+ for each wind speed
for phase = 1:4

    for s = 1:length(wind_speeds)
        wind_speed = wind_speeds{s};
        % Loop through waves
        for w = 1:length(waves)
    
            % Get numbers
            wave = waves{w};
            caze = strcat(wind_speed, '_WV', wave, '_AG0');
            disp(caze)
            u_plus = LL_Waves.(caze).phase(phase).u_plus;
            y_plus = LL_Waves.(caze).phase(phase).y_plus;
           
            % % Law of the wall
            % u_plus = u_profile / friction_velocity;
            % y_plus = (y * friction_velocity) / nu;
            % 
    
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


% Phase \Delta U^+ vs wave phase
clc; close all
figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,7]); %#ok<*NASGU>
tiledlayout(1,1,'padding', 'tight')
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

    % Shaded region
    hFill = patch( ...
    [1:5, fliplr(1:5)], ...
    [tmp_max, fliplr(tmp_min)], ...
    hex2rgb(wave_colors{w}), ...
    'FaceAlpha', 0.15, ...        % transparency
    'EdgeColor', 'none', ...      % no outline
    'HandleVisibility', 'off');   % keep legend clean
    uistack(hFill, 'bottom')
end

% Mark were smooth wall is
yline(0, 'linestyle', '--', 'Label', 'Smooth Wall', 'interpreter', 'latex', ...
      'fontsize', 6, 'HandleVisibility', 'off', 'linewidth', 1)

% Plot representtive wave profile
fake_x = 0:0.01:6;
amplitude = 1.5;
offset = -4;
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
    plot(nan, nan, 'linewidth', lw, 'color', wave_colors{w}, ...
         'DisplayName', label)
end

leg = legend('interpreter', 'latex', 'box', 'off', 'orientation', 'horizontal', 'fontsize', legendFontSize);
leg.Layout.Tile = 'north';
leg.IconColumnWidth = 19;

hold off
xticks(1:4)
xlim([0.5, 6])
ylim([-7 18])
yticks(0:5:15)
ylabel('$\overline{\left(\Delta U^+ \right)}_{u_{\infty}}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xticks([1 2 3 4 5])
xticklabels([phase_names(reordered_phases), phase_names(1)])
tmp = gca;
tmp.FontSize = labelFontSize;
set(tmp, 'TickLabelInterpreter', 'latex')

% Save figure
% pause(3)
% figure_name = 'LogLaw_Cartesian_RoughnessFunction_WavePhase.pdf';
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
ylabel('$\overline{\left(\Delta U^+\right)}_{\varphi}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlabel('$c \mathbin{/} u_{\infty}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
tmp = gca;
tmp.FontSize = labelFontSize;
set(tmp, 'TickLabelInterpreter', 'latex')


% Save figure
% pause(3)
% figure_name = 'LogLaw_Cartesian_RoughnessFunction_WaveAge.pdf';
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


        % avg_friction_velocity = mean(LL_Waves.(caze).phase(:).u_star);

        
    
        % Get data across all 4 phases
        friction_velocities = nan(1,4);
        tmpY = nan(1,4);
        for phase = 1:4
            tmpY(phase) = wave_DeltaUs(phase).(caze);
            friction_velocities(phase) = LL_Waves.(caze).phase(phase).u_star;
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
ylabel('$\overline{\left(\Delta U^+\right)}_{\varphi}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlabel('$c \mathbin{/} \overline{ \left( u^* \right)}_{\varphi}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
tmp = gca;
tmp.FontSize = labelFontSize;
set(tmp, 'TickLabelInterpreter', 'latex')


% Save figure
% pause(3)
% figure_name = 'LogLaw_Cartesian_RoughnessFunction_FrictionWaveAge.pdf';
% exportgraphics(ax, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 600, 'ContentType', 'image'); 
% close all













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

