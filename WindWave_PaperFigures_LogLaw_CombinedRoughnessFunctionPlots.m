%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMBINED CARTESIAN AND CURVILINEAR LOG LAW COMPARISON
% Generates side-by-side comparison plots for:
%   1. Delta U+ vs wave phase
%   2. Delta U+ vs wave age (c/u_inf)
%   3. Delta U+ vs friction wave age (c/u*)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% PATHS AND SETUP
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions/');
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')

experiment_log  = readcell('Offshore_Inflow_Log.xlsx');
wave_parameters = readcell('Offshore_Waves.xlsx');

s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";

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

% Physical parameters
karman_constant = 0.41;
c_plus = 5;
nu = 1.46E-5;

% Case identifiers
WTs = {'4', '6', '8'};
WVs = {'A', 'B', 'C', 'D'};
wind_speeds = {'WT4', 'WT6', 'WT8'};
waves = {'A', 'B', 'C', 'D'};

% Wave properties
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

% Plot settings
wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
wind_speed_markers = {'^', 'square', 'o'};
phase_names = {'$\varphi = 0$', '$\varphi = \lambda / 4$', '$\varphi = \lambda / 2$', '$\varphi = 3 \lambda / 4$'};
reordered_phases = [1, 4, 3, 2];


%% =====================================================
% SECTION 1: LOAD CARTESIAN DATA
% ======================================================
fprintf('=== Loading Cartesian Data ===\n')

BL_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/boundary_layer';
LL_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/log_law_fixed/';
MN_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means/';

% Load cases with waves
for s_idx = 1:length(WTs)
    speed = WTs{s_idx};
    for i = 1:length(WVs)
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
for s_idx = 1:length(WTs)
    no_wave_case = strcat('WT', num2str(WTs{s_idx}), '_WV0_AGP');
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

fprintf('\n')


%% =====================================================
% SECTION 2: LOAD CURVILINEAR DATA
% ======================================================
fprintf('=== Loading Curvilinear Data ===\n')

project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
means_path = fullfile(project_path, 'means');
curvilinear_path = fullfile(project_path, 'curvilinear_new');

for s_idx = 1:length(wind_speeds)
    wind_speed = wind_speeds{s_idx};
    
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
        fprintf('Loading Case: %s...\n', caze)

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

            % Clean up U profile
            u_profile = StatisticalGradientFilter(u_profile, stat_grad_tol, 9);
            u_profile(find(~isnan(u_profile), 1, 'first')) = nan;
        
            % Get Friction Velocity
            u_star_curv.(caze)(phase).raw = max(sqrt(-uv_profile), [], 'all', 'omitnan');
            u_star_curv.(caze)(phase).u_profile = u_profile;
            u_star_curv.(caze)(phase).u_profile_normalized = u_profile / (u_inf);
            u_star_curv.(caze)(phase).uv_profile = uv_profile;
            u_star_curv.(caze)(phase).uv_profile_normalized = uv_profile / (u_inf^2);
            u_star_curv.(caze)(phase).y = y_profile;
        end
    end
    fprintf('\n')
end

% Also include the no-wave cases for curvilinear
for s_idx = 1:length(wind_speeds)
    wind_speed = wind_speeds{s_idx};
    caze = strcat(wind_speed, '_WV0_AGP');
    fprintf('Loading Case: %s...\n', caze)

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

    u_profile(find(~isnan(u_profile), 1, 'first')) = nan;

    % Get Friction Velocity
    u_star_curv.(caze).raw = max(sqrt(-uv_profile), [], 'all', 'omitnan');
    u_star_curv.(caze).u_profile = u_profile;
    u_star_curv.(caze).u_profile_normalized = u_profile / (u_inf);
    u_star_curv.(caze).uv_profile = uv_profile;
    u_star_curv.(caze).uv_profile_normalized = uv_profile / (u_inf^2);
    u_star_curv.(caze).y = y_profile;
end

fprintf('\n')


%% =====================================================
% SECTION 3: COMPUTE DELTA U+ FOR CARTESIAN
% ======================================================
fprintf('=== Computing Delta U+ (Cartesian) ===\n')

for phase = 1:4
    for s_idx = 1:length(wind_speeds)
        wind_speed = wind_speeds{s_idx};
        for w = 1:length(waves)
            wave = waves{w};
            caze = strcat(wind_speed, '_WV', wave, '_AG0');
            
            u_plus = LL_Waves.(caze).phase(phase).u_plus;
            y_plus = LL_Waves.(caze).phase(phase).y_plus;

            % Find constant slope region (log region)
            [bounds, ~] = find_constant_slope_region(y_plus, u_plus, ...
                                                              'plateau_tol', 0.2, ...
                                                              'smooth_window', 7);
        
            % Compute Delta U+
            in_region = (y_plus >= bounds(2)) & (y_plus <= bounds(1));
            y_fit = y_plus(in_region);
            u_fit = u_plus(in_region);
            
            % Smooth wall reference
            u_smooth = (1/karman_constant) * log(y_fit) + c_plus;
            
            % Delta U+
            deltaU = mean(u_smooth - u_fit, 'omitnan');

            % Save Delta U+
            wave_DeltaUs_Cart(phase).(caze) = deltaU;
        end
    end
end
fprintf('Done.\n\n')


%% =====================================================
% SECTION 4: COMPUTE DELTA U+ FOR CURVILINEAR
% ======================================================
fprintf('=== Computing Delta U+ (Curvilinear) ===\n')

for phase = 1:4
    for s_idx = 1:length(wind_speeds)
        wind_speed = wind_speeds{s_idx};
        for w = 1:length(waves)
            wave = waves{w};
            caze = strcat(wind_speed, '_WV', wave, '_AG0');

            % Get data
            y = u_star_curv.(caze)(phase).y;
            u_profile = u_star_curv.(caze)(phase).u_profile;
            friction_velocity = u_star_curv.(caze)(phase).raw;
           
            % Law of the wall
            u_plus = u_profile / friction_velocity;
            y_plus = (y * friction_velocity) / nu;

            % Find constant slope region (log region)
            [bounds, ~] = find_constant_slope_region(y_plus, u_plus, ...
                                                              'plateau_tol', 0.2, ...
                                                              'smooth_window', 7);
        
            % Compute Delta U+
            in_region = (y_plus >= bounds(2)) & (y_plus <= bounds(1));
            y_fit = y_plus(in_region);
            u_fit = u_plus(in_region);
            
            % Smooth wall reference
            u_smooth = (1/karman_constant) * log(y_fit) + c_plus;
            
            % Delta U+
            deltaU = mean(u_smooth - u_fit, 'omitnan');

            % Save Delta U+
            wave_DeltaUs_Curv(phase).(caze) = deltaU;
        end
    end
end
fprintf('Done.\n\n')


%% =====================================================
% FIGURE 1: DELTA U+ VS WAVE PHASE (SIDE-BY-SIDE)
% ======================================================
fprintf('=== Generating Figure 1: Delta U+ vs Wave Phase ===\n')

wave_transparency = 0.25;

% Font sizes
tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;
% titleFontSize = 12;

% Plot sizes
sz = 20;
lw = 1.5;

clc; close all
fig1 = figure('color', 'white', 'units', 'centimeters', 'position', [10, 10, 12, 9]);
tiledlayout(2, 1, 'padding', 'loose', 'TileSpacing', 'compact')

% -------------------- LEFT: CARTESIAN --------------------
ax1 = nexttile;
hold on
% title('$\overline{\left(\Delta U^+ \right)}_{u_{\infty}}$', 'interpreter', 'latex', 'fontsize', titleFontSize)
set(ax1, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% Compute mean/max/min for Cartesian
for p = 1:4
    phase = reordered_phases(p);
    for w = 1:length(waves)
        wave = waves{w};
        tmp = nan(1,3);
        for s_idx = 1:length(wind_speeds)
            wind_speed = wind_speeds{s_idx};
            caze = strcat(wind_speed, '_WV', wave, '_AG0');
            tmp(s_idx) = wave_DeltaUs_Cart(phase).(caze);
        end
        
        scatter(p, mean(tmp, 'all', 'omitnan'), sz, 'filled', ...
                'markerfacecolor', wave_colors{w}, 'HandleVisibility', 'off')
        if p == 1
             scatter(5, mean(tmp, 'all', 'omitnan'), sz, 'filled', ...
                'markerfacecolor', wave_colors{w}, 'HandleVisibility', 'off')
        end
        mean_DeltaUs_Cart(phase).(wave) = mean(tmp, 'all', 'omitnan');
        max_DeltaUs_Cart(phase).(wave) = max(tmp, [], 'all', 'omitnan');
        min_DeltaUs_Cart(phase).(wave) = min(tmp, [], 'all', 'omitnan');
    end
end

% Plot lines and shading for Cartesian
for w = 1:length(waves)
    wave = waves{w};
    tmp = nan(1,5); tmp_max = nan(1,5); tmp_min = nan(1,5);
    for p = 1:4
        phase = reordered_phases(p);
        tmp(p) = mean_DeltaUs_Cart(phase).(wave);
        tmp_max(p) = max_DeltaUs_Cart(phase).(wave);
        tmp_min(p) = min_DeltaUs_Cart(phase).(wave);
    end
    tmp(5) = tmp(1); tmp_max(5) = tmp_max(1); tmp_min(5) = tmp_min(1);
    
    plot(1:5, tmp, 'linewidth', lw, 'color', wave_colors{w}, 'HandleVisibility', 'off')
    % hFill = patch([1:5, fliplr(1:5)], [tmp_max, fliplr(tmp_min)], ...
    %     hex2rgb(wave_colors{w}), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    % uistack(hFill, 'bottom')
end

% Smooth wall line
yline(0, 'linestyle', '--', 'Label', 'Smooth Wall', 'interpreter', 'latex', ...
      'fontsize', 6, 'HandleVisibility', 'off', 'linewidth', 1)

% Representative wave profile
fake_x = 0:0.01:6;
amplitude = 2.5;
offset = -5;
reference_wave = amplitude * cos(((2*pi) / 4) * (fake_x - 1)) + offset;
plot(fake_x, reference_wave, 'color', 'black', 'linewidth', 1.5, 'HandleVisibility', 'off')
scatter(1:5, [offset + amplitude, offset, offset - amplitude, offset, offset + amplitude], sz, 'filled', ...
        'MarkerFaceColor', 'black', 'HandleVisibility', 'off')

X = [fake_x(:); flipud(fake_x(:))];
Y = [reference_wave(:); (-10)*ones(size(reference_wave(:)))];
patch(X, Y, 'k', 'FaceAlpha', wave_transparency, 'EdgeColor', 'none', 'handlevisibility', 'off');

hold off
% xticks([1 2 3 4 5])
% xticklabels([phase_names(reordered_phases), phase_names(1)])
set(ax1,'xtick',[])
set(ax1,'xticklabel',[])

xlim([0.5, 6])
ylim([-7 18])
yticks(0:5:15)
ylabel('$\overline{\left(\Delta U^+ \right)}_{u_{\infty}}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
set(ax1, 'FontSize', tickFontSize, 'TickLabelInterpreter', 'latex')


% -------------------- RIGHT: CURVILINEAR --------------------
ax2 = nexttile;
hold on
% title('$\overline{\left(\Delta U_{\xi}^+ \right)}_{u_{\infty}}$', 'interpreter', 'latex', 'fontsize', titleFontSize)
set(ax2, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% Compute mean/max/min for Curvilinear
for p = 1:4
    phase = reordered_phases(p);
    for w = 1:length(waves)
        wave = waves{w};
        tmp = nan(1,3);
        for s_idx = 1:length(wind_speeds)
            wind_speed = wind_speeds{s_idx};
            caze = strcat(wind_speed, '_WV', wave, '_AG0');
            tmp(s_idx) = wave_DeltaUs_Curv(phase).(caze);
        end
        
        scatter(p, mean(tmp, 'all', 'omitnan'), sz, 'filled', ...
                'markerfacecolor', wave_colors{w}, 'HandleVisibility', 'off')
        if p == 1
             scatter(5, mean(tmp, 'all', 'omitnan'), sz, 'filled', ...
                'markerfacecolor', wave_colors{w}, 'HandleVisibility', 'off')
        end
        mean_DeltaUs_Curv(phase).(wave) = mean(tmp, 'all', 'omitnan');
        max_DeltaUs_Curv(phase).(wave) = max(tmp, [], 'all', 'omitnan');
        min_DeltaUs_Curv(phase).(wave) = min(tmp, [], 'all', 'omitnan');
    end
end

% Plot lines and shading for Curvilinear
for w = 1:length(waves)
    wave = waves{w};
    tmp = nan(1,5); tmp_max = nan(1,5); tmp_min = nan(1,5);
    for p = 1:4
        phase = reordered_phases(p);
        tmp(p) = mean_DeltaUs_Curv(phase).(wave);
        tmp_max(p) = max_DeltaUs_Curv(phase).(wave);
        tmp_min(p) = min_DeltaUs_Curv(phase).(wave);
    end
    tmp(5) = tmp(1); tmp_max(5) = tmp_max(1); tmp_min(5) = tmp_min(1);
    
    plot(1:5, tmp, 'linewidth', lw, 'color', wave_colors{w}, 'HandleVisibility', 'off')
    % hFill = patch([1:5, fliplr(1:5)], [tmp_max, fliplr(tmp_min)], ...
    %     hex2rgb(wave_colors{w}), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    % uistack(hFill, 'bottom')
end

% Smooth wall line
yline(0, 'linestyle', '--', 'Label', 'Smooth Wall', 'interpreter', 'latex', ...
      'fontsize', 6, 'HandleVisibility', 'off', 'linewidth', 1)

% Representative wave profile
plot(fake_x, reference_wave, 'color', 'black', 'linewidth', 1.5, 'HandleVisibility', 'off')
scatter(1:5, [offset + amplitude, offset, offset - amplitude, offset, offset + amplitude], sz, 'filled', ...
        'MarkerFaceColor', 'black', 'HandleVisibility', 'off')
patch(X, Y, 'k', 'FaceAlpha', wave_transparency, 'EdgeColor', 'none', 'handlevisibility', 'off');

% Add legend for wavelengths
for w = 1:length(waves)
    wave = waves{w};
    label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses_names.(wave));
    plot(nan, nan, 'linewidth', lw, 'color', wave_colors{w}, 'DisplayName', label)
end

leg = legend('interpreter', 'latex', 'box', 'off', 'orientation', 'horizontal', 'fontsize', legendFontSize);
leg.Layout.Tile = 'north';
leg.IconColumnWidth = 19;

hold off
xticks([1 2 3 4 5])
xticklabels([phase_names(reordered_phases), phase_names(1)])
yticks(0:5:15)
ylabel('$\overline{\left(\Delta U_{\xi}^+ \right)}_{u_{\infty}}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
set(ax2, 'FontSize', tickFontSize, 'TickLabelInterpreter', 'latex')

% Link axes
linkaxes([ax1, ax2], 'xy')
xlim([0.5, 6])
ylim([-10 15])

% Add a,b labels
addPanelLabels([ax1, ax2], {'a', 'b'}, 'FontSize', 10, 'Offset', [-0.12,1.2])

% Uncomment to save
pause(3)
figure_name = 'LogLaw_Combined_RoughnessFunction_WavePhase.pdf';
exportgraphics(fig1, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 600, 'ContentType', 'image');
close all



%% =====================================================
% FIGURE 2: DELTA U+ VS WAVE AGE (c/u_inf) (SIDE-BY-SIDE)
% ======================================================
fprintf('=== Generating Figure 2: Delta U+ vs Wave Age (c/u_inf) ===\n')

% Font sizes
tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;
titleFontSize = 12;

% Plot line sizes
sz = 20;
lw = 1.5;

clc; close all
fig2 = figure('color', 'white', 'units', 'centimeters', 'position', [10, 10, 13, 12]);
tiledlayout(2, 1, 'padding', 'compact', 'TileSpacing', 'compact')

% -------------------- LEFT: CARTESIAN --------------------
ax1 = nexttile;
hold on
set(ax1, 'FontSize', tickFontSize, 'TickLabelInterpreter', 'latex')
% title('Cartesian', 'interpreter', 'latex', 'fontsize', titleFontSize)

for w = 1:length(waves)
    wave = waves{w};
    wavelength = wavelengths.(wave);
    frequency = frequencies.(wave);
    wave_speed = (wavelength * 1E-3 * frequency);

    xplot = nan(1,3); yplot = nan(1,3);
    maxs = nan(1,3); mins = nan(1,3);

    for s_idx = 1:length(wind_speeds)
        wind_speed = wind_speeds{s_idx};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        freestream = freestreams.(wind_speed);
        wave_age = wave_speed / freestream;
    
        tmpY = nan(1,4);
        for phase = 1:4
            tmpY(phase) = wave_DeltaUs_Cart(phase).(caze);
        end

        scatter(wave_age, mean(tmpY, 'all', 'omitnan'), sz, 'filled', wind_speed_markers{s_idx}, ...
                'markerfacecolor', wave_colors{w}, 'HandleVisibility', 'off')

        xplot(s_idx) = wave_age;
        yplot(s_idx) = mean(tmpY, 'all', 'omitnan');
        maxs(s_idx) = max(tmpY, [], 'all', 'omitnan');
        mins(s_idx) = min(tmpY, [], 'all', 'omitnan');
    end
    plot(xplot, yplot, 'linewidth', lw, 'color', wave_colors{w}, 'HandleVisibility', 'off');
    hFill = patch([xplot, fliplr(xplot)], [maxs, fliplr(mins)], ...
        hex2rgb(wave_colors{w}), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    uistack(hFill, 'bottom')
end

yline(0, 'linestyle', '--', 'Label', 'Smooth Wall', 'interpreter', 'latex', ...
      'fontsize', 6, 'HandleVisibility', 'off', 'linewidth', 1, 'labelhorizontalAlignment', 'left')
hold off
xlim([0.05, 0.35])
ylim([-4, 17])
yticks(-4:4:16)
ylabel('$\overline{\left(\Delta U^+\right)}_{\varphi}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% xlabel('$c / u_{\infty}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
set(ax1, 'FontSize', labelFontSize, 'TickLabelInterpreter', 'latex')
% set(ax1,'xtick',[])
set(ax1,'xticklabel',[])


% -------------------- RIGHT: CURVILINEAR --------------------
ax2 = nexttile;
set(ax2, 'FontSize', tickFontSize, 'TickLabelInterpreter', 'latex')
hold on
% title('Curvilinear', 'interpreter', 'latex', 'fontsize', titleFontSize)

for w = 1:length(waves)
    wave = waves{w};
    wavelength = wavelengths.(wave);
    frequency = frequencies.(wave);
    wave_speed = (wavelength * 1E-3 * frequency);

    xplot = nan(1,3); yplot = nan(1,3);
    maxs = nan(1,3); mins = nan(1,3);

    for s_idx = 1:length(wind_speeds)
        wind_speed = wind_speeds{s_idx};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        freestream = freestreams.(wind_speed);
        wave_age = wave_speed / freestream;
    
        tmpY = nan(1,4);
        for phase = 1:4
            tmpY(phase) = wave_DeltaUs_Curv(phase).(caze);
        end

        scatter(wave_age, mean(tmpY, 'all', 'omitnan'), sz, 'filled', wind_speed_markers{s_idx}, ...
                'markerfacecolor', wave_colors{w}, 'HandleVisibility', 'off')

        xplot(s_idx) = wave_age;
        yplot(s_idx) = mean(tmpY, 'all', 'omitnan');
        maxs(s_idx) = max(tmpY, [], 'all', 'omitnan');
        mins(s_idx) = min(tmpY, [], 'all', 'omitnan');
    end
    plot(xplot, yplot, 'linewidth', lw, 'color', wave_colors{w}, 'HandleVisibility', 'off');
    hFill = patch([xplot, fliplr(xplot)], [maxs, fliplr(mins)], ...
        hex2rgb(wave_colors{w}), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    uistack(hFill, 'bottom')
end

% Add legends
for w = 1:length(waves)
    wave = waves{w};
    label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses_names.(wave));
    plot(nan, nan, 'linewidth', 2, 'color', wave_colors{w}, 'DisplayName', label)
end
plot(nan, nan, 'color', 'white', 'displayname', ' ')
for s_idx = 1:length(wind_speeds)
    wind_speed = wind_speeds{s_idx};
    freestream = freestreams.(wind_speed);
    label = sprintf('$u_{\\infty} = %1.2f$ m/s', freestream);
    plot(nan, nan, wind_speed_markers{s_idx}, 'MarkerFaceColor', 'black', ...
        'MarkerEdgeColor', 'black', 'MarkerSize', 4, 'LineWidth', 1, ...
        'LineStyle', 'none', 'DisplayName', label);
end

leg = legend('interpreter', 'latex', 'box', 'off', 'fontsize', legendFontSize);
leg.Layout.Tile = 'east';
leg.IconColumnWidth = 19;

yline(0, 'linestyle', '--', 'Label', 'Smooth Wall', 'interpreter', 'latex', ...
      'fontsize', 6, 'HandleVisibility', 'off', 'linewidth', 1, 'labelhorizontalAlignment', 'left')
hold off
xlim([0.08, 0.35])
ylim([-4, 17])
yticks(-4:4:16)
ylabel('$\overline{\left(\Delta U_{\xi}^+\right)}_{\varphi}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlabel('$c / u_{\infty}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
set(gca, 'FontSize', labelFontSize, 'TickLabelInterpreter', 'latex')

linkaxes([ax1, ax2], 'xy')
xlim([0.08, 0.35])


% Uncomment to save
% pause(3)
% figure_name = 'LogLaw_Combined_RoughnessFunction_WaveAge.pdf';
% exportgraphics(fig2, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 600, 'ContentType', 'image');
% close all



%% =====================================================
% FIGURE 3: DELTA U+ VS FRICTION WAVE AGE (c/u*) (SIDE-BY-SIDE)
% ======================================================
fprintf('=== Generating Figure 3: Delta U+ vs Friction Wave Age (c/u*) ===\n')

% Font sizes
tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;
titleFontSize = 12;

% Plot line sizes
sz = 20;
lw = 1.5;

fig3 = figure('color', 'white', 'units', 'centimeters', 'position', [10, 10, 13, 12]);
tiledlayout(2, 1, 'padding', 'compact', 'TileSpacing', 'compact')

% -------------------- LEFT: CARTESIAN --------------------
ax1 = nexttile;
hold on
set(ax1, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% title('Cartesian', 'interpreter', 'latex', 'fontsize', titleFontSize)

for w = 1:length(waves)
    wave = waves{w};
    wavelength = wavelengths.(wave);
    frequency = frequencies.(wave);
    wave_speed = (wavelength * 1E-3 * frequency);

    xplot = nan(1,3); yplot = nan(1,3);
    maxs = nan(1,3); mins = nan(1,3);

    for s_idx = 1:length(wind_speeds)
        wind_speed = wind_speeds{s_idx};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
    
        friction_velocities = nan(1,4);
        tmpY = nan(1,4);
        for phase = 1:4
            tmpY(phase) = wave_DeltaUs_Cart(phase).(caze);
            friction_velocities(phase) = LL_Waves.(caze).phase(phase).u_star;
        end

        avg_friction_velocity = mean(friction_velocities, 'all', 'omitnan');
        wave_age = wave_speed / avg_friction_velocity;

        scatter(wave_age, mean(tmpY, 'all', 'omitnan'), sz, 'filled', wind_speed_markers{s_idx}, ...
                'markerfacecolor', wave_colors{w}, 'HandleVisibility', 'off')

        xplot(s_idx) = wave_age;
        yplot(s_idx) = mean(tmpY, 'all', 'omitnan');
        maxs(s_idx) = max(tmpY, [], 'all', 'omitnan');
        mins(s_idx) = min(tmpY, [], 'all', 'omitnan');
    end
    plot(xplot, yplot, 'linewidth', lw, 'color', wave_colors{w}, 'HandleVisibility', 'off');
    hFill = patch([xplot, fliplr(xplot)], [maxs, fliplr(mins)], ...
        hex2rgb(wave_colors{w}), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    uistack(hFill, 'bottom')
end

yline(0, 'linestyle', '--', 'Label', 'Smooth Wall', 'interpreter', 'latex', ...
      'fontsize', 6, 'HandleVisibility', 'off', 'linewidth', 1, 'labelhorizontalAlignment', 'left')
hold off
xlim([1, 7])
ylim([-4, 17])
yticks(-4:4:16)
ylabel('$\overline{\left(\Delta U^+\right)}_{\varphi}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlabel('$c \, / \, \overline{ \left( u^* \right)}_{\varphi}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
set(ax1, 'FontSize', labelFontSize, 'TickLabelInterpreter', 'latex')


% -------------------- RIGHT: CURVILINEAR --------------------
ax2 = nexttile;
set(ax2, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
hold on
% title('Curvilinear', 'interpreter', 'latex', 'fontsize', titleFontSize)

for w = 1:length(waves)
    wave = waves{w};
    wavelength = wavelengths.(wave);
    frequency = frequencies.(wave);
    wave_speed = (wavelength * 1E-3 * frequency);

    xplot = nan(1,3); yplot = nan(1,3);
    maxs = nan(1,3); mins = nan(1,3);

    for s_idx = 1:length(wind_speeds)
        wind_speed = wind_speeds{s_idx};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
    
        friction_velocities = nan(1,4);
        tmpY = nan(1,4);
        for phase = 1:4
            tmpY(phase) = wave_DeltaUs_Curv(phase).(caze);
            friction_velocities(phase) = u_star_curv.(caze)(phase).raw;
        end

        avg_friction_velocity = mean(friction_velocities, 'all', 'omitnan');
        wave_age = wave_speed / avg_friction_velocity;

        scatter(wave_age, mean(tmpY, 'all', 'omitnan'), sz, 'filled', wind_speed_markers{s_idx}, ...
                'markerfacecolor', wave_colors{w}, 'HandleVisibility', 'off')

        xplot(s_idx) = wave_age;
        yplot(s_idx) = mean(tmpY, 'all', 'omitnan');
        maxs(s_idx) = max(tmpY, [], 'all', 'omitnan');
        mins(s_idx) = min(tmpY, [], 'all', 'omitnan');
    end
    plot(xplot, yplot, 'linewidth', lw, 'color', wave_colors{w}, 'HandleVisibility', 'off');
    hFill = patch([xplot, fliplr(xplot)], [maxs, fliplr(mins)], ...
        hex2rgb(wave_colors{w}), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    uistack(hFill, 'bottom')
end

% Add legends
for w = 1:length(waves)
    wave = waves{w};
    label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses_names.(wave));
    plot(nan, nan, 'linewidth', 2, 'color', wave_colors{w}, 'DisplayName', label)
end
plot(nan, nan, 'color', 'white', 'displayname', ' ')
for s_idx = 1:length(wind_speeds)
    wind_speed = wind_speeds{s_idx};
    freestream = freestreams.(wind_speed);
    label = sprintf('$u_{\\infty} = %1.2f$ m/s', freestream);
    plot(nan, nan, wind_speed_markers{s_idx}, 'MarkerFaceColor', 'black', ...
        'MarkerEdgeColor', 'black', 'MarkerSize', 4, 'LineWidth', 1, ...
        'LineStyle', 'none', 'DisplayName', label);
end

leg = legend('interpreter', 'latex', 'box', 'off', 'fontsize', legendFontSize);
leg.Layout.Tile = 'east';
leg.IconColumnWidth = 19;

yline(0, 'linestyle', '--', 'Label', 'Smooth Wall', 'interpreter', 'latex', ...
      'fontsize', 6, 'HandleVisibility', 'off', 'linewidth', 1, 'labelhorizontalAlignment', 'left')
hold off
xlim([1, 7])
ylim([-4, 17])
yticks(-4:4:16)
ylabel('$\overline{\left(\Delta U_{\xi}^+\right)}_{\varphi}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlabel('$c \, / \, \overline{ \left( u_{\xi}^* \right)}_{\varphi}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
set(ax2, 'FontSize', labelFontSize, 'TickLabelInterpreter', 'latex')

linkaxes([ax1, ax2], 'xy')

% save figure
% pause(3)
% figure_name = 'LogLaw_Combined_RoughnessFunction_FrictionWaveAge.pdf';
% exportgraphics(fig3, fullfile(figure_folder, 'LogLaw', figure_name), 'Resolution', 600, 'ContentType', 'image');
% close all



%% =====================================================
% HELPER FUNCTIONS
% ======================================================

function addPanelLabels(ax, labels, varargin)
% addPanelLabels(ax, labels) adds (a),(b),... just OUTSIDE top-left of each axes.
% ax     : array of axes handles (e.g., from tiledlayout / findall)
% labels : cellstr like {'a','b','c'} or string array ["a" "b" "c"]
%
% Optional name-value:
% 'Offset'   : [dx dy] in normalized axes units (default [-0.10 1.02])
% 'FontSize' : default 12
% 'FontName' : default 'Times New Roman'

p = inputParser;
addParameter(p,'Offset',[-0.10 1.1]);
addParameter(p,'FontSize', 10);
addParameter(p,'FontName','Times New Roman');
parse(p,varargin{:});
off = p.Results.Offset;

labels = string(labels);
for k = 1:numel(ax)
    if ~isgraphics(ax(k),'axes'), continue; end

    % Plain parentheses + italic letter:
    % TeX interpreter: \it turns italic ON, \rm returns to roman.
    s = sprintf('(\\ita\\rm)');              % placeholder
    s = sprintf('(\\it%s\\rm)', labels(k));  % actual label

    text(ax(k), off(1), off(2), s, ...
        'Units','normalized', ...
        'Interpreter','tex', ...           % keeps italics control simple
        'FontName',p.Results.FontName, ...
        'FontSize',p.Results.FontSize, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'Clipping','off');                 % critical: allow outside axes
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
    q = RemoveIslands(q, max_island_width);
    output = q;
end


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