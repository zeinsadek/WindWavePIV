%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions');

% Cases with waves
wind_speeds = {'WT4', 'WT6', 'WT8'};
wind_speed_colors = {'#0075F2', '#FF8C42', '#D30C7B'};

% Paths
MN_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means/';
WV_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/wave/';
LL_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/log_law_fixed/';

% Wind speeds
freestreams.('WT4') = 2.4181;
freestreams.('WT6') = 3.8709;
freestreams.('WT8') = 5.4289;

% % Log Law Shifts
% DeltaUPluses.('WT4') = -2.35;
% DeltaUPluses.('WT6') = -2.02;
% DeltaUPluses.('WT8') = 2.54;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load cases with waves
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
   
    % Case
    caze = strcat(wind_speed, '_WV0_AGP');
    fprintf('Loading Case: %s...\n', caze)
    MN_path = fullfile(MN_folder, strcat(caze, '_MEANS.mat'));
    WV_path = fullfile(WV_folder, strcat(caze, '_WAVE.mat'));
    LL_path = fullfile(LL_folder, strcat(caze, '_LL_Fixed.mat'));
   
    % Save MN to structure
    means_temp = load(MN_path);
    means_temp = means_temp.output;
    means.(caze) = means_temp;

    % Save WV to structure
    waves_tmp = load(WV_path);
    waves_tmp = waves_tmp.output;
    waves.(caze) = waves_tmp;

    % Save LL to structure
    LL_temp = load(LL_path);
    LL_temp = LL_temp.output;
    log_law.(caze) = LL_temp;

end


% Compute log law profiles
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
    data = means.(caze);

    % Get friction velocity and cf
    u = data.ensemble.u;
    uv = data.ensemble.uv;

    % Get profiles
    idx = round(length(uv)/2);
    y_profile = data.Y(:, 1) * 1E-3; 
    uv_profile = uv(:, idx); 
    u_profile = u(:, idx);

    % CLEAN UP U PROFILE
    u_profile(find(~isnan(u_profile), 1, 'first')) = nan;

    % Save values
    friction_profiles.(caze).raw = max(sqrt(-uv_profile), [], 'all', 'omitnan');
    friction_profiles.(caze).u_profile = u_profile;
    friction_profiles.(caze).u_profile_normalized = u_profile / (u_inf);
    friction_profiles.(caze).uv_profile = uv_profile;
    friction_profiles.(caze).uv_profile_normalized = uv_profile / (u_inf^2);
    friction_profiles.(caze).y = y_profile;
end

% Clear temporary variables
clear BL_path LL_path MN_path
clear BL_temp LL_temp means_temp
clear caze i no_wave_case s speed
clear MN_folder LL_folder s WV_folder WV_path waves_tmp
clear caze wind_speed u uv data idx y_profile uv_profile u_profile u_inf


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOOK AT HIGH-RES SURFACE WAVES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% How mamy waves to plot
n_waves = 100;
skip = 5;
linewidth = 1;

% Plot
clc; close all
figure('color', 'white')
t = tiledlayout(3 ,1, 'TileSpacing', 'tight', 'padding', 'compact');
sgtitle("No Forced Waves: Free Surface Profile Fluctuations $\eta_{i}'$", 'interpreter', 'latex')

% Loop through wind speeds
for s = 1:length(wind_speeds)

    % Get case
    wind_speed = wind_speeds{s};
    caze = strcat(wind_speed, '_WV0_AGP');
    freestream = freestreams.(wind_speed);
    disp(caze)
    
    % Get wave information
    x = waves.(caze).x;
    test_waves = waves.(caze).wave_profiles;
    time_mean_waves = mean(test_waves, 1, 'omitnan');
    wave_fluctuations = test_waves - time_mean_waves;

    % Get friction velocity
    u_star = log_law.(caze).ensemble.u_star.u_star;

    % Plot tile
    h(s) = nexttile;
    title(sprintf('$u_{\\infty} = %1.2f, \\ u^* = %1.3f$ m/s', freestream, u_star), 'Interpreter', 'latex')
    hold on
    plot(x, wave_fluctuations(1:skip:n_waves,:), 'linewidth', linewidth)
    hold off
    
    % Hide x axis for top two plots
    if s ~= 3
        ax = gca;
        ax.XTickLabel = [];
        ax.XAxis.Visible = 'off';
    end

    % Limits
    axis equal
    xlim([-100, 100])
    ylim([-5,5])
    xticks(-100:20:100)
end

xlabel(t, '$x$ [mm]', 'interpreter', 'latex')
ylabel(t, '$y$ [mm]', 'interpreter', 'latex')

clear n_waves skip t s wind_speed caze freestream x test_waves time_mean_waves wave_fluctuations ax h linewidth


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE RMS VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kinematic viscosity
nu = 1.46E-5;

% Make arrays to save values
rms_elevations = nan(1,3);
rms_slopes = nan(1,3);
sublayers = nan(1,3);
Re_roughs = nan(1,3);

% Loop through wind speeds
clc; close all
for s = 1:length(wind_speeds)

    % Get wind speed
    wind_speed = wind_speeds{s};
    caze = strcat(wind_speed, '_WV0_AGP');
    % freestream = freestreams.(wind_speed);
    disp(caze)
    
    % Get wave data
    x = waves.(caze).x;
    dx = mean(diff(x), 'all');
    test_waves = waves.(caze).wave_profiles;

    % Get friction velocity
    u_star = log_law.(caze).ensemble.u_star.u_star;

    % Get average wave
    time_mean_waves = mean(test_waves, 1, 'omitnan');

    % Get Wave fluctuations
    wave_fluctuations = test_waves - time_mean_waves;

    % Compute instantaneous wave slopes
    [eta_dx, ~] = gradient(wave_fluctuations, dx);

    % Comepute RMS of wave elevation [m]
    mean_square_elevation = mean(wave_fluctuations.^2, 'all', 'omitnan');
    rms_elevation = sqrt(mean_square_elevation) * 1E-3;

    % Compute RMS of wave slope across everything
    mean_square_slope = mean(eta_dx.^2, 'all', 'omitnan');
    rms_slope = sqrt(mean_square_slope);

    % Compute roughness reynolds number
    Re_rough = (rms_elevation * u_star) / nu;

    % Just check what the actual viscous sublayer height would be
    sublayer = ((5 * nu) / u_star) * 1E3;

    % Save values
    rms_elevations(s) = rms_elevation;
    rms_slopes(s) = rms_slope;
    Re_roughs(s) = Re_rough;
    sublayers(s) = sublayer;
    
    % Prints
    fprintf('RMS Elevation: %2.3f mm\n', rms_elevation * 1E3)
    fprintf('RMS Slope: %2.3f\n',rms_slope)
    fprintf('u* = %2.3f m/s\n', u_star)
    fprintf('Re_rough = %2.3f\n', Re_rough)
    fprintf('Sublayer = %3.2f mm\n\n', sublayer)

end

clear caze dx eta_dx freestream mean_square_elevation mean_square_slope rms_elevation rms_slope Re_rough s 
clear sublayer test_waves time_mean_waves u_star wave_fluctuations wind_speed x 



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING LOG LAW + ROUGHNESS REYNOLDS NUMBER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Plot details
linewidth = 2;
marker_size = 30;
wind_speed_markers = {'^', 'square', 'o'};

% Plot fontsizes
tickFontSize = 14;
labelFontSize = 16;
legendFontSize = 12;
titleFontSize = 18;

% Smooth wall reference
c_plus = 5;
karman_constant = 0.41;
smooth_y_plus = 0:1:1E4;
smooth_u_plus = (1/karman_constant) * log(smooth_y_plus) + c_plus;

% Make array to hold onto log-region bounds
no_wave_log_regions = nan(3,2);



% Plot
clc; close all;
ax = figure('color', 'white', 'units', 'centimeters', 'position', [2,4,30,10]);
tiledlayout(1, 3, 'Padding', 'tight')

h(1) = nexttile([1,2]);
hold on
% Smooth wall reference
plot(smooth_y_plus, smooth_u_plus, 'linestyle', '--', 'HandleVisibility', 'off', 'color', 'black', 'linewidth', 2)

% Loop through wind speeds
for s = 1:length(wind_speeds)

    % Get wind speed
    wind_speed = wind_speeds{s};
    u_inf = freestreams.(wind_speed);

    % No-Wave
    no_wave_caze = strcat(wind_speed, '_WV0_AGP');
    y = friction_profiles.(no_wave_caze).y;
    u_profile = friction_profiles.(no_wave_caze).u_profile;
    friction_velocity = friction_profiles.(no_wave_caze).raw;
   
    % Law of the wall
    u_plus = u_profile / friction_velocity;
    y_plus = (y * friction_velocity) / nu;

    %%% NEW METHOD
    % Find constant slope region (log region)
    [bounds, ~] = find_constant_slope_region(y_plus, u_plus, ...
                                                      'plateau_tol', 0.2, ...
                                                      'smooth_window', 7);

    % Save bounds
    no_wave_log_regions(s,:) = bounds;

    % Compute \Delta U^+
    in_region = (y_plus >= bounds(2)) & (y_plus <= bounds(1));
    y_fit = y_plus(in_region);
    u_fit = u_plus(in_region);
    
    % Smooth wall reference
    u_smooth = (1/karman_constant) * log(y_fit) + c_plus;
    
    % Delta U+
    deltaU = mean(u_smooth - u_fit, 'omitnan');

    % Prints
    disp(wind_speed)
    fprintf('Log Region: %3.1f - %3.1f\n', bounds(2), bounds(1))
    fprintf('Delta U = %1.3f\n\n', deltaU)

    % Plot
    spacing = 1;
    label = sprintf('$u_{\\infty} = %1.2f$ m/s, $\\Delta U^+ = %1.2f$',  u_inf, deltaU);

    % Plot lines
    plot(y_plus(1:spacing:end), u_plus(1:spacing:end) + deltaU, ...
         'linewidth', linewidth, 'displayname', 'No Wave', ...
         'color', wind_speed_colors{s}, 'HandleVisibility', 'off')

    % Plot points
    scatter(y_plus(1:spacing:end), u_plus(1:spacing:end) + deltaU, marker_size, wind_speed_markers{s}, 'filled', ...
        'MarkerFaceColor', wind_speed_colors{s}, ...
        'displayname', label)
end
% Plot average log-region
% l = xline(mean(no_wave_log_regions, 1, 'omitnan'), 'color', 'black', ...
%           'linewidth', 1, 'HandleVisibility', 'off');
% uistack(l, 'bottom')
hold off

grid on
xscale('log')
xlim([5, 5E3])
xlabel('$y^+$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$u^+ + \Delta U^+$', 'interpreter', 'latex', 'fontsize', labelFontSize)
legend('location', 'northwest', 'interpreter', 'latex', 'fontsize', legendFontSize, 'box', 'off')



% Plotting Roughness Reynolds number vs u_infty (or u*)
h(2) = nexttile;
hold on
for s = 1:length(wind_speeds)

    % Get wind speed
    wind_speed = wind_speeds{s};
    caze = strcat(wind_speed, '_WV0_AGP');
    freestream = freestreams.(wind_speed);
    disp(caze)

    % Get friction velocity
    u_star = log_law.(caze).ensemble.u_star.u_star;

    % Plot
    scatter(freestream, Re_roughs(s), 4 * marker_size, 'filled', wind_speed_markers{s}, ...
            'markerfacecolor', wind_speed_colors{s})
    tmpX(s) = freestream;
    tmpY(s) = Re_roughs(s);
end
p = plot(tmpX, tmpY, 'linewidth', linewidth, 'color', 'black');
uistack(p, 'bottom')
hold off

yline(5, 'linestyle', '--', 'color', 'black', ...
      'label', '$k_{s}^+ = 5$', ...
      'LabelHorizontalAlignment', 'Left', ...
      'interpreter', 'latex', ...
      'fontsize', 8)

ylabel('$k_{s}^+ = \frac{\eta_{rms} u^*}{\nu}$', 'Interpreter', 'latex', 'fontsize', labelFontSize)
xlabel("$u_{\infty}$ [m/s]", 'interpreter', 'latex', 'fontsize', labelFontSize)
ylim([0, 14])

clear bounds c_plus deltaU freestream friction_velocity in_region karman_constant l label no_wave_cazeu_star wind_speed
clear no_wave_log_regions p slope_info smooth_u_plus smooth_y_plus spacing u_fit u_plus u_profile u_smooth h wind s
clear y y_fit y_plus u_inf tmpX tmpY caze


%% No-waves Log Law












%%

% % Plotting \Delta U^+ against Roughness Reynolds number
% h(2) = nexttile;
% hold on
% for s = 1:3
%     speed = WTs{s};
%     caze = strcat('WT', num2str(speed), '_WV0_AGP');
%     freestream = freestreams.(['WT', num2str(speed)]);
%     disp(caze)
% 
%     % Get friction velocity
%     u_star = log_law.(caze).ensemble.u_star.u_star;
% 
%     % Get \Delta U^+
%     DeltaUPlus = DeltaUPluses.(['WT', num2str(speed)]);
% 
%     % Plot
%     scatter(Re_roughs(s), DeltaUPlus, 100, 'filled', 'o', ...
%             'markerfacecolor', wind_speed_colors{s})
%     tmpX(s) = Re_roughs(s);
%     tmpY(s) = DeltaUPlus;
% end
% p = plot(tmpX, tmpY, 'linewidth', 2, 'color', 'black');
% uistack(p, 'bottom')
% hold off
% xline(5, 'linestyle', '--', 'color', 'black', ...
%       'label', 'Viscous Sublayer $k_{s}^+ = 5$', ...
%       'LabelHorizontalAlignment', 'Left', ...
%       'interpreter', 'latex', ...
%       'fontsize', 8)
% 
% yline(0, 'linestyle', '--', 'color', 'black', ...
%       'label', '$\Delta U^+ = 0$', ...
%       'LabelHorizontalAlignment', 'Left', ...
%       'interpreter', 'latex', ...
%       'fontsize', 8)
% 
% xlabel('$k_{s}^+ = \frac{\eta_{rms} u^*}{\nu}$', 'Interpreter', 'latex')
% ylabel("$\Delta U^+$", 'interpreter', 'latex')




%% Testing: Plot u* againt freestream

figure('color', 'white')
hold on
for s = 1:3
    speed = WTs{s};
    caze = strcat('WT', num2str(speed), '_WV0_AGP');
    freestream = freestreams.(['WT', num2str(speed)]);
    disp(caze)

    % Get friction velocity
    u_star = log_law.(caze).ensemble.u_star.u_star;

    % Get \Delta U^+
    DeltaUPlus = DeltaUPluses.(['WT', num2str(speed)]);

    % Plot
    scatter(freestream, u_star, 100, 'filled', 'o', ...
            'markerfacecolor', wind_speed_colors{s})
    tmpX(s) = freestream;
    tmpY(s) = u_star;
end
p = plot(tmpX, tmpY, 'linewidth', 2, 'color', 'black');
uistack(p, 'bottom')
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

