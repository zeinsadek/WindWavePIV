%%% Wind + Wave paper wave-deviations test

% WindWave paper figures: Curvilinear Phase Averages

clc; clear; close all
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions/Inpaint_nans/Inpaint_nans')
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps/slanCM')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
figure_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/pdf_test7';
curvilinear_path = fullfile(project_path, 'curvilinear_new');

% Cases
wind_speeds = {'WT4', 'WT6', 'WT8'};

% Cases manuel wants
waves = {'A', 'C'};
wave_colors = {'#FE6202', '#775EEF'};

% Freestream velocities
freestreams.('WT4') = 2.4181;
freestreams.('WT6') = 3.8709;
freestreams.('WT8') = 5.4289;

% Approximate wavelengths in mm for labeling plots
wavelength_names.A = '410';
wavelength_names.B = '313';
wavelength_names.C = '189';
wavelength_names.D = '124';

steepnesses.A = '0.180';
steepnesses.B = '0.211';
steepnesses.C = '0.305';
steepnesses.D = '0.267';

% Load data
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        disp(caze)
        tmp = load(fullfile(curvilinear_path, strcat(caze, '_CURVILINEAR.mat')));
        tmp = tmp.output;
        curvilinear.(caze) = tmp;
    end
end

clc; fprintf('All cases loaded\n')
clear caze tmp w no_wave_caze cartesian_path data_path
clear wind_speed wave s w caze tmp project_path


% Wave parameters
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


%% Compute u* per phase

idx = 86;

% Loop through phases
for phase = 1:4

    % Loop through wind speeds
    for s = 1:length(wind_speeds)
        wind_speed = wind_speeds{s};
        u_inf = freestreams.(wind_speed);
    
        % Loop through waves 
        for w = 1:length(waves)
            wave = waves{w};
            frequency = frequencies.(wave);
            wavelength = wavelengths.(wave);
            wave_speed = wavelength * frequency;
            caze = strcat(wind_speed, '_WV', wave, '_AG0');
    
            % Calculate u*
            uv_stress_profile = curvilinear.(caze).phase(phase).uv(:, idx);
            u_star = sqrt(max(-uv_stress_profile, [], 'omitnan'));
    
            % Save u*
            friction_velocities(phase).(caze) = u_star;
        end
    end
end


% Compute mean friction velocity per case
% Loop through wind speeds
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    u_inf = freestreams.(wind_speed);

    % Loop through waves 
    for w = 1:length(waves)
        wave = waves{w};
        frequency = frequencies.(wave);
        wavelength = wavelengths.(wave);
        wave_speed = wavelength * 1E-3 * frequency;
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        % Collect phase averaged u*
        tmp = nan(1,4);
        for phase = 1:4
            tmp(phase) = friction_velocities(phase).(caze);
        end

        % Average across phases
        mean_friction_velocities.(caze) = mean(tmp, 'all', 'omitnan');

        % Compute mean friction wave age to compare to paper
        friction_wave_ages.(caze) = wave_speed / mean_friction_velocities.(caze);

    end
end



%% Compute average of phase averages for each case

% How much to chop off edges [mm]
edge_buffer = 1;
top_buffer = 1;
bottom_buffer = 2;

% Loop through wind speeds
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};

    % Loop through waves
    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        disp(caze)

        X = curvilinear.(caze).X;
        Y = curvilinear.(caze).Y;

        u_tmp  = nan([size(X), 4]);
        v_tmp  = nan([size(X), 4]);
        hl_tmp = nan([size(X), 4]);
        vl_tmp = nan([size(X), 4]);

        % Loop through phases
        for phase = 1:4
            u  = curvilinear.(caze).phase(phase).u;
            v  = -1 * curvilinear.(caze).phase(phase).v;
            hl = curvilinear.(caze).phase(phase).horizontal_lines;
            vl = curvilinear.(caze).phase(phase).vertical_lines;
        
            % Crop the sides
            u(vl < min(X, [], 'all') + edge_buffer) = nan;
            u(vl > max(X, [], 'all') - edge_buffer) = nan;
            v(vl < min(X, [], 'all') + edge_buffer) = nan;
            v(vl > max(X, [], 'all') - edge_buffer) = nan;

            % Crop the top
            u(hl > max(Y, [], 'all') - top_buffer) = nan;            
            v(hl < max(Y, [], 'all') - top_buffer) = nan;

            % Trimp edges
            % X = curvilinear.(caze).X;
            % Y = curvilinear.(caze).Y;
            % tmp(vl < min(X, [], 'all') + edge_buffer) = nan;
            % tmp(vl > max(X, [], 'all') - edge_buffer) = nan;
            % tmp(hl > max(Y, [], 'all') + 2 * edge_buffer) = nan;
        
            % Trim first rows of real data right above surface
            is_all_nans = all(isnan(u), 2);
            is_not_all_nans = ~is_all_nans;
            last_non_nan_row_index = find(is_not_all_nans, 1, 'last');
            first_non_nan_row_index = find(is_not_all_nans, 1, 'first');
        
            u(last_non_nan_row_index - bottom_buffer:end, :) = nan;
            u(1:first_non_nan_row_index + bottom_buffer, :) = nan;

            % Save
            u_tmp(:,:,phase)  = u;
            v_tmp(:,:,phase)  = v;
            hl_tmp(:,:,phase) = hl;
            vl_tmp(:,:,phase) = vl;
        end

        % Compute mean of phase
        phase_mean.(caze).u = mean(u_tmp, 3, 'omitnan');
        phase_mean.(caze).v = mean(v_tmp, 3, 'omitnan');
        phase_mean.(caze).horizontal_lines = mean(hl_tmp, 3, 'omitnan');
        phase_mean.(caze).vertical_lines = mean(vl_tmp, 3, 'omitnan');

    end
end

clear edge_buffer s wind_speed w wave caze X Y 
clear phase u_tmp v_tmp hl_tmp vl_tmp u v hl vl



% Plot: Compare the mean of phases
wind_speed = 'WT6';
freestream = freestreams.(wind_speed);

cmin = 0.5;
cmax = 1;

% % Plot
% clc; close all
% figure('color', 'white')
% tiledlayout(1, length(waves))
% 
% % Loop through waves
% for w = 1:length(waves)
%     wave = waves{w};
%     caze = strcat(wind_speed, '_WV', wave, '_AG0');
%     disp(caze)
% 
%     h(w) = nexttile;
%     contourf(phase_mean.(caze).vertical_lines, ...
%              phase_mean.(caze).horizontal_lines, ...
%              phase_mean.(caze).u / freestream, ...
%             100, 'linestyle', 'none')
%     yline(0, 'linewidth', 2)
%     axis equal
%     clim([cmin, cmax])
%     label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses.(wave));
%     title(label, 'Interpreter', 'latex')
%     if w == 4
%         colorbar()
%     end
% end
% 
% linkaxes(h, 'xy')
% ylim([-5, 200])
% 
% clear wind_speed freestream cmin cmax w wave caze h 


%% Comparing log law + roughness functions using two different methods

LL_cases = {'WT4_WVA_AG0', 'WT6_WVC_AG0'};
nu = 1.46E-5;
karman_constant = 0.41;
c_plus = 5;

% Method we use in the paper (average \Delta U+)
for c = 1:length(LL_cases)
    caze = LL_cases{c};

    tmp = nan(1,4);
    for phase = 1:4
        % Get friction velocities
        friction_velocity = friction_velocities(phase).(caze);
    
        % Get profiles
        Uc = curvilinear.(caze).phase(phase).u;
        Yc = curvilinear.(caze).phase(phase).horizontal_lines;
        u_profile = Uc(:, idx);
        y_profile = Yc(:, idx) * 1E-3;

        u_plus = u_profile / friction_velocity;
        y_plus = (y_profile * friction_velocity) / nu;
    
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

        % Phase Delta U+
        phase_DeltaU(phase).(caze) = deltaU;
        tmp(phase) = deltaU;
    end

    % Mean delta U+
    mean_DeltaU.(caze) = mean(tmp, 'all', 'omitnan');

end



% Other method, using mean velocity profile (Manuel)
for c = 1:length(LL_cases)
    caze = LL_cases{c};

    friction_velocity = mean_friction_velocities.(caze);
    % Get mean profiles
    Xc = phase_mean.(caze).vertical_lines;
    Yc = phase_mean.(caze).horizontal_lines;
    Uc = phase_mean.(caze).u;

    u_profile = Uc(:, idx);
    y_profile = Yc(:, idx) * 1E-3;

    u_plus = u_profile / friction_velocity;
    y_plus = (y_profile * friction_velocity) / nu;

    % Compute \Delta U^+
    in_region = (y_plus >= bounds(2)) & (y_plus <= bounds(1));
    y_fit = y_plus(in_region);
    u_fit = u_plus(in_region);
    
    % Smooth wall reference
    u_smooth = (1/karman_constant) * log(y_fit) + c_plus;
    
    % Delta U+
    deltaU = mean(u_smooth - u_fit, 'omitnan');

    % Save 
    meanprofile_deltaU.(caze) = deltaU;
end


















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


















