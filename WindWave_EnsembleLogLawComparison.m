%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions/');

project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
means_path = fullfile(project_path, 'means');
curvilinear_path = fullfile(project_path, 'curvilinear');
wave_parameters = readcell("Offshore_Waves.xlsx");

top_bound_value   = 205;       % relative to Y centered at still water
left_bound_value  = -121;      % relative to X centered at DaVis default
right_bound_value = 115;       % relative to X centered at DaVis default
range = abs(left_bound_value) + abs(right_bound_value);

figure_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new';

% Approximate wavelengths in mm for labeling plots
wavelength_names.A = '410';
wavelength_names.B = '313';
wavelength_names.C = '189';
wavelength_names.D = '124';


%% Method used in curvilinear log-law code for cartesian no-wave


nu = 1.46E-5;

wind_speed = 'WT6';
caze = strcat(wind_speed, '_WV0_AGP');

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
curvilinear_u_profile = u(:, idx);

% CLEAN UP U PROFILE
% curvilinear_u_profile = StatisticalGradientFilter(curvilinear_u_profile, 2.4, 10);
curvilinear_u_profile(find(~isnan(curvilinear_u_profile), 1, 'first')) = nan;
curvilinear_u_star = max(sqrt(-uv_profile), [], 'all', 'omitnan');
disp(curvilinear_u_star)

% Log law
curvilinear_y_plus = (y_profile * curvilinear_u_star) / nu;
curvilinear_u_plus = curvilinear_u_profile / curvilinear_u_star;

% Plot
clc;
figure()
plot(curvilinear_y_plus, curvilinear_u_plus)
xscale('log')



%% Older method used in cartesian log law


LL_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/log_law/';
wind_speed = 'WT6';

% Case
no_wave_case = strcat(wind_speed, '_WV0_AGP');
LL_path = fullfile(LL_folder, strcat(no_wave_case, '_LL.mat'));
cartesian_LogLaw = load(LL_path);
cartesian_LogLaw = cartesian_LogLaw.output;

cartesian_LogLaw_u_profile = cartesian_LogLaw.ensemble.u_profile;
cartesian_LogLaw_u_star = cartesian_LogLaw.ensemble.u_star.u_star;
cartesian_LogLaw_u_plus = cartesian_LogLaw.ensemble.u_plus;
cartesian_LogLaw_y_plus = cartesian_LogLaw.ensemble.y_plus;
disp(cartesian_LogLaw_u_star)


figure()
plot(cartesian_LogLaw_y_plus, cartesian_LogLaw_u_plus)
xscale('log')

%% Plot u profiles together

close all; clc
figure()
hold on
plot(curvilinear_u_profile, 'color', 'black')
plot(cartesian_LogLaw_u_profile, 'color', 'red')
hold off
xscale('log')

%% Plot log laws together

close all; clc
figure()
hold on
plot(curvilinear_y_plus, curvilinear_u_plus, 'color', 'black')
plot(cartesian_LogLaw_y_plus, cartesian_LogLaw_u_plus, 'color', 'red')
% plot((y_profile * cartesian_LowLaw_u_star) / nu, cartesian_LogLaw_u_plus, 'color', 'red')
hold off
xscale('log')












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


