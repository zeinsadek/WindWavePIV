%%% Wind + Wave paper curvilinear coordinates figure for demonstration

% WindWave paper figures: Curvilinear Phase Averages

clc; clear; close all
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions/Inpaint_nans/Inpaint_nans')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
fig_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/pdf_test6';
curvilinear_path = fullfile(project_path, 'curvilinear_new');
cartesian_path = fullfile(project_path, 'means');
data_path = fullfile(project_path, 'data');

% Cases
wind_speed = 'WT6';
wave = 'D';
caze = strcat(wind_speed, '_WV', wave, '_AG0');

% Wind speeds
freestreams.('WT4') = 2.4181;
freestreams.('WT6') = 3.8709;
freestreams.('WT8') = 5.4289;
u_inf = freestreams.(wind_speed);

% Load data
instantaneous = load(fullfile(data_path, strcat(caze, '_DATA.mat')));
cartesian_mean = load(fullfile(cartesian_path, strcat(caze, '_MEANS.mat')));
cartesian_mean = cartesian_mean.output;
curvilinear_mean = load(fullfile(curvilinear_path, strcat(caze, '_CURVILINEAR.mat')));
curvilinear_mean = curvilinear_mean.output;


clc; fprintf('All %s cases loaded\n', wind_speed)
clear caze tmp w no_wave_caze cartesian_path curvilinear_path data_path


%% Set up curvilinear information

% Plot settings
phase = 1;
profile_idx = 86;

% Cartesian Coordinates
X = cartesian_mean.X;
x = X(1,:);
Y = cartesian_mean.Y;
y = Y(:,1);

% Curvilinear Coordinates
vertical_lines = curvilinear_mean.phase(phase).vertical_lines;
horizontal_lines = curvilinear_mean.phase(phase).horizontal_lines;
max_wave_profile = curvilinear_mean.phase(phase).max_wave_profile;
reference_wave = curvilinear_mean.phase(phase).wave_profile;

% Compute used u* value (cartesian)
curvilinear_uv_profile = curvilinear_mean.phase(phase).uv(:, profile_idx);
curvilinear_star = max(sqrt(-curvilinear_uv_profile), [], 'all', 'omitnan');
clc; fprintf('Full Curvilinear u* (phase %1.0f) = %1.3f m/s\n\n', phase, curvilinear_star)



% Get angles for curvilinear transformation 
% Angles
dx = mean(diff(x), 'all');
dy = mean(diff(y), 'all');
[dz_dx, dz_dy] = gradient(horizontal_lines, dx, dy);
dz_dx(horizontal_lines < reference_wave) = nan;
dz_dy(horizontal_lines < reference_wave) = nan;
angles = atan2(dz_dx, dz_dy);
angles(horizontal_lines < reference_wave) = nan;

% Compute unit vectors
% Tangent
xi_hat_x = cos(angles);    
xi_hat_y = sin(angles);

% Normal
zeta_hat_x = -sin(angles);   
zeta_hat_y =  cos(angles);


%% Bootstrapping

% Image batches to average over
phase_image_indicies = cartesian_mean.phase(phase).idxs;
n_min = 750;
n_step = 50;
n_max = length(phase_image_indicies);
total_number_images = length(phase_image_indicies);

% How many times to repeat the random draw
B = 500;                         
Nvec = n_min:n_step:n_max;       

% Bootstrap estimates for each N
u_star_boot = nan(B, numel(Nvec)); 

rng(1,'twister')  
for i = 1:numel(Nvec)
    N = Nvec(i);
    disp(N)

    for b = 1:B
        % --- choose indices ---
        % Bootstrap (with replacement):
        disp(b)
        idx = randi(total_number_images, [1 N]);

        % Load data
        u = instantaneous.output.U(:,:,phase_image_indicies(idx));
        v = instantaneous.output.V(:,:,phase_image_indicies(idx));

        % --- compute uv ---
        mean_u = mean(u, 3, 'omitnan');
        mean_v = mean(v, 3, 'omitnan');
        fluc_u = u - mean_u;
        fluc_v = v - mean_v;

        uv = mean(fluc_u .* fluc_v, 3, 'omitnan');
        uu = mean(fluc_u .^2, 3, 'omitnan');
        vv = mean(fluc_v .^2, 3, 'omitnan');




        % --- convert into Curvilinear ---
        % Fill NaNs in original data just in case
        uu(isnan(uu)) = 0;
        vv(isnan(vv)) = 0;
        uv(isnan(uv)) = 0;
        
        % Interpolate
        uu_interp = interp2(X, Y, uu, vertical_lines, horizontal_lines, 'cubic', nan);
        vv_interp = interp2(X, Y, vv, vertical_lines, horizontal_lines, 'cubic', nan);
        uv_interp = interp2(X, Y, uv, vertical_lines, horizontal_lines, 'cubic', nan);
        
        % Return nans
        uu_interp(uu_interp == 0) = nan;
        vv_interp(vv_interp == 0) = nan;
        uv_interp(vv_interp == 0) = nan;
        
        % Crop below wave
        uu_interp(horizontal_lines < max_wave_profile) = nan;
        vv_interp(horizontal_lines < max_wave_profile) = nan;
        uv_interp(horizontal_lines < max_wave_profile) = nan;

        % Projected shear stress (symmetric in turbulence, so either direction is fine)
        curvilinear_uv = uu_interp .* xi_hat_x .* zeta_hat_x + uv_interp .* (xi_hat_x .* zeta_hat_y + xi_hat_y .* zeta_hat_x) + vv_interp .* xi_hat_y .* zeta_hat_y;




        % --- extract u* from belly max ---
        uv_profile = curvilinear_uv(:, profile_idx);

        % Remove probelmatic values
        uv_profile = -1 * uv_profile;
        uv_profile(uv_profile < 0) = nan;

        % Crop a little near the surface
        buffer = 2;
        % i_nan = find(isnan(uv_profile), 1, 'last');   % last NaN in the vector
        i_surface = find(~isnan(uv_profile(1:end-1)) & isnan(uv_profile(2:end)), 1, 'last');
        uv_profile = uv_profile(1:(i_surface - buffer));

        % Smooth profile
        u_star_boot(b,i) = real(max(sqrt(uv_profile), [], 'all', 'omitnan'));
    end
    clc;
end


% Compute variance from converged region
kc = ceil(0.7 * total_number_images);        % e.g. last 30% of N values
isPlateau = Nvec >= kc;
u_pool = u_star_boot(:, isPlateau);          % B x (#plateau N values)
u_pool = u_pool(:);                          % vectorize
lo = prctile(u_pool, 2.5);
hi = prctile(u_pool, 97.5);

% Compare variation to mean of cenverged region
curvilinear_star = mean(u_pool, 'all');

Xminus = (curvilinear_star - lo)/curvilinear_star;
Xplus  = (hi - curvilinear_star)/curvilinear_star;
Xsym   = max(Xminus, Xplus);

fprintf('u* = %.4f m/s, 95%% interval [%0.4f, %0.4f] => ±%.2f%% (conservative)\n', ...
        curvilinear_star, lo, hi, 100*Xsym);
%%
% Histrogram
figure; histogram(u_pool, 60); xline(curvilinear_star,'k--');

% Plot
medN = median(u_star_boot, 1, 'omitnan');
loN  = prctile(u_star_boot, 2.5, 1);
hiN  = prctile(u_star_boot, 97.5, 1);

figure('color','white'); hold on
fill([Nvec fliplr(Nvec)], [loN fliplr(hiN)], 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
plot(Nvec, medN, 'r-', 'LineWidth', 2);
yline(curvilinear_star, 'r--', 'LineWidth', 2, 'Alpha', 0.5);
xline(kc, 'k--', 'LineWidth', 1.5);
xlabel('Number of Images $N$', 'Interpreter','latex')
ylabel('Friction Velocity $u^*$ [m/s]', 'Interpreter','latex')
% title(caze, 'Interpreter','none')
ylim([0, 0.5])