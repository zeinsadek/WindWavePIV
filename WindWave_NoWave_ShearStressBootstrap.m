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
INST_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/data/';
MN_folder   = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means/';
% WV_folder   = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/wave/';
% LL_folder   = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/log_law_fixed/';

% Wind speeds
freestreams.('WT4') = 2.4181;
freestreams.('WT6') = 3.8709;
freestreams.('WT8') = 5.4289;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load cases with waves
% for s = 1:length(wind_speeds)
for s = 3
    wind_speed = wind_speeds{s};
   
    % Case
    caze = strcat(wind_speed, '_WV0_AGP');
    fprintf('Loading Case: %s...\n', caze)
    MN_path = fullfile(MN_folder, strcat(caze, '_MEANS.mat'));
    inst_path = fullfile(INST_folder, strcat(caze, '_DATA.mat'));
    % WV_path = fullfile(WV_folder, strcat(caze, '_WAVE.mat'));
    % LL_path = fullfile(LL_folder, strcat(caze, '_LL_Fixed.mat'));
   
    % Save MN to structure
    means_temp = load(MN_path);
    means_temp = means_temp.output;
    means.(caze) = means_temp;

    % Save instantaneous images
    inst_tmp = load(inst_path);
    inst_tmp = inst_tmp.output;
    inst.(caze) = inst_tmp;

end

clear caze inst_tmp means_temp s wind_speed


%% Loop through different numbers of images and compute u* from u'v' max profile
% Take both a single slice (in middle of plane) and spanwise average


% Which case to look at
caze = 'WT8_WV0_AGP';
profile_idx = 86;
X = means.(caze).X;
Y = means.(caze).Y;

% Full image u* value
uv_full = means.(caze).ensemble.uv;
uv_profile_full = uv_full(:, profile_idx);
u_star_profile_full = max(sqrt(-uv_profile_full), [], 'all', 'omitnan');
u_star_profile_average = max(sqrt(-mean(uv_full, 2, 'omitnan')), [], 'all', 'omitnan');

% Image batches to average over
n_min = 1400;
n_step = 200;
n_max = inst.(caze).D;
total_number_images = inst.(caze).D;



% How many times to repeat the random draw
B = 1000;                         
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

        % If you prefer subsampling (without replacement), use:
        % idx = randperm(total_number_images, N);

        % --- snapshots ---
        u = inst.(caze).U(:,:,idx);
        v = inst.(caze).V(:,:,idx);

        % --- compute uv ---
        mean_u = mean(u, 3, 'omitnan');
        mean_v = mean(v, 3, 'omitnan');

        fluc_u = u - mean_u;
        fluc_v = v - mean_v;

        uv = mean(fluc_u .* fluc_v, 3, 'omitnan');

        % --- extract u* from belly max ---
        uv_profile = uv(:, profile_idx);
        u_star_boot(b,i) = max(sqrt(-uv_profile), [], 'all', 'omitnan');
    end
    clc;
end


% Compute variance from converged region
kc = ceil(0.7 * total_number_images);  % e.g. last 30% of N values
isPlateau = Nvec >= kc;
u_pool = u_star_boot(:, isPlateau);          % B x (#plateau N values)
u_pool = u_pool(:);                          % vectorize
lo = prctile(u_pool, 2.5);
hi = prctile(u_pool, 97.5);

Xminus = (u_star_profile_full - lo)/u_star_profile_full;
Xplus  = (hi - u_star_profile_full)/u_star_profile_full;
Xsym   = max(Xminus, Xplus);

fprintf('u* = %.4f m/s, 95%% interval [%0.4f, %0.4f] => ±%.2f%% (conservative)\n', ...
        u_star_profile_full, lo, hi, 100*Xsym);

%% Plot

medN = median(u_star_boot, 1, 'omitnan');
loN  = prctile(u_star_boot, 2.5, 1);
hiN  = prctile(u_star_boot, 97.5, 1);

figure('color','white'); hold on
fill([Nvec fliplr(Nvec)], [loN fliplr(hiN)], 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
plot(Nvec, medN, 'r-', 'LineWidth', 2);
yline(u_star_profile_full, 'k--', 'LineWidth', 2, 'Alpha', 0.5);
xline(kc, 'k--', 'LineWidth', 1.5);
xlabel('Number of Images $N$', 'Interpreter','latex')
ylabel('Friction Velocity $u^*$ [m/s]', 'Interpreter','latex')
title(caze, 'Interpreter','none')
ylim([0, 0.3])

