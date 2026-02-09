%% Log Law Testing
% Zein Sadek

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions/');
wave_parameters = readcell('Offshore_Waves.xlsx');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data
experiment_name_waves = 'WT4_WVA_AG0';
experiment_name_no_waves = 'WT4_WV0_AGP';

% Save paths
results_path   = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/';
means_waves    = fullfile(results_path, 'means', strcat(experiment_name_waves, '_MEANS.mat'));
means_no_waves = fullfile(results_path, 'means', strcat(experiment_name_no_waves, '_MEANS.mat'));

% Load 
means_waves = load(means_waves);
means_waves = means_waves.output;

means_no_waves = load(means_no_waves);
means_no_waves = means_no_waves.output;

clear experiment_name_no_waves experiment_name_waves

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REPAIR CAT-SCRATCH IN DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = means_waves.X;
Y = means_waves.Y;
x = X(1,:);
y = Y(:,1);

components = {'u', 'v', 'uu', 'vv', 'uv'};

fprintf("Phase\n")
% Fixing ensemble average for case without waves
for c = 1:length(components)
    component = components{c};
    for p = 1:4

        % Monitor
        fprintf("%s: Phase %1.0f\n", component, p)

        % Load component
        tmp = means_waves.phase(p).(component);

        % Shift RHS
        LHS_values = tmp(:,86);
        RHS_values = tmp(:,87);
        differences = RHS_values - LHS_values;
        differences(isnan(differences)) = 0;
        corrected = tmp;
        corrected(:,87:end) = corrected(:,87:end) - differences;

        % Interpolate across strip
        [~, left] = min(abs(x - 115));
        [~, right] = min(abs(x - 120));
        corrected(:, left:right) = nan;
        corrected = fillmissing(corrected, 'makima', 2);
        corrected(Y < means_waves.phase(p).max_wave_profile) = nan;

        % Resave data
        means_waves.phase(p).(component) = corrected;
        clear corrected
    end
    fprintf("\n")
end


fprintf("Ensemble\n")
% Fixing phase average for case with waves
clc;
for c = 1:length(components)
    component = components{c};

    % Monitor
    fprintf("%s\n", component)

    % Load component
    tmp = means_no_waves.ensemble.(component);

    % Shift RHS
    LHS_values = tmp(:,86);
    RHS_values = tmp(:,87);
    differences = RHS_values - LHS_values;
    differences(isnan(differences)) = 0;
    corrected = tmp;
    corrected(:,87:end) = corrected(:,87:end) - differences;

    % Interpolate across strip
    [~, left] = min(abs(x - 115));
    [~, right] = min(abs(x - 120));
    corrected(:, left:right) = nan;
    corrected = fillmissing(corrected, 'makima', 2);
    corrected(Y < means_no_waves.ensemble.max_wave_profile) = nan;

    % Resave data
    means_no_waves.ensemble.(component) = corrected;
    fprintf("\n")
    clear corrected
end
clc;

clear LHS_values RHS_values differences left right
clear p c component tmp


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT PROFILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = 86;

phase = 1;

profile_wave = means_waves.phase(phase).u(:,idx);
profile_no_wave = means_no_waves.ensemble.u(:,idx);

figure()
hold on
plot(profile_no_wave, y)
plot(profile_wave, y)
hold off

figure()
hold on
plot(y * 1E-3, profile_no_wave)
plot(y * 1E-3, profile_wave)
hold off
xscale("log")




