%% WindWave Testing with phase average tolerance

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions/');
wave_parameters = readcell('Offshore_Waves.xlsx');

% Mat file paths
project_path = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV';
wave_path = fullfile(project_path, 'wave');
data_path = fullfile(project_path, 'data');

%% Load data

caze = 'WT4_WVA_AG0';

% Load wave parameters
wave_type      = caze(strfind(caze, 'WV') + 2);
wave_length    = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 2};
wave_amplitude = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 3};
wave_frequency = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 4};
phase_offset   = [0, wave_length/4, wave_length/2, 3*wave_length/4];

%% Phase average with different tolerance

clc; close all;

% Specift phase tolerance as a function of wavelength
percent_phase_tolerance = 0.005;
phase_tolerance         =  percent_phase_tolerance * wave_length;

% Constants Structure
left_bound_value  = -121;
right_bound_value = 115;

constants.recording_name = caze;
constants.left_bound_value        = left_bound_value;
constants.right_bound_value       = right_bound_value;
constants.percent_phase_tolerance = percent_phase_tolerance;
constants.phase_tolerance         = phase_tolerance; 
constants.phase_offset            = phase_offset;
constants.wave_type          = wave_type;
constants.wave_amplitude     = wave_amplitude;
constants.wave_length        = wave_length;

phases = phaseaverage(constants, fullfile(wave_path, strcat(caze, '_WAVE.mat')), 'xxx');
close all;

%% Compute phase + ensemble averages

% Load data
data = load(fullfile(data_path, strcat(caze, '_DATA.mat')));
data = data.output;

% Load waves
waves = load(fullfile(wave_path, strcat(caze, '_WAVE.mat')));
waves = waves.output;

clc;
means = data2means2D('xxx', data, waves, phases, constants, data.D);

%% Plot phase average quantities

lw = 3;

close all;
figure()
tiledlayout(2,2,'TileSpacing','tight')
for p = 1:4
    nexttile
    hold on
    contourf(means.X, means.Y, means.phase(p).uv, 100, 'linestyle', 'none')
    plot(means.X(1,:), means.phase(p).max_wave_profile, ...
         'linewidth', lw, 'color', 'red')
    plot(means.X(1,:), means.phase(p).reference_wave, ...
         'LineWidth', lw, 'color', 'black')
    hold off
    axis equal
    colorbar()
    clim([-0.02, 0])
end
sgtitle(strcat(caze, ' uv phase average', {' '}, num2str(percent_phase_tolerance * 100), '% tol'), ...
        'interpreter', 'none')



