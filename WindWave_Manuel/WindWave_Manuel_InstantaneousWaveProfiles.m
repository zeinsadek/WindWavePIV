% Wrap up the curvilinear-grid Cartesian velocities instantaneous 
% and phase averages for Cerrina

clc; close all; clear

% Get case to process
wind_speed = 'WT4';
wave_type = 'C';
caze = [wind_speed, '_WV', wave_type, '_AG0'];

% Load all useful data
% data   = load(strcat('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/data/', caze, '_DATA.mat'));
waves  = load(strcat('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/wave/', caze, '_WAVE.mat'));
means  = load(strcat('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means/', caze, '_MEANS.mat'));
phases = load(strcat('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/phase/', caze, '_PHASE.mat'));

% Rename structure
% data   = data.output;
waves  = waves.output;
means  = means.output;
phases = phases.output;
clc; fprintf('Data loaded\n\n')

% Get wave parameters
wave_parameters = readcell("/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/Offshore_Waves.xlsx");
wavelength      = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 2}; %#ok<*FNDSB>
amplitude       = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 3};
frequency       = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 4};
phase_offset    = [0, wavelength/4, wavelength/2, 3*wavelength/4];
wavenumber      = (2 * pi) / wavelength;

% Get freestream
freestreams.('WT4') = 2.4181;
freestreams.('WT6') = 3.8709;
freestreams.('WT8') = 5.4289;

clear wave_type


%% Load relevent data

% Coordinates
x = waves.x;

% Bin instantaneous snapshots by phase
for phase = 1:4

    % Indicies of images (out of the 6000 frames) that belong to 'phase'
    indicies = find(phases.phase_average_idx == phase);
    fprintf('Phase %1.0f: %4.0f Images\n\n', phase, length(indicies))

    PhaseBinnedWaves(phase).waves = waves.wave_profiles(indicies,:); %#ok<*SAGROW>

    % Save reference wave profiles
    PhaseBinnedWaves(phase).reference = imresize(means.phase(phase).reference_wave.', [1, length(x)]);

end

clear phases wave_phases_indicies phase indicies


%% Plot one profile per phase

figure('color', 'white')
title(caze, 'interpreter', 'none')
hold on
for phase = 1:4
    plot(x, PhaseBinnedWaves(phase).waves(1,:), 'linewidth', 1)
    % plot(x, PhaseBinnedWaves(phase).reference, 'linewidth', 1, 'color', 'black')

    % Save profiles for manuel
    output(phase).InstantaneousProfile = PhaseBinnedWaves(phase).waves(1,:);
    output(phase).ReferenceProfile = PhaseBinnedWaves(phase).reference;
    output(phase).x = x;
end
hold off
axis equal

% save_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/manuel';
% file_name = 
% save()


