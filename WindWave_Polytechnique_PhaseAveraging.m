%% Wind-Wave PIV Phase Averaging for Polytechnique Montreal
% Zein Sadek, Portland State University
% sadek@pdx.edu
% 11/2/2025


%% Paths and names for all imported files

clc; clear; close all;

% Name of the matfile being opened
caze = 'WT6_WVD_AG0_PolyMTL_Instantaneous.mat';

% Name of the folder where the matfile is located
caze_folder = '/Users/zeinsadek/Downloads';

% Path to wave spreadsheet to automate getting wave parameters
spreadsheet_path = '/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/Offshore_Waves.xlsx';

% Path to the phase averaging function
function_folder = '/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions';
addpath(function_folder)

% Clear variables we wont nead anymore
clear wave_parameters function_folder


%% Get wave parameters

% Open spreadsheet
wave_parameters = readcell(spreadsheet_path);

% Get wave type from caze name
wave_type = caze(7);

% Read spreadsheet to get wave parameters
wave_row = find(strcmp(wave_parameters, wave_type) == 1);
wavelength = wave_parameters{wave_row, 2};
amplitude = wave_parameters{wave_row, 3};
frequency = wave_parameters{wave_row, 4};

clc;
fprintf('Wave %s\nWavelength = %3.2f mm\nAmplitude = %3.2f mm\nFrequency = %1.2f Hz\n\n', ...
        wave_type, wavelength, amplitude, frequency);

% Clear variables we wont nead anymore
clear wave_parameters spreadsheet_path frequency wave_row


%% Load PIV data

% Import matfile
data = load(fullfile(caze_folder, caze));
data = data.output;

% Load coordinates
X = data.X;
Y = data.Y;
x = X(1,:);
y = Y(:,1);

% Load instantaneous velocities
u = data.U;
v = data.V;

% Load waves
waves = data.waves;

% Load phase offsets fitted to instantaneous wave profiles
phases = data.phases;


%% Plot a single snapshot to check

% Which snapshot to plot
frame = 1;

% Plot a single snapshot to check
figure('color', 'white')
tiledlayout(1,2)
sgtitle(sprintf('Frame = %1.0f', frame))

% Plot streamwise velocity
h(1) = nexttile;
hold on
contourf(X, Y, u(:,:,frame), 50, 'linestyle', 'none')

% Plot instantaneous wave surface 
plot(x, waves(frame, :), 'linewidth', 2, 'color', 'red')
hold off
axis equal
xlabel('x [mm]')
ylabel('y [mm]')
title('Streamwise Velocity')


% Plot vertical velocity
h(2) = nexttile;
hold on
contourf(X, Y, v(:,:,frame), 50, 'linestyle', 'none')

% Plot instantaneous wave surface 
plot(x, waves(frame, :), 'linewidth', 2, 'color', 'red')
hold off
axis equal
xlabel('x [mm]')
ylabel('y [mm]')
title('Vertical Velocity')

% Make both tiles be the same size
linkaxes(h, 'xy')

% Clear variables we wont nead anymore
clear h frame


%% Phase-Averaging binning

% Constant inputs to phase-averaging function: no need to change these
constants.recording_name = caze;
constants.amplitude = amplitude;
constants.wavelength = wavelength;

% These are the phases that we will be binning into. You can explore
% changing these to see how the velocity field changes based on the
% 'frozen' wave position. They should span from 0 up to 1 full wavelength
phase_offsets = [0, wavelength/4, wavelength/2, 3 * wavelength/4];
constants.phase_offset = phase_offsets;

% This is the tolerance of the binning, expressed as percentage of the
% wavelength. You can play with this to see exapnding the tolerance will
% give you more images per phase, leading to better statistics but the
% waves will be more spread out leading to less 'crisp' wave alignment
constants.percent_phase_tolerance = 0.1;

% This will run the function and output an array with values 1 through the
% number of phases you are binning into. We will use this to select the
% corresponding PIV images and compute the phase average. Entries of 0 are
% images that did not fit into any of the phases
phaseaverage_indicies = phaseaverage_Polytechnique(constants, data);


%% Phase-Average quantities

% This part will be looped to phase average overall of the different phase
% offsets, and will be saved to a single structure. First we collect all
% the instantaneous images for u and v, per phase
for phase = 1:length(phase_offsets)
    indicies = find(phaseaverage_indicies == phase);
    
    % Grab the instantaneous images that belong to a specific phase
    instantaneous_u = u(:,:,indicies);
    instantaneous_v = v(:,:,indicies);
    instantaneous_waves = waves(indicies, :);

    % Now we will compute the phase-averaged u and v velocities
    phaseaverage_u = mean(instantaneous_u, 3, 'omitnan');
    phaseaverage_v = mean(instantaneous_v, 3, 'omitnan');
    max_wave_profile = max(instantaneous_waves, [], 1, 'omitnan');

    % We need to crop below the local, maximum wave height
    phaseaverage_u(Y < max_wave_profile) = nan;
    phaseaverage_v(Y < max_wave_profile) = nan;


    % Compute Reynolds stresses be getting the velocity fluctuations
    u_fluctuations = instantaneous_u - phaseaverage_u;
    v_fluctuations = instantaneous_v - phaseaverage_v;

    % Covariance
    phaseaverage_uu = mean(u_fluctuations .* u_fluctuations, 3, 'omitnan');
    phaseaverage_vv = mean(v_fluctuations .* v_fluctuations, 3, 'omitnan');
    phaseaverage_uv = mean(u_fluctuations .* v_fluctuations, 3, 'omitnan');

    % We need to crop below the local, maximum wave height
    phaseaverage_uu(Y < max_wave_profile) = nan;
    phaseaverage_vv(Y < max_wave_profile) = nan;
    phaseaverage_uv(Y < max_wave_profile) = nan;


    % Save to array for later use
    % Velocities
    phaseaverage_mean(phase).u = phaseaverage_u;
    phaseaverage_mean(phase).v = phaseaverage_v;

    % Stresses
    phaseaverage_mean(phase).uu = phaseaverage_uu;
    phaseaverage_mean(phase).vv = phaseaverage_vv;
    phaseaverage_mean(phase).uv = phaseaverage_uv;
    phaseaverage_mean(phase).max_wave = max_wave_profile;

    % Number of images used
    phaseaverage_mean(phase).D = length(indicies);

    % Clear variables we wont nead anymore
    clear phase indicies phaseaverage_u phaseaverage_v instantaneous_u instantaneous_v
    clear phaseaverage_uu phaseaverage_vv phaseaverage_uv
end


%% Plotting the phase-averaged quantites

% Which component to look at: u, v, uu, vv, uv
component = 'uv';

figure('color', 'white')
t = tiledlayout(1, length(phase_offsets));
for phase = 1:length(phase_offsets)

    h(phase) = nexttile;
    hold on
    contourf(X, Y, phaseaverage_mean(phase).(component), 50, 'linestyle', 'none')
    plot(x, phaseaverage_mean(phase).max_wave, 'color', 'red', 'linewidth', 2)
    hold off
    colorbar()
    axis equal

    % Clear variables we wont nead anymore
    clear phase
end

xlabel(t, 'x [mm]')
ylabel(t, 'y [mm]')
linkaxes(h, 'xy')

% Handy bit of code to automatically sync all the color bar ranges
ax = findall(gcf, 'Type', 'axes');
clims = cell2mat(get(ax, 'CLim'));
global_clim = [min(clims(:,1)), max(clims(:,2))];
set(ax, 'CLim', global_clim);

% Clear variables we wont nead anymore
clear t h





