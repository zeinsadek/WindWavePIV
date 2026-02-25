% Demo code computing the wave-coherent deviations from the 
% curvilinear-projected instantaneous and phase averages for
% Montreal Polytechnique


clc; close all; clear

% Get case to process
wind_speed = 'WT6';
wave_type = 'C';
caze = [wind_speed, '_WV', wave_type, '_AG0'];

matfile_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/curvilinear_montreal';
matfile_name   = [caze, '_CurvilinearInstantaneous_Polytechnique.mat'];

% Load data
output = load(fullfile(matfile_folder, matfile_name));
output = output.output;


%% Breakdown constants

% Wave constants
wavelength = output.constants.wavelength_mm;
amplitude  = output.constants.amplitude_mm;
frequency  = output.constants.frequency_hz;
steepness  = output.constants.steepness;
wave_speed = output.constants.wave_speed;

% Freestream velocity
u_inf = output.constants.freestream;


%% Break down curvililinear grid

curvilinear_grid = output.curvilinear_grid;

% curvilinear_grid(phase).vertical_lines ~ curvilinear 'x'
% curvilinear_grid(phase).horizontal_lines ~ curvilinear 'y'

% curvilinear_grid(phase).wave_profile ~ reference wave profiles
% curvilinear_grid(phase).max_wave_profile ~ max wave profile per phase bin


%% Load Cartesian coordinates

X = output.Cartesian_grid.X;
Y = output.Cartesian_grid.Y;
x = output.Cartesian_grid.x;
y = output.Cartesian_grid.y;


%% Pull out curvilinear means/instantaneous

curvilinear_means = output.curvilinear_means;
curvilinear_instantaneous = output.curvilinear_instantaneous;


%% Compute average over the 4 different phases

% Empty arrays to hold values
u_tmp  = nan([size(X), 4]);
v_tmp  = nan([size(X), 4]);
vl_tmp = nan([size(X), 4]);
hl_tmp = nan([size(X), 4]);

% How much to chop off edges [mm]
edge_buffer = 3;
top_crop = 5;

 % Loop through phases
for phase = 1:4

    % Load data
    u  = curvilinear_means(phase).u;
    v  = curvilinear_means(phase).v;
    vl = curvilinear_grid(phase).vertical_lines;
    hl = curvilinear_grid(phase).horizontal_lines;
    
    %%% Crop data/artifacts
    % Crop sides
    u(vl < min(X, [], 'all') + edge_buffer) = nan;
    u(vl > max(X, [], 'all') - edge_buffer) = nan;
    v(vl < min(X, [], 'all') + edge_buffer) = nan;
    v(vl > max(X, [], 'all') - edge_buffer) = nan;

    % Crop top
    u(hl > max(Y, [], 'all') - top_crop) = nan;
    v(hl > max(Y, [], 'all') - top_crop) = nan;

    % Crop near surface
    

    u_tmp(:,:,phase)  = u;
    v_tmp(:,:,phase)  = v;
    vl_tmp(:,:,phase) = vl;
    hl_tmp(:,:,phase) = hl;
end

% Compute mean of phase
phase_mean.u = mean(u_tmp, 3, 'omitnan');
phase_mean.v = mean(v_tmp, 3, 'omitnan');
phase_mean.vertical_lines = mean(vl_tmp, 3, 'omitnan');
phase_mean.horizontal_lines = mean(hl_tmp, 3, 'omitnan');

% Plot
figure('color', 'white')
contourf(phase_mean.vertical_lines, ...
         phase_mean.horizontal_lines, ...
         phase_mean.u, ...
         100, 'LineStyle', 'none')
axis equal

clear edge_crop top_crop phase u v vl hl u_tmp v_tmp vl_tmp hl_tmp


%% Compute wave-coherent velocities
% (wave coherent velocities) = (phase average) - (average of phases)

% Loop through phases
 for phase = 1:4
    wave_coherent(phase).u = curvilinear_means(phase).u - phase_mean.u;
    wave_coherent(phase).v = curvilinear_means(phase).v - phase_mean.v;
 end

 % Plot 
figure('color', 'white')
tiledlayout(1,4)

% Loop through phases
for phase = 1:4
    h(phase) = nexttile;
    hold on

    % Plot quantity
    contourf(curvilinear_grid(phase).vertical_lines, ...
             curvilinear_grid(phase).horizontal_lines, ...
             wave_coherent(phase).u, ...
             100, 'linestyle', 'none')

    % Plot reference wave profile
    plot(x, curvilinear_grid(phase).wave_profile, 'linewidth', 2, 'color', 'black')

    % Plot max wave profile
    plot(x, curvilinear_grid(phase).max_wave_profile, 'linewidth', 2, 'color', 'red')
    hold off
    axis equal
    clim([-0.5,0.5])
    title(sprintf('Phase %1.0f', phase))
    
end


