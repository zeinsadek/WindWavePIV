%% 2D Offshore PIV Processing: Ondrej Fercak, Zein Sadek, 1/2023
% Converts DaVis vector data (.vc7) files into a Matlab file

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Functions/');
experiment_log  = readcell('Offshore_Inflow_Log.xlsx');
wave_parameters = readcell('Offshore_Waves.xlsx');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data paths
project_path   = '/Volumes/WT4_WVA_AG0/WT4_WVA_AG0';
recording_name = 'WT4_WVA_AG0';
processing     = '/PIV_MP(2x24x24_50%ov)/PostProc';
still_name     = strcat(experiment_log{find(strcmp(experiment_log, recording_name) == 1), 2}, '_waterlevel');
inpt_name      = recording_name;
% inpt_name      = strcat(recording_name, '_test2');

% Experiment Specifics
tunnel_freq    = recording_name(strfind(recording_name, 'WT') + 2);
wave_type      = recording_name(strfind(recording_name, 'WV') + 2);
% active_grid    = recording_name(strfind(recording_name, 'AG') + 2: strfind(recording_name, 'AG') + 3);

% Image paths
still_path     = strcat('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Offshore/Water_Level/', still_name, '/Correction');
raw_path       = strcat(project_path, '/', recording_name, '/Correction');
piv_path       = strcat(project_path, '/', recording_name, processing);

% Save paths
results_path = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/';
mtlb_file    = strcat(results_path, 'data'   , '/', inpt_name, '_DATA.mat');
% mean_file    = strcat(results_path, 'means'  , '/', inpt_name, '_MEANS.mat');
% phase_file   = strcat(results_path, 'phase'  , '/', inpt_name, '_PHASE.mat');
wave_file    = strcat(results_path, 'wave'   , '/', inpt_name, '_WAVE.mat');
figure_file  = strcat(results_path, 'figures', '/', recording_name);

% Make specific folder for figures of an experiment
if ~exist(figure_file, 'dir')
    mkdir(figure_file)
end

% Edge Detection Constants
std_tol            = 3;
background_removal = 1000;
canny_lower        = 0.1;
canny_upper        = 0.4;
grad_tol           = 0.4;
nan_dist           = 5;

% Filtering Parameters
top_bound_value   = 205;       % relative to Y centered at still water
left_bound_value  = -121;      % relative to X centered at DaVis default
right_bound_value = 115;       % relative to X centered at DaVis default
wave_buffer       = 0.5;
cmap_perc         = 0.05;
% phase_tolerance    = 10;
% phase_tolerance    = 33;


% Wave Parameters
if wave_type == '0'
    wave_length    = 0;
    wave_amplitude = 0;
    phase_offset   = [0, wave_length/4, wave_length/2, 3*wave_length/4];
else
    wave_length    = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 2};
    wave_amplitude = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 3};
    wave_frequency = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 4};
    phase_offset   = [0, wave_length/4, wave_length/2, 3*wave_length/4];
end

% Specift phase tolerance as a function of wavelength
percent_phase_tolerance = 0.05;
phase_tolerance         =  percent_phase_tolerance * wave_length;

% Video Settings
num_frames = 50;
fps        = 1;

% Constants Structure
constants.recording_name     = inpt_name;
constants.figure_file        = figure_file;
constants.num_frames         = num_frames;
constants.fps                = fps;
constants.cmap_perc          = cmap_perc;

constants.background_removal = background_removal;
constants.std_tol            = std_tol;
constants.canny_lower        = canny_lower;
constants.canny_upper        = canny_upper;
constants.grad_tol           = grad_tol;
constants.nan_dist           = nan_dist;

constants.left_bound_value        = left_bound_value;
constants.right_bound_value       = right_bound_value;
constants.top_bound_value         = top_bound_value;
constants.percent_phase_tolerance = percent_phase_tolerance;
constants.phase_tolerance         = phase_tolerance; 
constants.phase_offset            = phase_offset;

constants.wave_type          = wave_type;
constants.wave_amplitude     = wave_amplitude;
constants.wave_length        = wave_length;
constants.vertical_offset    = waterlevel(still_path, still_name, constants);
constants.wave_buffer        = wave_buffer;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIDE BY SIDE COMPARISON OF TERMS USING DIFFERENT TOL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

percents = [1, 2, 5, 10];


for i = 1:length(percents)
    stress_name = strcat('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means/WT4_WVA_AG0_', num2str(percents(i)), '_percent_tol_MEANS.mat');
    temp    = load(stress_name);
    temp    = temp.output;
    data(i) = temp;
end

%%

range = abs(left_bound_value) + abs(right_bound_value);
components = fields(data(1).ensemble);

% Reference Wave Function
cos_fit = @(b, v) wave_amplitude * cos(2 * pi * (v - (range/2) - b(1)) / wave_length);

t = tiledlayout(1, 4, 'Padding', 'compact', 'TileSpacing', 'none');
% title(t, 'uv stress')

for i = 1:length(percents)
    
    u_inf = mean(data(1).ensemble.u(:, 1:3), 'all', 'omitnan');
    
    X = data(i).X;
    Y = data(i).Y;
    x = unique(X);
    

    term  = data(i).phase(1).uv;
    waves = data(i).phase(1).waves;
    
    % Crop below wave
    ref = cos_fit(phase_offset(1), x);
    term(Y < ref.') = nan;
    
    num_waves = size(waves);
    num_waves = num_waves(1);
% 
%     nexttile
%     hold on
%     contourf(X, Y, term/(u_inf^2) , 500, 'LineStyle', 'none')
%     plot(x, ref, 'Color', 'red', 'LineWidth', 2)
%     hold off
%     axis equal
%     xlim([0, range])
%     ylim([-20, top_bound_value])
%     title(strcat('+/- ', num2str(percents(i)), '% of Wavelength:', {' '}, num2str(num_waves), ' Images'))
%     caxis([-5E-3, 5E-3])
    %colorbar()
    
    nexttile(i)
    num_waves = size(waves);
    num_waves = num_waves(1);
    
    hold on
    for j = 1:num_waves
        p = plot(x, imresize(waves(j,:), [1, length(x)]) - constants.vertical_offset, 'Color', 'black');
        p.Color(4) = 0.1;
    end
%     plot(x, ref, 'Color', 'red', 'LineWidth', 2)
    hold off
    axis equal
    xlim([0, range])
    ylim([-2.5 * wave_amplitude, 2.5 * wave_amplitude])
end

% cb = colorbar;
% cb.Layout.Tile = 'east';















