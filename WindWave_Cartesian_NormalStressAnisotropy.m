%%% Wind + Wave paper: checking that u'u' > v'v' for cartesian data

% WindWave paper figures: Curvilinear Phase Averages

clc; clear; close all
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions/Inpaint_nans/Inpaint_nans')
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps/slanCM')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
figure_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/pdf_test6';
cartesian_path = fullfile(project_path, 'means');

% Cases
wind_speeds = {'WT4', 'WT6', 'WT8'};
waves = {'A', 'B', 'C', 'D'};

freestreams.('WT4') = 2.4181;
freestreams.('WT6') = 3.8709;
freestreams.('WT8') = 5.4289;


% Load data
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    for w = 1:length(waves)
        wave = waves{w};

        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        disp(caze)
        tmp = load(fullfile(cartesian_path, strcat(caze, '_MEANS.mat')));
        cartesian.(caze) = tmp.output;
    end
end

clear caze tmp w no_wave_caze cartesian_path curvilinear_path data_path
clear s  wave wind_speed


%% Loop through all wind+waves and phases and check if u'u' > v'v'

c = 1;
clc;
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    for w = 1:length(waves)
        wave = waves{w};

        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        disp(caze)

        data = cartesian.(caze);

        for phase = 1:4
            uu = data.phase(phase).uu;
            vv = data.phase(phase).vv;

            anisotropy_diff = uu - vv;
            anisotropy_ratio = uu ./ vv;
            percent_off = sum(anisotropy_diff < 0, 'all', 'omitnan') / numel(uu);
            % disp(percent_off)
            disp(mean(anisotropy_ratio, 'all', 'omitnan'))

            tmp(c) = mean(anisotropy_ratio, 'all', 'omitnan');
            c = c + 1;
        end
    end
    fprintf('\n')
end