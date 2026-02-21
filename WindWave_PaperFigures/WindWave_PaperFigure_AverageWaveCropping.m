% WindWave paper figures: Curvilinear Phase Averages

clc; clear; close all
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
curvilinear_path = fullfile(project_path, 'curvilinear');

% Figure Folder
figure_folder = 'pdf_test3';

% Cases
% wind_speed = 'WT6';
wind_speeds = {'WT4', 'WT6', 'WT8'};
waves = {'A', 'B', 'C', 'D'};

% Approximate wavelengths in mm for labeling plots
% wavelengths.A = '410';
% wavelengths.B = '313';
% wavelengths.C = '189';
% wavelengths.D = '124';
% 
% steepnesses.A = '0.180';
% steepnesses.B = '0.211';
% steepnesses.C = '0.305';
% steepnesses.D = '0.267';

% if ismember(wind_speed(end), {'4'})
%     u_inf = 2.4181;
% elseif ismember(wind_speed(end), {'6'})
%     u_inf = 3.8709;
% elseif ismember(wind_speed(end), {'8'})
%     u_inf = 5.4289;
% end


% Load all wave cases
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        disp(caze)
    
        % Load data
        tmp = load(fullfile(curvilinear_path, strcat(caze, '_CURVILINEAR.mat')));
        tmp = tmp.output;
    
        % Save data
        data.(caze) = tmp;
    end
end
cazes = fields(data);
clear caze tmp w wave no_wave_caze wind_speed


%%

caze = 'WT6_WVA_AG0';

phase = 3;
x = data.(caze).X(1,:);
ref_wave = data.(caze).phase(phase).wave_profile;
crop_wave = data.(caze).phase(phase).max_wave_profile;

avg_gap = mean(crop_wave - ref_wave, 'all', 'omitnan');

% figure('color', 'white')
% hold on
% plot(x, ref_wave, 'color', 'black')
% plot(x, crop_wave, 'color', 'red')
% plot(x, ref_wave + avg_gap, 'color', 'red', 'linestyle', '--')
% hold off
% axis equal



% Loop through all cases and find average wave gap
c = 1;
for c = 1:length(cazes)
    caze = cazes{c};
    for phase = 1:4
        ref_wave = data.(caze).phase(phase).wave_profile;
        crop_wave = data.(caze).phase(phase).max_wave_profile;
        avg_gap = mean(crop_wave - ref_wave, 'all', 'omitnan');

        tmp(c) = avg_gap;
        c = c + 1;
    end
end

avg_wave_gap = mean(tmp, 'all', 'omitnan');
clc;
fprintf('Average phase average vertical crop = %1.3f mm\n\n', avg_wave_gap)