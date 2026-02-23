clc; close all;

caze = 'WT6_WVC_AG0';
data = load(strcat('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/data/', caze, '_DATA.mat'));
waves = load(strcat('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/wave/', caze, '_WAVE.mat'));

data = data.output;
waves = waves.output;

clc; fprintf('Data loaded\n\n')

%% Save a subset of frames and waves

num_frames = 1000;

clear output
output.X = data.X;
output.Y = data.Y;
output.D = num_frames;

for f = 1:num_frames
    output.U(:,:,f) = data.U(:,:,f);
    output.V(:,:,f) = data.V(:,:,f);
    output.waves(f, :) = imresize(waves.wave_profiles(f,:), [1, size(data.X, 2)]);
    output.phases(f) = waves.fitted_phases(f);
    clear f
end

clear num_frames

%% Test plot

frame = 5;

figure('color', 'white')
hold on
contourf(output.X, output.Y, output.U(:,:,frame), 10, 'linestyle', 'none')
plot(output.X(1,:), output.waves(frame,:), 'color', 'red', 'linewidth', 2)
hold off

axis equal

clear frame

%% Save package

save_dir = '/Users/zeinsadek/Downloads';
save(fullfile(save_dir, strcat(caze, '_PolyMTL_Instantaneous.mat')), 'output');

%% Test open

clear; clc;

caze = 'WT6_WVC_AG0';
save_dir = '/Users/zeinsadek/Downloads';
test = load(fullfile(save_dir, strcat(caze, '_PolyMTL_Instantaneous.mat')));
test = test.output;


frame = 5;

figure('color', 'white')
hold on
contourf(test.X, test.Y, test.U(:,:,frame), 10, 'linestyle', 'none')
plot(test.X(1,:), test.waves(frame,:), 'color', 'red', 'linewidth', 2)
hold off

axis equal

clear frame