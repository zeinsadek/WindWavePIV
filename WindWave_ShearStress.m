%%% Plotting log law profiles w/out ensemble average data with waves
% Zein Sadek

clear; close all; clc;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions/');
wave_parameters = readcell('Offshore_Waves.xlsx');
fprintf('ALl paths loaded\n\n')

%% Load means

% Inputs
WT = "8";
waves = {"A", "B", "C", "D"};
Means_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means";

% Paths
MN_no_waves = load(fullfile(Means_path, strcat("WT", WT, "_WV0_AGP_MEANS.mat")));
MN_no_waves = MN_no_waves.output;

% Load waves cases
for w = 1:length(waves)
    wave = waves{w};
    tmp = load(fullfile(Means_path, strcat("WT", WT, "_WV", wave, "_AG0_MEANS.mat")));
    MN_waves.(wave) = tmp.output;
end

clear tmp w wave Means_path

%% Coordinates

X = MN_no_waves.X;
Y = MN_no_waves.Y;

x = X(1,:);
y = Y(:,1);

%%

phase = 4;
lw = 3;

clc; close all;
figure()
tiledlayout(1,4)
for w = 1:length(waves)
    wave = waves{w};
    wavelength = wave_parameters{w + 1, 2};
    
    X_norm = (X - (range(x) / 2)) / wavelength;
    data = MN_waves.(wave).phase(phase).uv;

    data(data < -0.18) = nan;

    nexttile
    hold on
    contourf(X_norm, Y / wavelength, data, 100, 'LineStyle', 'none')
    % plot(X_norm(1,:), MN_waves.(wave).phase(phase).max_wave_profile / wavelength, ...
         % 'color', 'red', 'linewidth', lw)
    plot(X_norm(1,:), MN_waves.(wave).phase(phase).reference_wave / wavelength, ...
         'color', 'black', 'linewidth', lw)
    hold off
    axis equal
    xlim([-0.25, 0.5])
    ylim([-0.1, 0.7])
    colorbar
    clim([-0.2, 0])
    title(strcat(wave, "\lambda = ", num2str(wavelength)))
end

clear w wave X_norm wavelength data

%% Plot phase contours with uv profile plotted on top

wave = 'B';
phase = 4;
idx = 86;
scale = 80;
lw = 2;

phase_4_data = MN_waves.(wave).phase(phase).uv;

figure()
hold on
% Contour
contourf(X, Y, phase_4_data, 100, 'linestyle', 'none');

% Wave crop
plot(x, MN_waves.(wave).phase(phase).max_wave_profile, ...
     'color', 'red', 'linewidth', lw)

% Reference Wave
plot(x, MN_waves.(wave).phase(phase).reference_wave, ...
     'color', 'black', 'linewidth', lw)

for i = 10:20:length(phase_4_data)
    % slice = imgaussfilt(phase_4_data(:,i), 0);
    slice = phase_4_data(:,i);
    offset = min(y(~isnan(slice)));
    xloc = X(1,i);
    disp(i)
    xline(xloc)
    plot(scale * slice + xloc, y, ...
         'color', 'k', 'linewidth', lw)
end
hold off
axis equal
ylim([ -15, 100])
xlim([0, 235])
title(wave)

%% Investigate phase 4 UV stress to see when the wake plume begins

wave = "A";
phase = 4;
lw = 3;

phase_4_data = MN_waves.(wave).phase(4).uv;
colors = parula(length(phase_4_data));

figure()
hold on
for i = 5:10:length(phase_4_data)
    % slice = imgaussfilt(phase_4_data(:,i), 0);
    slice = phase_4_data(:,i);
    offset = min(y(~isnan(slice)));
    xloc = X(1,i);
    disp(i)

    plot(slice + 0.01 * i, y - offset, ...
         'displayname', num2str(xloc), 'color', colors(i,:), 'linewidth', lw)
    scatter(slice + 0.01 * i, y - offset, 40, 'filled', ...
         'markerfacecolor', colors(i,:), 'linewidth', lw)
end
hold off
ylim([ -10, 100])
% legend()



%% Plot phase contours with uv profile plotted on top + tracking plume

close all;
wave = 'D';
phase = 1;
idx = 86;
scale = 80;
lw = 2;

phase_4_data = MN_waves.(wave).phase(phase).uv;

figure()
hold on
% Contour
contourf(X, Y, phase_4_data, 100, 'linestyle', 'none');

% Wave crop
plot(x, MN_waves.(wave).phase(phase).max_wave_profile, ...
     'color', 'red', 'linewidth', lw)

% Reference Wave
plot(x, MN_waves.(wave).phase(phase).reference_wave, ...
     'color', 'black', 'linewidth', lw)

for i = 10:5:length(phase_4_data)
    slice = phase_4_data(:,i);
    [~, idx] = max(abs(slice));
    offset = min(y(~isnan(slice)));
    xloc = X(1,i);
    % xline(xloc)
    % plot(scale * slice + xloc, y, ...
    %      'color', 'k', 'linewidth', lw)
    scatter(xloc, y(idx), 50, 'filled', 'MarkerFaceColor', 'magenta')
end
hold off
axis equal
ylim([ -15, 100])
xlim([0, 235])
title(wave)

%% Plot phase contours with u profile plotted on top

wave = 'A';
phase = 4;
idx = 86;
scale = 2;
lw = 2;

phase_4_data = MN_waves.(wave).phase(4).u;

figure()
hold on
% Contour
contourf(X, Y, phase_4_data, 100, 'linestyle', 'none');

% Wave crop
plot(x, MN_waves.(wave).phase(phase).max_wave_profile, ...
     'color', 'red', 'linewidth', lw)

% Reference Wave
plot(x, MN_waves.(wave).phase(phase).reference_wave, ...
     'color', 'black', 'linewidth', lw)

for i = 10:20:length(phase_4_data)
    % slice = imgaussfilt(phase_4_data(:,i), 0);
    slice = phase_4_data(:,i);
    offset = min(y(~isnan(slice)));
    xloc = X(1,i);
    disp(i)
    xline(xloc)
    plot(scale * slice + xloc, y, ...
         'color', 'k', 'linewidth', lw)
end
hold off
axis equal
ylim([ -15, 100])
xlim([0, 235])
title(wave)  



%% Investigate phase 4 UV stress to see when the wake plume begins

wave = "A";
phase = 4;
lw = 3;

phase_4_data = MN_waves.(wave).phase(4).u;
colors = parula(length(phase_4_data));

figure()
hold on
for i = 5:20:length(phase_4_data)
    % slice = imgaussfilt(phase_4_data(:,i), 0);
    slice = phase_4_data(:,i);
    offset = min(y(~isnan(slice)));
    xloc = X(1,i);
    disp(i)

    plot(slice, y, ...
         'displayname', num2str(xloc), 'color', colors(i,:), 'linewidth', lw)
    scatter(slice, y, 40, 'filled', ...
         'markerfacecolor', colors(i,:), 'linewidth', lw)
end
hold off
ylim([ -10, 200])
% legend()