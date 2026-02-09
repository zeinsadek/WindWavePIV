%%% Plotting log law profiles w/out ensemble average data with waves
% Zein Sadek

clear; close all; clc;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions/');

%% Specify windspeed and import all different wave cases
% WV0, WVA, WVB, WVC, WVD

% Inputs
WT = "8";
waves = {"A", "B", "C", "D"};
LL_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/log_law";
Means_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means";

% Load no wave case
LL_no_waves = load(fullfile(LL_path, strcat("WT", WT, "_WV0_AGP_LL.mat")));
LL_no_waves = LL_no_waves.output;

% Load waves cases
for w = 1:length(waves)
    wave = waves{w};
    tmp = load(fullfile(LL_path, strcat("WT", WT, "_WV", wave, "_AG0_LL.mat")));
    LL_waves.(wave) = tmp.output;
end

clear tmp w wave LL_path

%% Load means

MN_no_waves = load(fullfile(Means_path, strcat("WT", WT, "_WV0_AGP_MEANS.mat")));
MN_no_waves = MN_no_waves.output;

% Load waves cases
for w = 1:length(waves)
    wave = waves{w};
    tmp = load(fullfile(Means_path, strcat("WT", WT, "_WV", wave, "_AG0_MEANS.mat")));
    MN_waves.(wave) = tmp.output;
end

clear tmp w wave Means_path

%% Plot phase contours with u profile plotted on top

X = MN_no_waves.X;
Y = MN_no_waves.Y;

x = X(1,:);
y = Y(:,1);

wave = 'A';
idx = 86;
lw = 3;
scale = 1;

figure()
tiledlayout(1,4)
for p = 1:4
    nexttile
    hold on
    % 2D contourf data
    contourf(X, Y, MN_waves.(wave).phase(p).u, 100, 'LineStyle', 'none')
    
    % Velocity profile at middle
    xline(X(1,idx))
    plot(X(1, idx) + MN_waves.(wave).phase(p).u(:,idx) * scale, y, ...
         'color', 'black', 'linewidth', 2)

    % Wave crop
    plot(x, MN_waves.(wave).phase(p).max_wave_profile, ...
         'color', 'red', 'linewidth', lw)

    % Reference Wave
    plot(x, MN_waves.(wave).phase(p).reference_wave, ...
         'color', 'black', 'linewidth', lw)

    hold off
    title(strcat("Phase: ", num2str(p)))
    axis equal
    xlim([0, 235])
    ylim([-20, 200])
end  


%% No wave u contour

figure()
hold on
contourf(X, Y, MN_no_waves.ensemble.u, 100, 'LineStyle', 'none')
plot(X(1,idx) + MN_no_waves.ensemble.u(:,idx), y, 'color', 'black', 'LineWidth', lw)
hold off
axis equal
xlim([0, 235])
ylim([-20, 200])


%% Plot u profiles: as is

wave = 'A';
idx = 86;
lw = 3;
buffer = 5;

figure()
hold on

% No waves
plot(mean(MN_no_waves.ensemble.u(:,idx - buffer:idx + buffer), 2, 'omitnan'), y, ...
     'LineWidth', lw, 'DisplayName', "No Waves")

for p = 1:4
    plot(mean(MN_waves.(wave).phase(p).u(:,idx - buffer:idx + buffer), 2, 'omitnan'), y, ...
         'LineWidth', lw, 'DisplayName', strcat("Phase: ", num2str(p)))
end

hold off
legend('location', 'northwest')
ylim([-20, 200])

%% Plot u profiles: shifted

wave = 'A';
idx = 86;
lw = 3;
buffer = 5;

figure()
hold on

% No waves
plot(mean(MN_no_waves.ensemble.u(:,idx - buffer:idx + buffer), 2, 'omitnan'), y, ...
     'LineWidth', lw, 'DisplayName', "No Waves")

for p = 1:4
    profile = mean(MN_waves.(wave).phase(p).u(:,idx - buffer:idx + buffer), 2, 'omitnan');
    virtual_origin = min(y(~isnan(profile)));
    plot(profile, y - virtual_origin, ...
         'LineWidth', lw, 'DisplayName', strcat("Phase: ", num2str(p)))
end

hold off
legend('location', 'northwest')
ylim([-20, 200])


%% Plot phase contours with uv profile plotted on top

wave = 'A';
idx = 86;
lw = 3;
scale = 100;

figure()
tiledlayout(1,4)
for p = 1:4
    nexttile
    hold on
    % 2D contourf data
    contourf(X, Y, MN_waves.(wave).phase(p).uv, 100, 'LineStyle', 'none')
    
    % Velocity profile at middle
    xline(X(1,idx))
    plot(X(1, idx) + MN_waves.(wave).phase(p).uv(:,idx) * scale, y, ...
         'color', 'black', 'linewidth', 2)

    % Wave crop
    plot(x, MN_waves.(wave).phase(p).max_wave_profile, ...
         'color', 'red', 'linewidth', lw)

    % Reference Wave
    plot(x, MN_waves.(wave).phase(p).reference_wave, ...
         'color', 'black', 'linewidth', lw)

    hold off
    title(strcat("Phase: ", num2str(p)))
    axis equal
    xlim([0, 235])
    ylim([-20, 200])
end  

%% No wave u contour

figure()
hold on
contourf(X, Y, MN_no_waves.ensemble.uv, 100, 'LineStyle', 'none')
plot(X(1,idx) + scale * MN_no_waves.ensemble.uv(:,idx), y, 'color', 'black', 'LineWidth', lw)
hold off
axis equal
xlim([0, 235])
ylim([-20, 200])


%% Plot uv profiles: averaging over small slice

wave = 'A';
idx = 86;
lw = 3;
buffer = 10;

figure()
hold on

% No waves
plot(mean(MN_no_waves.ensemble.uv(:,idx - buffer:idx + buffer), 2, 'omitnan'), y, ...
     'LineWidth', lw, 'DisplayName', "No Waves")

for p = 1:4
    plot(mean(MN_waves.(wave).phase(p).uv(:,idx - buffer:idx + buffer), 2, 'omitnan'), y, ...
         'LineWidth', lw, 'DisplayName', strcat("Phase: ", num2str(p)))
end

hold off
legend('location', 'northwest')
ylim([-20, 200])

%% Plot uv profiles: shifted

wave = 'A';
idx = 86;
lw = 3;
buffer = 5;

figure()
hold on

% No waves
plot(mean(MN_no_waves.ensemble.uv(:,idx - buffer:idx + buffer), 2, 'omitnan'), y, ...
     'LineWidth', lw, 'DisplayName', "No Waves")

for p = 1:4
    profile = mean(MN_waves.(wave).phase(p).uv(:,idx - buffer:idx + buffer), 2, 'omitnan');
    virtual_origin = min(y(~isnan(profile)));
    plot(profile, y - virtual_origin, ...
         'LineWidth', lw, 'DisplayName', strcat("Phase: ", num2str(p)))
end

hold off
legend('location', 'northwest')
ylim([-20, 200])


%% Plotting log law: Aligning at the bottom
% Recalculate y+ with shifts

nu = 1.48E-5;  

X = MN_no_waves.X;
Y = MN_no_waves.Y;

x = X(1,:);
y = Y(:,1) * 1E-3;

figure()
hold on
markersize = 40;
lw = 3;

% no waves
no_waves_profile = LL_no_waves.ensemble.u_profile;
no_waves_u_star = LL_no_waves.ensemble.u_star.u_star;
no_waves_y_shift = min(y(~isnan(no_waves_profile)));

no_waves_u_plus = no_waves_profile / no_waves_u_star;
no_waves_y_plus = (y - no_waves_y_shift) * (no_waves_u_star / nu);

plot(no_waves_y_plus, no_waves_u_plus, 'linewidth', lw, 'color', 'black', 'DisplayName', 'No Waves')
scatter(no_waves_y_plus, no_waves_u_plus, markersize, 'filled', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'black', ...
        'HandleVisibility', 'off')

colors = {"#52AA5E", "#FFC15E", "#FF9F1C", "#BCB6FF"};

% Phase average
for p = 1:4

    profile = LL_waves.(wave).phase(p).u_profile;
    u_star = LL_waves.(wave).phase(p).u_star;
    y_shift = min(y(~isnan(profile)));

    u_plus = profile / u_star;
    y_plus = (y - y_shift) * (u_star / nu);

    u_shift = max(no_waves_u_plus) - max(u_plus);

    plot(y_plus, u_plus + u_shift, 'linewidth', lw, 'color', colors{p}, 'DisplayName', strcat("Phase: ", num2str(p)))
    scatter(y_plus, u_plus + u_shift, markersize, 'filled', ...
            'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors{p}, ...
            'HandleVisibility', 'off')
end

hold off
xscale log
legend('location', 'northwest', 'fontsize', 18)
xlabel("$ \left( y + \epsilon \right)^+$", 'Interpreter', 'latex', 'fontsize', 18)
ylabel("$ u^+ + B$", 'Interpreter', 'latex', 'FontSize', 18)
title(strcat("WT", num2str(WT), " WV ", wave))

%% Investigate phase 4 UV stress to see when the wake plume begins

wave = "C";
phase = 4;

phase_4_data = MN_waves.(wave).phase(4).uv;

figure()
hold on
for i = 5:10:length(phase_4_data)
    % slice = imgaussfilt(phase_4_data(:,i), 0);
    slice = phase_4_data(:,i);
    offset = min(y(~isnan(slice)));
    xloc = X(1,i);
    disp(i)
    plot(slice + 0.01 * i, y - offset, 'displayname', num2str(xloc))
end
hold off
ylim([-0.01, 0.1])
legend()

%%

figure()
contourf(X, Y, phase_4_data, 100, 'LineStyle', 'none')
axis equal
xline(X(1,115))
xline(X(1,86))

%%

clc; close all;
figure()
hold on
contourf(X, Y, phase_4_data, 100, 'LineStyle', 'none')

scale = 100;

for i = 5:10:length(phase_4_data)
    slice = imgaussfilt(phase_4_data(:,i), 0.1);
    disp(i)
    % slice = phase_4_data(:,i);
    xloc = X(1,i);
    xline(xloc)
    plot(scale * slice + xloc, y * 1E3, 'linewidth', 2, 'color', 'black')
end
plot(x, MN_waves.(wave).phase(phase).max_wave_profile, ...
     'color', 'red', 'linewidth', 5)
hold off
axis equal
xlim([0, 230])
ylim([-20, 200])


%%

clc; close all;
test = phase_4_data;
test(test > -0.11) = nan;

figure();
hold on
contourf(X, Y, test, 100, 'LineStyle', 'none');
plot(x, MN_waves.(wave).phase(phase).max_wave_profile, ...
     'color', 'red', 'linewidth', 5)

plot(x, MN_waves.(wave).phase(phase).reference_wave, ...
     'color', 'k', 'linewidth', 3)
hold off
axis equal




