%%% Plotting log law profiles w/out ensemble average data with waves
% Zein Sadek

clear; close all; clc;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions/');
LL_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/log_law";
Means_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means";

%% Specify windspeed and import all different wave cases
% WV0, WVA, WVB, WVC, WVD

WT = "8";

% Load no wave case
LL_no_waves = load(fullfile(LL_path, strcat("WT", WT, "_WV0_AGP_LL.mat")));
LL_no_waves = LL_no_waves.output;

% Load waves cases
waves = {"A", "B", "C", "D"};
for w = 1:length(waves)
    wave = waves{w};
    tmp = load(fullfile(LL_path, strcat("WT", WT, "_WV", wave, "_AG0_LL.mat")));
    LL_waves.(wave) = tmp.output;
end

clear tmp w wave 

%% Plot log law

linewidth = 3;

figure()
tiledlayout(1,4)

for p = 1:4
    ax(p) = nexttile;
    hold on
    % Plot LL no waves, ensemble
    plot(LL_no_waves.ensemble.y_plus, LL_no_waves.ensemble.u_plus, ...
        'LineWidth', linewidth, 'DisplayName', 'WV0')
    
    % Plot LL w/ waves, phase
    for w = 1:length(waves)
        wave = waves{w};
        plot(LL_waves.(wave).phase(p).y_plus, LL_waves.(wave).phase(p).u_plus, ...
            'LineWidth', linewidth, 'DisplayName', strcat("WV", wave))
    end
    hold off
    xscale log
    xlim([1, 5E3])
    title(strcat("Phase ", num2str(p)))
    
end
leg = legend('Orientation', 'Horizontal');
leg.Layout.Tile = 'north';
linkaxes(ax, 'xy')

%% Plot uv profiles that give us u*

% Get coordinates
means_no_waves = load(fullfile(Means_path, strcat("WT", WT, "_WV0_AGP_MEANS.mat")));
means_no_waves = means_no_waves.output;

y = means_no_waves.Y(:,1);

idx = LL_no_waves.profile_idx;

% Plot
linewidth = 3;

figure()
tiledlayout(1,4)

for p = 1:4
    ax(p) = nexttile;
    hold on
    % Plot LL no waves, ensemble
    plot(-means_no_waves.ensemble.uv(:,idx), y, ...
        'LineWidth', linewidth, 'DisplayName', 'WV0')
    
    % Plot LL w/ waves, phase
    for w = 1:length(waves)
        wave = waves{w};
        plot(LL_waves.(wave).phase(p).uv_profile, y, ...
            'LineWidth', linewidth, 'DisplayName', strcat("WV", wave))
    end
    hold off
    % xscale log
    ylim([-10, 100])
    title(strcat("Phase ", num2str(p)))
    
end
leg = legend('Orientation', 'Horizontal');
leg.Layout.Tile = 'north';
linkaxes(ax, 'xy')


%%

% up = LL_no_waves.ensemble.u_plus;
% yp = LL_no_waves.ensemble.y_plus;

up = LL_waves.A.phase(1).u_plus;
yp = LL_waves.A.phase(1).y_plus;
 
func = 1 * log(yp) + 12;

figure();
hold on
plot(yp, up, 'linewidth', 3)
plot(yp, func)
hold off
xscale log


%% shifting phase profiles by virtual origin
nu = 1.48E-5;   

clc; close all;
figure()
hold on
% no waves

u_plus_no_waves = LL_no_waves.ensemble.u_plus;
y_plus_no_waves = LL_no_waves.ensemble.y_plus;

plot(y_plus_no_waves, u_plus_no_waves, ...
        'LineWidth', linewidth, 'DisplayName', 'WV0')
% w/ waves
for p = 1:4
    profile = LL_waves.(wave).phase(p).u_profile;
    u_star = LL_waves.(wave).phase(p).u_star;
    virtual_origin = min(y(~isnan(profile)));

    y_plus = ((y - virtual_origin) * 1E-3) * (u_star / nu);
    u_plus = profile / u_star;

    delta_u = max(u_plus_no_waves) - max(u_plus);

    plot(y_plus, u_plus + delta_u)
end
hold off
xscale log

%% Comapre friction velocities

phase_positions = {"Phase 1: Peak", "Phase 2: Ascending", "Phase 3: Trough", "Phase 4: Descending"};
wave = 'A';

figure()
hold on
for p = 1:4
    u_star = LL_waves.(wave).phase(p).u_star;
    bar(phase_positions{p}, u_star)
end
bar("No Waves", LL_no_waves.ensemble.u_star.u_star);
hold off
ylabel("Friction Velocity")
title(strcat("WT", num2str(WT), " Wave: ", wave))




