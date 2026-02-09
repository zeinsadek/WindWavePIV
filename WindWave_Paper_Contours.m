%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions');
experiment_log  = readcell('Offshore_Inflow_Log.xlsx');
wave_parameters = readcell('Offshore_Waves.xlsx');

% Offsets and shifts
top_bound_value   = 205;       % relative to Y centered at still water
left_bound_value  = -121;      % relative to X centered at DaVis default
right_bound_value = 115;       % relative to X centered at DaVis default
range = abs(left_bound_value) + abs(right_bound_value);

% Cases with waves
WTs = {'4', '6', '8'};
WVs = {'A', 'B', 'C', 'D'};

% Paths
BL_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/boundary_layer';
LL_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/log_law/';
MN_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means/';
figure_path = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures';

% s = settings;
% s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load cases with waves
for s = 1:length(WTs)
    speed = WTs{s};
    for i = 1:length(WVs)
        % Case
        caze = strcat('WT', num2str(speed), '_WV', WVs{i}, '_AG0');
        fprintf('Loading Case: %s...\n', caze)
    
        % Paths
        BL_path = fullfile(BL_folder, strcat(caze, '_BL.mat'));
        LL_path = fullfile(LL_folder, strcat(caze, '_LL.mat'));
        MN_path = fullfile(MN_folder, strcat(caze, '_MEANS.mat'));
        
        % Save LL to Structure
        LL_temp = load(LL_path);
        LL_temp = LL_temp.output;
        LL_Waves.(caze) = LL_temp;

        % Save BL to Structure
        BL_temp = load(BL_path);
        BL_temp = BL_temp.output;
        BL_Waves.(caze) = BL_temp;
       
        % Save MN to Structure
        means_temp = load(MN_path);
        means_temp = means_temp.output;
        MN_Waves.(caze) = means_temp;
    end
    fprintf('\n')
end

% Load cases without waves
for s = 1:length(WTs)
    % Case
    no_wave_case = strcat('WT', num2str(WTs{s}), '_WV0_AGP');
    fprintf('Loading Case: %s...\n', no_wave_case)
    
    % Paths
    BL_path = fullfile(BL_folder, strcat(no_wave_case, '_BL.mat'));
    LL_path = fullfile(LL_folder, strcat(no_wave_case, '_LL.mat'));
    MN_path = fullfile(MN_folder, strcat(no_wave_case, '_MEANS.mat'));
    
    % Save LL to Structure
    LL_temp = load(LL_path);
    LL_temp = LL_temp.output;
    LL_No_Waves.(no_wave_case) = LL_temp;

    % Save BL to Structure
    BL_temp = load(BL_path);
    BL_temp = BL_temp.output;
    BL_No_Waves.(no_wave_case) = BL_temp;
    
    % Save MN to Structure
    means_temp = load(MN_path);
    means_temp = means_temp.output;
    MN_No_Waves.(no_wave_case) = means_temp;
end

% Get Case Names
wave_cases    = fields(MN_Waves);
no_wave_cases = fields(MN_No_Waves);

% Clear temporary variables
clear BL_path LL_path MN_path
clear BL_temp LL_temp means_temp
clear caze i no_wave_case s speed


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET FREESTREAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
for s = 1:length(WTs)
    
    c = 1;
    all_u_infs = zeros(1,5);
    WT = WTs{s};

    % Waves
    for w = 1:length(WVs)
        WV = WVs{w};
        caze = strcat('WT', WT, '_WV', WV, '_AG0');
        u_inf = BL_Waves.(caze).u_inf;
        all_u_infs(c) = u_inf;
        c = c + 1;
    end

    % No Waves
    caze = strcat('WT', WT, '_WV0_AGP');
    u_inf = BL_No_Waves.(caze).u_inf;
    all_u_infs(5) = u_inf;
    avg_u_inf = mean(all_u_infs, 'all');
    fprintf("WT %s = %2.2f m/s\n", WT, avg_u_inf);
    u_infs.(strcat("WT", WT)) = avg_u_inf;
end

% Clear temp variables
clear c s w all_u_infs WT WV caze u_inf  avg_u_inf


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE AVERAGE WAVE PROFILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

caze = "WT4_WVC_AG0";

ax = figure();
tiledlayout(4, 1, "TileSpacing", "tight", "Padding", "tight")

for p = 1:4
    
    % Load data
    waves = MN_Waves.(caze).phase(p).waves;
    reference = MN_Waves.(caze).phase(p).reference_wave;
    max_wave = MN_Waves.(caze).phase(p).max_wave_profile;
    x = MN_Waves.(caze).X(1,:);
    [num_waves, ~] = size(waves);

    % Plot
    nexttile()
    hold on

    % Plot all instantaneous wave profiles
    for i = 1:num_waves
        plt = plot(x, imresize(waves(i,:), size(x)), 'color', 'black', 'linewidth', 0.1, "HandleVisibility", "off");
        plt.Color(4) = 0.1;
    end

    % Plot reference wave profile
    plot(x, reference, 'color', 'red', 'linewidth', 3, "DisplayName", "Reference Wave")

    % Plot max wave profile
    plot(x, max_wave, 'color', 'green', 'linewidth', 3, "DisplayName", "Max Height")
    hold off
    
    axis equal
    xlim([0, range])
    ylim([-20, 20])

    if p == 2
        ylabel('y [mm]')
    end
    
    if p == 4
        xlabel('x [mm]')
    end
end

leg = legend('Orientation', 'Horizontal');
leg.Layout.Tile = 'north';
leg.Box = "off";

% exportgraphics(ax, fullfile(figure_path, strcat(caze, "_PhaseAverageWaves.png")), "Resolution", 300);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE AVERAGE CONTOURS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

caze = "WT4_WVC_AG0";

lw = 3;
X = MN_Waves.(caze).X;
Y = MN_Waves.(caze).Y;

ax = figure('Position', [200,200,1000,300]);
t = tiledlayout(1,4,'padding', 'tight');
for p = 1:4
    nexttile()
    hold on
    contourf(X, Y, MN_Waves.(caze).phase(p).u / u_infs.("WT4"), 100, 'linestyle', 'none');
    plot(unique(X), MN_Waves.(caze).phase(p).max_wave_profile, 'color', 'red', 'linewidth', lw)
    hold off
    axis equal
    ylim([min(Y, [], 'all'), 200])
    xlim([min(X, [], 'all') + 2, max(X, [], 'all') - 2])
    xlabel('x [mm]', 'interpreter', 'latex')
    if p == 1
        ylabel('y [mm]', 'interpreter', 'latex')
    end
    clim([0.3,1.0])
    tit = strcat("Phase", {' '}, num2str(p));
    title(tit{1})
    
    if p == 4
        c = colorbar();
        c.Label.String = "$u / u_{\infty}$";
        c.Label.Interpreter = "Latex";
        c.Label.FontSize = 14;
    end
end

% exportgraphics(ax, fullfile(figure_path, strcat(caze, "_phases.png")), resolution=300)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENSEMBLE SHEAR STRESS PROFILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;

WT = "4";
lw = 3;
pad = 12;
lower_limit = 18;
WV_colors = cool(length(WVs));
ax = figure('Position', [200,200,500,500]);

hold on

% Waves
for w = 1:4
    caze = "WT" + WT + "_WV" + WVs{w} + "_AG0";
    disp(caze)
    [~, cols] = size(MN_Waves.(caze).ensemble.uv);
    uv   = MN_Waves.(caze).ensemble.uv(:,floor(cols/2) - pad:floor(cols/2) + pad);
    uv   = mean(uv, 2, 'omitnan');
    uv   = uv / (u_infs.("WT" + WT)^2);
    y    = flipud(unique(MN_Waves.(caze).Y));
    % uv(y < lower_limit) = nan;
    plot(-uv, y, linewidth = lw, displayname = "WV" + WVs{w}, color = WV_colors(w,:))

end

% No Waves
caze = "WT" + WT + "_WV0_AGP";
disp(caze)
[~, cols] = size(MN_No_Waves.(caze).ensemble.uv);
uv   = MN_No_Waves.(caze).ensemble.uv(:,floor(cols/2) - pad:floor(cols/2) + pad);
uv   = mean(uv, 2, 'omitnan');
uv   = uv / (u_infs.("WT" + WT)^2);
y    = flipud(unique(MN_No_Waves.(caze).Y));
% uv(y < lower_limit) = nan;
plot(-uv, y, linewidth = lw, displayname = "WV0", color = 'black')

hold off
ylim([0, 200])
legend('location', 'northeast')
ylabel('y [mm]', 'interpreter', 'latex')
xlabel("$-{\overline{u'v'}}^2 / {u_{\infty}}^2$", 'Interpreter', 'latex', 'fontsize', 16)

% exportgraphics(ax, fullfile(figure_path, strcat("WT", WT, "_Ensemble_uv_Profiles.png")), resolution=300)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE SHEAR STRESS PROFILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;

WT = "4";

lw = 3;
pad = 8;
fontsize = 16;
WV_colors = cool(length(WVs));
ax = figure('Position', [200,200,1000,500]);
tiledlayout(1,4)

% Waves
for p = 1:4
    nexttile()
    hold on
    for w = 1:4
        caze = "WT" + WT + "_WV" + WVs{w} + "_AG0";
        disp(caze)

        [~, cols] = size(MN_Waves.(caze).phase(p).uv);
        
        % Take a slice in the cente of the frame and average in the
        % x-direction
        uv   = MN_Waves.(caze).phase(p).uv(:,floor(cols/2) - pad:floor(cols/2) + pad);
        uv   = mean(uv, 2, 'omitnan');

        % Normalize
        uv   = uv / (u_infs.("WT" + WT)^2);
        y    = flipud(unique(MN_Waves.(caze).Y));
        wavelength = wave_parameters{w + 1, 2};

        % Plot
        plot(-uv, y / wavelength, linewidth=lw, displayname="WV" + WVs{w}, color=WV_colors(w,:))
    end
    hold off
    title(strcat("Phase: ", num2str(p)))
    ylim([-0.1,1])
    xlim([-1E-5, 8E-3])

    if p == 1
        ylabel('$y / \lambda$ [mm]', 'interpreter', 'latex', 'fontsize', fontsize)
    end

        xlabel("$-{\overline{u'v'}}^2 / {u_{\infty}}^2$", 'Interpreter', 'latex', 'fontsize', fontsize)
end

leg = legend('Orientation', 'Horizontal');
leg.Layout.Tile = 'north';
leg.Box = "off";


% exportgraphics(ax, fullfile(figure_path, strcat("WT", WT, "_Phase_uv_Profiles.png")), resolution=300)






%% Colored Wavelength

A = 2;
x = linspace(-0.6,0.6,100);
wave = A * cos(2 * pi  * x);
LW = 4;

ax = figure('Position', [300,300,1200,200]);
hA = axes(ax);
hold on

for j = 1:3
    if j == 1
        colors = summer(length(x));
        os = 0;
    end
    if j == 2
        colors = autumn(length(x));
        os = 1;
    end
    if j == 3
        colors = winter(length(x));
        os = 2;
    end

    for i = 1:(length(x) - 1)
        plot(x(i:i+1), wave(i:i+1) + os, 'color', colors(i,:), 'linewidth', LW)
    end
end
xline(0, 'color', 'black', 'linestyle', '--')
hold off
set(hA, 'box', 'off', 'YTick',[])
set(get(hA, 'YAxis'), 'Visible', 'off');
xlim([-0.5, 0.5])
ylim([-5,5])
xlabel('$x / \lambda$', 'interpreter', 'Latex')

% exportgraphics(ax, fullfile(figure_path, "Wave_Path_Legend.png"), resolution=300)




























