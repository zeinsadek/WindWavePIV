%%% Wind + Wave paper curvilinear coordinates figure for demonstration

% WindWave paper figures: Curvilinear Phase Averages

clc; clear; close all
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions/Inpaint_nans/Inpaint_nans')
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps/slanCM')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
figure_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/pdf_test7';
curvilinear_path = fullfile(project_path, 'curvilinear_new');
cartesian_path = fullfile(project_path, 'means');

% Cases
wind_speeds = {'WT4', 'WT6', 'WT8'};
waves = {'A', 'B', 'C', 'D'};
wave_colors = {'#FE6202', '#DC2680', '#775EEF', '#648FFF'};

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

        cartesian_tmp = load(fullfile(cartesian_path, strcat(caze, '_MEANS.mat')));
        curvilinear_tmp = load(fullfile(curvilinear_path, strcat(caze, '_CURVILINEAR.mat')));

        cartesian.(caze) = cartesian_tmp.output;
        curvilinear.(caze) = curvilinear_tmp.output;

        clear cartesian_tmp curvilinear_tmp
    end
end


% Also load no-wave data
clc;
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    caze = strcat(wind_speed, '_WV0_AGP');
    disp(caze)

    cartesian_tmp = load(fullfile(cartesian_path, strcat(caze, '_MEANS.mat')));
    cartesian.(caze) = cartesian_tmp.output;
end


clear caze tmp w no_wave_caze cartesian_path curvilinear_path data_path
clear s wind_speed w wave


% Approximate wavelengths in mm for labeling plots
wavelength_names.A = '410';
wavelength_names.B = '313';
wavelength_names.C = '189';
wavelength_names.D = '124';

steepnesses.A = '0.180';
steepnesses.B = '0.211';
steepnesses.C = '0.305';
steepnesses.D = '0.267';

% Wave parameters
wavelengths.A = 410.8 * 1E-3;
wavelengths.B = 313.3 * 1E-3;
wavelengths.C = 189.6 * 1E-3;
wavelengths.D = 124.3 * 1E-3;

amplitudes.A = 11.78 * 1E-3;
amplitudes.B = 10.53 * 1E-3;
amplitudes.C = 9.21 * 1E-3;
amplitudes.D = 5.29 * 1E-3;

frequencies.A = 1.96;
frequencies.B = 2.27;
frequencies.C = 3;
frequencies.D = 3.93;

% wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};


%% Plot no wave uv profiles

clc; close all
figure('color', 'white')
tiledlayout(1,3)
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    caze = strcat(wind_speed, '_WV0_AGP');
    u_inf = freestreams.(wind_speed);

    x_avg_uv = mean(-cartesian.(caze).ensemble.uv / (u_inf^2), 2, 'omitnan');

    h(s) = nexttile;
    hold on
    for i = 1:5:171
        plot(-cartesian.(caze).ensemble.uv(:,i) / (u_inf^2), cartesian.(caze).Y(:,1), 'color', 'black')
    end
    plot(x_avg_uv, cartesian.(caze).Y(:,1), 'color', 'red', 'linewidth', 3)
    hold off
    axis square

    % Save no-wave max shear stress
    no_wave_max_stress.(wind_speed) = max(mean(-cartesian.(caze).ensemble.uv, 2, 'omitnan'), [], 'all', 'omitnan');

    ustar = sqrt(max(mean(-cartesian.(caze).ensemble.uv, 2, 'omitnan'), [], 'all', 'omitnan'));
    disp(ustar)

end

linkaxes(h, 'xy')
ylim([0, 60])

clear h x_avg_uv

%% Compare shear profiles


idx = 86;

figure('color', 'white')
tiledlayout(2,4)

% Cartesian
for phase = 1:4
    h(phase) = nexttile;
    hold on
    for s = 1:length(wind_speeds)
        wind_speed = wind_speeds{s};
        u_inf = freestreams.(wind_speed);

        % Plot wave profiles
        for w = 1:length(waves)
            wave = waves{w};

            caze = strcat(wind_speed, '_WV', wave, '_AG0');
            disp(caze)

            local_wave_surface = cartesian.(caze).phase(phase).reference_wave(idx);
            plot(-cartesian.(caze).phase(phase).uv(:, idx) / (u_inf^2), cartesian.(caze).Y(:,1) - local_wave_surface, ...
                 'color', wave_colors{w})
            title(sprintf('Cartesian Phase %1.0f', phase))
        end

        % Plot no wave profiles
        caze = strcat(wind_speed, '_WV0_AGP');
        plot(-cartesian.(caze).ensemble.uv(:, idx) / (u_inf^2), cartesian.(caze).Y(:,1), ...
                 'color', 'black')
    end
    hold off
end


% Curvilinear
for phase = 1:4
    h(phase + 4) = nexttile;
    hold on
    for s = 1:length(wind_speeds)
        wind_speed = wind_speeds{s};
        u_inf = freestreams.(wind_speed);

        % Plot wave cases
        for w = 1:length(waves)
            wave = waves{w};

            caze = strcat(wind_speed, '_WV', wave, '_AG0');
            disp(caze)

            buffer = 2;
            profile = -curvilinear.(caze).phase(phase).uv(:, idx) / (u_inf^2);
            y_profile = curvilinear.(caze).phase(phase).horizontal_lines(:, idx);
            y_profile = y_profile(1:end - buffer);
            profile = profile(1:end - buffer);

            plot(profile, y_profile, 'color', wave_colors{w})
            title(sprintf('Curvilinear Phase %1.0f', phase))
        end

        % Plot no wave profiles
        caze = strcat(wind_speed, '_WV0_AGP');
        plot(-cartesian.(caze).ensemble.uv(:, idx) / (u_inf^2), cartesian.(caze).Y(:,1), ...
                 'color', 'black')
    end    
    hold off
end

linkaxes(h, 'xy')
ylim([-10, 200])


%% Plot curvilinear u^* for a specific phase against wave age

idx = 86;

figure('color', 'white')
tiledlayout(1,4)

for phase = 1:4
    h(phase) = nexttile;
    title(sprintf('Phase %1.0f', phase))
    hold on
    for s = 1:length(wind_speeds)
        wind_speed = wind_speeds{s};
        u_inf = freestreams.(wind_speed);
    
        for w = 1:length(waves)
            wave = waves{w};
            frequency = frequencies.(wave);
            wavelength = wavelengths.(wave);
            wave_speed = wavelength * frequency;
            wave_age = wave_speed / u_inf;
            caze = strcat(wind_speed, '_WV', wave, '_AG0');
    
            % Calculate u*
            uv_stress_profile = curvilinear.(caze).phase(phase).uv(:, idx);
            u_star = sqrt(max(-uv_stress_profile, [], 'omitnan'));
    
            % Convert to friction velocity
            % Cf = 2 * (u_star / u_inf)^2;
            Cf = u_star / wave_speed;
            wave_age = u_inf / wave_speed;

            % Plot
            scatter(wave_age, Cf, 75, 'filled', 'markerfacecolor', wave_colors{w})
    
        end
    end
    hold off
    axis square
end

linkaxes(h, 'xy')
% ylim([0, 0.55])
% xlim([0.05, 0.35])




%% Fit power law to all data + per phase

% All stacked
clc;
phase_markers = {'o', 'square', 'diamond', '^'};
phase_names = {'\phi = 0', '\phi = \lambda/4', '\phi = \lambda/2', '\phi = 3\lambda/4'};

% First, collect all data for fitting
all_wave_age = [];
all_u_star = [];
phase_wave_age = cell(1, 4);
phase_u_star = cell(1, 4);

for phase = 1:4
    phase_wave_age{phase} = [];
    phase_u_star{phase} = [];
    
    for s = 1:length(wind_speeds)
        wind_speed = wind_speeds{s};
        u_inf = freestreams.(wind_speed);
    
        for w = 1:length(waves)
            wave = waves{w};
            frequency = frequencies.(wave);
            wavelength = wavelengths.(wave);
            wave_speed = wavelength * frequency;
            wave_age = wave_speed / u_inf;
            caze = strcat(wind_speed, '_WV', wave, '_AG0');
    
            % Calculate u*
            uv_stress_profile = curvilinear.(caze).phase(phase).uv(:, idx);
            u_star = sqrt(max(-uv_stress_profile, [], 'omitnan'));

            % Convert to friction velocity
            % Cf = 2 * (u_star / u_inf)^2;
            % Cf = u_star / u_inf;
            % Cf = u_star / wave_speed;
            Cf = u_star / wave_speed;
            wave_age = u_inf / wave_speed;
    
            % Store for fitting
            all_wave_age(end+1) = wave_age;
            all_u_star(end+1) = Cf;
            phase_wave_age{phase}(end+1) = wave_age;
            phase_u_star{phase}(end+1) = Cf;
        end
    end
end

% Fit power law: u* = A * (c/U_inf)^B
% In log space: log(u*) = log(A) + B*log(c/U_inf)

% Fit all data combined
log_wave_age_all = log(all_wave_age);
log_u_star_all = log(all_u_star);
p_all = polyfit(log_wave_age_all, log_u_star_all, 1);
B_all = p_all(1);
A_all = exp(p_all(2));

% Calculate R^2 for all data
u_star_pred_all = A_all * all_wave_age.^B_all;
SS_res_all = sum((all_u_star - u_star_pred_all).^2);
SS_tot_all = sum((all_u_star - mean(all_u_star)).^2);
R2_all = 1 - SS_res_all / SS_tot_all;

fprintf('\n=== POWER LAW FIT: ALL PHASES ===\n')
fprintf('u* = %.4f * (c/U_inf)^(%.4f)\n', A_all, B_all)
fprintf('R^2 = %.4f\n\n', R2_all)

% Fit per phase
A_phase = zeros(1, 4);
B_phase = zeros(1, 4);
R2_phase = zeros(1, 4);

fprintf('=== POWER LAW FIT: PER PHASE ===\n')
for phase = 1:4
    log_wa = log(phase_wave_age{phase});
    log_us = log(phase_u_star{phase});
    p = polyfit(log_wa, log_us, 1);
    B_phase(phase) = p(1);
    A_phase(phase) = exp(p(2));
    
    % R^2
    u_star_pred = A_phase(phase) * phase_wave_age{phase}.^B_phase(phase);
    SS_res = sum((phase_u_star{phase} - u_star_pred).^2);
    SS_tot = sum((phase_u_star{phase} - mean(phase_u_star{phase})).^2);
    R2_phase(phase) = 1 - SS_res / SS_tot;
    
    fprintf('Phase %d: u* = %.4f * (c/U_inf)^(%.4f), R^2 = %.4f\n', ...
            phase, A_phase(phase), B_phase(phase), R2_phase(phase))
end

% Generate fit curves
wave_age_fit = linspace(0, 20, 100);
u_star_fit_all = A_all * wave_age_fit.^B_all;

% Phase-specific fit curves
u_star_fit_phase = zeros(4, length(wave_age_fit));
for phase = 1:4
    u_star_fit_phase(phase, :) = A_phase(phase) * wave_age_fit.^B_phase(phase);
end

% Colors for phase fit lines
phase_fit_colors = [0.5, 0.5, 0.5;    % Gray for phase 1
                    0.7, 0.4, 0.1;    % Brown for phase 2
                    0.1, 0.5, 0.5;    % Teal for phase 3
                    0.5, 0.1, 0.5];   % Purple for phase 4






%% Make test for paper figure with all phases plotted together and sub-tiles per phase

clear h

reordered_phases = [1,4,3,2];
% reordered_phases = [1,2,3,4];

phase_labels = {'$\varphi = 0$', '$\varphi = \lambda / 4$', ...
                '$\varphi = \lambda / 2$', '$\varphi = 3 \lambda / 4$'};

% Fontsizes
tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;

% Marker sizes
fit_lw = 1;
sz = 20;


close all;
fig = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,5]);
t = tiledlayout(1, 4, 'TileSpacing', 'compact', 'padding', 'compact');




% Plot per-phase tiles
for p = 1:4    

    phase = reordered_phases(p);
    h(p) = nexttile;
    set(h(p), 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
    hold on
    title(phase_labels{p}, 'interpreter', 'latex', 'fontsize', labelFontSize)
    % Plot data points
    for s = 1:length(wind_speeds)
        wind_speed = wind_speeds{s};
        u_inf = freestreams.(wind_speed);
    
        for w = 1:length(waves)
            wave = waves{w};
            frequency = frequencies.(wave);
            wavelength = wavelengths.(wave);
            wave_speed = wavelength * frequency;
            wave_age = wave_speed / u_inf;
            caze = strcat(wind_speed, '_WV', wave, '_AG0');
    
            % Calculate u*
            uv_stress_profile = curvilinear.(caze).phase(phase).uv(:, idx);
            u_star = sqrt(max(-uv_stress_profile, [], 'omitnan'));

            % Convert to friction velocity
            Cf = u_star / wave_speed;
            wave_age = u_inf / wave_speed;
    
            % Plot
            scatter(wave_age, Cf, sz, phase_markers{phase}, 'filled', ...
                    'markerfacecolor', wave_colors{w}, 'handlevisibility', 'off')
        end
    end

    % if mod(phase, 2) == 0
    if p > 1
        set(h(p), 'yticklabel',[])
    end

    % Fit line per-phase
    P = plot(wave_age_fit, u_star_fit_phase(phase, :), '-', 'linewidth', fit_lw, 'color', 'black', 'handlevisibility', 'off');
    uistack(P, 'bottom')


    % Add legend
    if p == 4
        for w = 1:length(waves)
            wave = waves{w};
            label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses.(wave));
            hLeg = plot(nan, nan, 'o', ...
            'MarkerFaceColor', wave_colors{w}, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerSize', 4, ...
            'LineStyle', 'none', ...
            'DisplayName', label);
        end
    end
    hold off
    % uistack([P, PP], 'bottom')
end

% Figure sizes
linkaxes(h, 'xy')
% ylim(h, [0, 0.55])
% xlim(h, [0.07, 0.35])
% xticks(h, 0.1:0.1:0.3)
xticks(h, 4:4:12)
% yticks(h, 0:0.1:0.5)
% yticks(h, 0:0.25:0.5)
xlim([2, 13])
ylim([0, 1])
xlabel(t, '$u_{\infty} \mathbin{/} c$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(t, '$u_{\xi}^*  \mathbin{/} c$', 'interpreter', 'latex', 'fontsize', labelFontSize)

% Legend
leg = legend('interpreter', 'latex', 'box', 'off', 'orientation', 'horizontal', 'fontsize', legendFontSize);
leg.Layout.Tile = 'north';
leg.ItemTokenSize(1) = 10;
leg.IconColumnWidth = 19;

% Add a,b,c,d,e panels
% addPanelLabelsFixed(fig, h, {'a', 'b', 'c', 'd', 'e'}, 'FontSize', 10, 'OffsetPts', [-20,-4])



% Save figure
pause(3)
fig_folder = fullfile(figure_folder, 'Friction');
fig_name = 'FrictionVelocity_Curvilinear_ScalingCollapse_New.pdf';
exportgraphics(fig, fullfile(fig_folder, fig_name), 'Resolution', 600, 'ContentType', 'image')
close all 








%% Functions

function hAnn = addPanelLabelsFixed(fig, ax, labels, varargin)
% Places (a),(b),... just outside top-left with FIXED physical offsets.
%
% fig    : figure handle (e.g., totalFigure)
% ax     : array of axes handles
% labels : {'a','b',...} or ["a","b",...]
%
% Name-value:
% 'OffsetPts' : [dx dy] in points (default [-10 6])  (left, up)
% 'FontSize'  : default 12
% 'FontName'  : default 'Times New Roman'

p = inputParser;
addParameter(p,'OffsetPts',[-10 6]);
addParameter(p,'FontSize',12);
addParameter(p,'FontName','Times New Roman');
parse(p,varargin{:});

labels = string(labels);
offPts = p.Results.OffsetPts;

ppi = get(0,'ScreenPixelsPerInch');
offPx = offPts/72 * ppi;   % points -> pixels

hAnn = gobjects(numel(ax),1);

for k = 1:numel(ax)
    if ~isgraphics(ax(k),'axes'), continue; end

    % Axes outer position in figure pixels
    op = getpixelposition(ax(k), true);  % [x y w h] in pixels relative to figure

    % Anchor point: outside top-left of axes
    x = op(1) + offPx(1);
    y = op(2) + op(4) + offPx(2);

    % TeX gives roman parentheses + italic letter
    str = sprintf('(\\it%s\\rm)', labels(k));

    hAnn(k) = annotation(fig,'textbox', ...
        'Units','pixels', ...
        'Position',[x y 40 20], ...   % small box; big enough for "(a)"
        'String',str, ...
        'Interpreter','tex', ...
        'FontName',p.Results.FontName, ...
        'FontSize',p.Results.FontSize, ...
        'LineStyle','none', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom');
    end
end

