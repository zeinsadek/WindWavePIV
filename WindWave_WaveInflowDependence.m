%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions');

% Cases with waves
WTs = {'4', '6', '8'};
WVs = {'A', 'B', 'C', 'D'};
wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};

% Paths
MN_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means/';

% Wind speeds
freestreams.('WT4') = 2.4181;
freestreams.('WT6') = 3.8709;
freestreams.('WT8') = 5.4289;

% Wave parameters
wavelengths.A = 410.8;
wavelengths.B = 313.3;
wavelengths.C = 189.6;
wavelengths.D = 124.3;

amplitudes.A = 11.78;
amplitudes.B = 10.53;
amplitudes.C = 9.21;
amplitudes.D = 5.29;

frequencies.A = 1.96;
frequencies.B = 2.27;
frequencies.C = 3;
frequencies.D = 3.93;

% Approximate wavelengths in mm for labeling plots
wavelength_names.A = '410';
wavelength_names.B = '313';
wavelength_names.C = '189';
wavelength_names.D = '124';

steepness_names.A = '0.180';
steepness_names.B = '0.211';
steepness_names.C = '0.305';
steepness_names.D = '0.267';

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
        MN_path = fullfile(MN_folder, strcat(caze, '_MEANS.mat'));
       
        % Save MN to Structure
        means_temp = load(MN_path);
        means_temp = means_temp.output;
        means.(caze) = means_temp;
    end
    fprintf('\n')
end

% Clear temporary variables
clear BL_path LL_path MN_path
clear BL_temp LL_temp means_temp
clear caze i no_wave_case s speed
clear MN_folder

%% Plot Phase-Average wave profile for each wave
% Fixed phase, different wind speeds

% Titles
phase_titles = {'$\varphi = 0$',...
                '$\varphi = \lambda / 4$',...
                '$\varphi = \lambda / 2$',...
                '$\varphi = 3 \lambda / 4$'};

phase = 4;
linestyles = {'-.', '--', '-'};
lw = 2;

tickFontSize = 14;
labelFontSize = 16;
legendFontSize = 14;

figure('color', 'white', 'units', 'centimeters', 'position', [2, 2, 16, 15]);
t = tiledlayout(4,1, 'TileSpacing', 'tight');
sgtitle(phase_titles{phase}, 'interpreter', 'latex', 'fontsize', labelFontSize)

for w = 1:length(WVs)

    h(w) = nexttile;
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
    hold on
    for s = 1:length(WTs)
        speed = WTs{s};

        % Get data
        caze = strcat('WT', num2str(speed), '_WV', WVs{w}, '_AG0');
        waves = means.(caze).phase(phase).waves;
        x = means.(caze).X(1,:);
        [num_waves, ~] = size(waves);
        waves = imresize(waves, [num_waves, length(x)]);

        % Plot
        plot(x, mean(waves, 1, 'omitnan'), 'color', wave_colors{w}, ...
             'linewidth', lw, 'linestyle', linestyles{s}, ...
             'HandleVisibility', 'off')

    end
    
    % Plot reference wave
    reference = means.(caze).phase(phase).reference_wave;
    P = plot(x, reference, 'color', 'black', 'linewidth', lw, 'HandleVisibility', 'off');
    P.Color(4) = 0.5;
    hold off

    axis equal
    xlim([0, 234])
    ylim([-25, 25])
    yticks(-15:15:15)

    if w ~= 4
        ax = gca;
        ax.XTickLabel = [];
        ax.XAxis.Visible = 'off';
    end
end

hold on
% Add legend for wavelength
for w = 1:length(WVs)
    wave = WVs{w};
    label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepness_names.(wave));
    plot(nan, nan, 'linewidth', lw, 'color', wave_colors{w}, ...
         'Displayname', label)
end

% Add white space
plot(nan, nan, 'color', 'white', 'DisplayName', ' ')

% Add linestyle legend

for s = 1:length(WTs)
    plot(nan, nan, 'linewidth', lw, 'Color', 'black', ...
         'LineStyle', linestyles{s}, ...
         'Displayname', sprintf('$u_{\\infty} = %1.2f$ m/s', freestreams.(['WT', WTs{s}])));
end
hold off
leg = legend('interpreter', 'latex', 'box', 'off');
leg.Layout.Tile = 'east';
xlabel(t, '$x$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(t, '$y$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)

linkaxes(h, 'xy')


%% Plot Difference of Phase-Average wave profile for each wave
% Fixed phase, different wind speeds

% Titles
phase_titles = {'$\varphi = 0$',...
                '$\varphi = \lambda / 4$',...
                '$\varphi = \lambda / 2$',...
                '$\varphi = 3 \lambda / 4$'};

phase = 4;
linestyles = {'-.', '--', '-'};
lw = 2;

tickFontSize = 14;
labelFontSize = 16;
legendFontSize = 14;

clc; close all
figure('color', 'white', 'units', 'centimeters', 'position', [2, 2, 16, 15]);
t = tiledlayout(4,1, 'TileSpacing', 'tight');
sgtitle(phase_titles{phase}, 'interpreter', 'latex', 'fontsize', labelFontSize)

for w = 1:length(WVs)

    h(w) = nexttile;
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
    hold on
    for s = 1:length(WTs)
        speed = WTs{s};

        % Get data
        caze = strcat('WT', num2str(speed), '_WV', WVs{w}, '_AG0');
        waves = means.(caze).phase(phase).waves;
        reference = means.(caze).phase(phase).reference_wave;
        x = means.(caze).X(1,:);
        [num_waves, ~] = size(waves);
        waves = imresize(waves, [num_waves, length(x)]);
        avg_wave = mean(waves, 1, 'omitnan');
        wave_diff = avg_wave(:) - reference(:);

        % Plot
        plot(x, wave_diff, 'color', wave_colors{w}, ...
             'linewidth', lw, 'linestyle', linestyles{s}, ...
             'HandleVisibility', 'off')

    end
    
    % Plot reference wave
    reference = means.(caze).phase(phase).reference_wave;
    P = plot(x, reference, 'color', 'black', 'linewidth', lw, 'HandleVisibility', 'off');
    P.Color(4) = 0.5;
    yline(0, 'color', 'black', 'alpha', 0.5, 'linewidth', lw, 'HandleVisibility', 'off');
    hold off

    axis equal
    xlim([0, 234])
    ylim([-25, 25])
    yticks(-15:15:15)

    if w ~= 4
        ax = gca;
        ax.XTickLabel = [];
        ax.XAxis.Visible = 'off';
    end
end

hold on
% Add legend for wavelength
for w = 1:length(WVs)
    wave = WVs{w};
    label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepness_names.(wave));
    plot(nan, nan, 'linewidth', lw, 'color', wave_colors{w}, ...
         'Displayname', label)
end

% Add white space
plot(nan, nan, 'color', 'white', 'DisplayName', ' ')

% Add linestyle legend

for s = 1:length(WTs)
    plot(nan, nan, 'linewidth', lw, 'Color', 'black', ...
         'LineStyle', linestyles{s}, ...
         'Displayname', sprintf('$u_{\\infty} = %1.2f$ m/s', freestreams.(['WT', WTs{s}])));
end
hold off
leg = legend('interpreter', 'latex', 'box', 'off');
leg.Layout.Tile = 'east';
xlabel(t, '$x$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(t, '$y$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)

linkaxes(h, 'xy')


%% Plot all phases of a single wavelength: to show calculation



%% Compute the max and min difference between reference and wave, for each wavelength and wind speed

% Make arrays to save values
% Row ~ wave
% Column ~ wind speed
avg_abs_differences = nan(4, 3);


clc; close all
figure('color', 'white')
hold on
for w = 1:length(WVs)
    for s = 1:length(WTs)
        speed = WTs{s};
        freestream = freestreams.(['WT', num2str(speed)]);

        % Get data
        caze = strcat('WT', num2str(speed), '_WV', WVs{w}, '_AG0');
        disp(caze)
        tmp_plot= nan(4,length(x));
        for phase = 1:4

            % Get waves
            waves = means.(caze).phase(phase).waves;
            x = means.(caze).X(1,:);
            [num_waves, ~] = size(waves);
            waves = imresize(waves, [num_waves, length(x)]);
            reference = means.(caze).phase(phase).reference_wave;
            avg_wave = mean(waves, 1, 'omitnan');

            % Difference between reference and wave
            wave_diff = avg_wave(:) - reference(:);
            abs_wave_diff = abs(wave_diff);
            tmp_plot(phase, :) = abs_wave_diff;
        end

        p05_mm = prctile(tmp_plot, 5, 'all');
        p95_mm = prctile(tmp_plot, 95, 'all');

        disp(p05_mm)
        disp(p95_mm)

        tmpMax(s) = p95_mm;
        tmpMin(s) = p05_mm;

        % Worst difference
        max_mm  = max(tmp_plot(:), [], 'all', 'omitnan');  

        % Avg difference
        error = mean(tmp_plot, 'all', 'omitnan');

        % Plot
        scatter(freestream, error, 50, 'filled', 'MarkerFaceColor', wave_colors{w})
        tmpX(s) = freestream;
        tmpY(s) = error;
    end   
    plot(tmpX, tmpY, 'linewidth', 2, 'color', wave_colors{w})

    avg_abs_differences(w,:) = tmpY;

    % Shaded region
    % hFill = patch( ...
    % [tmpX, fliplr(tmpX)], ...
    % [tmpMax, fliplr(tmpMin)], ...
    % hex2rgb(wave_colors{w}), ...
    % 'FaceAlpha', 0.15, ...        % transparency
    % 'EdgeColor', 'none', ...      % no outline
    % 'HandleVisibility', 'off');   % keep legend clean
    % uistack(hFill, 'bottom')

end
hold off
ylim([0,5])
% yline(7.5, 'linestyle', '--')
xlabel('$u_{\infty}$ [m/s]', 'interpreter', 'latex', 'Fontsize', labelFontSize)
ylabel('$\overline{ \left| \Delta \eta \right|_{\varphi, x}}$ [mm]', 'Interpreter', 'latex', 'Fontsize', labelFontSize)



%% Same plot but normalized by the first wind-speed, to show variation across wind speeds

% Normalize by slowest wind speed
normalized_avg_abs_differences = avg_abs_differences ./ avg_abs_differences(:,1);

figure('color', 'white')
hold on
for w = 1:length(WVs)
    disp(w)
    scatter(tmpX, normalized_avg_abs_differences(w,:) - 1 , 50, 'filled', 'markerfacecolor', wave_colors{w})
    plot(tmpX, normalized_avg_abs_differences(w,:) - 1, 'linewidth', 2, 'color', wave_colors{w})
end
hold off
ylim([-0.5,0.5])