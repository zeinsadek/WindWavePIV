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
MN_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means/';

% Save folder
fig_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/pdf_test6';


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
        MN_path = fullfile(MN_folder, strcat(caze, '_MEANS.mat'));
       
        % Save MN to Structure
        means_temp = load(MN_path);
        means_temp = means_temp.output;
        MN_Waves.(caze) = means_temp;
    end
    fprintf('\n')
end

% Clear temporary variables
clear BL_path LL_path MN_path
clear BL_temp LL_temp means_temp
clear caze i no_wave_case s speed


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE AVERAGE WAVE PROFILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

caze = 'WT4_WVC_AG0';

linewidth = 1;

% Colors
reference_color = 'green';
average_color = 'green';
max_color = 'red';

% Plot fontsizes
tickFontSize = 8;
labelFontSize = 10;
titleFontSize = 10;
legendFontSize = 8;
annotationFontSize = 6;


% Titles
phase_titles = {'$\varphi = 0$',...
                '$\varphi = \lambda / 4$',...
                '$\varphi = \lambda / 2$',...
                '$\varphi = 3 \lambda / 4$'};


% Plot
clc; close all
totalFigure = figure('color', 'white', 'units', 'centimeters', 'position', [10, 10, 13, 6.0]);
t = tiledlayout(2, 2, "TileSpacing", "tight", "Padding", "loose");

for p = 1:4
    
    % Load data
    waves = MN_Waves.(caze).phase(p).waves;
    reference = MN_Waves.(caze).phase(p).reference_wave;
    max_wave = MN_Waves.(caze).phase(p).max_wave_profile;
    x = MN_Waves.(caze).X(1,:);
    [num_waves, ~] = size(waves);

    fprintf('Phase %1.0f\n', p)
    fprintf('# Images = %4.0f\n\n', num_waves)

    % Plot
    h(p) = nexttile();
    ax = gca;
    set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)

    hold on
    % Plot all instantaneous wave profiles
    for i = 1:num_waves
        plt = plot(x, imresize(waves(i,:), size(x)), 'color', 'black', 'linewidth', 0.1, "HandleVisibility", "off");
        plt.Color(4) = 0.1;
    end

    % Plot reference wave profile
    plot(x, reference, 'color', reference_color, 'linewidth', linewidth, "HandleVisibility", "off")

    % Plot max wave profile
    plot(x, smoothdata(max_wave, 'gaussian', 15), 'color', max_color, 'linewidth', linewidth, "HandleVisibility", 'off')


    % Make legend in specific order
    % Add legend entry for instantaneous profiles
    plot(nan, nan, 'color', 'black', 'linewidth', 1, 'DisplayName', '$\eta_{i}$')
    % Add legend entry for reference profiles
    plot(nan, nan, 'color', reference_color, 'linewidth', 1, "DisplayName", "$\eta$")
    % Add legend entry for max profiles
    plot(nan, nan, 'color', max_color, 'linewidth', 1, "DisplayName", "$\eta_{max}$")
    hold off
    
    axis equal
    xlim([0, range])
    ylim([-30, 30])
    yticks(-20:20:20)
    title(phase_titles{p}, 'interpreter', 'latex', 'FontSize', titleFontSize)


    % Add image count text box
    alf = 0.5;
    text(1, 1, sprintf('$%d$ Images', num_waves), ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'right', ...
        'Interpreter', 'latex', ...
        'FontSize', annotationFontSize, ...
        'Color', alf .* [1,1,1])

    if mod(p,2) == 0
        ax.YTickLabel = [];
        ax.YAxis.Visible = 'off';
    end

end

leg = legend('Orientation', 'horizontal', 'interpreter', 'latex', ...
             'FontSize', legendFontSize, 'box', 'off');
leg.Layout.Tile = 'north';
leg.Box = "off";
leg.IconColumnWidth = 19;
ylabel(t, '$y$ [mm]', 'Interpreter', 'latex', 'FontSize', labelFontSize)
xlabel(t, '$x$ [mm]', 'Interpreter', 'latex', 'FontSize', labelFontSize)


% Save figure
pause(3)
fig_name = strcat(caze, '_PhaseAverageWaves.pdf');
exportgraphics(totalFigure, fullfile(fig_folder, fig_name), 'Resolution', 600, 'ContentType', 'image')
close all 




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of images used in phase averages


clc; close all

c = 1;
for s = 1:length(WTs)
    wind_speed = WTs{s};

    for w = 1:length(WVs)
        wave = WVs{w};
        caze = strcat('WT', wind_speed, '_WV', wave, '_AG0');
        fprintf('%s \n', caze)

        for p = 1:4
            
            % Load data
            waves = MN_Waves.(caze).phase(p).waves;
            % reference = MN_Waves.(caze).phase(p).reference_wave;
            % max_wave = MN_Waves.(caze).phase(p).max_wave_profile;
            x = MN_Waves.(caze).X(1,:);
            [num_waves, ~] = size(waves);
        
            fprintf('Phase %1.0\nf', p)
            fprintf('# Images = %4.0f\n\n', num_waves)

            % Save number
            phase_avg_images(c) = num_waves;
            c = c + 1;
        end
    end
end


clc;
% Min
fprintf('Min number of images: %3.1f\n', min(phase_avg_images, [], 'all'))

% Max
fprintf('Max number of images: %3.1f\n', max(phase_avg_images, [], 'all'))

% Mean
fprintf('Mean number of images: %3.1f\n', mean(phase_avg_images, 'all'))



