%%% Wind + Wave paper wave-deviations test

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

% Cases
wind_speeds = {'WT4', 'WT6', 'WT8'};
waves = {'A', 'B', 'C', 'D'};

% Freestream velocities
freestreams.('WT4') = 2.4181;
freestreams.('WT6') = 3.8709;
freestreams.('WT8') = 5.4289;

% Approximate wavelengths in mm for labeling plots
wavelength_names.A = '410';
wavelength_names.B = '313';
wavelength_names.C = '189';
wavelength_names.D = '124';

steepnesses.A = '0.180';
steepnesses.B = '0.211';
steepnesses.C = '0.305';
steepnesses.D = '0.267';

% Load data
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        disp(caze)
        tmp = load(fullfile(curvilinear_path, strcat(caze, '_CURVILINEAR.mat')));
        tmp = tmp.output;
        curvilinear.(caze) = tmp;
    end
end

clc; fprintf('All cases loaded\n')
clear caze tmp w no_wave_caze cartesian_path curvilinear_path data_path
clear wind_speed wave s w caze tmp project_path


%% Compute average of phase averages for each case

% How much to chop off edges [mm]
edge_buffer = 3;

% Loop through wind speeds
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};

    % Loop through waves
    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        disp(caze)

        X = curvilinear.(caze).X;
        Y = curvilinear.(caze).Y;

        u_tmp  = nan([size(X), 4]);
        v_tmp  = nan([size(X), 4]);
        hl_tmp = nan([size(X), 4]);
        vl_tmp = nan([size(X), 4]);

        % Loop through phases
        for phase = 1:4
            u  = curvilinear.(caze).phase(phase).u;
            v  = -1 * curvilinear.(caze).phase(phase).v;
            hl = curvilinear.(caze).phase(phase).horizontal_lines;
            vl = curvilinear.(caze).phase(phase).vertical_lines;
        
            u(vl < min(X, [], 'all') + edge_buffer) = nan;
            u(vl > max(X, [], 'all') - edge_buffer) = nan;
            u(hl > max(Y, [], 'all') - 3 * edge_buffer) = nan;
            v(vl < min(X, [], 'all') + edge_buffer) = nan;
            v(vl > max(X, [], 'all') - edge_buffer) = nan;
            v(hl < max(Y, [], 'all') - edge_buffer) = nan;
        
            u_tmp(:,:,phase)  = u;
            v_tmp(:,:,phase)  = v;
            hl_tmp(:,:,phase) = hl;
            vl_tmp(:,:,phase) = vl;
        end

        % Compute mean of phase
        phase_mean.(caze).u = mean(u_tmp, 3, 'omitnan');
        phase_mean.(caze).v = mean(v_tmp, 3, 'omitnan');
        phase_mean.(caze).horizontal_lines = mean(hl_tmp, 3, 'omitnan');
        phase_mean.(caze).vertical_lines = mean(vl_tmp, 3, 'omitnan');

    end
end

clear edge_buffer s wind_speed w wave caze X Y 
clear phase u_tmp v_tmp hl_tmp vl_tmp u v hl vl


%% Example plot of the phase averages, and average of phases

% wind_speed = 'WT6';
% wave = 'C';
% caze = strcat(wind_speed, '_WV', wave, '_AG0');
% freestream = freestreams.(wind_speed);
% 
% cmin = 0.3;
% cmax = 1;
% 
% figure('color', 'white')
% tiledlayout(1,5)
% for phase = 1:4
%     h(phase) = nexttile;
%     hold on
%     contourf(curvilinear.(caze).phase(phase).vertical_lines, ...
%              curvilinear.(caze).phase(phase).horizontal_lines, ...
%              curvilinear.(caze).phase(phase).u / freestream, ...
%              100, 'linestyle', 'none');
%     plot(curvilinear.(caze).X(1,:), curvilinear.(caze).phase(phase).wave_profile, ...
%          'linewidth', 2, 'color', 'black')
%     hold off
%     axis equal
%     clim([cmin, cmax])
%     title(sprintf('Phase %1.0f', phase))
% end
% 
% h(5) = nexttile;
% box off
% contourf(phase_mean.(caze).vertical_lines, ...
%          phase_mean.(caze).horizontal_lines, ...
%          phase_mean.(caze).u / freestream, ...
%          100, 'linestyle', 'none');
% axis equal
% clim([cmin, cmax])
% C = colorbar();
% C.Label.String = '$u / u_{\infty}$';
% C.Label.Interpreter = 'latex';
% C.Label.FontSize = 16;
% title('Average of Phase Averages')
% 
% linkaxes(h, 'xy')
% xlim([-20, 250])
% ylim([-20, 200])
% 
% clear wind_speed wave caze freestream phase h

%% Plot: Compare the mean of phases

% wind_speed = 'WT6';
% freestream = freestreams.(wind_speed);
% 
% cmin = 0.5;
% cmax = 1;
% 
% % Plot
% figure('color', 'white')
% tiledlayout(1, length(waves))
% 
% % Loop through waves
% for w = 1:length(waves)
%     wave = waves{w};
%     caze = strcat(wind_speed, '_WV', wave, '_AG0');
%     disp(caze)
% 
%     h(w) = nexttile;
%     contourf(phase_mean.(caze).vertical_lines, ...
%              phase_mean.(caze).horizontal_lines, ...
%              phase_mean.(caze).u / freestream, ...
%             100, 'linestyle', 'none')
%     yline(0, 'linewidth', 2)
%     axis equal
%     clim([cmin, cmax])
%     label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses.(wave));
%     title(label, 'Interpreter', 'latex')
%     if w == 4
%         colorbar()
%     end
% end
% 
% linkaxes(h, 'xy')
% ylim([-5, 210])
% 
% clear wind_speed freestream cmin cmax w wave caze h 


%% Compute the wave-coherent fluctuations per phase

% Loop through wind speeds
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};

    % Loop through waves
    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        disp(caze)

        for phase = 1:4
            wave_coherent.(caze)(phase).u = curvilinear.(caze).phase(phase).u - phase_mean.(caze).u;
            wave_coherent.(caze)(phase).v = (-1 * curvilinear.(caze).phase(phase).v) - phase_mean.(caze).v;
        end
    end
end


%% Plot for a fixed wind speed and wave across all wind speeds


wave_transparency = 0.25;


% Fontsizes
levels = 100;
linewidth = 1;

% colorbar_fontsize = 16;
% tickFontSize = 14;
% labelFontSize = 16;
% legendFontSize = 14;
% titleFontSize = labelFontSize;

colorbar_fontsize = 10;
tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;
titleFontSize = labelFontSize;

% How much to chop off edges [mm]
edge_buffer = 3;

% What to plot
wind_speed = 'WT6';
freestream = freestreams.(wind_speed);
phase = 2;

% Plot
clc; close all
figure('color', 'white', 'Units', 'centimeters', 'position', [10,10,13,4])
t = tiledlayout(1, 4, 'Padding', 'tight', 'TileSpacing', 'tight');

% Loop through waves
for w = 1:length(waves)
    wave = waves{w};
    caze = strcat(wind_speed, '_WV', wave, '_AG0');
    disp(caze)

    % Load data
    tmp = wave_coherent.(caze)(phase).u / freestream;
    vl = curvilinear.(caze).phase(phase).vertical_lines;
    hl = curvilinear.(caze).phase(phase).horizontal_lines;

    % Trimp edges
    X = curvilinear.(caze).X;
    Y = curvilinear.(caze).Y;
    tmp(vl < min(X, [], 'all') + edge_buffer) = nan;
    tmp(vl > max(X, [], 'all') - edge_buffer) = nan;
    tmp(hl > max(Y, [], 'all') + 2 * edge_buffer) = nan;

    % Trim first rows of real data right above surface
    is_all_nans = all(isnan(tmp), 2);
    is_not_all_nans = ~is_all_nans;
    last_non_nan_row_index = find(is_not_all_nans, 1, 'last');
    first_non_nan_row_index = find(is_not_all_nans, 1, 'first');

    buff = 2;
    tmp(last_non_nan_row_index - buff:end, :) = nan;
    tmp(1:first_non_nan_row_index + buff, :) = nan;

    % Plot tile
    h(w) = nexttile;
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
    hold on
    contourf(vl, hl, tmp, levels, 'linestyle', 'none')
    plot(curvilinear.(caze).X(1,:), curvilinear.(caze).phase(phase).wave_profile, 'linewidth', linewidth, 'color', 'black')
    
    % Shaded region below wave
    hFill = patch( ...
    [curvilinear.(caze).X(1,:), fliplr(curvilinear.(caze).X(1,:))], ...
    [curvilinear.(caze).phase(phase).wave_profile, -100 * ones(size(curvilinear.(caze).phase(phase).wave_profile))], ...
    'k', ...
    'FaceAlpha', wave_transparency, ...      
    'EdgeColor', 'none', ...    
    'HandleVisibility', 'off'); 
    uistack(hFill, 'bottom')

    hold off

    axis equal
    clim([-0.1, 0.1])
    colormap(slanCM('bwr'))
    label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses.(wave));
    title(label, 'Interpreter', 'latex', 'fontsize', titleFontSize)
    xticks(0:100:200)
    yticks(0:100:200)

    if w ~= 1
        ax = gca;
        ax.YTickLabel = [];
    end
    if w == 4
        C = colorbar;
        C.Label.String = '$\tilde{u}_{\xi} \mathbin{/} u_{\infty}$';
        C.FontSize = tickFontSize;
        C.Label.Interpreter = 'latex';
        C.Label.FontSize = labelFontSize;
        C.TickLabelInterpreter = 'latex';
        C.Ticks = -0.1:0.1:0.1;
    end
end

linkaxes(h, 'xy')
ylim([-20, 200])

% xlim([-10, 250])
xlim([0, max(curvilinear.(caze).X(1,:))])
xlabel(t, '$x$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(t, '$y$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)



% Save figure
pause(3)
figure_name = sprintf('%s_Phase%1.0f_u_WaveCoherentDeviations.pdf', wind_speed, phase);
exportgraphics(t, fullfile(figure_folder, 'Phase', figure_name), 'Resolution', 600, 'ContentType', 'image');
close all
clc; fprintf('Generated figure: %s\n\n', figure_name)