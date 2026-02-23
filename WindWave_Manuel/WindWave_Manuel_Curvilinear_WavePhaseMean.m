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
% waves = {'A', 'B', 'C', 'D'};
% wave_colors = {'#FE6202', '#DC2680', '#775EEF', '#648FFF'};

% Cases manuel wants
waves = {'A', 'C'};
wave_colors = {'#FE6202', '#775EEF'};

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


%% Compute average of phase averages for each case

% How much to chop off edges [mm]
edge_buffer = 1;
top_buffer = 1;
bottom_buffer = 2;

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
        
            % Crop the sides
            u(vl < min(X, [], 'all') + edge_buffer) = nan;
            u(vl > max(X, [], 'all') - edge_buffer) = nan;
            v(vl < min(X, [], 'all') + edge_buffer) = nan;
            v(vl > max(X, [], 'all') - edge_buffer) = nan;

            % Crop the top
            u(hl > max(Y, [], 'all') - top_buffer) = nan;            
            v(hl < max(Y, [], 'all') - top_buffer) = nan;

            % Trimp edges
            % X = curvilinear.(caze).X;
            % Y = curvilinear.(caze).Y;
            % tmp(vl < min(X, [], 'all') + edge_buffer) = nan;
            % tmp(vl > max(X, [], 'all') - edge_buffer) = nan;
            % tmp(hl > max(Y, [], 'all') + 2 * edge_buffer) = nan;
        
            % Trim first rows of real data right above surface
            is_all_nans = all(isnan(u), 2);
            is_not_all_nans = ~is_all_nans;
            last_non_nan_row_index = find(is_not_all_nans, 1, 'last');
            first_non_nan_row_index = find(is_not_all_nans, 1, 'first');
        
            u(last_non_nan_row_index - bottom_buffer:end, :) = nan;
            u(1:first_non_nan_row_index + bottom_buffer, :) = nan;

            % Save
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



% Plot: Compare the mean of phases
wind_speed = 'WT6';
freestream = freestreams.(wind_speed);

cmin = 0.5;
cmax = 1;

% Plot
clc; close all
figure('color', 'white')
tiledlayout(1, length(waves))

% Loop through waves
for w = 1:length(waves)
    wave = waves{w};
    caze = strcat(wind_speed, '_WV', wave, '_AG0');
    disp(caze)

    h(w) = nexttile;
    contourf(phase_mean.(caze).vertical_lines, ...
             phase_mean.(caze).horizontal_lines, ...
             phase_mean.(caze).u / freestream, ...
            100, 'linestyle', 'none')
    yline(0, 'linewidth', 2)
    axis equal
    clim([cmin, cmax])
    label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses.(wave));
    title(label, 'Interpreter', 'latex')
    if w == 4
        colorbar()
    end
end

linkaxes(h, 'xy')
ylim([-5, 200])

clear wind_speed freestream cmin cmax w wave caze h 


%% Plot profiles (y non-normalized)

wind_speed_linestyles = {':', '--', '-'};
labelFontSize = 14;

clc; close all
figure('color', 'white')
set(gca, 'ticklabelinterpreter', 'latex')
hold on

% Loop through wind speeds
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    freestream = freestreams.(wind_speed);

    % Loop through waves
    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        wavelength = wavelengths.(wave);
        disp(caze)

        % Load averaged u velocity
        Xc = phase_mean.(caze).vertical_lines;
        Yc = phase_mean.(caze).horizontal_lines;
        Uc = phase_mean.(caze).u / freestream;

        % Get profiles
        xc = Xc(1,:);
        idx = round(length(xc)/2);

        u_profile = Uc(:, idx);
        y_profile = Yc(:, idx);

        % Plot
        plot(u_profile, y_profile, 'color', wave_colors{w}, 'linewidth', 2, ...
             'linestyle', wind_speed_linestyles{s}, 'HandleVisibility', 'off')


        % Save profiles to sent to Manuel
        output.(caze).mean.u_normalized = u_profile;
        output.(caze).mean.y_unnonormalized = y_profile;
        output.(caze).mean.y_normalized = y_profile / wavelength;

    end
end
hold off

% Add legend for waves
hold on
for w = 1:length(waves)
    wave = waves{w};
    label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses.(wave));
    plot(nan, nan, 'color', wave_colors{w}, 'linewidth', 2, ...
         'displayname', label)
end

% Add legend for windspeeds
plot(nan, nan, 'color', 'white', 'displayname', ' ')
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    freestream = freestreams.(wind_speed);

    label = sprintf('$u_{\\infty} = %1.2f~m/s$', freestream);
    plot(nan, nan, 'color', 'black', 'linewidth', 2, ...
         'linestyle', wind_speed_linestyles{s}, 'displayname', label)
end
hold off

leg = legend('box', 'off', 'interpreter', 'latex', 'location', 'northwest');
xlabel('$\overline{\langle u_{\xi} \rangle}_{\varphi} \mathbin{/} u_{\infty}$', ...
       'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$y$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)

ylim([0, 200])
xlim([0.4, 1.05])


%% Plot profiles (y normalized by wavelength)

wind_speed_linestyles = {':', '--', '-'};
labelFontSize = 14;

clc; close all
figure('color', 'white')
set(gca, 'ticklabelinterpreter', 'latex')
hold on

% Loop through wind speeds
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    freestream = freestreams.(wind_speed);

    % Loop through waves
    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        wavelength = wavelengths.(wave);
        disp(caze)

        % Load averaged u velocity
        Xc = phase_mean.(caze).vertical_lines;
        Yc = phase_mean.(caze).horizontal_lines;
        Uc = phase_mean.(caze).u / freestream;

        % Get profiles
        xc = Xc(1,:);
        idx = round(length(xc)/2);

        u_profile = Uc(:, idx);
        y_profile =  Yc(:, idx) / wavelength;

        % Plot
        plot(u_profile, y_profile, 'color', wave_colors{w}, 'linewidth', 2, ...
             'linestyle', wind_speed_linestyles{s}, 'HandleVisibility', 'off')

    end
end
hold off

% Add legend for waves
hold on
for w = 1:length(waves)
    wave = waves{w};
    label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses.(wave));
    plot(nan, nan, 'color', wave_colors{w}, 'linewidth', 2, ...
         'displayname', label)
end

% Add legend for windspeeds
plot(nan, nan, 'color', 'white', 'displayname', ' ')
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    freestream = freestreams.(wind_speed);

    label = sprintf('$u_{\\infty} = %1.2f~m/s$', freestream);
    plot(nan, nan, 'color', 'black', 'linewidth', 2, ...
         'linestyle', wind_speed_linestyles{s}, 'displayname', label)
end
hold off
leg = legend('box', 'off', 'interpreter', 'latex', 'location', 'northwest');

xlabel('$\overline{\langle u_{\xi} \rangle}_{\varphi} \mathbin{/} u_{\infty}$', ...
       'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$y \mathbin{/} \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)

ylim([0, 1.1])
xlim([0.4, 1.05])







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


%% Profiles of wave coherent deviations at peak and trough
% (y non-normalized)

clear h 
wind_speed_linestyles = {':', '--', '-'};
labelFontSize = 14;

clc; close all
figure('color', 'white')
tile = tiledlayout(1, 2, 'TileSpacing', 'compact');

% Loop through peak and trough phases
for phase = 1:2:4
    fprintf('Phase %1.0f\n', phase)

    h(phase) = nexttile;
    hold on
    % Loop through waves
    for w = 1:length(waves)
        wave = waves{w};
        wavelength = wavelengths.(wave);
        amplitude = amplitudes.(wave);
    
        % Loop through wind speeds
        for s = 1:length(wind_speeds)
            wind_speed = wind_speeds{s};
            freestream = freestreams.(wind_speed);
        
            caze = strcat(wind_speed, '_WV', wave, '_AG0');
            disp(caze)
        
            % Load data
            Uc = wave_coherent.(caze)(phase).u / freestream;
            Xc = curvilinear.(caze).phase(phase).vertical_lines;
            Yc = curvilinear.(caze).phase(phase).horizontal_lines;

            % Trim first rows of real data right above surface (only for
            % trough)
            if phase == 3
                is_all_nans = all(isnan(Uc), 2);
                is_not_all_nans = ~is_all_nans;
                last_non_nan_row_index = find(is_not_all_nans, 1, 'last');
                first_non_nan_row_index = find(is_not_all_nans, 1, 'first');
    
                buff = 0;
                Uc(last_non_nan_row_index - buff:end, :) = nan;
                Uc(1:first_non_nan_row_index + buff, :) = nan;
            end
    
            % Take vertical profile
            xc = Xc(1,:);
            idx = round(length(xc)/2);

            u_profile = Uc(:, idx);
            y_profile = Yc(:, idx);

            % Zero to wave surface
            if phase == 1
                offset = amplitude;
            elseif phase == 3
                offset = -amplitude;
            end

            % Plot
            plot(u_profile, (y_profile - offset), 'color', wave_colors{w}, 'linewidth', 2, ...
                 'linestyle', wind_speed_linestyles{s}, 'HandleVisibility', 'off')

            % x limits
            if phase == 1
                xlim([-0.01, 0.3])
            elseif phase == 3
                xlim([-0.3, 0.01])
                set(gca, 'YTickLabel', []);
            end


            % Save profiles to sent to Manuel
            if phase == 1
                field = 'peak';
            elseif phase == 3
                field = 'trough';
            end
            output.(caze).coherent.(field).u_normalized = u_profile;
            output.(caze).coherent.(field).y_unnormalized = (y_profile - offset);
            output.(caze).coherent.(field).y_normalized = (y_profile - offset) / wavelength;

        end
    end
    hold off
    
    % Titles
    if phase == 1
        title('Wave Peak')
    elseif phase == 3
        title('Wave Trough')
    end

end

% Add legend for waves
hold on
for w = 1:length(waves)
    wave = waves{w};
    label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses.(wave));
    plot(nan, nan, 'color', wave_colors{w}, 'linewidth', 2, ...
         'displayname', label)
end

% Add legend for windspeeds
plot(nan, nan, 'color', 'white', 'displayname', ' ')
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    freestream = freestreams.(wind_speed);

    label = sprintf('$u_{\\infty} = %1.2f~m/s$', freestream);
    plot(nan, nan, 'color', 'black', 'linewidth', 2, ...
         'linestyle', wind_speed_linestyles{s}, 'displayname', label)
end
hold off
leg = legend('box', 'off', 'interpreter', 'latex', 'location', 'northwest');
leg.Layout.Tile = 'east';

linkaxes(h, 'y')
ylim([0, 200])
xlabel(tile, '$\tilde{u}_{\xi} \mathbin{/} u_{\infty}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(tile, '$y$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)


%% Profiles of wave coherent deviations at peak and trough
% (y non-normalized)

clear h 
wind_speed_linestyles = {':', '--', '-'};
labelFontSize = 14;

clc; close all
figure('color', 'white')
tile = tiledlayout(1, 2, 'TileSpacing', 'compact');

% Loop through peak and trough phases
for phase = 1:2:4
    fprintf('Phase %1.0f\n', phase)

    h(phase) = nexttile;
    hold on
    % Loop through waves
    for w = 1:length(waves)
        wave = waves{w};
        wavelength = wavelengths.(wave);
        amplitude = amplitudes.(wave);
    
        % Loop through wind speeds
        for s = 1:length(wind_speeds)
            wind_speed = wind_speeds{s};
            freestream = freestreams.(wind_speed);
        
            caze = strcat(wind_speed, '_WV', wave, '_AG0');
            disp(caze)
        
            % Load data
            Uc = wave_coherent.(caze)(phase).u / freestream;
            Xc = curvilinear.(caze).phase(phase).vertical_lines;
            Yc = curvilinear.(caze).phase(phase).horizontal_lines;

            % Trim first rows of real data right above surface (only for
            % trough)
            if phase == 3
                is_all_nans = all(isnan(Uc), 2);
                is_not_all_nans = ~is_all_nans;
                last_non_nan_row_index = find(is_not_all_nans, 1, 'last');
                first_non_nan_row_index = find(is_not_all_nans, 1, 'first');
    
                buff = 0;
                Uc(last_non_nan_row_index - buff:end, :) = nan;
                Uc(1:first_non_nan_row_index + buff, :) = nan;
            end
    
            % Take vertical profile
            xc = Xc(1,:);
            idx = round(length(xc)/2);

            u_profile = Uc(:, idx);
            y_profile = Yc(:, idx);

            % Zero to wave surface
            if phase == 1
                offset = amplitude;
            elseif phase == 3
                offset = -amplitude;
            end

            % Plot
            plot(u_profile, (y_profile - offset) / wavelength, 'color', wave_colors{w}, 'linewidth', 2, ...
                 'linestyle', wind_speed_linestyles{s}, 'HandleVisibility', 'off')

            % x limits
            if phase == 1
                xlim([-0.01, 0.3])
            elseif phase == 3
                xlim([-0.3, 0.01])
                set(gca, 'YTickLabel', []);
            end

        end
    end
    hold off
    
    % Titles
    if phase == 1
        title('Wave Peak')
    elseif phase == 3
        title('Wave Trough')
    end
end

% Add legend for waves
hold on
for w = 1:length(waves)
    wave = waves{w};
    label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses.(wave));
    plot(nan, nan, 'color', wave_colors{w}, 'linewidth', 2, ...
         'displayname', label)
end

% Add legend for windspeeds
plot(nan, nan, 'color', 'white', 'displayname', ' ')
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    freestream = freestreams.(wind_speed);

    label = sprintf('$u_{\\infty} = %1.2f~m/s$', freestream);
    plot(nan, nan, 'color', 'black', 'linewidth', 2, ...
         'linestyle', wind_speed_linestyles{s}, 'displayname', label)
end
hold off
leg = legend('box', 'off', 'interpreter', 'latex', 'location', 'northwest');
leg.Layout.Tile = 'east';


linkaxes(h, 'y')
ylim([0, 1.1])
xlabel(tile, '$\tilde{u}_{\xi} \mathbin{/} u_{\infty}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(tile, '$y \mathbin{/} \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)




%% Save matfile to sent to manuel

% Add wave parameters and freestream velocities
output.wavelengths_mm = wavelengths;
output.amplitudes_mm = amplitudes;
output.frequencies_hz = frequencies;
output.freestreams = freestreams;


save_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/manuel';
filename = 'WindWave_ManuelProfiles.mat';
fprintf('Saving Matfile...\n')
pause(3)
save(fullfile(save_folder, filename), 'output');
clc; fprintf('Matfile saved!\n')


%% Plot for a fixed wind speed and wave across all wind speeds


% wave_transparency = 0.25;
% 
% % Fontsizes
% levels = 100;
% linewidth = 1;
% 
% colorbar_fontsize = 10;
% tickFontSize = 8;
% labelFontSize = 10;
% legendFontSize = 8;
% titleFontSize = labelFontSize;
% 
% % How much to chop off edges [mm]
% edge_buffer = 3;
% 
% % What to plot
% wind_speed = 'WT6';
% freestream = freestreams.(wind_speed);
% phase = 3;
% 
% % Plot
% clc; close all
% figure('color', 'white', 'Units', 'centimeters', 'position', [10,10,13,4])
% t = tiledlayout(1, length(waves), 'Padding', 'tight', 'TileSpacing', 'tight');
% 
% % Loop through waves
% for w = 1:length(waves)
%     wave = waves{w};
%     caze = strcat(wind_speed, '_WV', wave, '_AG0');
%     disp(caze)
% 
%     % Load data
%     tmp = wave_coherent.(caze)(phase).u / freestream;
%     vl = curvilinear.(caze).phase(phase).vertical_lines;
%     hl = curvilinear.(caze).phase(phase).horizontal_lines;
% 
%     % Trimp edges
%     X = curvilinear.(caze).X;
%     Y = curvilinear.(caze).Y;
%     tmp(vl < min(X, [], 'all') + edge_buffer) = nan;
%     tmp(vl > max(X, [], 'all') - edge_buffer) = nan;
%     tmp(hl > max(Y, [], 'all') + 2 * edge_buffer) = nan;
% 
%     % Trim first rows of real data right above surface
%     is_all_nans = all(isnan(tmp), 2);
%     is_not_all_nans = ~is_all_nans;
%     last_non_nan_row_index = find(is_not_all_nans, 1, 'last');
%     first_non_nan_row_index = find(is_not_all_nans, 1, 'first');
% 
%     buff = 1;
%     tmp(last_non_nan_row_index - buff:end, :) = nan;
%     tmp(1:first_non_nan_row_index + buff, :) = nan;
% 
%     % Plot tile
%     h(w) = nexttile;
%     set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
%     hold on
%     contourf(vl, hl, tmp, levels, 'linestyle', 'none')
%     plot(curvilinear.(caze).X(1,:), curvilinear.(caze).phase(phase).wave_profile, 'linewidth', linewidth, 'color', 'black')
% 
%     % Shaded region below wave
%     hFill = patch( ...
%     [curvilinear.(caze).X(1,:), fliplr(curvilinear.(caze).X(1,:))], ...
%     [curvilinear.(caze).phase(phase).wave_profile, -100 * ones(size(curvilinear.(caze).phase(phase).wave_profile))], ...
%     'k', ...
%     'FaceAlpha', wave_transparency, ...      
%     'EdgeColor', 'none', ...    
%     'HandleVisibility', 'off'); 
%     uistack(hFill, 'bottom')
% 
%     hold off
% 
%     axis equal
%     clim([-0.1, 0.1])
%     colormap(slanCM('bwr'))
%     label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses.(wave));
%     title(label, 'Interpreter', 'latex', 'fontsize', titleFontSize)
%     xticks(0:100:200)
%     yticks(0:100:200)
% 
%     if w ~= 1
%         ax = gca;
%         ax.YTickLabel = [];
%     end
%     if w == 4
%         C = colorbar;
%         C.Label.String = '$\tilde{u}_{\xi} \mathbin{/} u_{\infty}$';
%         C.FontSize = tickFontSize;
%         C.Label.Interpreter = 'latex';
%         C.Label.FontSize = labelFontSize;
%         C.TickLabelInterpreter = 'latex';
%         C.Ticks = -0.1:0.1:0.1;
%     end
% end
% 
% linkaxes(h, 'xy')
% ylim([-20, 200])
% 
% % xlim([-10, 250])
% xlim([0, max(curvilinear.(caze).X(1,:))])
% xlabel(t, '$x$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylabel(t, '$y$ [mm]', 'interpreter', 'latex', 'fontsize', labelFontSize)

