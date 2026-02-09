%%% Curvilinear Friction Velocity

addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
means_path = fullfile(project_path, 'means');
curvilinear_path = fullfile(project_path, 'curvilinear_new');
waves = {'A', 'B', 'C', 'D'};
wind_speeds = {'WT4', 'WT6', 'WT8'};

wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};

figure_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/pdf_test6';

% Approximate wavelengths in mm for labeling plots
wavelength_names.A = '410';
wavelength_names.B = '313';
wavelength_names.C = '189';
wavelength_names.D = '124';

steepnesses.A = '0.180';
steepnesses.B = '0.211';
steepnesses.C = '0.305';
steepnesses.D = '0.267';

% Wind speed to consider
for s = 1:length(wind_speeds)

    % Loop through wind speeds
    wind_speed = wind_speeds{s};
 
    if ismember(wind_speed(end), {'4'})
        u_inf = 2.4181;
    elseif ismember(wind_speed(end), {'6'})
        u_inf = 3.8709;
    elseif ismember(wind_speed(end), {'8'})
        u_inf = 5.4289;
    end
    
    % Loop through waves
    for w = 1:length(waves)
    
        % Define case
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        disp(caze)
    
        % Load data
        curvilinear = load(fullfile(curvilinear_path, strcat(caze, '_CURVILINEAR.mat')));
        curvilinear = curvilinear.output;
    
        % Wave parameters
        wave_type       = caze(strfind(caze, 'WV') + 2);
        wave_parameters = readcell("Offshore_Waves.xlsx");
        wavelength      = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 2};
        amplitude       = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 3};

        wavelengths.(wave) = wavelength * 1E-3;
        amplitudes.(wave) = amplitude * 1E-3;
    
        for phase = 1:4
        
            % Get components
            uv = curvilinear.phase(phase).uv;
            % X (xi) in mm
            vertical_lines = curvilinear.phase(phase).vertical_lines * 1E-3;
            % Y (zeta) in mm
            horizontal_lines = curvilinear.phase(phase).horizontal_lines * 1E-3;
            % Waves
            wave_profile = curvilinear.phase(phase).wave_profile * 1E-3;
            max_wave_profile = curvilinear.phase(phase).max_wave_profile * 1E-3;

            %%% Save wave_profiles
            wave_profiles.(wave)(phase).profile = wave_profile;
            wave_profiles.(wave)(phase).x = curvilinear.X(1,:) * 1E-3;
            
            % Get profiles
            idx = round(length(uv)/2);
            uv_profile = uv(:, idx); 

            % Save profiles
            normalized_uv_profiles.(caze)(phase).uv = uv_profile / (u_inf^2);
            normalized_uv_profiles.(caze)(phase).y = horizontal_lines(:, idx) - wave_profile(idx);

        
            %%% Get Friction Velocity
            friction_velocity.(caze)(phase) = max(sqrt(-uv_profile), [], 'all', 'omitnan');
            skin_friction.(caze)(phase) = 2 * (friction_velocity.(caze)(phase) / u_inf)^2;

            %%% Get Skin Friciton profiles
            for i = 1:size(uv,2)
                column = uv(:,i);
                if all(isnan(column))
                    continue;
                end
        
                [M, I] = max(sqrt(-column), [], 'omitnan');
                cfs_tmp(i) = 2 * (M / u_inf).^2;
                x_tmp(i) = vertical_lines(I, i);
            end
        
        
            % Save values
            skin_friction_profiles.(caze)(phase).x = x_tmp;
            skin_friction_profiles.(caze)(phase).cf = cfs_tmp;
        
        end
    end
end

% Also include the no-wave cases
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    caze = strcat(wind_speed, '_WV0_AGP');

    % Loop through wind speeds
    wind_speed = wind_speeds{s};

    if ismember(wind_speed(end), {'4'})
        u_inf = 2.4181;
    elseif ismember(wind_speed(end), {'6'})
        u_inf = 3.8709;
    elseif ismember(wind_speed(end), {'8'})
        u_inf = 5.4289;
    end

    % Load data
    means = load(fullfile(means_path, strcat(caze, '_MEANS.mat')));
    means = means.output;

    % Get friction velocity and cf
    uv = means.ensemble.uv;

    % Get profiles
    idx = round(length(uv)/2);
    uv_profile = uv(:, idx); 

    % Save profiles
    normalized_uv_profiles.(caze).uv = uv_profile / (u_inf^2);
    normalized_uv_profiles.(caze).y = means.Y(:,1) * 1E-3;

    %%% Get Friction Velocity
    friction_velocity.(caze) = max(sqrt(-uv_profile), [], 'all', 'omitnan');
    skin_friction.(caze) = 2 * (friction_velocity.(caze) / u_inf)^2;
end

clc;
clear caze curvilinear horizontal_lines idx phase uv uv_profile 
clear vertical_lines w wave wave_profile z_profile wave_type s wind_speed 
clear curvilinear_path max_wave_profile means means_path project_path u_inf
clear wave_parameters wavelength amplitude cfs_tmp column i I M x_tmp


%% Plot shear stress profiles used to get u*
% rows ~ fixed wind speeds
% cols ~ fixed phases

lw = 2; 
ax = figure();
tiledlayout(length(wind_speeds), 4)

c = 1;
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    for phase = 1:4
        h(c) = nexttile;
        hold on
        title(sprintf("%s Phase %1.0f", wind_speed, phase))
        for w = 1:length(waves)
            wave = waves{w};
            caze = strcat(wind_speed, '_WV', wave, '_AG0');

            % Phase
            plot(normalized_uv_profiles.(caze)(phase).uv, normalized_uv_profiles.(caze)(phase).y, ...
                 'linewidth', lw, 'DisplayName', wave)
        end

        % No-Wave
        no_wave_caze = strcat(wind_speed, '_WV0_AGP');
        plot(normalized_uv_profiles.(no_wave_caze).uv, normalized_uv_profiles.(no_wave_caze).y, ...
             'color', 'black', 'linewidth', lw, 'DisplayName', 'No Waves')

        c = c + 1;
        hold off
    end
end

linkaxes(h, 'xy')
leg = legend('Orientation', 'Horizontal');
leg.Layout.Tile = 'north';

clear c s phase w wave wind_speed h ax caze no_wave_caze leg lw 





%% Plot Cf profiles

tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;
titleFontSize = 10;

phase = 2;
wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
linestyles = {'-.', '--', '-'};
lw = 1;
buff = 5;
scale = 1.25;

smoothin_kernel = 1;

scale = 0.55 * 1E3;
offset = -0.002 * 1E3;
alpha = 0.25;

close all; 
ax = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,6]);
t = tiledlayout(1,1, 'Padding','loose');
ax = nexttile;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)

hold on
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};

    if ismember(wind_speed(end), {'4'})
        u_inf = 2.4181;
    elseif ismember(wind_speed(end), {'6'})
        u_inf = 3.8709;
    elseif ismember(wind_speed(end), {'8'})
        u_inf = 5.4289;
    end

    % Legend shit
    if s == 3
        vis = 'on';
    else
        vis = 'off';
    end

    for w = 1:length(waves)
        wave = waves{w};

        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        disp(caze)

        % Crop ends of data
        data = 1E3 * skin_friction_profiles.(caze)(phase).cf(buff:end-buff);

        % Smooth data
        data = smoothdata(data, 'movmean', smoothin_kernel);
        data = filloutliers(data, 'pchip', 'gesd');

        % Sgolay filt
        % data = sgolayfilt(data, 2, 7);

        wavelength = wavelengths.(wave); 
        amplitude = amplitudes.(wave); 
        steepness = (2 * amplitude) / wavelength;

        x = skin_friction_profiles.(caze)(phase).x;
        label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses.(wave));

        plot((x(buff:end-buff) - mean(x, 'all', 'omitnan')) / wavelengths.(wave), data, ...
             'color', wave_colors{w}, 'linestyle', linestyles{s}, 'linewidth', lw, 'DisplayName', label, 'HandleVisibility', vis)

        % Plot phase for reference
        if wave == 'D' && s == 1
            testX = (wave_profiles.(wave)(phase).x - mean(wave_profiles.(wave)(phase).x)) / wavelengths.(wave);
            testY = scale *  wave_profiles.(wave)(phase).profile + offset;
            extrapX = -1:0.001:1;
            extrap_testY = interp1(testX, testY, extrapX, 'spline', 'extrap');
            plot(extrapX, extrap_testY, ...
                 'color', 'black', 'linewidth', lw, 'linestyle', '-', 'HandleVisibility', 'off');
            tmpX = extrapX;
            tmpY = extrap_testY;
            X = [tmpX(:); flipud(tmpX(:))];
            Y = [tmpY(:); (-1000)*ones(size(tmpY(:)))];
            
            patch(X, Y, 'k', ...
                  'FaceAlpha', alpha, ...
                  'EdgeColor', 'none', ...
                  'handlevisibility', 'off');
        end
    end
end

% % Dummy line to add a break
% plot(nan, nan, 'color', 'white', 'displayname', '')
% 
% % Dummy lines to get wind-speed legend
% for s = 1:length(wind_speeds)
%     wind_speed = wind_speeds{s};
% 
%     if ismember(wind_speed(end), {'4'})
%         u_inf = 2.4181;
%     elseif ismember(wind_speed(end), {'6'})
%         u_inf = 3.8709;
%     elseif ismember(wind_speed(end), {'8'})
%         u_inf = 5.4289;
%     end
%     P = plot(nan, nan, 'linestyle', linestyles{s}, ...
%          'linewidth', lw, 'color', 'black', 'DisplayName', sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf));
% end

hold off

% Make y axis easier to read
% ax = gca;
% ax.YAxis.Exponent = -3;          % shows ×10^{-6} once
% ax.YAxis.TickLabelFormat = '%2.0f';  % two digits on ticks

xticks(-1:0.25:1)
xlim([-1, 1])
% ylim([0, 0.02])
ylim([-0.007 * 1E3, 0.02 * 1E3])
yticks(0:0.005 * 1E3 :0.02 * 1E3)
leg = legend('interpreter', 'latex', 'fontsize', legendFontSize, 'orientation',  'horizontal', 'box', 'off');
leg.Layout.Tile = 'north';
leg.IconColumnWidth = 19;


xlabel('$\xi \mathbin{/} \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$C_f \times 10^{3}$', 'interpreter', 'latex', 'fontsize', labelFontSize)

clc;


% Save figure
% pause(3)
% figure_name = strcat('SkinFriction_Curvilinear_Phase', num2str(phase), '.pdf');
% exportgraphics(ax, fullfile(figure_folder, 'Friction', figure_name), 'Resolution', 600, 'ContentType', 'image'); 
% close all





%% Plot average of phases Cf against wave steepness

% buff = 5;
% smoothin_kernel = 4;
% 
% tickFontSize = 8;
% labelFontSize = 16;
% legendFontSize = 12;
% 
% wind_speed_markers = {'^', 'square', 'o'};
% close all; clc;
% 
% ax = figure('color', 'white', 'units', 'centimeters', 'position', [2,4,30,10]);
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% hold on
% 
% % Add wave markers
% for w = 1:4
%     wave = waves{w};
%     label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelength_names.(wave), steepnesses.(wave));
%     scatter(nan, nan, 100, 'filled', 'MarkerFaceColor', wave_colors{w}, 'displayname', label)
% end
% 
% % Legend space
% plot(nan, nan, 'color', 'white', 'displayname', '')
% 
% for s = 1:length(wind_speeds)
%     wind_speed = wind_speeds{s};
% 
%     if ismember(wind_speed(end), {'4'})
%         u_inf = 2.4181;
%     elseif ismember(wind_speed(end), {'6'})
%         u_inf = 3.8709;
%     elseif ismember(wind_speed(end), {'8'})
%         u_inf = 5.4289;
%     end
% 
%     for w = 1:length(waves)
% 
%         fprintf('%s\n', wind_speed)
%         wave = waves{w};
%         wavelength = wavelengths.(wave); 
%         amplitude = amplitudes.(wave); 
%         steepness = (2 * pi * amplitude) / wavelength;
%         disp(steepness)
% 
% 
%         % Save steepnesses for line plot
%         steepnessez(w) = steepness;
% 
%         % Plot max and min Cf values from phase 2
%         caze = strcat(wind_speed, '_WV', wave, '_AG0');
%         disp(caze)
% 
%         % Crop ends of data and smooth
%         Cf_profile = skin_friction_profiles.(caze)(2).cf(buff:end-buff);
%         Cf_profile = smoothdata(Cf_profile, 'movmean', smoothin_kernel);
% 
% 
%         scatter(steepness, max(Cf_profile, [], 'all', 'omitnan'), 50, wind_speed_markers{s}, 'filled', ...
%                 'MarkerFaceColor', wave_colors{w}, 'handlevisibility', 'off');
% 
%         scatter(steepness, min(Cf_profile, [], 'all', 'omitnan'), 50, wind_speed_markers{s}, 'filled', ...
%                 'MarkerFaceColor', wave_colors{w}, 'handlevisibility', 'off');
% 
% 
% 
%         % Compute average of phases quantity
%         for phase = 1:4
%             caze = strcat(wind_speed, '_WV', wave, '_AG0');
%             u_star_data(phase) = friction_velocity.(caze)(phase);
%             cf_data(phase) = skin_friction.(caze)(phase);
% 
%             % scatter(steepness, skin_friction.(caze)(phase), 50, 'filled', 'markerfacecolor', wave_colors{w})
%         end
% 
%         average_of_phases_Cf = mean(cf_data, 'all', 'omitnan');
% 
%         % Save average of phases for line plot
%         if w == 1
%             vis = 'on';
%         else
%             vis = 'off';
%         end
% 
%         label = sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf);
% 
%         pseudo_ensemble(w) = average_of_phases_Cf;
% 
%         PhaseMean = scatter(steepness, average_of_phases_Cf, 100, wind_speed_markers{s}, 'filled', ...
%                 'markerfacecolor', 'black', 'displayname', label, 'handlevisibility', vis);
% 
%         uistack(PhaseMean, 'top');
% 
%         clear u_star_data cf_data
%     end
% 
%     [sorted_steepness, sort_indices] = sort(steepnessez);
%     plot(sorted_steepness, pseudo_ensemble(sort_indices), 'color', 'black', 'linewidth', 3, 'handlevisibility', 'off')
% 
% 
% end
% 
% hold off
% clear u_star_data cf_data wind_speed b caze i phase temp_uv temp_y w wave
% 
% legend('interpreter', 'latex', 'fontsize', legendFontSize, 'location',  'eastoutside', 'box', 'off')
% xlabel('$ak$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylabel('$\overline{{C_f}_{\phi}}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% 
% xlim([0.17, 0.32])
% ylim([0, 16E-3])
% yticks(0:0.005:0.02)
% 
% 
% % Name figure
% pause(3)
% figure_name = 'SkinFriction_CurvilinearCf_vs_Steepness.pdf';
% exportgraphics(ax, fullfile(figure_folder, 'pdf_test2', 'Friction', figure_name), 'Resolution', 600, 'ContentType', 'image'); 
% close all












%% Plot average of phases Cf against wave steepness: all wind speeds

% wind_speed_shapes = {'o', '^', 'square'};
% 
% close all; figure('color', 'white')
% hold on
% 
% for s = 1:length(wind_speeds)
%     wind_speed = wind_speeds{s};
% 
%     for w = 1:length(waves)
% 
%         wave = waves{w};
%         wavelength = wavelengths.(wave); 
%         amplitude = amplitudes.(wave); 
%         steepness = (2 * amplitude) / wavelength;
%         disp(steepness)
% 
%         for phase = 1:4
%             caze = strcat(wind_speed, '_WV', wave, '_AG0');
%             u_star_data(phase) = friction_velocity.(caze)(phase);
%             cf_data(phase) = skin_friction.(caze)(phase);
% 
%             scatter(steepness, skin_friction.(caze)(phase), 100, 'filled', 'markerfacecolor', wave_colors{w})
%         end
% 
%         scatter(steepness, mean(cf_data, 'all', 'omitnan'), 200, 'filled', wind_speed_shapes{s}, 'markerfacecolor', 'black')
% 
%         clear u_star_data cf_data
%     end
% end
% hold off




%% Functions

function output = StatisticalGradientFilter(q, tol, max_island_width)
    % Take gradient
    grad = gradient(q);
    % Get standard deviation
    sd = std(grad, 0, 'omitnan');
    % Get average
    avg = mean(grad, 'all', 'omitnan');
    % Logical indexing
    q(grad > avg + tol * sd | grad < avg - tol * sd) = nan;
    % Remove Islands
    q = RemoveIslands(q,max_island_width);
    output = q;
end


% Chat GPT
function cleaned = RemoveIslands(signal, max_island_width)

    % Create logical mask of NaNs
    is_nan = isnan(signal);

    % Find islands of valid data (i.e., ~isnan) between NaNs
    cleaned = signal; 
    n = length(signal);
    i = 1;
    while i <= n
        % Skip NaNs
        if is_nan(i)
            i = i + 1;
            continue;
        end

        % Start of a non-NaN segment
        start_idx = i;
        while i <= n && ~is_nan(i)
            i = i + 1;
        end
        end_idx = i - 1;

        % Check if island is surrounded by NaNs
        is_surrounded = false;
        if start_idx > 1 && end_idx < n
            if is_nan(start_idx - 1) && is_nan(end_idx + 1)
                is_surrounded = true;
            end
        end

        % Remove if surrounded and short
        if is_surrounded && (end_idx - start_idx + 1) <= max_island_width
            cleaned(start_idx:end_idx) = NaN;
        end
    end
end
