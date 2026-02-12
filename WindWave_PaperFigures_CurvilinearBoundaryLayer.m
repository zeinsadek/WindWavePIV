%%% Curvilinear Integral Boundary Layer Parameters

addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
curvilinear_path = fullfile(project_path, 'curvilinear_new');
wave_parameters = readcell("Offshore_Waves.xlsx");

figure_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/pdf_test6';


% Boundary layer detection parameters
boundary_layer_percent = 0.96;
left_mask = 4E-3;
right_mask = 234E-3;

% Approximate wavelengths in mm for labeling plots
wavelengths.A = '410';
wavelengths.B = '313';
wavelengths.C = '189';
wavelengths.D = '124';

steepnesses.A = '0.180';
steepnesses.B = '0.211';
steepnesses.C = '0.305';
steepnesses.D = '0.267';

frequencies.A = 1.96;
frequencies.B = 2.27;
frequencies.C = 3;
frequencies.D = 3.93;

% Edge killing values
% OG was 0.008
threshold.thickness = 0.005;
threshold.displacement = 0.005;
threshold.momentum = 0.008;

% Data filtering and smoothing
smooth_kernel = 9;
gradient_tolerance = 1;
island_length = 5;
sgolay_length = 35;
n = 9;


% Loop through all waves and wind speeds
wind_speeds = {'WT4', 'WT6', 'WT8'};
waves = {'A', 'B', 'C', 'D'};

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


    for w = 1:length(waves)
    
        % Define case
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        disp(caze)
    
        % Load data
        curvilinear = load(fullfile(curvilinear_path, strcat(caze, '_CURVILINEAR.mat')));
        curvilinear = curvilinear.output;
        
        % Wave Parameters
        wave_type  = wave;
        wavelength = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 2};
        amplitude  = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 3};
    
        % Save wavelength and amplitude
        integral.(caze).wavelength = wavelength * 1E-3;
        integral.(caze).amplitude = amplitude * 1E-3;
    
    
        %%% Boundary Layer Edge Detection
        for phase = 1:4
    
            fprintf("Phase %1.0f\n", phase)
        
            % Curvilinear coordinates
            vertical_lines = curvilinear.phase(phase).vertical_lines * 1E-3;
            horizontal_lines = curvilinear.phase(phase).horizontal_lines * 1E-3;
            
            % Waves
            % wave_profile = curvilinear.phase(phase).wave_profile * 1E-3;
            % max_wave_profile = curvilinear.phase(phase).max_wave_profile * 1E-3;
            
            % Find boundary layer edge
            u = curvilinear.phase(phase).u;
            u = juliaan_smooth(u, smooth_kernel);
        
            % Mask the data properly
            u(vertical_lines < left_mask | vertical_lines > right_mask) = nan;
            vertical_lines(vertical_lines < left_mask | vertical_lines > right_mask) = nan;
            horizontal_lines(vertical_lines < left_mask | vertical_lines > right_mask) = nan;
            
            % Find BL
            u(u > u_inf * boundary_layer_percent) = nan;
           
            %%% Boundary layer thickness
            thickness = nan(1, size(u,2));
            x_BL = nan(1, size(u,2));
        
            for i = 1:size(u,2)
                column = u(:,i);
                if all(isnan(column))
                    continue;
                end
        
                [~, idx] = max(column, [], 'omitnan');
                thickness(i) = horizontal_lines(idx, i);
                x_BL(i) = vertical_lines(idx, i);
            end
        
            % Kil bad edge values
            thickness_edge_mask = gradient(thickness) > threshold.thickness | gradient(thickness) < -threshold.thickness;
            thickness(thickness_edge_mask) = nan;
            x_BL(thickness_edge_mask) = nan;
        
        
    
            %%% Momentum and displacement thicknesses
            displacement = nan(1, size(u,2));
            momentum = nan(1, size(u,2));
            
            for i = 1:size(u,2)
            
                u_column = u(:,i);
                zeta_column = horizontal_lines(:,i);
            
                % Remove NaNs
                nan_mask = ~isnan(u_column) & ~isnan(zeta_column);
                u_column = u_column(nan_mask);
                zeta_column = zeta_column(nan_mask);
            
                % not enough points to integrate
                if length(u_column) < 2
                    continue; 
                end
            
                % Normalize velocity
                u_normalized = u_column / u_inf;
                    
                % Decreasing ζ
                if zeta_column(end) < zeta_column(1) 
                    zeta_column = flipud(zeta_column);
                    u_normalized = flipud(u_normalized);
                end
            
                % Integrands
                displacement_integrand = (1 - u_normalized);
                momentum_integrand = u_normalized .* (1 - u_normalized);
        
                % Integrals
                displacement_cum = trapz(zeta_column, displacement_integrand);
                momentum_cum = trapz(zeta_column, momentum_integrand);
            
                % Final value = full integral over zeta
                displacement(i) = displacement_cum;
                momentum(i) = momentum_cum;
        
            end
        
            % Kill bad edge values
            displacement_edge_mask = gradient(displacement) > threshold.displacement | gradient(displacement) < -threshold.displacement;
            momentum_edge_mask = gradient(momentum) > threshold.momentum | gradient(momentum) < -threshold.momentum;
        
            displacement(displacement_edge_mask) = nan;
            momentum(momentum_edge_mask) = nan;
        
    
    
            %%% Save signals
            % Save reference wave profile
            integral.(caze)(phase).wave.wave_profile = curvilinear.phase(phase).wave_profile * 1E-3;
            integral.(caze)(phase).wave.x = curvilinear.X(1,:) * 1E-3;
             
            % Save raw signals
            integral.(caze)(phase).raw.thickness = thickness;
            integral.(caze)(phase).raw.displacement = displacement;
            integral.(caze)(phase).raw.momentum = momentum;
            integral.(caze)(phase).raw.shape = displacement ./ momentum;
            integral.(caze)(phase).x = x_BL;
        
            % Save filtered
            integral.(caze)(phase).filtered.thickness = FilterData(x_BL, thickness, 10, 10, sgolay_length);
            integral.(caze)(phase).filtered.displacement = FilterData(x_BL, displacement, gradient_tolerance, island_length, sgolay_length);
            integral.(caze)(phase).filtered.momentum = FilterData(x_BL, momentum, gradient_tolerance, island_length, sgolay_length);
            integral.(caze)(phase).filtered.shape = FilterData(x_BL, displacement ./ momentum, gradient_tolerance, island_length, sgolay_length);
        
            % Save polynomial fit
            integral.(caze)(phase).fitted.thickness = PolynomialFit(x_BL, integral.(caze)(phase).filtered.thickness, n);
            integral.(caze)(phase).fitted.displacement = PolynomialFit(x_BL, integral.(caze)(phase).filtered.displacement, n);
            integral.(caze)(phase).fitted.momentum = PolynomialFit(x_BL, integral.(caze)(phase).filtered.momentum, n);
            integral.(caze)(phase).fitted.shape = PolynomialFit(x_BL, integral.(caze)(phase).filtered.shape, n);
        
        end
    end
end

clear amplitude caze column curvilinear displacement displacement_cum displacement_edge_mask displacement_integrand
clear horizontal_lines i idx island_length left_mask max_wave_profile momentum momentum_cum momentum_edge_mask momentum_integrand
clear n nan_mask phase right_mask s sgolay_length smooth_kernerl thickness thickness_edge_mask threshold u u_column u_inf u_normalized
clear vertical_lines w wave wave_profile wave_parameters wave_type wave_length wind_speed x_BL zeta_column
clear gradient_tolerance project_path smooth_kernel wavelength curvilinear_path
clc;



%%%
%%% Plot all quantitues by phase and wind speed
%%%

% Quantities to be looped
parameters = {'thickness', 'displacement', 'momentum', 'shape'};
wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
% versions = {'filtered', 'fitted'};
versions = {'filtered'};

% version_colors = {'red', 'blue'};
lw = 2;

wind_speed = 'WT4';
ylim_max = [0.18, 0.035, 0.025, 2.0];
% ylim_min = [0.075, 0, 0, 1];
alpha = 0.5;
scales = [2, 1, 1, 20];


clc; close all
figure('color', 'white')
tiledlayout(length(waves), length(parameters));
sgtitle(wind_speed)

count = 1;
% Interate through waves
for p = 1:4
    
    % Iterate through integral parameters
    for i = 1:length(parameters)
        h(count) = nexttile;
        hold on
        
        % Iterate through different filters
        for w = 1:length(waves)
            wave = waves{w};
            caze = strcat(wind_speed, '_WV', wave, '_AG0');
            disp(caze)

            % Plot data
            for v = 1:length(versions)
                data = integral.(caze)(p).(versions{v}).(parameters{i});
                H(w) = plot((integral.(caze)(p).x - mean(integral.(caze)(p).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, data, ...
                            'color', wave_colors{w}, 'linewidth', lw, ...
                            'displayname', sprintf("Wave %s", wave));
            end

            % Plot phase for reference
            if wave == 'D'
                disp(scales(i))
                reference = plot((integral.(caze)(p).wave.x  - mean(integral.(caze)(p).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, (scales(i) * integral.(caze)(p).wave.wave_profile) + mean(data, 'all', 'omitnan'), ...
                           'color', 'black', 'linewidth', lw, 'linestyle', '--');
                reference.Color(4) = alpha;
            end
       end
       count = count + 1;
       hold off
       ylim([0, ylim_max(i)])
       title(sprintf("%s: Phase %1.0f", parameters{i}, p))

    end
end

leg = legend(H, 'Orientation', 'Horizontal');
leg.Layout.Tile = 'north';
set(leg,'Box','off')

% Share axes for same parameter across different plots
for i = 1:4
    linkaxes(h(i:4:end), 'xy')
end
clc;

clear alpha caze count data h H i leg lw p v reference scales c version_colors verrions w wave wind_speed ylim_max
clc;





%% Plot specific parameter for a specific phase for all waves and wind speeds

phase = 1;
parameter = 'thickness';
version = 'filtered';

% Plot specs
linestyles = {'-.', '--', '-'};
lw = 1;
alpha = 0.25;
tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;

% Scale reference wave based on parameter
if strcmp(parameter, 'thickness')
    scale = 1.5;
elseif strcmp(parameter, 'displacement')
    scale = 1;
elseif strcmp(parameter, 'momentum')
    scale = 1;
elseif strcmp(parameter, 'shape')
    scale = 15;
end

% Offset based on parameter
if strcmp(parameter, 'thickness')
    offset = 0.072;
elseif strcmp(parameter, 'displacement')
    offset = -0.005;
elseif strcmp(parameter, 'momentum')
    offset = -0.005;
elseif strcmp(parameter, 'shape')
    offset = 0.95;
end

% Plot
clc; close all; 
ax = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,6]);
t = tiledlayout(1,1,'Padding','loose');
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

    % Legend
    if s == 3
        vis = 'on';
    else
        vis = 'off';
    end

    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        wavelength = integral.(caze).wavelength; 
        amplitude = integral.(caze).amplitude;
        steepness = (2 * amplitude) / wavelength;
        
        label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelengths.(wave), steepnesses.(wave));
        data = integral.(caze)(phase).(version).(parameter);
        H(w) = plot((integral.(caze)(phase).x - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, data, ...
                    'color', wave_colors{w}, 'linewidth', lw, 'linestyle', linestyles{s}, ...
                    'displayname', label, 'HandleVisibility', vis);

        % Plot phase for reference
        if wave == 'D' && s == 1
            testX = (integral.(caze)(phase).wave.x  - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength;integral.(caze)(phase).wave.wave_profile;
            testY = (scale * integral.(caze)(phase).wave.wave_profile) + offset;
            extrapX = -1:0.001:1;
            extrap_testY = interp1(testX, testY, extrapX, 'spline', 'extrap');
            plot(extrapX, extrap_testY, ...
                 'color', 'black', 'linewidth', lw, 'linestyle', '-', 'HandleVisibility', 'off');
            tmpX = extrapX;
            tmpY = extrap_testY;
            X = [tmpX(:); flipud(tmpX(:))];
            Y = [tmpY(:); (-1)*ones(size(tmpY(:)))];
            
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


% Axes labels
if strcmp(parameter, 'thickness')
    vertLabel = '$\delta$ [m]'; %#ok<*UNRCH>
    ylim([0.055, 0.14])
    yticks(0.08:0.02:0.16)
elseif strcmp(parameter, 'displacement')
    vertLabel = '$\delta^*$ [m]';
    ylim([-0.015, 0.04])
elseif strcmp(parameter, 'momentum')
    vertLabel = '$\theta$ [m]';
    ylim([-0.015, 0.04])
elseif strcmp(parameter, 'shape')
    vertLabel = '$H$';
    ylim([0.8, 1.7])
end

% Legend + Labels
lgd = legend('interpreter', 'latex', 'fontsize', legendFontSize, 'orientation',  'horizontal', 'box', 'off');
lgd.Layout.Tile = 'north';
lgd.IconColumnWidth = 19;
xlabel('$\xi \mathbin{/} \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(vertLabel, 'interpreter', 'latex', 'fontsize', labelFontSize)
xticks(-1:0.25:1)


% Save figure
% pause(3)
% figure_name = strcat('BoundaryLayer_Curvilinear', [upper(parameter(1)), lower(parameter(2:end))], '_Phase', num2str(phase), '.pdf');
% exportgraphics(ax, fullfile(figure_folder, 'BoundaryLayer', figure_name), 'Resolution', 600, 'ContentType', 'image');
% close all
% clc; fprintf('Generated figure: %s\n\n', figure_name)

clear alpha limits caze data H linestyles lw parameter phase reference s scale w wave wind_speed version ans vertLabel; clc




%% Plot displacement thickness + momentum thickness + shape factor + collapse parameter


phase = 1;
version = 'filtered';
linestyles = {'-.', '--', '-'};
lw = 1;
alpha = 0.25;

tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;

clc; close all; 
ax = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,8]);
tile = tiledlayout(2,2, 'tilespacing', 'loose', 'padding', 'loose');

% Displacement
parameter = 'displacement';

% Scale reference wave based on parameter
if strcmp(parameter, 'thickness')
    scale = 2;
elseif strcmp(parameter, 'displacement')
    scale = 1;
elseif strcmp(parameter, 'momentum')
    scale = 1;
elseif strcmp(parameter, 'shape')
    scale = 15;
end

% Offset based on parameter
if strcmp(parameter, 'thickness')
    offset = 0.06;
elseif strcmp(parameter, 'displacement')
    offset = -0.005;
elseif strcmp(parameter, 'momentum')
    offset = -0.005;
elseif strcmp(parameter, 'shape')
    offset = 0.95;
end

% Plot displacement thickness
h(1) = nexttile;
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

    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        wavelength = integral.(caze).wavelength; 
        amplitude = integral.(caze).amplitude;
        steepness = (2 * amplitude) / wavelength;
        
        label = sprintf('$\\lambda_{%s}, u_{\\infty} = %1.2f$', wavelengths.(wave), u_inf);
        
        data = integral.(caze)(phase).(version).(parameter);
        H(w) = plot((integral.(caze)(phase).x - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, data, ...
                    'color', wave_colors{w}, 'linewidth', lw, 'linestyle', linestyles{s}, ...
                    'displayname', label);

        % Plot phase for reference
        if wave == 'D' && s == 1
            testX = (integral.(caze)(phase).wave.x  - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength;integral.(caze)(phase).wave.wave_profile;
            testY = (scale * integral.(caze)(phase).wave.wave_profile) + offset;
            extrapX = -1:0.001:1;
            extrap_testY = interp1(testX, testY, extrapX, 'spline', 'extrap');
            plot(extrapX, extrap_testY, ...
                 'color', 'black', 'linewidth', lw, 'linestyle', '-', 'HandleVisibility', 'off');
            tmpX = extrapX;
            tmpY = extrap_testY;
            X = [tmpX(:); flipud(tmpX(:))];
            Y = [tmpY(:); (-1)*ones(size(tmpY(:)))];
            
            patch(X, Y, 'k', ...
                  'FaceAlpha', alpha, ...
                  'EdgeColor', 'none', ...
                  'handlevisibility', 'off');
        end
    end
end
hold off



vertLabel = '$\delta^*$ [m]';
ylim([-0.015, 0.035])

% xlabel('$\xi / \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(vertLabel, 'interpreter', 'latex', 'fontsize', labelFontSize)
xticks(-1:0.5:1)
yticks(0:0.01:0.04)







%%% Momentum
parameter = 'momentum';

% Scale reference wave based on parameter
if strcmp(parameter, 'thickness')
    scale = 2;
elseif strcmp(parameter, 'displacement')
    scale = 1;
elseif strcmp(parameter, 'momentum')
    scale = 1;
elseif strcmp(parameter, 'shape')
    scale = 15;
end

% Offset based on parameter
if strcmp(parameter, 'thickness')
    offset = 0.06;
elseif strcmp(parameter, 'displacement')
    offset = -0.005;
elseif strcmp(parameter, 'momentum')
    offset = -0.005;
elseif strcmp(parameter, 'shape')
    offset = 0.95;
end

h(2) = nexttile;
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

    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        wavelength = integral.(caze).wavelength; 
        amplitude = integral.(caze).amplitude;
        steepness = (2 * amplitude) / wavelength;
        
        label = sprintf('$\\lambda_{%s}, u_{\\infty} = %1.2f$', wavelengths.(wave), u_inf);
        data = integral.(caze)(phase).(version).(parameter);
        H(w) = plot((integral.(caze)(phase).x - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, data, ...
                    'color', wave_colors{w}, 'linewidth', lw, 'linestyle', linestyles{s}, ...
                    'displayname', label);

        % Plot phase for reference
        if wave == 'D' && s == 1
            testX = (integral.(caze)(phase).wave.x  - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength;integral.(caze)(phase).wave.wave_profile;
            testY = (scale * integral.(caze)(phase).wave.wave_profile) + offset;
            extrapX = -1:0.001:1;
            extrap_testY = interp1(testX, testY, extrapX, 'spline', 'extrap');
            plot(extrapX, extrap_testY, ...
                 'color', 'black', 'linewidth', lw, 'linestyle', '-', 'HandleVisibility', 'off');
            tmpX = extrapX;
            tmpY = extrap_testY;
            X = [tmpX(:); flipud(tmpX(:))];
            Y = [tmpY(:); (-1)*ones(size(tmpY(:)))];
            
            patch(X, Y, 'k', ...
                  'FaceAlpha', alpha, ...
                  'EdgeColor', 'none', ...
                  'handlevisibility', 'off');
        end
    end
end
hold off

% ax = gca;
% ax.YAxis.Exponent = -3;
% ax.YAxis.TickLabelFormat = '%2.0f';


vertLabel = '$\theta$ [m]';
ylim([-0.015, 0.035])


% xlabel('$\xi / \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(vertLabel, 'interpreter', 'latex', 'fontsize', labelFontSize)
xticks(-1:0.5:1)
yticks(0:0.01:0.04)
linkaxes(h, 'xy')






%%% Shape 
parameter = 'shape';

% Scale reference wave based on parameter
if strcmp(parameter, 'thickness')
    scale = 2;
elseif strcmp(parameter, 'displacement')
    scale = 1;
elseif strcmp(parameter, 'momentum')
    scale = 1;
elseif strcmp(parameter, 'shape')
    scale = 15;
end

% Offset based on parameter
if strcmp(parameter, 'thickness')
    offset = 0.06;
elseif strcmp(parameter, 'displacement')
    offset = -0.005;
elseif strcmp(parameter, 'momentum')
    offset = -0.005;
elseif strcmp(parameter, 'shape')
    offset = 0.92;
end

% Plot displacement thickness
h(3) = nexttile;
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

    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        wavelength = integral.(caze).wavelength; 
        amplitude = integral.(caze).amplitude;
        steepness = (2 * amplitude) / wavelength;
        
        label = sprintf('$\\lambda_{%s}, u_{\\infty} = %1.2f$', wavelengths.(wave), u_inf);
        
        data = integral.(caze)(phase).(version).(parameter);
        H(w) = plot((integral.(caze)(phase).x - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, data, ...
                    'color', wave_colors{w}, 'linewidth', lw, 'linestyle', linestyles{s}, ...
                    'displayname', label);

        % Plot phase for reference
        if wave == 'D' && s == 1
            testX = (integral.(caze)(phase).wave.x  - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength;integral.(caze)(phase).wave.wave_profile;
            testY = (scale * integral.(caze)(phase).wave.wave_profile) + offset;
            extrapX = -1:0.001:1;
            extrap_testY = interp1(testX, testY, extrapX, 'spline', 'extrap');
            plot(extrapX, extrap_testY, ...
                 'color', 'black', 'linewidth', lw, 'linestyle', '-', 'HandleVisibility', 'off');
            tmpX = extrapX;
            tmpY = extrap_testY;
            X = [tmpX(:); flipud(tmpX(:))];
            Y = [tmpY(:); (-1)*ones(size(tmpY(:)))];
            
            patch(X, Y, 'k', ...
                  'FaceAlpha', alpha, ...
                  'EdgeColor', 'none', ...
                  'handlevisibility', 'off');
        end
    end
end
hold off

% ax = gca;
% ax.YAxis.Exponent = -3;
% ax.YAxis.TickLabelFormat = '%2.0f';


vertLabel = '$H$';
ylim([0.74, 1.7])


% xlabel('$\xi / \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(vertLabel, 'interpreter', 'latex', 'fontsize', labelFontSize)
xticks(-1:0.5:1)
yticks(1:0.2:1.8)






%%% Collapse parameter: displacement / thickness
offset = -0.05;
scale = 9;

% Plot displacement thickness
h(4) = nexttile;
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

    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        % Legend
        if s == 3
            vis = 'on';
        else
            vis = 'off';
        end

        wavelength = integral.(caze).wavelength; 
        amplitude = integral.(caze).amplitude;
        steepness = (2 * amplitude) / wavelength;
        
        label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelengths.(wave), steepnesses.(wave));
        
        thickness_data = integral.(caze)(phase).(version).('thickness');
        displacement_data = integral.(caze)(phase).(version).('displacement');
        collapse_parameter = displacement_data ./ thickness_data;

        H(w) = plot((integral.(caze)(phase).x - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, collapse_parameter, ...
                    'color', wave_colors{w}, 'linewidth', lw, 'linestyle', linestyles{s}, ...
                    'displayname', label, 'HandleVisibility', vis);

        % Plot phase for reference
        if wave == 'D' && s == 1
            testX = (integral.(caze)(phase).wave.x  - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength;integral.(caze)(phase).wave.wave_profile;
            testY = (scale * integral.(caze)(phase).wave.wave_profile) + offset;
            extrapX = -1:0.001:1;
            extrap_testY = interp1(testX, testY, extrapX, 'spline', 'extrap');
            plot(extrapX, extrap_testY, ...
                 'color', 'black', 'linewidth', lw, 'linestyle', '-', 'HandleVisibility', 'off');
            tmpX = extrapX;
            tmpY = extrap_testY;
            X = [tmpX(:); flipud(tmpX(:))];
            Y = [tmpY(:); (-1)*ones(size(tmpY(:)))];
            
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
% % Add legend for linestyles
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

vertLabel = '$\delta^* \mathbin{/} \delta$';
ylim([-0.15, 0.35])
xlabel(tile, '$\xi \mathbin{/} \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(vertLabel, 'interpreter', 'latex', 'fontsize', labelFontSize)
xticks(-1:0.5:1)
yticks(0:0.1:0.4)



addPanelLabels(h, {'a', 'b', 'c', 'd'}, 'Offset', [-0.2, 1.18])
% [-0.25, 1.18]

% Legend
leg = legend('interpreter', 'latex', 'box', 'off', 'fontsize', legendFontSize, 'orientation', 'horizontal');
leg.Layout.Tile = 'north';
leg.IconColumnWidth = 19;

% Save figure
% pause(3)
% figure_name = strcat('BoundaryLayer_Curvilinear_Displacement_Momentum_Shape_Collapse_Phase', num2str(phase), '.pdf');
% exportgraphics(tile, fullfile(figure_folder, 'BoundaryLayer', figure_name), 'Resolution', 300, 'Resolution', 600, 'ContentType', 'image');
% close all
% fprintf('Generated figure: %s\n\n', figure_name)








%% 3 Panel plot of \delta^* / \delta, \theta / \delta and H


clear h

phase = 1;
version = 'filtered';
linestyles = {'-.', '--', '-'};
lw = 1;
alpha = 0.25;

tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;

clc; close all; 
ax = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,5]);
tile = tiledlayout(1,3, 'tilespacing', 'loose', 'padding', 'tight');



%%% Collapse parameter: displacement / thickness
offset = -0.027;
scale = 5;

% Plot displacement thickness
h(1) = nexttile;
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

    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        % Legend
        if s == 3
            vis = 'on';
        else
            vis = 'off';
        end

        wavelength = integral.(caze).wavelength; 
        amplitude = integral.(caze).amplitude;
        steepness = (2 * amplitude) / wavelength;
        
        label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelengths.(wave), steepnesses.(wave));
        
        thickness_data = integral.(caze)(phase).(version).('thickness');
        displacement_data = integral.(caze)(phase).(version).('displacement');
        collapse_parameter = displacement_data ./ thickness_data;

        H(w) = plot((integral.(caze)(phase).x - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, collapse_parameter, ...
                    'color', wave_colors{w}, 'linewidth', lw, 'linestyle', linestyles{s}, ...
                    'displayname', label, 'HandleVisibility', vis);

        % Plot phase for reference
        if wave == 'D' && s == 1
            testX = (integral.(caze)(phase).wave.x  - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength;integral.(caze)(phase).wave.wave_profile;
            testY = (scale * integral.(caze)(phase).wave.wave_profile) + offset;
            extrapX = -1:0.001:1;
            extrap_testY = interp1(testX, testY, extrapX, 'spline', 'extrap');
            plot(extrapX, extrap_testY, ...
                 'color', 'black', 'linewidth', lw, 'linestyle', '-', 'HandleVisibility', 'off');
            tmpX = extrapX;
            tmpY = extrap_testY;
            X = [tmpX(:); flipud(tmpX(:))];
            Y = [tmpY(:); (-1)*ones(size(tmpY(:)))];
            
            patch(X, Y, 'k', ...
                  'FaceAlpha', alpha, ...
                  'EdgeColor', 'none', ...
                  'handlevisibility', 'off');
        end
    end
end


hold off
vertLabel = '$\delta^* \mathbin{/} \delta$';
xlabel(tile, '$\xi \mathbin{/} \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(vertLabel, 'interpreter', 'latex', 'fontsize', labelFontSize)
xticks(-1:1:1)
yticks(0:0.1:0.4)
ylim([-0.08, 0.3001])
axis square



% Legend
leg = legend('interpreter', 'latex', 'box', 'off', 'fontsize', legendFontSize, 'orientation', 'horizontal');
leg.Layout.Tile = 'north';
leg.IconColumnWidth = 19;


% set(gca, 'XTick', []); 





%%% Momentum
parameter = 'momentum';

offset = -0.027;
scale = 5;

h(2) = nexttile;
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

    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        wavelength = integral.(caze).wavelength; 
        amplitude = integral.(caze).amplitude;
        steepness = (2 * amplitude) / wavelength;
        
        label = sprintf('$\\lambda_{%s}, u_{\\infty} = %1.2f$', wavelengths.(wave), u_inf);

        thickness_data = integral.(caze)(phase).(version).('thickness');
        data = integral.(caze)(phase).(version).(parameter);
        
        data = data ./ thickness_data;
        
        H(w) = plot((integral.(caze)(phase).x - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, data, ...
                    'color', wave_colors{w}, 'linewidth', lw, 'linestyle', linestyles{s}, ...
                    'displayname', label);

        % Plot phase for reference
        if wave == 'D' && s == 1
            testX = (integral.(caze)(phase).wave.x  - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength;integral.(caze)(phase).wave.wave_profile;
            testY = (scale * integral.(caze)(phase).wave.wave_profile) + offset;
            extrapX = -1:0.001:1;
            extrap_testY = interp1(testX, testY, extrapX, 'spline', 'extrap');
            plot(extrapX, extrap_testY, ...
                 'color', 'black', 'linewidth', lw, 'linestyle', '-', 'HandleVisibility', 'off');
            tmpX = extrapX;
            tmpY = extrap_testY;
            X = [tmpX(:); flipud(tmpX(:))];
            Y = [tmpY(:); (-1)*ones(size(tmpY(:)))];
            
            patch(X, Y, 'k', ...
                  'FaceAlpha', alpha, ...
                  'EdgeColor', 'none', ...
                  'handlevisibility', 'off');
        end
    end
end
hold off


vertLabel = '$\theta \mathbin{/} \delta$';
ylim([-0.08, 0.30001])
ylabel(vertLabel, 'interpreter', 'latex', 'fontsize', labelFontSize)
yticks(0:0.1:0.4)
axis square

% set(gca, 'XTick', []); 




%%% Shape 
parameter = 'shape';


scale = 10;
offset = 0.95;

% Plot displacement thickness
h(3) = nexttile;
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

    for w = 1:length(waves)
        wave = waves{w};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');

        wavelength = integral.(caze).wavelength; 
        amplitude = integral.(caze).amplitude;
        steepness = (2 * amplitude) / wavelength;
        
        label = sprintf('$\\lambda_{%s}, u_{\\infty} = %1.2f$', wavelengths.(wave), u_inf);
        
        data = integral.(caze)(phase).(version).(parameter);
        H(w) = plot((integral.(caze)(phase).x - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, data, ...
                    'color', wave_colors{w}, 'linewidth', lw, 'linestyle', linestyles{s}, ...
                    'displayname', label);

        % Plot phase for reference
        if wave == 'D' && s == 1
            testX = (integral.(caze)(phase).wave.x  - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength;integral.(caze)(phase).wave.wave_profile;
            testY = (scale * integral.(caze)(phase).wave.wave_profile) + offset;
            extrapX = -1:0.001:1;
            extrap_testY = interp1(testX, testY, extrapX, 'spline', 'extrap');
            plot(extrapX, extrap_testY, ...
                 'color', 'black', 'linewidth', lw, 'linestyle', '-', 'HandleVisibility', 'off');
            tmpX = extrapX;
            tmpY = extrap_testY;
            X = [tmpX(:); flipud(tmpX(:))];
            Y = [tmpY(:); (-1)*ones(size(tmpY(:)))];
            
            patch(X, Y, 'k', ...
                  'FaceAlpha', alpha, ...
                  'EdgeColor', 'none', ...
                  'handlevisibility', 'off');
        end
    end
end
hold off
vertLabel = '$H$';
ylim([0.84, 1.6])


% xlabel('$\xi / \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(vertLabel, 'interpreter', 'latex', 'fontsize', labelFontSize)
xticks(-1:1:1)
yticks(1:0.2:1.8)
axis square

addPanelLabels(h, {'a', 'b', 'c'}, 'Offset', [-0.3, 1.18])






% Save figure
% pause(3)
% figure_name = strcat('BoundaryLayer_Curvilinear_Displacement_Momentum_Shape_Collapse_Phase', num2str(phase), '_New.pdf');
% exportgraphics(tile, fullfile(figure_folder, 'BoundaryLayer', figure_name), 'Resolution', 300, 'Resolution', 600, 'ContentType', 'image');
% close all
% fprintf('Generated figure: %s\n\n', figure_name)











%% Functions


function addPanelLabels(ax, labels, varargin)
% addPanelLabels(ax, labels) adds (a),(b),... just OUTSIDE top-left of each axes.
% ax     : array of axes handles (e.g., from tiledlayout / findall)
% labels : cellstr like {'a','b','c'} or string array ["a" "b" "c"]
%
% Optional name-value:
% 'Offset'   : [dx dy] in normalized axes units (default [-0.10 1.02])
% 'FontSize' : default 12
% 'FontName' : default 'Times New Roman'

p = inputParser;
addParameter(p,'Offset',[-0.10 1.1]);
addParameter(p,'FontSize', 10);
addParameter(p,'FontName','Times New Roman');
parse(p,varargin{:});
off = p.Results.Offset;

labels = string(labels);
for k = 1:numel(ax)
    if ~isgraphics(ax(k),'axes'), continue; end

    % Plain parentheses + italic letter:
    % TeX interpreter: \it turns italic ON, \rm returns to roman.
    s = sprintf('(\\ita\\rm)');              % placeholder
    s = sprintf('(\\it%s\\rm)', labels(k));  % actual label

    text(ax(k), off(1), off(2), s, ...
        'Units','normalized', ...
        'Interpreter','tex', ...           % keeps italics control simple
        'FontName',p.Results.FontName, ...
        'FontSize',p.Results.FontSize, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'Clipping','off');                 % critical: allow outside axes
    end
end




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


function interpolated = InterpolateInterior(x, y)
    % x: 1D array of independent variable (may contain NaNs)
    % y: 1D array of dependent variable (may contain NaNs)

    % Step 1: Find joint validity mask
    valid_mask = ~isnan(x) & ~isnan(y);

    % Step 2: Identify first and last *jointly valid* indices
    first_idx = find(valid_mask, 1, 'first');
    last_idx  = find(valid_mask, 1, 'last');

    if isempty(first_idx) || isempty(last_idx) || first_idx == last_idx
        % Nothing valid to interpolate
        interpolated = y;
        return
    end

    % Step 3: Trim x and y to internal valid range
    x_trim = x(first_idx:last_idx);
    y_trim = y(first_idx:last_idx);

    % Step 4: Mask internal NaNs in y, only where x is valid
    internal_valid = ~isnan(x_trim);
    nan_mask = isnan(y_trim) & internal_valid;

    % Step 5: Interpolate
    y_interp = y_trim;
    y_interp(nan_mask) = interp1(x_trim(~isnan(y_trim) & ~isnan(x_trim)), ...
                                 y_trim(~isnan(y_trim) & ~isnan(x_trim)), ...
                                 x_trim(nan_mask), 'spline');

    % Step 6: Rebuild full output with original NaN pads
    interpolated = y;
    interpolated(first_idx:last_idx) = y_interp;
end


function output = FilterData(x, q, tol, max_island_length, sgolay_length)
    output = smoothdata(InterpolateInterior(x, StatisticalGradientFilter(q, tol, max_island_length)), 'sgolay', sgolay_length);
end


function output = PolynomialFit(x, y, n)
    mask = ~isnan(x) & ~isnan(y);
    p = polyfit(x(mask), y(mask), n);
    output = polyval(p, x);
    output(isnan(y)) = nan;
end





