%%% Curvilinear Integral Boundary Layer Parameters

addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
curvilinear_path = fullfile(project_path, 'curvilinear_new');
wave_parameters = readcell("Offshore_Waves.xlsx");

% Where to save figures
figure_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/pdf_test6';


% Boundary layer detection parameters
boundary_layer_percent = 0.95;
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



% Quantities to be looped
parameters = {'thickness', 'displacement', 'momentum', 'shape'};
wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};




%% Plot specific parameter for a specific phase for all waves and wind speeds
% 
% phase = 1;
% parameter = 'thickness';
% version = 'filtered';
% linestyles = {'-.', '--', '-'};
% lw = 2;
% alpha = 0.5;
% 
% % Scale reference wave based on parameter
% if strcmp(parameter, 'thickness')
%     scale = 6;
% elseif strcmp(parameter, 'displacement')
%     scale = 3;
% elseif strcmp(parameter, 'momentum')
%     scale = 1;
% elseif strcmp(parameter, 'shape')
%     scale = 35;
% end
% 
% 
% tickFontSize = 8;
% labelFontSize = 16;
% legendFontSize = 12;
% 
% close all; 
% ax = figure('color', 'white', 'units', 'centimeters', 'position', [2,4,30,10]);
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% 
% sgtitle
% hold on
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
%     % Legend shit
%     if s == 3
%         vis = 'on';
%     else
%         vis = 'off';
%     end
% 
%     for w = 1:length(waves)
%         wave = waves{w};
%         caze = strcat(wind_speed, '_WV', wave, '_AG0');
% 
%         wavelength = integral.(caze).wavelength; 
%         amplitude = integral.(caze).amplitude;
%         steepness = (2 * amplitude) / wavelength;
% 
%         % label = sprintf('$\\lambda_{%s}, \\hspace{1mm} u_{\\infty} = %1.2f$ m/s', wavelengths.(wave), u_inf);
%         label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelengths.(wave), steepnesses.(wave));
% 
%         data = integral.(caze)(phase).(version).(parameter);
%         H(w) = plot((integral.(caze)(phase).x - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, data, ...
%                     'color', wave_colors{w}, 'linewidth', lw, 'linestyle', linestyles{s}, ...
%                     'displayname', label, 'HandleVisibility', vis);
% 
%         % Plot phase for reference
%         if wave == 'D' && s == 1
%             reference = plot((integral.(caze)(phase).wave.x  - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, (scale * integral.(caze)(phase).wave.wave_profile) + mean(data, 'all', 'omitnan'), ...
%                        'color', 'black', 'linewidth', lw, 'linestyle', ':', 'HandleVisibility', 'off');
%             reference.Color(4) = alpha;
%         end
%     end
% end
% 
% % Dummy line to add a break
% plot(nan, nan, 'color', 'white', 'displayname', '')
% % plot(nan, nan, 'color', 'white', 'displayname', '')
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
%     P.Color(4) = 0.5;
% end
% 
% % Vertical limits
% if strcmp(parameter, 'thickness')
%     ylim([0.06, 0.16])
% end
% 
% hold off
% 
% 
% % Legend
% lgd = legend('interpreter', 'latex', 'fontsize', 10, 'location',  'eastoutside', 'box', 'off');
% 
% 
% % Axes labels
% if strcmp(parameter, 'thickness')
%     vertLabel = '$\delta$ [m]';
% elseif strcmp(parameter, 'displacement')
%     vertLabel = '$\delta^*$ [m]';
% elseif strcmp(parameter, 'momentum')
%     vertLabel = '$\theta$ [m]';
% elseif strcmp(parameter, 'shape')
%     vertLabel = '$H$';
%     limits = [1, 1.7];
% end
% 
% xlabel('$\xi / \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylabel(vertLabel, 'interpreter', 'latex', 'fontsize', labelFontSize)
% xticks(-1:0.25:1)


% Save figure
% pause(1)
% figure_name = strcat('BoundaryLayer_Curvilinear', [upper(parameter(1)), lower(parameter(2:end))], '_Phase', num2str(phase), '.pdf');
% exportgraphics(ax, fullfile(figure_folder, 'BoundaryLayer', figure_name), 'Resolution', 600, 'ContentType', 'image');
% close all
% clc; fprintf('Generated figure: %s\n\n', figure_name)

% clear alpha limits caze data H linestyles lw parameter phase reference s scale w wave wind_speed version ans vertLabel; clc





%% Boundary layer thickness at different phases versus wave age


% phase = 2;
% parameter = 'thickness';
% version = 'filtered';
% linestyles = {'-.', '--', '-'};
% lw = 2;
% alpha = 0.5;
% 
% wind_speed_markers = {'^', 'square', 'o'};
% 
% tickFontSize = 8;
% labelFontSize = 16;
% legendFontSize = 12;
% 
% 
% clear tmp tmp2
% clc; close all; 
% ax = figure('color', 'white');
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% title(sprintf('Phase %1.0f', phase))
% hold on
% 
% for w = 1:length(waves)
%     wave = waves{w};
% 
% 
%     for s = 1:length(wind_speeds)
%         wind_speed = wind_speeds{s};
%         caze = strcat(wind_speed, '_WV', wave, '_AG0');
% 
%         if ismember(wind_speed(end), {'4'})
%             u_inf = 2.4181;
%         elseif ismember(wind_speed(end), {'6'})
%             u_inf = 3.8709;
%         elseif ismember(wind_speed(end), {'8'})
%             u_inf = 5.4289;
%         end
% 
% 
% 
% 
%         wavelength = integral.(caze).wavelength; 
%         amplitude = integral.(caze).amplitude;
%         frequency = frequencies.(wave);
%         steepness = (2 * pi * amplitude) / wavelength;
%         wave_speed = frequency * wavelength;
%         wave_age = wave_speed / u_inf;
% 
%         fprintf('%s: Wave speed = %1.3f\n', wave, wave_speed)
%         fprintf('Wave_speed / u_inf = %1.3f\n\n', wave_age)
% 
% 
%         % label = sprintf('$\\lambda_{%s}, \\hspace{1mm} u_{\\infty} = %1.2f$ m/s', wavelengths.(wave), u_inf);
%         label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelengths.(wave), steepnesses.(wave));
% 
%         data = integral.(caze)(phase).(version).(parameter);
%         idx = round(length(data) / 2);
% 
%         % H(w) = plot((integral.(caze)(phase).x - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, data, ...
%         %             'color', wave_colors{w}, 'linewidth', lw, 'linestyle', linestyles{s}, ...
%         %             'displayname', label, 'HandleVisibility', vis);
% 
% 
%         scatter(wave_age, data(idx), 100,  wind_speed_markers{s}, 'filled', 'MarkerFaceColor', wave_colors{w})
% 
%         tmp(s) = wave_age;
%         tmp2(s) = data(idx);
% 
%     end
% 
%     plot(tmp, tmp2, 'linewidth', 2, 'color', wave_colors{w})
% end
% 
% hold off
% ylim([0.05, 0.15])
% xlabel('Wave Age: $c / u_{\infty}$', 'interpreter', 'latex')
% ylabel('Boundary Layer Thickness $\delta$', 'interpreter', 'latex')



%% Boundary layer thickness at all phases + mean, versus wave age: Single wave

% clear tmp
% wv = 1;
% parameter = 'thickness';
% version = 'filtered';
% linestyles = {'-.', '--', '-'};
% lw = 2;
% alpha = 0.5;
% 
% phase_trans = [0.25, 0.5, 0.75, 1];
% 
% wind_speed_markers = {'^', 'square', 'o'};
% 
% tickFontSize = 8;
% labelFontSize = 16;
% legendFontSize = 12;
% 
% clc; close all; 
% ax = figure('color', 'white');
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% 
% hold on
% for phase = 1:4
%     for w = wv
%         wave = waves{w};
% 
% 
%         % Loop through wind speeds
%         for s = 1:length(wind_speeds)
%             wind_speed = wind_speeds{s};
%             caze = strcat(wind_speed, '_WV', wave, '_AG0');
% 
%             if ismember(wind_speed(end), {'4'})
%                 u_inf = 2.4181;
%             elseif ismember(wind_speed(end), {'6'})
%                 u_inf = 3.8709;
%             elseif ismember(wind_speed(end), {'8'})
%                 u_inf = 5.4289;
%             end
% 
%             % Get wave parameters
%             wavelength = integral.(caze).wavelength; 
%             amplitude = integral.(caze).amplitude;
%             frequency = frequencies.(wave);
%             steepness = (2 * pi * amplitude) / wavelength;
%             wave_speed = frequency * wavelength;
%             wave_age = wave_speed / u_inf;
% 
%             fprintf('%s: Wave speed = %1.3f\n', wave, wave_speed)
%             fprintf('Wave_speed / u_inf = %1.3f\n\n', wave_age)
% 
% 
%             label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelengths.(wave), steepnesses.(wave));
% 
% 
%             data = integral.(caze)(phase).(version).(parameter);
%             idx = round(length(data) / 2);
% 
%             scatter(wave_age, data(idx), 100,  wind_speed_markers{s}, 'filled', ...
%                     'MarkerFaceColor', wave_colors{w}, ...
%                     'markerfacealpha', phase_trans(phase), ...
%                     'HandleVisibility', 'off')
% 
%             tmp(s) = wave_age;
%             tmp2(s) = data(idx);
% 
%         end
%         phase_label = sprintf('Phase %1.0f', phase);
%         P = plot(tmp, tmp2, 'linewidth', 2, 'color', wave_colors{w}, 'HandleVisibility', 'on', 'displayname', phase_label);
%         P.Color(4) = phase_trans(phase);
%     end
% 
%     xlabel('Wave Age: $c / u_{\infty}$', 'interpreter', 'latex')
%     ylabel('Boundary Layer Thickness $\delta$', 'interpreter', 'latex')
% end
% 
% % Loop through a given wavelength and average of all 4 phases
% for w = wv
%     wave = waves{w};
%     tmpX = nan(1,3);
%     tmpY = nan(1,3);
%     for s = 1:length(wind_speeds)
%         wind_speed = wind_speeds{s};
% 
%         if ismember(wind_speed(end), {'4'})
%             u_inf = 2.4181;
%         elseif ismember(wind_speed(end), {'6'})
%             u_inf = 3.8709;
%         elseif ismember(wind_speed(end), {'8'})
%             u_inf = 5.4289;
%         end
% 
%         caze = strcat(wind_speed, '_WV', wave, '_AG0');
%         wavelength = integral.(caze).wavelength; 
%         amplitude = integral.(caze).amplitude;
%         frequency = frequencies.(wave);
%         steepness = (2 * pi * amplitude) / wavelength;
%         wave_speed = frequency * wavelength;
%         wave_age = wave_speed / u_inf;
% 
% 
% 
%         tmp = nan(1,4);
%         for phase = 1:4
%             data = integral.(caze)(phase).(version).(parameter);
%             idx = round(length(data) / 2);
%             tmp(phase) = data(idx);
%         end
%         disp(mean(tmp, 'all', 'omitnan'))
%         disp(wave_age)
%         scatter(wave_age, mean(tmp, 'all', 'omitnan'), 100,  wind_speed_markers{s}, 'filled', 'MarkerFaceColor', 'black', 'HandleVisibility', 'off')
%         tmpX(s) = wave_age;
%         tmpY(s) = mean(tmp, 'all', 'omitnan');
%     end
%     plot(tmpX, tmpY, 'linewidth', 2, 'color', 'black', 'displayname', 'Mean of Phases')
% end
% hold off
% 
% 
% legend('interpreter', 'latex', 'box', 'off', 'location', 'northeast')
% xlim([0.05, 0.35])
% ylim([0.05, 0.15])



%% Boundary layer thickness at crest versus wave age all phases


% clear tmp
% 
% % phase = 1;
% parameter = 'thickness';
% version = 'filtered';
% linestyles = {'-.', '--', '-'};
% lw = 2;
% alpha = 0.5;
% 
% wind_speed_markers = {'^', 'square', 'o'};
% 
% tickFontSize = 8;
% labelFontSize = 16;
% legendFontSize = 12;
% 
% clc; close all; 
% ax = figure('color', 'white');
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
% 
% hold on
% 
% for phase = 1:4
%     for w = 1:length(waves)
%         wave = waves{w};
% 
% 
%         for s = 1:length(wind_speeds)
%             wind_speed = wind_speeds{s};
%             caze = strcat(wind_speed, '_WV', wave, '_AG0');
% 
%             if ismember(wind_speed(end), {'4'})
%                 u_inf = 2.4181;
%             elseif ismember(wind_speed(end), {'6'})
%                 u_inf = 3.8709;
%             elseif ismember(wind_speed(end), {'8'})
%                 u_inf = 5.4289;
%             end
% 
% 
% 
%             wavelength = integral.(caze).wavelength; 
%             amplitude = integral.(caze).amplitude;
%             frequency = frequencies.(wave);
%             steepness = (2 * pi * amplitude) / wavelength;
%             wave_speed = frequency * wavelength;
%             wave_age = wave_speed / u_inf;
% 
%             fprintf('%s: Wave speed = %1.3f\n', wave, wave_speed)
%             fprintf('Wave_speed / u_inf = %1.3f\n\n', wave_age)
% 
% 
%             % label = sprintf('$\\lambda_{%s}, \\hspace{1mm} u_{\\infty} = %1.2f$ m/s', wavelengths.(wave), u_inf);
%             label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelengths.(wave), steepnesses.(wave));
% 
%             data = integral.(caze)(phase).(version).(parameter);
%             idx = round(length(data) / 2);
% 
%             scatter(wave_age, data(idx), 100,  wind_speed_markers{s}, 'filled', 'MarkerFaceColor', wave_colors{w})
% 
%             tmp(s) = wave_age;
%             tmp2(s) = data(idx);
% 
%         end
% 
%         plot(tmp, tmp2, 'linewidth', 2, 'color', wave_colors{w})
%     end
% end
% 
% hold off
% ylim([0.05, 0.15])
% 
% xlabel('Wave Age: $c / u_{\infty}$', 'interpreter', 'latex')
% ylabel('Boundary Layer Thickness $\delta$', 'interpreter', 'latex')


%% Plotting the mean of phases with phase curves in the background
% Max/min of BL across the 4 phases is shaded
% MAIN SECTION HERE

tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;

sz = 20;
lw = 1;
wind_speed_markers = {'^', 'square', 'o'};

clear tmpX tmpY 
clear tmpY_max tmpY_min

% Loop through a given wavelength and average of all 4 phases
clc; close all
figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,6]);
tiledlayout(1,1, 'padding', 'loose')
ax = nexttile;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
hold on

% Loop through wavelengths
for w = 1:length(waves)
    wave = waves{w};
    
    % Temporary structures
    tmpX = nan(1,3);
    tmpY = nan(1,3);
    tmpY_max = nan(1,3);
    tmpY_min = nan(1,3);

    % Loop through wind speeds
    for s = 1:length(wind_speeds)
        wind_speed = wind_speeds{s};

        if ismember(wind_speed(end), {'4'})
            u_inf = 2.4181;
        elseif ismember(wind_speed(end), {'6'})
            u_inf = 3.8709;
        elseif ismember(wind_speed(end), {'8'})
            u_inf = 5.4289;
        end

        % Get wave parameters
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        wavelength = integral.(caze).wavelength; 
        amplitude = integral.(caze).amplitude;
        frequency = frequencies.(wave);
        steepness = (2 * pi * amplitude) / wavelength;
        wave_speed = frequency * wavelength;
        wave_age = wave_speed / u_inf;
    
        % Loop through phases
        tmp = nan(1,4);
        for phase = 1:4
            data = integral.(caze)(phase).('fitted').('thickness');
            idx = round(length(data) / 2);
            tmp(phase) = data(idx);
        end
        
        % Plot
        scatter(wave_age, mean(tmp, 'all', 'omitnan'), sz,  wind_speed_markers{s}, 'filled', ...
                'MarkerFaceColor', wave_colors{w}, 'HandleVisibility', 'off')

        % Save to plot line
        tmpX(s) = wave_age;
        tmpY(s) = mean(tmp, 'all', 'omitnan');

        tmpY_max(s) = max(tmp);
        tmpY_min(s) = min(tmp);
    end

    % Plot line
    label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelengths.(wave), steepnesses.(wave));
    plot(tmpX, tmpY, 'linewidth', lw, 'color', wave_colors{w}, 'displayname', label)

    % Shaded region
    % hFill = patch( ...
    % [tmpX, fliplr(tmpX)], ...
    % [tmpY_max, fliplr(tmpY_min)], ...
    % hex2rgb(wave_colors{w}), ...
    % 'FaceAlpha', 0.15, ...        % transparency
    % 'EdgeColor', 'none', ...      % no outline
    % 'HandleVisibility', 'off');   % keep legend clean
    % uistack(hFill, 'bottom')
end

% Add legend for markers
plot(nan, nan, 'color', 'white', 'displayname', '')

% Dummy lines to get wind-speed legend
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};

    if ismember(wind_speed(end), {'4'})
        u_inf = 2.4181;
    elseif ismember(wind_speed(end), {'6'})
        u_inf = 3.8709;
    elseif ismember(wind_speed(end), {'8'})
        u_inf = 5.4289;
    end

    hLeg = plot(nan, nan, wind_speed_markers{s}, ...
    'MarkerFaceColor','black', ...
    'MarkerEdgeColor','black', ...
    'MarkerSize', 4, ...        % <-- CONTROL LEGEND SIZE HERE
    'LineWidth', 1, ...
    'LineStyle','none', ...
    'DisplayName', sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf));

end

leg = legend('interpreter', 'latex', 'location', 'eastoutside', 'box', 'off', 'fontsize', legendFontSize);
leg.IconColumnWidth = 19;
hold off
axis square
% ylim([0.085, 0.115])
xlim([0.05, 0.35])
xticks(0.1:0.1:0.3)
% yticks(0.08:0.01:0.12)
xlabel('$c \mathbin{/} u_{\infty}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$\overline{\delta}_{\varphi}$ [m]', 'interpreter', 'latex', 'fontsize', labelFontSize)



% Save figure
pause(3)
figure_name = 'BoundaryLayer_Curvilinear_WaveAge.pdf';
exportgraphics(ax, fullfile(figure_folder, 'BoundaryLayer', figure_name), 'Resolution', 600, 'ContentType', 'image');
close all
clc; fprintf('Generated figure: %s\n\n', figure_name)



%% Plot with power law fits

tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;
sz = 40;
lw = 1;
wind_speed_markers = {'^', 'square', 'o'};
clear tmpX tmpY
clear tmpY_max tmpY_min

% Storage for fitting
all_X = cell(1, length(waves));  % wave age for each wavelength
all_Y = cell(1, length(waves));  % delta for each wavelength

% Loop through a given wavelength and average of all 4 phases
clc; close all
figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,7]);
tiledlayout(1,1, 'padding', 'tight')
ax = nexttile;
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)
hold on

% Loop through wavelengths
for w = 1:length(waves)
    wave = waves{w};
    % Temporary structures
    tmpX = nan(1,3);
    tmpY = nan(1,3);
    tmpY_max = nan(1,3);
    tmpY_min = nan(1,3);
    
    % Loop through wind speeds
    for s = 1:length(wind_speeds)
        wind_speed = wind_speeds{s};
        if ismember(wind_speed(end), {'4'})
            u_inf = 2.4181;
        elseif ismember(wind_speed(end), {'6'})
            u_inf = 3.8709;
        elseif ismember(wind_speed(end), {'8'})
            u_inf = 5.4289;
        end
        
        % Get wave parameters
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
        wavelength = integral.(caze).wavelength;
        amplitude = integral.(caze).amplitude;
        frequency = frequencies.(wave);
        steepness = (2 * pi * amplitude) / wavelength;
        wave_speed = frequency * wavelength;
        wave_age = wave_speed / u_inf;
        
        % Loop through phases
        tmp = nan(1,4);
        for phase = 1:4
            data = integral.(caze)(phase).('fitted').('thickness');
            idx = round(length(data) / 2);
            tmp(phase) = data(idx);
        end
        
        % Plot
        scatter(wave_age, mean(tmp, 'all', 'omitnan'), sz,  wind_speed_markers{s}, 'filled', ...
            'MarkerFaceColor', wave_colors{w}, 'HandleVisibility', 'off')
        
        % Save to plot line
        tmpX(s) = wave_age;
        tmpY(s) = mean(tmp, 'all', 'omitnan');
        tmpY_max(s) = max(tmp);
        tmpY_min(s) = min(tmp);
    end
    
    % Store for fitting
    all_X{w} = tmpX;
    all_Y{w} = tmpY;
    
    % Plot line (will be replaced by fit lines, so comment out or set HandleVisibility off)
    % label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelengths.(wave), steepnesses.(wave));
    % plot(tmpX, tmpY, 'linewidth', lw, 'color', wave_colors{w}, 'displayname', label)
end

% =====================
% Shared-exponent power law fit
% Model: delta = A_w * (c/u_inf)^b, where b is shared, A_w varies per wavelength
% =====================

% Objective function: sum of squared residuals in log space
% Parameters: [b, log(A_1), log(A_2), log(A_3), log(A_4)]
n_waves = length(waves);

% Build arrays for fitting
X_all = [];
Y_all = [];
wave_idx = [];
for w = 1:n_waves
    X_all = [X_all, all_X{w}];
    Y_all = [Y_all, all_Y{w}];
    wave_idx = [wave_idx, w*ones(1, length(all_X{w}))];
end

% Remove NaN
valid = isfinite(X_all) & isfinite(Y_all);
X_all = X_all(valid);
Y_all = Y_all(valid);
wave_idx = wave_idx(valid);

% Log transform
log_X = log(X_all);
log_Y = log(Y_all);

% Objective function
obj_fun = @(params) compute_residuals(params, log_X, log_Y, wave_idx, n_waves);

% Initial guess: b = -0.3, A_w ≈ mean(Y) for each wave
b0 = -0.3;
logA0 = zeros(1, n_waves);
for w = 1:n_waves
    logA0(w) = log(mean(all_Y{w}, 'omitnan'));
end
params0 = [b0, logA0];

% Optimize
options = optimset('Display', 'off', 'TolFun', 1e-10, 'TolX', 1e-10);
params_fit = fminsearch(obj_fun, params0, options);

% Extract fitted parameters
b_fit = params_fit(1);
A_fit = exp(params_fit(2:end));

% Calculate R^2
y_pred = zeros(size(log_Y));
for i = 1:length(log_Y)
    y_pred(i) = log(A_fit(wave_idx(i))) + b_fit * log_X(i);
end
SS_res = sum((log_Y - y_pred).^2);
SS_tot = sum((log_Y - mean(log_Y)).^2);
R2 = 1 - SS_res / SS_tot;

% Print results
fprintf('Shared exponent fit: delta = A * (c/u_inf)^(%.2f)\n', b_fit);
fprintf('R^2 = %.3f\n', R2);
for w = 1:n_waves
    fprintf('  Wave %s: A = %.4f m\n', waves{w}, A_fit(w));
end

% Plot fitted curves
x_fit = linspace(0.05, 0.35, 100);
for w = 1:n_waves
    y_fit = A_fit(w) * x_fit.^b_fit;
    label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelengths.(waves{w}), steepnesses.(waves{w}));
    plot(x_fit, y_fit, 'linewidth', lw, 'color', wave_colors{w}, 'displayname', label)
end

% Add fit annotation
fit_text = sprintf('$\\bar{\\delta}_\\varphi = A_\\lambda \\cdot (c/u_\\infty)^{%.2f}$\n$R^2 = %.2f$', b_fit, R2);
text(0.22, 0.118, fit_text, 'interpreter', 'latex', 'fontsize', 8)

% Add legend for markers
plot(nan, nan, 'color', 'white', 'displayname', '')
% Dummy lines to get wind-speed legend
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};
    if ismember(wind_speed(end), {'4'})
        u_inf = 2.4181;
    elseif ismember(wind_speed(end), {'6'})
        u_inf = 3.8709;
    elseif ismember(wind_speed(end), {'8'})
        u_inf = 5.4289;
    end
    hLeg = plot(nan, nan, wind_speed_markers{s}, ...
        'MarkerFaceColor','black', ...
        'MarkerEdgeColor','black', ...
        'MarkerSize', 4, ...
        'LineWidth', 1, ...
        'LineStyle','none', ...
        'DisplayName', sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf));
end
leg = legend('interpreter', 'latex', 'location', 'eastoutside', 'box', 'off', 'fontsize', legendFontSize);
leg.IconColumnWidth = 19;
hold off
axis square
% ylim([0.08, 0.12])
yticks(0.08:0.01:0.13)
xlim([0.05, 0.35])
xticks(0.1:0.1:0.3)
xlabel('$c \mathbin{/} u_{\infty}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$\overline{\delta}_{\varphi}$ [m]', 'interpreter', 'latex', 'fontsize', labelFontSize)

% =====================
% Helper function for residuals
% =====================
function res = compute_residuals(params, log_X, log_Y, wave_idx, n_waves)
    b = params(1);
    logA = params(2:end);
    
    pred = zeros(size(log_Y));
    for i = 1:length(log_Y)
        pred(i) = logA(wave_idx(i)) + b * log_X(i);
    end
    
    res = sum((log_Y - pred).^2);
end



%% Plotting the mean of phases with phase curves in the background: As a tiledlayout
% Max/min of BL across the 4 phases is shaded


% clear tmpX tmpY 
% clear tmpY_max tmpY_min
% 
% % Loop through a given wavelength and average of all 4 phases
% figure('color', 'white')
% t = tiledlayout(1,4);
% sgtitle('Mean of Phases')
% 
% % Loop through wavelengths
% for w = 1:length(waves)
%     wave = waves{w};
% 
%     h(w) = nexttile;
%     hold on
%     % Temporary structures
%     tmpX = nan(1,3);
%     tmpY = nan(1,3);
%     tmpY_max = nan(1,3);
%     tmpY_min = nan(1,3);
% 
%     % Loop through wind speeds
%     for s = 1:length(wind_speeds)
%         wind_speed = wind_speeds{s};
% 
%         if ismember(wind_speed(end), {'4'})
%             u_inf = 2.4181;
%         elseif ismember(wind_speed(end), {'6'})
%             u_inf = 3.8709;
%         elseif ismember(wind_speed(end), {'8'})
%             u_inf = 5.4289;
%         end
% 
%         % Get wave parameters
%         caze = strcat(wind_speed, '_WV', wave, '_AG0');
%         wavelength = integral.(caze).wavelength; 
%         amplitude = integral.(caze).amplitude;
%         frequency = frequencies.(wave);
%         steepness = (2 * pi * amplitude) / wavelength;
%         wave_speed = frequency * wavelength;
%         wave_age = wave_speed / u_inf;
% 
%         % Loop through phases
%         tmp = nan(1,4);
%         for phase = 1:4
%             data = integral.(caze)(phase).(version).(parameter);
%             idx = round(length(data) / 2);
%             tmp(phase) = data(idx);
%         end
% 
%         % Plot
%         scatter(wave_age, mean(tmp, 'all', 'omitnan'), 100,  wind_speed_markers{s}, 'filled', ...
%                 'MarkerFaceColor', wave_colors{w}, 'HandleVisibility', 'off')
% 
%         % Save to plot line
%         tmpX(s) = wave_age;
%         tmpY(s) = mean(tmp, 'all', 'omitnan');
% 
%         tmpY_max(s) = max(tmp);
%         tmpY_min(s) = min(tmp);
%     end
% 
%     % Plot line
%     label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelengths.(wave), steepnesses.(wave));
%     % title(sprintf('$\\lambda_{%s}$', wavelengths.(wave)), 'interpreter', 'latex');
%     title(label, 'interpreter', 'latex');
%     plot(tmpX, tmpY, 'linewidth', 2, 'color', wave_colors{w}, 'displayname', label)
% 
%     % Fill in max/min range
%     hFill = patch( ...
%     [tmpX, fliplr(tmpX)], ...
%     [tmpY_max, fliplr(tmpY_min)], ...
%     hex2rgb(wave_colors{w}), ...
%     'FaceAlpha', 0.25, ...        % transparency
%     'EdgeColor', 'none', ...      % no outline
%     'HandleVisibility', 'off');   % keep legend clean
% 
%     uistack(hFill, 'bottom')
%     hold off
%     axis square
% end
% 
% linkaxes(h, 'xy')
% ylim([0.06, 0.13])
% xlim([0.05, 0.35])
% xlabel(t, 'Wave Age: $\lambda f / u_{\infty}$', 'interpreter', 'latex')
% ylabel(t, 'Boundary Layer Thickness $\delta$ [m]', 'interpreter', 'latex')



%% Boundary layer scaling approach, curve fitting to all 12 points

% % Max/min of BL across the 4 phases is shaded
% clear tmpX tmpY
% clear tmpY_max tmpY_min
% 
% % === Collect ALL data for fitting ===
% all_delta = [];      % Boundary layer thickness
% all_u_inf = [];      % Freestream velocity
% all_c = [];          % Wave phase speed
% all_wave_age = [];   % c/u_inf
% 
% % Loop through a given wavelength and average of all 4 phases
% clc; close all
% figure('color', 'white')
% title('Mean of Phases')
% hold on
% 
% % Loop through wavelengths
% for w = 1:length(waves)
%     wave = waves{w};
% 
%     % Temporary structures
%     tmpX = nan(1,3);
%     tmpY = nan(1,3);
%     tmpY_max = nan(1,3);
%     tmpY_min = nan(1,3);
% 
%     % Loop through wind speeds
%     for s = 1:length(wind_speeds)
%         wind_speed = wind_speeds{s};
% 
%         if ismember(wind_speed(end), {'4'})
%             u_inf = 2.4181;
%         elseif ismember(wind_speed(end), {'6'})
%             u_inf = 3.8709;
%         elseif ismember(wind_speed(end), {'8'})
%             u_inf = 5.4289;
%         end
% 
%         % Get wave parameters
%         caze = strcat(wind_speed, '_WV', wave, '_AG0');
%         wavelength = integral.(caze).wavelength;
%         amplitude = integral.(caze).amplitude;
%         frequency = frequencies.(wave);
%         steepness = (2 * pi * amplitude) / wavelength;
%         wave_speed = frequency * wavelength;
%         wave_age = wave_speed / u_inf;
% 
%         % Loop through phases
%         tmp = nan(1,4);
%         for phase = 1:4
%             data = integral.(caze)(phase).(version).(parameter);
%             idx = round(length(data) / 2);
%             tmp(phase) = data(idx);
%         end
% 
%         % Plot
%         scatter(wave_age, mean(tmp, 'all', 'omitnan'), 100, wind_speed_markers{s}, 'filled', ...
%             'MarkerFaceColor', wave_colors{w}, 'HandleVisibility', 'off')
% 
%         % Save to plot line
%         tmpX(s) = wave_age;
%         tmpY(s) = mean(tmp, 'all', 'omitnan');
%         tmpY_max(s) = max(tmp);
%         tmpY_min(s) = min(tmp);
% 
%         % === Store for fitting ===
%         all_delta = [all_delta, mean(tmp, 'all', 'omitnan')];
%         all_u_inf = [all_u_inf, u_inf];
%         all_c = [all_c, wave_speed];
%         all_wave_age = [all_wave_age, wave_age];
%     end
% 
%     % Plot line
%     label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelengths.(wave), steepnesses.(wave));
%     plot(tmpX, tmpY, 'linewidth', 2, 'color', wave_colors{w}, 'displayname', label)
% 
%     % Shaded region
%     hFill = patch( ...
%         [tmpX, fliplr(tmpX)], ...
%         [tmpY_max, fliplr(tmpY_min)], ...
%         hex2rgb(wave_colors{w}), ...
%         'FaceAlpha', 0.15, ...
%         'EdgeColor', 'none', ...
%         'HandleVisibility', 'off');
%     uistack(hFill, 'bottom')
% end
% 
% % === Fit: delta = B * (u_inf - c)^n ===
% % Take log of both sides: log(delta) = log(B) + n*log(u_inf - c)
% % Linear regression in log-space
% 
% % velocity_deficit = all_u_inf - all_c;  % (u_inf - c)
% 
% % === Fit: delta = A * (c/u_inf)^B ===
% % Take log of both sides: log(delta) = log(A) + B*log(c/u_inf)
% % Linear regression in log-space
% 
% % Log-transform for linear regression
% log_wave_age = log(all_wave_age);
% log_delta = log(all_delta);
% 
% % Linear fit: log(delta) = B * log(c/u_inf) + log(A)
% p = polyfit(log_wave_age, log_delta, 1);
% B_fit = p(1);
% A_fit = exp(p(2));
% 
% % Display results
% fprintf('\n=== Wave Age Power Law Fit ===\n')
% fprintf('delta = A * (c/u_inf)^B\n')
% fprintf('Fitted exponent B = %.4f\n', B_fit)
% fprintf('Fitted coefficient A = %.4f\n', A_fit)
% 
% % Compute R^2
% delta_predicted = A_fit * all_wave_age.^B_fit;
% SS_res = sum((all_delta - delta_predicted).^2);
% SS_tot = sum((all_delta - mean(all_delta)).^2);
% R_squared = 1 - SS_res / SS_tot;
% fprintf('R^2 = %.4f\n\n', R_squared)
% 
% % === Plot the fit curve ===
% wave_age_fit = linspace(0.05, 0.35, 100);
% delta_fit = A_fit * wave_age_fit.^B_fit;
% 
% % Plot fit curve
% plot(wave_age_fit, delta_fit, 'k--', 'linewidth', 2, ...
%     'displayname', sprintf('Fit: $\\delta = %.3f (c/u_\\infty)^{%.2f}$', A_fit, B_fit))

% % Log-transform for linear regression
% log_deficit = log(velocity_deficit);
% log_delta = log(all_delta);
% 
% % Linear fit: log(delta) = n * log(u_inf - c) + log(B)
% p = polyfit(log_deficit, log_delta, 1);
% n_fit = p(1);
% B_fit = exp(p(2));
% 
% % Display results
% fprintf('\n=== Turbulent BL Scaling Fit ===\n')
% fprintf('delta = B * (u_inf - c)^n\n')
% fprintf('Fitted exponent n = %.4f\n', n_fit)
% fprintf('Fitted coefficient B = %.4f\n', B_fit)
% fprintf('(Flat-plate turbulent BL predicts n = -0.2)\n\n')
% 
% % Compute R^2
% delta_predicted = B_fit * velocity_deficit.^n_fit;
% SS_res = sum((all_delta - delta_predicted).^2);
% SS_tot = sum((all_delta - mean(all_delta)).^2);
% R_squared = 1 - SS_res / SS_tot;
% fprintf('R^2 = %.4f\n\n', R_squared)
% 
% % === Plot the fit curve ===
% % Need to convert back to wave age for plotting
% % delta = B * (u_inf - c)^n = B * u_inf^n * (1 - c/u_inf)^n
% 
% % Generate smooth curve for each representative u_inf (or use mean)
% u_inf_mean = mean(all_u_inf);
% wave_age_fit = linspace(0.05, 0.35, 100);
% c_fit = wave_age_fit * u_inf_mean;
% deficit_fit = u_inf_mean - c_fit;
% delta_fit = B_fit * deficit_fit.^n_fit;
% 
% % Plot fit curve
% plot(wave_age_fit, delta_fit, 'k--', 'linewidth', 2, ...
%     'displayname', sprintf('Fit: $\\delta = %.3f (u_\\infty - c)^{%.2f}$', B_fit, n_fit))
% 
% legend('interpreter', 'latex', 'location', 'northeast', 'box', 'off')
% hold off
% ylim([0.06, 0.13])
% xlim([0.05, 0.35])
% xlabel('$c / u_{\infty}$', 'interpreter', 'latex')
% ylabel('$\delta$ [m]', 'interpreter', 'latex')
% 
% % === Additional: Compare with fixed n = -0.2 (flat-plate prediction) ===
% fprintf('=== Fixed Exponent Fit (n = -0.2) ===\n')
% n_fixed = -0.2;
% % Solve for B: delta = B * (u_inf - c)^(-0.2)
% % B = delta / (u_inf - c)^(-0.2)
% B_fixed = mean(all_delta ./ (velocity_deficit.^n_fixed));
% fprintf('With n = -0.2 fixed, B = %.4f\n', B_fixed)
% 
% delta_predicted_fixed = B_fixed * velocity_deficit.^n_fixed;
% SS_res_fixed = sum((all_delta - delta_predicted_fixed).^2);
% R_squared_fixed = 1 - SS_res_fixed / SS_tot;
% fprintf('R^2 (fixed n) = %.4f\n', R_squared_fixed)





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





