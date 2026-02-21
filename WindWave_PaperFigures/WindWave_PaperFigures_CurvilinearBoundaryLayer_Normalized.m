%%% Curvilinear Integral Boundary Layer Parameters

addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
curvilinear_path = fullfile(project_path, 'curvilinear');
wave_parameters = readcell("Offshore_Waves.xlsx");

figure_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/pdf_test3';


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


%% No waves

BL_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/boundary_layer/';
% Load cases without waves
for s = 1:length(wind_speeds)
    % Case
    no_wave_case = strcat(num2str(wind_speeds{s}), '_WV0_AGP');
    fprintf('Loading Case: %s...\n', no_wave_case)
    
    % Paths
    BL_path = fullfile(BL_folder, strcat(no_wave_case, '_BL.mat'));
    
    % Save BL to Structure
    BL_temp = load(BL_path);
    BL_temp = BL_temp.output;
    integral_no_waves.(no_wave_case) = BL_temp;
end


%%

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
parameter = 'momentum';
version = 'filtered';
linestyles = {'-.', '--', '-'};
lw = 2;
alpha = 0.5;

% Scale reference wave based on parameter
if strcmp(parameter, 'thickness')
    scale = 50;
    vshift = 1.2;
elseif strcmp(parameter, 'displacement')
    scale = 80;
    vshift = 1;
elseif strcmp(parameter, 'momentum')
    scale = 80;
    vshift = 1;
elseif strcmp(parameter, 'shape')
    scale = 80;
    vshift = 1;
end


tickFontSize = 8;
labelFontSize = 16;
legendFontSize = 12;

close all; 
ax = figure('color', 'white', 'units', 'centimeters', 'position', [2,4,30,10]);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)

sgtitle
hold on
% yline(1, 'linewidth', 2, 'HandleVisibility', 'off', 'Color', 'black', 'Alpha', 0.5);
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

        wavelength = integral.(caze).wavelength; 
        amplitude = integral.(caze).amplitude;
         
        label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelengths.(wave), steepnesses.(wave));
        
        data = integral.(caze)(phase).(version).(parameter);
        no_wave_mean_thickness = mean(integral_no_waves.(strcat(wind_speed, '_WV0_AGP')).ensemble.(parameter), 'all', 'omitnan') * 1E-3;

        H(w) = plot((integral.(caze)(phase).x - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, data ./ no_wave_mean_thickness, ...
                    'color', wave_colors{w}, 'linewidth', lw, 'linestyle', linestyles{s}, ...
                    'displayname', label, 'HandleVisibility', vis);

        % Plot phase for reference
        if wave == 'D' && s == 1
            reference = plot((integral.(caze)(phase).wave.x  - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, (scale * integral.(caze)(phase).wave.wave_profile) + vshift, ...
                       'color', 'black', 'linewidth', lw, 'linestyle', ':', 'HandleVisibility', 'off');
            reference.Color(4) = alpha;
        end
    end
end




%%% FOR LEGEND
% Dummy line to add a break
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
    P = plot(nan, nan, 'linestyle', linestyles{s}, ...
         'linewidth', lw, 'color', 'black', 'DisplayName', sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf));
    P.Color(4) = 0.5;
end

% Vertical limits
if strcmp(parameter, 'thickness')
    ylim([0.7, 1.7])
elseif strcmp(parameter, 'displacement')
    ylim([0, 2.5])
elseif strcmp(parameter, 'momentum')
    ylim([0, 2])
end

hold off


% Legend
lgd = legend('interpreter', 'latex', 'fontsize', 10, 'location',  'eastoutside', 'box', 'off');


% Axes labels
if strcmp(parameter, 'thickness')
    vertLabel = '$\delta / \delta_{0}$';
elseif strcmp(parameter, 'displacement')
    vertLabel = '$\delta^* / \delta^*_{0}$';
elseif strcmp(parameter, 'momentum')
    vertLabel = '$\theta / \theta_{0}$';
elseif strcmp(parameter, 'shape')
    vertLabel = '$H$';
    limits = [1, 1.7];
end

xlabel('$\xi / \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(vertLabel, 'interpreter', 'latex', 'fontsize', labelFontSize)
xticks(-1:0.25:1)


% Save figure
% pause(3)
% figure_name = strcat('BoundaryLayer_Curvilinear', [upper(parameter(1)), lower(parameter(2:end))], '_Phase', num2str(phase), '_Normalized.pdf');
% exportgraphics(ax, fullfile(figure_folder, 'BoundaryLayer', figure_name), 'Resolution', 600, 'ContentType', 'image');
% close all
% clc; fprintf('Generated figure: %s\n\n', figure_name)

clear alpha limits caze data H linestyles lw parameter phase reference s scale w wave wind_speed version ans vertLabel; clc








%% Plot displacement thickness next to momentum thickness


 
phase = 1;
version = 'filtered';
linestyles = {'-.', '--', '-'};
lw = 2;
alpha = 0.5;

tickFontSize = 8;
labelFontSize = 16;
legendFontSize = 12;


close all; 
ax = figure('color', 'white', 'units', 'centimeters', 'position', [2,4,30,7]);
tile = tiledlayout(1,2, 'tilespacing', 'compact');

% Displacement
parameter = 'displacement';

% Scale reference wave based on parameter
if strcmp(parameter, 'thickness')
    scale = 50;
    vshift = 1.2;
elseif strcmp(parameter, 'displacement')
    scale = 80;
    vshift = 1;
elseif strcmp(parameter, 'momentum')
    scale = 80;
    vshift = 1;
elseif strcmp(parameter, 'shape')
    scale = 80;
    vshift = 1;
end

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
        no_wave_mean_thickness = mean(integral_no_waves.(strcat(wind_speed, '_WV0_AGP')).ensemble.(parameter), 'all', 'omitnan') * 1E-3;

        H(w) = plot((integral.(caze)(phase).x - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, data / no_wave_mean_thickness, ...
                    'color', wave_colors{w}, 'linewidth', lw, 'linestyle', linestyles{s}, ...
                    'displayname', label);

        % Plot phase for reference
        if wave == 'D' && s == 1
            reference = plot((integral.(caze)(phase).wave.x  - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, (scale * integral.(caze)(phase).wave.wave_profile) + vshift, ...
                       'color', 'black', 'linewidth', lw, 'linestyle', ':', 'handlevisibility', 'off');
            reference.Color(4) = alpha;
        end
    end
end


% Vertical limits
if strcmp(parameter, 'thickness')
    ylim([0.7, 1.7])
elseif strcmp(parameter, 'displacement')
    ylim([0, 2.5])
elseif strcmp(parameter, 'momentum')
    ylim([0, 2])
end

hold off

% Axes labels
if strcmp(parameter, 'thickness')
    vertLabel = '$\delta / \delta_{0}$';
elseif strcmp(parameter, 'displacement')
    vertLabel = '$\delta^* / \delta^*_{0}$';
elseif strcmp(parameter, 'momentum')
    vertLabel = '$\theta / \theta_{0}$';
elseif strcmp(parameter, 'shape')
    vertLabel = '$H$';
    limits = [1, 1.7];
end

xlabel('$\xi / \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(vertLabel, 'interpreter', 'latex', 'fontsize', labelFontSize)
xticks(-1:0.25:1)
% ylim([0, 0.035])



%%% Momentum
parameter = 'momentum';

% Scale reference wave based on parameter
if strcmp(parameter, 'thickness')
    scale = 50;
    vshift = 1.2;
elseif strcmp(parameter, 'displacement')
    scale = 80;
    vshift = 1;
elseif strcmp(parameter, 'momentum')
    scale = 80;
    vshift = 1;
elseif strcmp(parameter, 'shape')
    scale = 80;
    vshift = 1;
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
        no_wave_mean_thickness = mean(integral_no_waves.(strcat(wind_speed, '_WV0_AGP')).ensemble.(parameter), 'all', 'omitnan') * 1E-3;
        H(w) = plot((integral.(caze)(phase).x - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, data / no_wave_mean_thickness, ...
                    'color', wave_colors{w}, 'linewidth', lw, 'linestyle', linestyles{s}, ...
                    'displayname', label);

        % Plot phase for reference
        if wave == 'D' && s == 1
            reference = plot((integral.(caze)(phase).wave.x  - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, (scale * integral.(caze)(phase).wave.wave_profile) + vshift, ...
                       'color', 'black', 'linewidth', lw, 'linestyle', ':', 'handlevisibility', 'off');
            reference.Color(4) = alpha;
        end

    end

end


% Vertical limits
if strcmp(parameter, 'thickness')
    ylim([0.7, 1.7])
elseif strcmp(parameter, 'displacement')
    ylim([0, 2.4])
elseif strcmp(parameter, 'momentum')
    ylim([0, 2])
end

hold off

% Axes labels
if strcmp(parameter, 'thickness')
    vertLabel = '$\delta / \delta_{0}$';
elseif strcmp(parameter, 'displacement')
    vertLabel = '$\delta^* / \delta^*_{0}$';
elseif strcmp(parameter, 'momentum')
    vertLabel = '$\theta / \theta_{0}$';
elseif strcmp(parameter, 'shape')
    vertLabel = '$H$';
    limits = [1, 1.7];
end

xlabel('$\xi / \lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(vertLabel, 'interpreter', 'latex', 'fontsize', labelFontSize)
xticks(-1:0.25:1)
linkaxes(h, 'xy')





% Save figure
figure_name = strcat('BoundaryLayer_Curvilinear_Displacement_vs_Momentum_Phase', num2str(phase), '_Normalized.pdf');
pause(3)
exportgraphics(tile, fullfile(figure_folder, 'BoundaryLayer', figure_name), 'Resolution', 300, 'Resolution', 600, 'ContentType', 'image');
close all
fprintf('Generated figure: %s\n\n', figure_name)




%% Boundary layer thickness at crest versus u_{wave} / u_{\infty}


phase = 1;
parameter = 'thickness';
version = 'filtered';
linestyles = {'-.', '--', '-'};
lw = 2;
alpha = 0.5;

wind_speed_markers = {'^', 'square', 'o'};

tickFontSize = 8;
labelFontSize = 16;
legendFontSize = 12;

clc; close all; 
% ax = figure('color', 'white', 'units', 'centimeters', 'position', [2,4,30,10]);
ax = figure('color', 'white');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)

sgtitle
hold on

for w = 1:length(waves)
    wave = waves{w};
   

    for s = 1:length(wind_speeds)
        wind_speed = wind_speeds{s};
        caze = strcat(wind_speed, '_WV', wave, '_AG0');
    
        if ismember(wind_speed(end), {'4'})
            u_inf = 2.4181;
        elseif ismember(wind_speed(end), {'6'})
            u_inf = 3.8709;
        elseif ismember(wind_speed(end), {'8'})
            u_inf = 5.4289;
        end


    

        wavelength = integral.(caze).wavelength; 
        amplitude = integral.(caze).amplitude;
        frequency = frequencies.(wave);
        steepness = (2 * pi * amplitude) / wavelength;
        wave_speed = frequency * wavelength;
        wave_age = wave_speed / u_inf;

        fprintf('%s: Wave speed = %1.3f\n', wave, wave_speed)
        fprintf('Wave_speed / u_inf = %1.3f\n\n', wave_age)

        
        % label = sprintf('$\\lambda_{%s}, \\hspace{1mm} u_{\\infty} = %1.2f$ m/s', wavelengths.(wave), u_inf);
        label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelengths.(wave), steepnesses.(wave));
        
        data = integral.(caze)(phase).(version).(parameter);
        idx = round(length(data) / 2);

        % H(w) = plot((integral.(caze)(phase).x - mean(integral.(caze)(phase).wave.x, 'all', 'omitnan')) / integral.(caze)(1).wavelength, data, ...
        %             'color', wave_colors{w}, 'linewidth', lw, 'linestyle', linestyles{s}, ...
        %             'displayname', label, 'HandleVisibility', vis);

        
        scatter(wave_age, data(idx), 100,  wind_speed_markers{s}, 'filled', 'MarkerFaceColor', wave_colors{w})

        tmp(s) = wave_age;
        tmp2(s) = data(idx);

    end

    plot(tmp, tmp2, 'linewidth', 2, 'color', wave_colors{w})
end

hold off
ylim([0.05, 0.15])
% xlim([1, 12])
% xscale('log')
% yscale('log')












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





