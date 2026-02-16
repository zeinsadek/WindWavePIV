% WindWave paper figures: Curvilinear Phase Averages

clc; clear; close all
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions')
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps/slanCM')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
curvilinear_path = fullfile(project_path, 'curvilinear_new');
wave_parameters = readcell("Offshore_Waves.xlsx");

% Figure Folder
figure_folder = 'pdf_test7';

% Boundary layer detection parameters
boundary_layer_percent = 0.96;
left_mask = 4E-3;
right_mask = 234E-3;
threshold.thickness = 0.005;
threshold.displacement = 0.005;
threshold.momentum = 0.008;
smooth_kernel = 9;
gradient_tolerance = 1;
island_length = 5;
sgolay_length = 35;
n = 9;


% Cases
wind_speed = 'WT6';
waves = {'A', 'B', 'C', 'D'};

% Approximate wavelengths in mm for labeling plots
wavelengths.A = '410';
wavelengths.B = '313';
wavelengths.C = '189';
wavelengths.D = '124';

% Wavelengths in mm
wavelength_values.A = 410.8;
wavelength_values.B = 313.3;
wavelength_values.C = 189.6;
wavelength_values.D = 124.3;

% Wave steepnesses 
steepnesses.A = '0.180';
steepnesses.B = '0.211';
steepnesses.C = '0.305';
steepnesses.D = '0.267';

if ismember(wind_speed(end), {'4'})
    u_inf = 2.4181;
elseif ismember(wind_speed(end), {'6'})
    u_inf = 3.8709;
elseif ismember(wind_speed(end), {'8'})
    u_inf = 5.4289;
end


% Load all wave cases
for w = 1:length(waves)
    wave = waves{w};
    caze = strcat(wind_speed, '_WV', wave, '_AG0');
    disp(caze)

    % Load data
    tmp = load(fullfile(curvilinear_path, strcat(caze, '_CURVILINEAR.mat')));
    tmp = tmp.output;

    % Save data
    data.(caze) = tmp;
end
cazes = fields(data);
clc; fprintf('All %s cases loaded\n', wind_speed)

% clear caze tmp w wave no_wave_caze


% Loop through all waves and wind speeds
% wind_speeds = {'WT4', 'WT6', 'WT8'};
waves = {'A', 'B', 'C', 'D'};

% for s = 1:length(wind_speeds)
% 
%     % Loop through wind speeds
%     wind_speed = wind_speeds{s};
% 
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


clear amplitude caze column curvilinear displacement displacement_cum displacement_edge_mask displacement_integrand
clear horizontal_lines i idx island_length left_mask max_wave_profile momentum momentum_cum momentum_edge_mask momentum_integrand
clear n nan_mask phase right_mask s sgolay_length smooth_kernerl thickness thickness_edge_mask threshold u u_column u_normalized
clear vertical_lines w wave wave_profile wave_parameters wave_type wave_length x_BL zeta_column
clear gradient_tolerance project_path smooth_kernel wavelength curvilinear_path
clc;

%% Plot all phases single case

levels = 100;
linewidth = 1;

colorbar_fontsize = 10;
tickFontSize = 8;
labelFontSize = 10;


reordered_phases = [1,4,3,2];

wave_transparency = 0.25;

wave = 'C';
component = 'u';
caze = strcat(wind_speed, '_WV', wave, '_AG0');
wave_length = 1;

phase_labels = {'$\varphi = 0$', '$\varphi = \lambda / 4$', ...
                '$\varphi = \lambda / 2$', '$\varphi = 3 \lambda / 4$'};

if strcmp(component, 'u')
    labelComponent = "u_{\xi}";
elseif strcmp(component, 'v')
    labelComponent = "u_{\zeta}";
elseif strcmp(component, 'uu')
    labelComponent = "u_{\xi}' u_{\xi}'";
elseif strcmp(component, 'vv')
    labelComponent = "u_{\zeta}' u_{\zeta}'";
elseif strcmp(component, 'uv')
    labelComponent = "u_{\xi}' u_{\zeta}'";
else
    labelComponent = component;
end


% Plot
clc; close all;
totalFigure = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,4]);
t = tiledlayout(1, 4, 'TileSpacing', 'tight', 'padding', 'tight');

% Normalize velocities vs stresses
if ismember(component, {'u', 'v'})
    norm = u_inf;
    colorbar_label = sprintf('$ \\langle %s \\rangle \\mathbin{/} u_{\\infty}$', labelComponent);
else
    norm = u_inf^2;
    colorbar_label = sprintf('$ \\langle %s \\rangle \\mathbin{/} u_{\\infty}^2$', labelComponent);
end

% Select appropriate colormap
if ismember(component, {'u', 'uu', 'vv'})
    cmap = 'plasma';
else
    cmap = 'coolwarm';
end

% Loop through all phases
for p = 1:4

    phase = reordered_phases(p);
    h(p) = nexttile;

    % Make tick marks look fancy
    ax = gca;
    set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)

    tmp = data.(caze).phase(phase).(component) / norm;
    tmp(data.(caze).phase(phase).vertical_lines < 0) = nan;
    tmp(data.(caze).phase(phase).vertical_lines > 236) = nan;

    
    index = round(size(tmp,2)/2);

    hold on
    % Plot contour
    contourf(data.(caze).phase(phase).vertical_lines / wave_length, ...
             data.(caze).phase(phase).horizontal_lines / wave_length, ...
             tmp, ...
             levels, 'linestyle', 'none')


    % Plot wave cropping
    % plot(data.(caze).X(1,:) / wave_length, data.(caze).phase(phase).max_wave_profile / wave_length, 'linewidth', linewidth, 'color', 'red')

    % Plot Reference wave profile
    plot(data.(caze).X(1,:) / wave_length, data.(caze).phase(phase).wave_profile / wave_length, 'linewidth', linewidth, 'color', 'black')

    hFill = patch( ...
    [data.(caze).X(1,:), fliplr(data.(caze).X(1,:))], ...
    [data.(caze).phase(phase).wave_profile / wave_length, -100 * ones(size(data.(caze).phase(phase).wave_profile / wave_length))], ...
    'k', ...
    'FaceAlpha', wave_transparency, ...      
    'EdgeColor', 'none', ...    
    'HandleVisibility', 'off'); 
    uistack(hFill, 'bottom')

    hold off

    axis equal
    ylim([-20 / wave_length, 200 / wave_length])
    xlim([0, 235])
    clim([0.4, 1])
    title(phase_labels{p}, 'interpreter', 'latex', 'FontSize', labelFontSize)

    if p == 1
        ylabel('$y$ [mm]', 'interpreter', 'latex', 'FontSize', labelFontSize)
    end

    if p ~= 1
        ax.YTickLabel = [];
    end

    if p == 4
        colormap(flipud(slanCM(cmap)))
        C = colorbar;
        C.Label.String = colorbar_label;
        C.FontSize = tickFontSize;
        C.Label.Interpreter = 'latex';
        C.Label.FontSize = colorbar_fontsize;
        C.TickLabelInterpreter = 'latex';
        C.Ticks = 0.2:0.2:1;
    end

    clim([0.4, 1])
    clear c
end

% Sync color ranges
ax = findall(gcf, 'Type', 'axes');
clims = cell2mat(get(ax, 'CLim'));
% global_clim = [min(clims(:,1)), max(clims(:,2))];
global_clim = [0.4, max(clims(:,2))];
set(ax, 'CLim', global_clim);

% Add a shared x-axis label to the entire layout
xlabel(t, '$x$ [mm]', 'interpreter', 'latex', 'FontSize', labelFontSize);
linkaxes(h,'xy')
% xlim([0, max(data.(caze).X(1,:))])

% Save figure
pause(3)
fig_folder = fullfile('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new', figure_folder, '/Phase');
fig_name = strcat(caze(1:end - 4), '_AllPhases_', component, '.pdf');
exportgraphics(totalFigure, fullfile(fig_folder, fig_name), 'Resolution', 600, 'ContentType', 'image')
close all 

clear ax C clims cmap colorbar_fontsize colorbar_label component global_clim h 
clear levels linewidth norm phase wave
clc


%% Plot single phase all cases

wave_transparency = 0.25;
% wave_colors = {'#FB3640', '#FFC324', '#09814A', '#1BE7FF'};
wave_colors = {'#FE6202', '#DC2680', '#775EEF', '#648FFF'};

levels = 100;
linewidth = 1;
linePlotLineWidth = 1.5;

colorbar_fontsize = 10;
tickFontSize = 8;
labelFontSize = 10;
legendFontSize = 8;

phase = 2;
component = 'vv';

% Convert to curvilinear notation
if strcmp(component, 'u')
    labelComponent = "u_{\xi}";
elseif strcmp(component, 'v')
    labelComponent = "u_{\zeta}";
elseif strcmp(component, 'uu')
    labelComponent = "u_{\xi}' u_{\xi}'";
elseif strcmp(component, 'vv')
    labelComponent = "u_{\zeta}' u_{\zeta}'";
elseif strcmp(component, 'uv')
    labelComponent = "u_{\xi}' u_{\zeta}'";
else
    labelComponent = component;
end


clc; close all;
totalFigure = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,8]);
tiledlayout(2, 4, 'TileSpacing', 'tight', 'padding', 'tight');

% Normalize velocities vs stresses
if ismember(component, {'u', 'v'})
    norm = u_inf;
    colorbar_label = sprintf('$\\langle %s \\rangle \\mathbin{/} u_{\\infty} $', labelComponent);
else
    norm = u_inf^2;
    colorbar_label = sprintf('$\\langle %s \\rangle \\mathbin{/} u_{\\infty}^2$', labelComponent);
    if strcmp(component, 'uv')
        colorbar_label = sprintf('$-\\langle %s \\rangle \\mathbin{/} u_{\\infty}^2$', labelComponent);
    end
end

% Select appropriate colormap
if ismember(component, {'u'})
    cmap = 'plasma';
else
    cmap = 'RdPu';
end

% Lopp through different waves
for c = 1:length(cazes)

    disp(c)
    caze = cazes{c};
    wave = caze(7);

    % wave_length = wavelength_values.(wave);
    wave_length = 1;

    % Make tick marks look fancy
    h(c) = nexttile;
    ax = gca;
    set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)

    tmp = data.(caze).phase(phase).(component) / norm;
    tmp(data.(caze).phase(phase).vertical_lines < 0) = nan;
    tmp(data.(caze).phase(phase).vertical_lines > 236) = nan;

    if strcmp(component, 'uv')
        tmp = -1 * tmp;
    end
    hold on

    % Plot contour
    contourf(data.(caze).phase(phase).vertical_lines / wave_length, ...
             data.(caze).phase(phase).horizontal_lines / wave_length, ...
             tmp, ...
             levels, 'linestyle', 'none')

    % Plot line of constant \xi
    index = round(size(tmp,2)/2);
    % P = plot(data.(caze).phase(phase).vertical_lines(:, index) / wave_length, data.(caze).Y(:,1) / wave_length, ...
    %          'color', 'black', 'linestyle', ':', 'linewidth', linewidth, 'HandleVisibility', 'off');
    % P.Color(4) = 0.5;

    % Plot wave crop
    % plot(data.(caze).X(1,:) / wave_length, data.(caze).phase(phase).max_wave_profile / wave_length, 'linewidth', linewidth, 'color', 'red')

    % Plot reference wave
    plot(data.(caze).X(1,:) / wave_length, data.(caze).phase(phase).wave_profile / wave_length, 'linewidth', linewidth, 'color', 'black')

    hFill = patch( ...
    [data.(caze).X(1,:), fliplr(data.(caze).X(1,:))], ...
    [data.(caze).phase(phase).wave_profile / wave_length, -100 * ones(size(data.(caze).phase(phase).wave_profile / wave_length))], ...
    'k', ...
    'FaceAlpha', wave_transparency, ...      
    'EdgeColor', 'none', ...    
    'HandleVisibility', 'off'); 
    uistack(hFill, 'bottom')
    hold off


    % Scale axes to wavelength but keep aspect ratio 1:1
    axis equal

    % absolute limits in [mm]
    % xmin = -10;
    % xmax = 250;
    xmin = 0;
    xmax = 237;
    ymin = -20;
    ymax = 200;

    x_range = (xmax) - (xmin);
    y_range = (ymax) - (ymin);

    aspect = x_range / y_range;

    xlim([xmin/wave_length, xmax/wave_length])
    ylim([ymin/wave_length, ymax/wave_length])
    xticks(0:100:250)
    yticks(0:100:250)

    % Titles
    if ismember(wave, waves)
        title_label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelengths.(wave), steepnesses.(wave));
    else
        title_label = sprintf('$\\lambda_{0}$');
    end
    title(title_label, 'interpreter', 'latex', 'fontsize', labelFontSize)

    if c == 1
        ylabel('$y$ [mm]', 'interpreter', 'latex', 'FontSize', labelFontSize)
    end

    if c == 2
        xlabel('     ', 'interpreter', 'latex', 'FontSize', labelFontSize)
    end

    if c ~= 1
        ax.YTickLabel = [];
    end

    if c == 4

        if ismember(component, {'u'})
            colormap(flipud(slanCM(cmap)))
        else
            colormap(slanCM(cmap))
        end

        % colormap(slanCM(cmap))

        C = colorbar;
        C.Label.String = colorbar_label;
        C.FontSize = tickFontSize;
        C.Label.Interpreter = 'latex';
        C.Label.FontSize = colorbar_fontsize;
        C.TickLabelInterpreter = 'latex';

        if ismember(component, {'uu', 'vv', 'uv'})
            C.Ruler.Exponent = -3;
            C.Ruler.TickLabelFormat = '%2.0f';
        end
    end

    if strcmp(component, 'u')
        clim([0.4, 1])
        C.Ticks = 0.2:0.2:1;
    elseif strcmp(component, 'uv')
        clim([0, 9E-3])
    elseif strcmp(component, 'vv')
        clim([0, 12E-3])
    end

    clear c
end

% Sync color ranges
ax = findall(gcf, 'Type', 'axes');
clims = cell2mat(get(ax, 'CLim'));
global_clim = [min(clims(:,1)), max(clims(:,2))];
set(ax, 'CLim', global_clim);

% Plot profiles
ax_profiles = nexttile([1,4]);

% Make tick marks look fancy
ax = gca;
set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', tickFontSize)

hold on
for c = 1:length(cazes)

    caze = cazes{c};
    wave = caze(7);

    tmp = data.(caze).phase(phase).(component) / norm;
    tmp(data.(caze).phase(phase).vertical_lines < 0) = nan;
    tmp(data.(caze).phase(phase).vertical_lines > 236) = nan;

    index = round(size(tmp,2)/2);

    y_profile = data.(caze).phase(phase).horizontal_lines(:, index);
    profile = tmp(:, index);

    valid = ~isnan(y_profile) & ~isnan(profile);

    y_profile = y_profile(valid);
    profile = profile(valid);

    buffer = 2;
    y_profile = y_profile(1:end - buffer);
    profile = profile(1:end - buffer);

    % Titles
    if ismember(wave, waves)
        label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelengths.(wave), steepnesses.(wave));
    else
        label = sprintf('$\\lambda_{0}$');
    end

    % Color profile based on wave
    if strcmp(wave, 'A')
        color = wave_colors{1};
    elseif strcmp(wave, 'B')
        color = wave_colors{2};
    elseif strcmp(wave, 'C')
        color = wave_colors{3};
    elseif strcmp(wave, 'D')
        color = wave_colors{4};
    else 
        color = 'black';
    end

    % Normalize 'y' scale by boundary layer thickness
    thickness = integral.(caze)(phase).filtered.thickness(1,index) * 1E3;

    if strcmp(component, 'uv')
        profile = -1 * profile;
    end

    % Plot a slice
    plot(y_profile / thickness, profile, ...
         'linewidth', linePlotLineWidth, 'displayname', label, 'color', color)

end
hold off
if strcmp(component, 'uv')
    ylim([0, 9E-3])
end

xlabel('$\zeta \mathbin{/} \delta$', 'interpreter', 'latex', 'FontSize', labelFontSize)

% make y axes nicer for small values
if ismember(component, {'uu', 'vv', 'uv'})
    ax = gca;
    ax.YAxis.Exponent = -3;
    ax.YAxis.TickLabelFormat = '%2.0f';
end
xlim([0, 2])
% xticks(0:0.25:2)
xticks(0:0.5:2.5)
ylabel(colorbar_label, 'interpreter', 'latex', 'FontSize', labelFontSize);



if strcmp(component, 'u')
    loc = 'southeast';
elseif strcmp(component, 'vv')
    loc = 'northeast';
elseif strcmp(component, 'uu')
    loc = 'northeast';
elseif strcmp(component, 'uv')
    loc = 'northeast';
end


lgd = legend('location', loc, 'interpreter', 'latex', 'box', 'off', 'fontsize', legendFontSize);
lgd.IconColumnWidth = 19;


% Add a centered xlabel across the top row
annotation(totalFigure, 'textbox', [0.21, 0.50, 0.5, 0.04], ...
    'String', '$x$ [mm]', ...
    'Interpreter', 'latex', ...
    'FontSize', labelFontSize, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'EdgeColor', 'none');


% Add (a) (b) annotations
addPanelLabelsFixed(totalFigure, [h(1), ax_profiles], {'a','b'}, ...
    'FontSize', 10, 'OffsetPts', [-30 3]);


% Save figure
pause(3)
fig_folder = fullfile('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures/new/', figure_folder, 'Phase');
fig_name = strcat(wind_speed, '_Phase', num2str(phase), '_', component, '_ThicknessScaled.pdf');
exportgraphics(totalFigure, fullfile(fig_folder, fig_name), 'Resolution', 600, 'ContentType', 'image')
close all 


clear ax C caze clims cmap colorbar_fontsize colorbar_label component global_clim h 
clear levels linewidth norm phase wave h






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






