%% 2D Offshore PIV Processing: Ondrej Fercak, Zein Sadek, 1/2023
% Converts DaVis vector data (.vc7) files into a Matlab file

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Functions/');
experiment_log  = readcell('Offshore_Inflow_Log.xlsx');
wave_parameters = readcell('Offshore_Waves.xlsx');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data paths
project_path   = '/Volumes/WT4_WVD_AG0/WT4_WVD_AG0';
recording_name = 'WT4_WVD_AG0';
processing     = '/PIV_MP(2x24x24_50%ov)/PostProc';
still_name     = strcat(experiment_log{find(strcmp(experiment_log, recording_name) == 1), 2}, '_waterlevel');
inpt_name      = recording_name;
% inpt_name      = strcat(recording_name, '_2');

% Experiment Specifics
tunnel_freq    = recording_name(strfind(recording_name, 'WT') + 2);
wave_type      = recording_name(strfind(recording_name, 'WV') + 2);
% active_grid    = recording_name(strfind(recording_name, 'AG') + 2: strfind(recording_name, 'AG') + 3);

% Image paths
still_path     = strcat('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Offshore/Water_Level/', still_name, '/Correction');
raw_path       = strcat(project_path, '/', recording_name, '/Correction');
piv_path       = strcat(project_path, '/', recording_name, processing);

% Save paths
results_path   = '/Users/zeinsadek/Desktop/Experiments/Offshore/Hopkins/Initial_Results/';
mtlb_file      = strcat(results_path, 'data'   , '/', inpt_name, '_DATA.mat');
mean_file      = strcat(results_path, 'means'  , '/', inpt_name, '_MEANS.mat');
phase_file     = strcat(results_path, 'phase'  , '/', inpt_name, '_PHASE.mat');
figure_file    = strcat('/Users/zeinsadek/Desktop/Experiments/Offshore/NSF/figures/', recording_name);

if ~exist(figure_file, 'dir')
    mkdir(figure_file)
end

% Edge Detection Constants
std_tol            = 2.5;
background_removal = 1000;
canny_lower        = 0.1;
canny_upper        = 0.4;
grad_tol           = 0.4;
nan_dist           = 5;

% Filtering Parameters
left_bound_value   = -115;
right_bound_value  = 115;
wave_buffer        = 0.5;
cmap_perc          = 0.05;
phase_tolerance    = 5;

% Wave Parameters
if wave_type == '0'
    wave_length    = 0;
    wave_amplitude = 0;
    phase_offset   = [0, wave_length/4, wave_length/2, 3*wave_length/4];
else
    wave_length    = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 2};
    wave_amplitude = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 3};
    wave_frequency = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 4};
    phase_offset   = [0, wave_length/4, wave_length/2, 3*wave_length/4];
end

% Video Settings
num_frames = 50;
fps        = 3;

% Constants Structure
constants.recording_name     = inpt_name;
constants.figure_file        = figure_file;
constants.num_frames         = num_frames;
constants.fps                = fps;
constants.cmap_perc          = cmap_perc;

constants.background_removal = background_removal;
constants.std_tol            = std_tol;
constants.canny_lower        = canny_lower;
constants.canny_upper        = canny_upper;
constants.grad_tol           = grad_tol;
constants.nan_dist           = nan_dist;

constants.left_bound_value   = left_bound_value;
constants.right_bound_value  = right_bound_value;
constants.phase_tolerance    = phase_tolerance; 
constants.phase_offset       = phase_offset;

constants.wave_type          = wave_type;
constants.wave_amplitude     = wave_amplitude;
constants.wave_length        = wave_length;
constants.vertical_offset    = waterlevel(still_path, still_name, constants);
constants.wave_buffer        = wave_buffer;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE AVERAGE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if constants.wave_type == '0'
    fprintf('* Skipping <phaseaverage>\n\n')
    phase_average = 0;
else
    if exist(phase_file, 'file')
        fprintf('* Loading PHASES from File\n')
        phase_average = load(phase_file);
        phase_average = phase_average.phaseaverage_output;
    else
        phase_average = phaseaverage(raw_path, constants, phase_file);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DAVIS TO MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist(mtlb_file, 'file')
    fprintf('* Loading DATA from File\n')
    data = load(mtlb_file);
    data = data.output;
else
     data = vector2matlab2D(raw_path, piv_path, constants, mtlb_file);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB DATA TO  ENSEMBLE/PHASE MEANS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist(mean_file, 'file')
    fprintf('* Loading MEANS from File\n')
    means = load(mean_file); 
    means = means.output;
else
     means = data2means2D(mean_file, data, phase_average);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOUNDARY LAYER THICKNESS (Ensemble)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Average over top row of image, 3 pixels down
u_inf = mean(means.ensemble.u(:, 1:3), 'all', 'omitnan');

% Remove lower part of mean 
max_wave_crop = max(data.waves, [], 'all') - constants.vertical_offset;

% Define boundary layer edge velocity
BL_edge = 0.99 * u_inf;
BL_crop = means.ensemble.u;

% Remove values above boundary layer
BL_crop(means.ensemble.u > BL_edge) = nan;

% Set Y zero at water level
Y = means.Y - (constants.vertical_offset + wave_buffer);
% y = unique(Y);
y = flipud(unique(Y));

% set X zero at left
X = means.X - (left_bound_value + 1);
x = unique(X);

% Remove lower part of image
BL_crop(Y < max_wave_crop) = nan;
means.ensemble.u(Y < max_wave_crop) = nan;

% BL Thickness
[M, ensemble_I] = max(BL_crop, [], 1, 'omitnan');
BL = Y(ensemble_I);

% Plot
ax = figure('Name', 'BL Contour');
hold on
contourf(X, Y, means.ensemble.u / u_inf, 500, 'Linestyle', 'none')
axis equal
yline(0, 'Color', 'red', 'LineWidth', 4)
hold off
C = colorbar();
caxis([0.55, 1]);
% C.Label.String = '${u}$ [m/s]';
% C.Label.Interpreter = 'Latex';
xlabel('x [mm]')
ylabel('z [mm]')
xlim([0, 2 * (right_bound_value - 1)])
ylim([-10, 210])
title('Ensemble')
% 
% figure_name = strcat('NSF_', recording_name, '_ensembleavg.png');
% exportgraphics(ax, strcat(figure_file, '/', figure_name))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLACEMENT THICKNESS (Ensemble)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DT = zeros(length(x), 1);

for i = 1:length(x)
    
    u         = BL_crop(:, i);
    mask      = ~isnan(u);
    u         = u(mask);
    y_mask    = y(mask);
    integrand = 1 - (u / u_inf);

    DT(i, 1)  = trapz(integrand, y_mask);
    
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOMENTUM THICKNESS (Ensemble)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MT = zeros(length(x), 1);

for i = 1:length(x)
    
    u         = BL_crop(:, i);
    mask      = ~isnan(u);
    u         = u(mask);
    y_mask    = y(mask);
    integrand = (u / u_inf) .* (1 - (u / u_inf));

    MT(i, 1)  = trapz(integrand, y_mask);
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENSEMBLE AVERAGE PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colors = colormap(jet(length(x)));

figure('Name', 'Integrands')
for i = 1:length(x)
    
    u         = BL_crop(:, i);
    mask      = ~isnan(u);
    u         = u(mask);
    y_mask    = y(mask);

    subplot(1,3,1)
    hold on
    plot(u, y_mask, 'Color', colors(i,:))
    xlabel('${u}$', 'Interpreter', 'Latex', 'FontSize', 18)
    ylabel('y [mm]')
    title('Streamwise Boundary Layer Velocity Profiles')
    ylim([0, 90])
    
    subplot(1,3,2)
    hold on
    plot(1 - (u/u_inf), y_mask, 'Color', colors(i,:))
    xlabel('${1 - \frac{u}{u_{\infty}}}$', 'Interpreter', 'Latex', 'FontSize', 18)
    title('Displacement Thickness Integrand')
    ylim([0, 90])
    
    subplot(1,3,3)
    hold on
    plot((u/u_inf) .* (1 - (u/u_inf)), y_mask, 'Color', colors(i,:))
    xlabel('$\frac{u}{u_{\infty}}\left(1 - \frac{u}{u_{\infty}}\right)$', 'Interpreter', 'Latex', 'FontSize', 18)
    title('Momentum Thickness Integrand')
    ylim([0, 90])
    
end
hold off
hold off
hold off

%%

figure('Name', 'Ens BL Lengths Contour');
hold on
contourf(X, Y, BL_crop, 500, 'LineStyle', 'none', 'DisplayName', 'u')
C = colorbar();
C.Label.String = '${u}$ [m/s]';
C.Label.Interpreter = 'Latex';
yline(0, 'Color', 'red', 'LineWidth', 3, 'DisplayName', 'Wave')
plot(x, BL, 'Color', 'black', 'LineWidth', 2, 'DisplayName', '\delta')
plot(x, DT, 'Color', 'black', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', '\delta^{*}')
plot(x, MT, 'Color', 'black', 'LineStyle', '-.', 'LineWidth', 2, 'DisplayName', '\theta')
hold off
axis equal
xlim([0, 2 * (right_bound_value - 1)])
ylim([-10, 210])
title(strcat(recording_name, ': Ensemble Avg Boundary Layer'), 'Interpreter', 'none')
xlabel('x [mm]', 'Interpreter', 'Latex')
ylabel('y [mm]', 'Interpreter', 'Latex')
legend()
% figure_name = strcat(recording_name, '_ensembleavg_BL.png');
% exportgraphics(ax, strcat(figure_file, '/', figure_name))


figure('Name', 'Ens BL Lengths Plot')
hold on
plot(x, BL, 'DisplayName', 'Boundary Layer Thickness')
plot(x, DT, 'DisplayName', 'Displacement Thickness')
plot(x, MT, 'DisplayName', 'Momentum Thickness')
hold off
axis equal
xlim([0, 220])
ylim([0, 80])
xlabel('x [mm]')
ylabel('z [mm]')
legend()
title('Ensemble Average Boundary Layer Thicknesses')

figure('Name', 'Ens BL Shape Factor')
plot(x, DT ./ MT)
xlim([0, 220])
xlabel('x [mm]')
ylabel('\delta^{*} / \theta')
title('Ensemble Average Shape Factor')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE AVERAGE CONTOURS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cos_fit = @(b, v) wave_amplitude * cos(2 * pi * (v - b(1)) / wave_length);


phase_BL = zeros(length(x), 4);
phase_DT = zeros(length(x), 4);
phase_MT = zeros(length(x), 4);


ax = figure('Name', 'Phase Contours');

for i = 1:4
    
    % Define boundary layer edge velocity
    BL_edge = 0.99 * u_inf;
    BL_crop = means.phase(i).u;
    % Remove values above boundary layer
    BL_crop(means.ensemble.u > BL_edge) = nan;

    % Remove lower part of image
    ref = cos_fit(phase_offset(i), unique(means.X)) + 1;
    BL_crop(Y < ref.') = nan;

    % BL Thickness
    [M, phase_I] = max(BL_crop, [], 1, 'omitnan');
    BL = Y(phase_I);
    phase_BL(:,i) = BL;
    
    % Displacement & Momentum Thickness
    for j = 1:length(x)
        
        u            = BL_crop(:, j);
        mask         = ~isnan(u);
        u            = u(mask);
        y_mask       = y(mask);
        
        DT_integrand = 1 - (u / u_inf);
        MT_integrand = (u / u_inf) .* (1 - (u / u_inf));
        
        phase_DT(j, i)  = trapz(DT_integrand, y_mask);
        phase_MT(j, i)  = trapz(MT_integrand, y_mask);
   
    end
    
    subplot(2,2,i)
    title(strcat('Phase', {' '}, num2str(i)))
    hold on
    contourf(X, Y, BL_crop / u_inf, 500, 'Linestyle', 'none', 'DisplayName', 'u')
    C = colorbar();
%     caxis([-0.02, 0.02])
%     C.Label.String = 'u [m/s]';
%     C.Label.Interpreter = 'Latex';
    axis equal
    plot(x, ref + 1, 'Color', 'red', 'LineWidth', 3, 'DisplayName', 'Wave')
    hold off
    xlim([0, 2 * (right_bound_value - 1)])
    ylim([-10, 210])
    xlabel('x [mm]')
    ylabel('z [mm]')
    
    title(strcat('Phase', {' '}, num2str(i)), 'Interpreter', 'none')
    
end

% figure_name = strcat('NSF_', recording_name, '_phaseavg__BL_contour_v.png');
% exportgraphics(ax, strcat(figure_file, '/', figure_name))


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE AVERAGE WAVEFORMS (FOR SANITY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(20)
hold on
for i = 1:4
    
    ref = cos_fit(phase_offset(i), unique(means.X));
    [x_peaks_troughs_value, peaks_troughs] = findpeaks(abs(ref));
   
    subplot(4,1,i)
    plot(x, ref, 'LineWidth', 2)

    xlim([0, 220])
    ylim([-8, 8])
    title(num2str(i))
end
hold off

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE AVERAGE PROFILES PEAKS VS TROUGHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1:4
    
    % Define boundary layer edge velocity
    BL_edge = 0.99 * u_inf;
    BL_crop = means.phase(j).u;
    % Remove values above boundary layer
    BL_crop(means.ensemble.u > BL_edge) = nan;

    % Remove lower part of image
    ref = cos_fit(phase_offset(j), unique(means.X)) + 1;
    BL_crop(Y < ref.') = nan;

    % BL Thickness
    [M, phase_I] = max(BL_crop, [], 1, 'omitnan');
    BL = Y(phase_I);
    
    % Find peaks/troughs
    [x_peaks_troughs_value, peaks_troughs] = findpeaks(abs(ref));

    figure('Name', strcat('BL Integrand Phase ', num2str(j)));
    for i = 1:length(peaks_troughs)
        
        u         = BL_crop(:, peaks_troughs(i));
        mask      = ~isnan(u);
        u         = u(mask);
        y_mask    = y(mask);
        
        DT_integrand = 1 - (u / u_inf);
        MT_integrand = (u / u_inf) .* (1 - (u / u_inf));

        hold on
        hold on
        if sign(ref(peaks_troughs(i))) == 1
            subplot(1,2,1)
            
%             plot(u, y_mask, 'Color', colors(i,:), 'LineWidth', 2)
            hold on
            %area(DT_integrand, y_mask, min(y_mask), 'FaceColor', 'green')
            plot(MT_integrand, y_mask, 'LineWidth', 4)
            hold off
            xlabel('${1 - \frac{u}{u_{\infty}}}$', 'Interpreter', 'Latex', 'FontSize', 18)
            ylabel('y [mm]', 'Interpreter', 'Latex')
            ylim([-10, 80])
            xlim([0, 0.6])
            yline(wave_amplitude, 'Color', 'red', 'LineWidth', 2)
            title('Peaks')

        elseif sign(ref(peaks_troughs(i))) == -1
            subplot(1,2,2)
%             plot(u, y_mask, 'Color', colors(i,:), 'LineWidth', 2)
            hold on
            %area(DT_integrand, y_mask, min(y_mask), 'FaceColor', 'green')
            plot(MT_integrand, y_mask, 'LineWidth', 4)
            hold off
            xlabel('${1 - \frac{u}{u_{\infty}}}$', 'Interpreter', 'Latex', 'FontSize', 18)
            title(strcat('Phase', {' '}, num2str(j)))
            ylim([-10, 80])
            xlim([0, 0.6])
            yline(-wave_amplitude, 'Color', 'red', 'LineWidth', 2)
            title('Troughs')
        end
    end

    hold off
    hold off
    sgtitle(strcat(recording_name, ' Phase', {' '},  num2str(j), ': Normalized Velocity Deficit'), 'Interpreter', 'none')
    
%     figure_name = strcat(recording_name, '_phaseavg_', num2str(j), '_DT_integrand_filled.png');
%     exportgraphics(ax, strcat(figure_file, '/', figure_name))
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE AVERAGE DISPLACEMENT/MOMENTUM THICKNESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax1 = figure(8);
title('Displacement Thickness')
hold on
for i = 1:4
    name = strcat('Phase', {' '}, num2str(i));
    plot(x, smooth(phase_DT(:,i)), 'LineWidth', 2, 'DisplayName', name{1})
end
%plot(x, smooth(DT), 'Color', 'black', 'Linestyle', '--', 'LineWidth', 2, 'DisplayName', 'Ensemble')
hold off
xlabel('x [mm]', 'Interpreter', 'Latex')
ylabel('Displacement Thickness $(\delta^{*})$ [mm]', 'Interpreter', 'Latex')
xlim([0, 220])
legend('Location', 'NorthWest')
% figure_name = strcat(recording_name, '_phaseavg_DT_smooth.png');
% exportgraphics(ax1, strcat(figure_file, '/', figure_name))


ax2 = figure(9);
title('Momentum Thickness')
hold on
for i = 1:4
    name = strcat('Phase', {' '}, num2str(i));
    plot(x, smooth(phase_MT(:,i)), 'LineWidth', 2, 'DisplayName', name{1})
end
%plot(x, smooth(MT), 'Color', 'black', 'Linestyle', '--', 'LineWidth', 2, 'DisplayName', 'Ensemble')
hold off
xlabel('x [mm]', 'Interpreter', 'Latex')
ylabel('Momentum Thickness $(\theta)$ [mm]', 'Interpreter', 'Latex')
xlim([0, 220])
legend('Location', 'NorthWest')
% figure_name = strcat(recording_name, '_phaseavg_MT_smooth.png');
% exportgraphics(ax2, strcat(figure_file, '/', figure_name))

%%

ax2 = figure(9);
title('Momentum Thickness')
hold on
for i = 1:1
    cos_fit = @(b, v) wave_amplitude * cos(2 * pi * (v - b(1)) / wave_length);
    
    ref = cos_fit(phase_offset(i), x);
    scatter(ref, smooth(phase_DT(:,i)))
end

hold off
xlabel('x [mm]', 'Interpreter', 'Latex')
ylabel('Momentum Thickness $(\theta)$ [mm]', 'Interpreter', 'Latex')



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE AVERAGE SHAPE FACTOR VS WAVE PROFILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ax = figure(11);
title('Displacement Thickness')
hold on
for i = 1:4
    subplot(2,2,i)
    hold on
    plot(x, smooth(phase_DT(:,i) - DT), 'LineWidth', 2, 'DisplayName', 'DT')
    plot(x, smooth(phase_MT(:,i) - MT), 'LineWidth', 2, 'DisplayName', 'MT')
    plot(x, smooth(phase_DT(:,i) ./ phase_MT(:,i)), 'LineWidth', 2, 'DisplayName', 'H')
    hold off
    ylim([-2, 2])
    legend()
end
%plot(x, smooth(DT), 'Color', 'black', 'Linestyle', '--', 'LineWidth', 2, 'DisplayName', 'Ensemble')
hold off
xlim([0, 220])
legend('Location', 'NorthWest')
% sgtitle(strcat(recording_name, ': Phase Average Shape Factor'), 'Interpreter', 'none')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE AVERAGE SHAPE FACTOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

new_cos_fit = @(b, v) b(1) * cos(2 * pi * (v - b(3)) / wave_length) + b(2);

ax = figure(11);
hold on
for i = 1:4
    
    subplot(2,2,i)
    H = phase_DT(:,i) ./ phase_MT(:,i);
    hold on
    e = plot(x, DT ./ MT, 'Color', 'black', 'LineWidth', 2, 'Displayname', 'Ensemble');
    plot(x, H, 'Color', 'blue', 'LineWidth', 2, 'DisplayName', 'Phase')
%     e.Color(4) = 0.25;
    ref = new_cos_fit([(max(H) - min(H)) / 2, mean(H, 'omitnan'), phase_offset(i)], unique(means.X));
    p = plot(x, ref, 'Color', 'red', 'LineWidth', 2,'DisplayName', 'Wave');
    p.Color(4) = 0.2;
    
    
    hold off
    ylim([1.2, 1.75])
    xlim([0, 220])
    xlabel('x [mm]')
    ylabel('H', 'Interpreter', 'Latex')
    title(strcat('Phase', {' '}, num2str(i)));
    legend()
end

hold off

% figure_name = strcat('NSF_', recording_name, '_phaseavg_shapefactor_subplot.png');
% exportgraphics(ax, strcat(figure_file, '/', figure_name))


%%

ax = figure(11);
hold on

e = plot(x, smooth(DT ./ MT), 'Color', 'black', 'LineWidth', 2, 'Displayname', 'Ensemble');

for i = 1:4
    H = phase_DT(:,i) ./ phase_MT(:,i);
    plot(x, smooth(H), 'LineWidth', 2, 'DisplayName', strcat('Phase', num2str(i)))
end

hold off
ylim([1.25, 1.6])
xlim([0, 220])
xlabel('x [mm]')
ylabel('H', 'Interpreter', 'Latex')
legend()

% figure_name = strcat('NSF_', recording_name, '_phaseavg_shapefactor.png');
% exportgraphics(ax, strcat(figure_file, '/', figure_name))

%%
figure(13)
hold on
for i = 1:4
    subplot(2,2,i)
    
    new_cos_fit = @(b, v) b(1) * cos(2 * pi * (v - b(3)) / wave_length) + b(2);
    H = phase_DT(:,i) ./ phase_MT(:,i);
    hold on
    plot(x, phase_DT(:,i), 'LineWidth', 2, 'DisplayName', '\delta^{*}') 
    plot(x, phase_MT(:,i), 'LineWidth', 2, 'DisplayName', '\theta') 
    plot(x, H, 'LineWidth', 2, 'Color', 'blue', 'DisplayName', 'H')
%     ref = new_cos_fit([1, 0, phase_offset(i)], unique(means.X));
%     p = plot(x, ref, 'Color', 'red', 'Linestyle', '--', 'DisplayName', 'Scaled Wave Profile');
%     p.Color(4) = 0.5;
    hold off
    ylim([0, 18])
    xlim([0, 220])
    xlabel('x [mm]')
    title(strcat('Phase', {' '}, num2str(i)));
    legend()
end
hold off





