%% 2D Offshore PIV Processing: Ondrej Fercak, Zein Sadek, 1/2023
% Converts DaVis vector data (.vc7) files into a Matlab file

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

% Paths for MacBook
% addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
% addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Offshore/Offshore_Functions');

% Paths for MME Desktop
addpath('G:/Other computers/Zein MacBook Pro/Offshore/Offshore_Functions/');
addpath('C:/Users/Zein/Documents/MATLAB/readimx-v2.1.8-win64/');

experiment_log  = readcell('Offshore_Inflow_Log.xlsx');
wave_parameters = readcell('Offshore_Waves.xlsx');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cases = {
         % 'WT4_WVD_AGP10', ...
         % 'WT5_WVD_AGP30', ...
         % 'WT7_WVD_AGS10', ...
         % 'WT7_WVD_AGS30', ...
         'WT8_WVD_AGP10', ...
         'WT9_WVD_AGP30', ...
         'WT15_WVD_AGS10', ...
         'WT15_WVD_AGS30'};

drives = {
          % 'K', ...
          % 'I', ...
          % 'P', ...
          % 'N', ...
          'J', ...
          'H', ...
          'O', ...
          'M'};



for iter = 1:length(cases)


    caze = cases{iter};
    
    % project_path   = strcat('Q:\Offshore_Inflow_No_Waves_Passive_Grid\', caze);
    % recording_name = caze;

    project_path   = strcat(drives{iter}, ':\', caze, '\', caze);
    recording_name = caze;
    
    
    % If PostProc not exist, use first vectors
    if ~exist(strcat(project_path, '\PIV_MP(2x24x24_50%ov)\PostProc'), 'dir')
        processing = '\PIV_MP(2x24x24_50%ov)';
    else
        processing = '\PIV_MP(2x24x24_50%ov)\PostProc';
    end
    
    still_name     = strcat(experiment_log{find(strcmp(experiment_log, recording_name) == 1), 2}, '_waterlevel');
    inpt_name      = recording_name;
    
    % Experiment Specifics
    tunnel_freq    = recording_name(strfind(recording_name, 'WT') + 2);
    wave_type      = recording_name(strfind(recording_name, 'WV') + 2);
    
    %%% Processing Active Grid No Waves
    % wave_type = '0';
    % still_name = '2_13_2023_waterlevel';
    
    % Image paths
    % Path for MacBook
    % still_path     = strcat('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Offshore/Water_Level/', still_name, '/Correction');
    
    % Path for MME Desktop
    still_path     = strcat('G:/Other computers/Zein MacBook Pro/Offshore/Water_Level/', still_name, '/Correction');
    raw_path       = strcat(project_path, '\Correction');
    piv_path       = strcat(project_path, processing);
    
    %%% Save paths
    % Path for MacBook
    % results_path = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/';
    
    % Path for MME Desktop
    results_path = 'G:/Other computers/Zein MacBook Pro/wind_wave_PIV/';
    
    mtlb_file    = strcat(results_path, 'data'   , '\', inpt_name, '_DATA.mat');
    mean_file    = strcat(results_path, 'means'  , '\', inpt_name, '_MEANS.mat');
    phase_file   = strcat(results_path, 'phase'  , '\', inpt_name, '_PHASE.mat');
    wave_file    = strcat(results_path, 'wave'   , '\', inpt_name, '_WAVE.mat');
    figure_file  = strcat(results_path, 'figures', '\', inpt_name);
    
    % Make specific folder for figures of an experiment
    if ~exist(figure_file, 'dir')
        mkdir(figure_file)
    end
    
    % Edge Detection Constants
    std_tol            = 3;
    background_removal = 1000;
    canny_lower        = 0.1;
    canny_upper        = 0.4;
    grad_tol           = 0.4;
    nan_dist           = 5;
    
    % Filtering Parameters
    top_bound_value   = 205;       % relative to Y centered at still water
    left_bound_value  = -121;      % relative to X centered at DaVis default
    right_bound_value = 115;       % relative to X centered at DaVis default
    wave_buffer       = 0.5;
    cmap_perc         = 0.05;
    
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
    
    % Specift phase tolerance as a function of wavelength
    percent_phase_tolerance = 0.10;
    phase_tolerance         =  percent_phase_tolerance * wave_length;
    
    % Video Settings
    num_frames = 50;
    fps        = 1;
    
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
    
    constants.left_bound_value        = left_bound_value;
    constants.right_bound_value       = right_bound_value;
    constants.top_bound_value         = top_bound_value;
    constants.percent_phase_tolerance = percent_phase_tolerance;
    constants.phase_tolerance         = phase_tolerance; 
    constants.phase_offset            = phase_offset;
    
    constants.wave_type          = wave_type;
    constants.wave_amplitude     = wave_amplitude;
    constants.wave_length        = wave_length;
    % constants.vertical_offset    = waterlevel(still_path, still_name, constants);
    constants.vertical_offset    = -100;

    constants.wave_buffer        = wave_buffer;
    
    constants.PIV_height_correction = 5.3944;

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % WAVE DETECTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % if exist(wave_file, 'file')
    %      fprintf('* Loading WAVES from File\n')
    %      waves = load(wave_file);
    %      waves = waves.output;
    % else
    %      waves = wavedetection(raw_path, constants, wave_file);
    % end
    
    waves = wavedetection(raw_path, constants, wave_file);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DAVIS TO MATLAB
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % if exist(mtlb_file, 'file')
    %     fprintf('* Loading DATA from File\n')
    %     data = load(mtlb_file);
    %     data = data.output;
    % else
    %     data = vector2matlab2D(waves, piv_path, constants, mtlb_file);
    % end
    
    data = vector2matlab2D(waves, piv_path, constants, mtlb_file);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % PHASE AVERAGE (only binning)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if constants.wave_type == '0'
        fprintf('* Skipping <phaseaverage>\n\n')
        phase_average = 0;
    else
        % if exist(phase_file, 'file')
        %     fprintf('* Loading PHASES from File\n')
        %     phase_average = load(phase_file);
        %     phase_average = phase_average.output;
        % else
        %     phase_average = phaseaverage(raw_path, constants, wave_file, phase_file);
        % end

        phase_average = phaseaverage(raw_path, constants, wave_file, phase_file);
    end
    
    %%
    
    % uq_Y = unique(data.Y);
    % for i = 1:1
    % 
    %     prof = waves.wave_profiles(i,:);
    %     prof = imresize(prof, [1, length(unique(data.X))]);
    % 
    %     figure(i)
    %     hold on
    %     contourf(data.X, data.Y, data.U(:,:,i), 100, 'Linestyle', 'none')
    %     plot(unique(data.X), prof, 'Color', 'red', 'LineWidth', 3)
    %     yline(0, 'linestyle', '--')
    % 
    %     for j = 1:length(uq_Y)
    %         yline(uq_Y(j))
    %     end
    % 
    %     hold off
    %     axis equal
    %     % ylim([-20, 20])
    %     % pause(2)
    %     % close all
    % end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MATLAB DATA TO ENSEMBLE/PHASE MEANS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % if exist(mean_file, 'file')
    %      fprintf('* Loading MEANS from File\n')
    %      means = load(mean_file); 
    %      means = means.output;
    % else
    %      means = data2means2D(mean_file, data, waves, phase_average, constants, waves.D);
    % end
    
    means = data2means2D(mean_file, data, waves, phase_average, constants, waves.D);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Dimensions
    X     = means.X;
    Y     = means.Y;
    x     = unique(X);
    y     = flipud(unique(Y));
    range = abs(left_bound_value) + abs(right_bound_value);
    
    
    %%% Wave Profiles
    wave_colors      = colormap(jet(means.D));
    max_wave_profile = max(waves.wave_profiles, [], 1);
    high_res_x       = imresize(x, [length(waves.wave_profiles(1,:)), 1]);
    
    ax = tiledlayout(1,1, 'Padding', 'compact');
    nexttile
    hold on
    for i = 1:5:means.D
        p = plot(high_res_x, waves.wave_profiles(i,:));
        p.Color(1:3) = wave_colors(i, :);
        p.Color(4)   = 0.1;
    end
    plot(high_res_x, max_wave_profile, 'Color', 'black') 
    hold off
    axis equal
    xlim([0, range])
    
    if wave_type ~= '0'
        ylim([-2 * wave_amplitude, 2 * wave_amplitude])
    else
        ylim([-7, 7])
    end
    
    xlabel('x [mm]')
    ylabel('z [mm]')
    title(strcat(inpt_name, {': '}, 'All Wave Profiles'), 'Interpreter', 'none')
    
    %Save figures
    figure_file_name = strcat('/', inpt_name, '_wave_profiles.png');
    exportgraphics(ax, strcat(figure_file, figure_file_name), 'Resolution', 300);
    close all
    
    %% Ensemble Averages
    
    components       = {'u', 'v', 'uu', 'vv', 'uv'};
    u_inf            = mean(means.ensemble.u(:, 1:3), 'all', 'omitnan');
    max_wave_profile = means.ensemble.max_wave_profile;
    
    for i = 1:length(components)
        component = components{i};
        dat = means.ensemble.(component);
    
        if length(component) == 1
            norm = u_inf;
            cbar_label = strcat('$', component, ' / u_{\infty}$');
        else
            norm = u_inf^2;
            cbar_label = strcat('$', component, ' / {u_{\infty}}^2$');
        end
    
        % if length(component) == 1
        %     norm = u_inf;
        %     cbar_label = strcat('$', component, '$');
        % else
        %     norm = u_inf^2;
        %     cbar_label = strcat('$', component, '$');
        % end
    
        % Plot
        ax = figure('Name', strcat('Ensemble ', component));
        hold on
        contourf(X, Y, dat/norm, 500, 'LineStyle', 'none')
        % contourf(X, Y, dat, 500, 'LineStyle', 'none')
        yline(0, 'Color', 'black', 'LineWidth', 3)
        plot(x, max_wave_profile, 'Color', 'red', 'LineWidth', 3)
        hold off
    
        axis equal
        colormap parula
        C = colorbar();
        C.Label.String = cbar_label;
        C.Label.Interpreter = 'latex';
    
        xlim([0, range])
        ylim([-20, top_bound_value])
        xlabel('x [mm]')
        ylabel('z [mm]')
        title(strcat(inpt_name, {' '}, component, ': Ensemble Average'), 'Interpreter', 'none')
    
        %%% Save figures
        figure_file_name = strcat('/', inpt_name, '_', component, '_ensemble_average.png');
        exportgraphics(ax, strcat(figure_file, figure_file_name), 'Resolution', 300);
    end
    close all


    %% Phase Averages
    if wave_type ~= '0'
        u_inf      = mean(means.ensemble.u(:, 1:3), 'all', 'omitnan');
        components = {'u', 'v', 'uu', 'vv', 'uv'};
        
        for j = 1:4
            for i = 1:length(components)
                component = components{i};
                dat = means.phase(j).(component);
        
        
                if length(component) == 1
                    norm = u_inf;
                    cbar_label = strcat('$', component, ' / u_{\infty}$');
                else
                    norm = u_inf^2;
                    cbar_label = strcat('$', component, ' / {u_{\infty}}^2$');
                end
        
                % if length(component) == 1
                %     norm = u_inf;
                %     cbar_label = strcat('$', component, '$');
                % else
                %     norm = u_inf^2;
                %     cbar_label = strcat('$', component, '$');
                % end
        
                % Plot
                ax = figure('Name', strcat('Phase ', num2str(j), component));
        
                hold on
                contourf(X, Y, dat/norm, 500, 'LineStyle', 'none')
                % contourf(X, Y, dat, 500, 'LineStyle', 'none')
                plot(x, means.phase(j).reference_wave, 'Color', 'black', 'LineWidth', 2)
                plot(x, means.phase(j).max_wave_profile, 'Color', 'red', 'LineWidth', 2)
                % yline(0, 'LineStyle', '--', 'Color', 'black')
                hold off
        
                axis equal
                colormap parula
                C = colorbar();
                C.Label.String = cbar_label;
                C.Label.Interpreter = 'latex';
        
                xlim([0, range])
                ylim([-20, top_bound_value])
                xlabel('x [mm]')
                ylabel('z [mm]')
                title(strcat(inpt_name, {' '}, component, ': Phase', {' '}, num2str(j),  ' Average'), 'Interpreter', 'none')
        
                %%% Save figures
                figure_file_name = strcat('/', inpt_name, '_', component, '_phase_', num2str(j), '_average.png');
                exportgraphics(ax, strcat(figure_file, figure_file_name), 'Resolution', 300);
        
            end
            close all
        end
    end
end

clc;
close all
fprintf('Batch Done!\n')

