%% 2D Offshore PIV Processing: Ondrej Fercak, Zein Sadek, 1/2023
% Converts DaVis vector data (.vc7) files into a Matlab file

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Functions/');
experiment_log  = readcell('Offshore_Inflow_Log.xlsx');
wave_parameters = readcell('Offshore_Waves.xlsx');

% Edge Detection Constants
std_tol            = 2;
background_removal = 1000;
canny_lower        = 0.1;
canny_upper        = 0.4;
grad_tol           = 0.4;
nan_dist           = 5;

% Filtering Parameters
left_bound_value   = -115;
right_bound_value  = 115;
phase_tolerance    = 5;
wave_buffer        = 1;
cmap_perc          = 0.05;

% Data paths
main_path     = '/Volumes/Zein_PIV_2/';
processing    = '/PIV_MP(2x24x24_50%ov)/PostProc';
project_names = {'Weitemyer_Alternating', 'Weitemyer_Parallel'};

for i = 1:length(project_names)
    project_path = strcat(main_path, project_names{i});
    
    cases = dir(strcat(main_path, project_names{i}));
    cases = cases([cases.isdir]);
    cases = {cases(4:end).name};
    
    results_path   = strcat(project_path, '_results/');
    if ~exist(results_path, 'dir')
        mkdir(results_path)
        mkdir(strcat(results_path, 'phase'))
        mkdir(strcat(results_path, 'data'))
        mkdir(strcat(results_path, 'means'))
        mkdir(strcat(results_path, 'figures'))
    end
    
    disp(project_names{i})
    
    
    for j = 1:length(cases)
        recording_name = cases{j};
        disp(recording_name)
       
        still_name     = strcat(experiment_log{find(strcmp(experiment_log, recording_name) == 1), 2}, '_waterlevel');
        inpt_name      = strcat(recording_name, '_2');
        
        % Image paths
        still_path     = strcat('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Offshore/Water_Level/', still_name, '/Correction');
        raw_path       = strcat(project_path, '/', recording_name, '/Correction');
        piv_path       = strcat(project_path, '/', recording_name, processing);
        
        % Save paths
        mtlb_file      = strcat(results_path, 'data'   , '/', inpt_name, '_DATA.mat');
        mean_file      = strcat(results_path, 'means'  , '/', inpt_name, '_MEANS.mat');
        phase_file     = strcat(results_path, 'phase'  , '/', inpt_name, '_PHASE.mat');
        figure_file    = strcat(results_path, 'figures', '/', recording_name);

        wave_type      = '0';
        wave_length    = 0;
        wave_amplitude = 0;
        phase_offset   = [0, wave_length/4, wave_length/2, 3*wave_length/4];
        
        % Constants Structure
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
        constants.vertical_offset    = waterlevel(still_path, piv_path, still_name, constants);
        constants.wave_buffer        = wave_buffer;
        
        %%% Phase Average
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
        
        %%% DaVis to Matlab
        if exist(mtlb_file, 'file')
            fprintf('* Loading DATA from File\n')
            data = load(mtlb_file);
            data = data.output;
        else
             data = vector2matlab2D(raw_path, piv_path, constants, mtlb_file);
        end

        %%% Means
        if exist(mean_file, 'file')
            fprintf('* Loading MEANS from File\n')
            means = load(mean_file); 
            means = means.output;
        else
             means = data2means2D(mean_file, data, phase_average);
        end

        %%% Figures
        if ~exist(figure_file, 'dir')
            mkdir(figure_file)
        end
    
        components = fields(means.ensemble);
        for k = 1:length(components)
            component = components{k};
            dat       = means.ensemble.(component);
            figure()
            contourf(means.X, means.Y - (constants.vertical_offset + wave_buffer), dat, 500, 'LineStyle', 'none')
            axis equal
            colorbar()
            colormap parula
            ylim([min(unique(means.Y)) - (constants.vertical_offset + wave_buffer), 205])
            xlim([left_bound_value, right_bound_value])
            yline(0, 'Color', 'red', 'LineWidth', 4)
            xlabel('x [mm]')
            ylabel('z [mm]')
            title(component)
            ax = gca;
            ax.CLim = [(1 + cmap_perc) * min(dat, [], 'all'), (1 - cmap_perc) * max(dat, [], 'all')];
            Tight = get(ax, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
            NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
            set(ax, 'Position', NewPos);
            figure_name = strcat(recording_name, '_ensemble_', component, '.png');
            saveas(ax, strcat(figure_file, '/', figure_name))
        end
        close all
    end
end


%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Functions/');


% Data paths
main_path     = '/Volumes/Zein_PIV_2/';
processing    = '/PIV_MP(2x24x24_50%ov)/PostProc';
project_names = {'Weitemyer_Alternating', 'Weitemyer_Parallel'};

for i = 1:length(project_names)
    project_path = strcat(main_path, project_names{i});
    results_path = strcat(project_path, '_results/means/');
    
    cases = dir(results_path);
    cases = {cases(3:end).name};
    
    figure(i)
    T = title(project_names{i});
    set(T, 'Interpreter', 'none')
    xlabel('Frequency [Hz]')
    ylabel('Turbulence Intensity')
    xlim([1, 11])
    ylim([0, 0.3])
    hold on
    for j = 1:length(cases)
        caze = cases{j};
        
        STD  = str2double(caze(strfind(caze, 'STD') + 3: strfind(caze, 'STD') + 4));
        freq = str2double(caze(strfind(caze, 'FQ') + 2: strfind(caze, 'FQ') + 3));
        
        data = load(strcat(results_path, caze));
        data = data.output;
        ti   = sqrt(data.ensemble.uu) ./ data.ensemble.u;
        ti   = mean(ti(1:20, :), 'all', 'omitnan');
        
        if STD == 10
            marker = 'o';
        elseif STD == 20
            marker = '^';
        elseif STD == 30
            marker = 'square';
        end
        
        scatter(freq, ti, 100, 'filled', marker, 'DisplayName', caze(1:end-10))
        
    end
    hold off
    legend('Interpreter', 'none')
end

%%



WT_speeds = {'4', '6', '8'};

% figure()
% hold on
for i = 1:length(WT_speeds)
    speed = WT_speeds{i};
    free_stream = load(strcat('/Volumes/No_Waves/Offshore_Inflow_No_Waves_Results/means/WT', speed, '_WV0_AG0_MEANS.mat'));
    free_stream = free_stream.output;
%     yline(mean(free_stream.ensemble.u(1:20, :), 'all', 'omitnan'))
    
    fprintf('WT %s Hz: %3.2f m/s\n', speed, mean(free_stream.ensemble.u(1:20, :), 'all', 'omitnan'))
end
% hold off


%%

% free_stream = load('/Volumes/No_Waves/Offshore_Inflow_No_Waves_Results/means/WT4_WV0_AG0_MEANS.mat');
% free_stream = free_stream.output;

% Data paths
main_path     = '/Volumes/Zein_PIV_2/';
processing    = '/PIV_MP(2x24x24_50%ov)/PostProc';
project_names = {'Weitemyer_Alternating', 'Weitemyer_Parallel'};


for i = 1:length(project_names)
    project_path = strcat(main_path, project_names{i});
    results_path = strcat(project_path, '_results/means/');
    
    cases = dir(results_path);
    cases = {cases(3:end).name};
    
    figure(i)
    yline(2.44, 'DisplayName', 'Passive Grid')
    T = title(project_names{i});
    set(T, 'Interpreter', 'none')
    xlabel('Frequency [Hz]')
    ylabel('Streamwise Velocity')
    xlim([1, 11])
    ylim([0, 3])
    hold on
    
    for j = 1:length(cases)
        caze = cases{j};
        
        STD  = str2double(caze(strfind(caze, 'STD') + 3: strfind(caze, 'STD') + 4));
        freq = str2double(caze(strfind(caze, 'FQ') + 2: strfind(caze, 'FQ') + 3));
        
        data = load(strcat(results_path, caze));
        data = data.output;
        u    = mean(data.ensemble.u(1:20, :), 'all', 'omitnan');
        
        if STD == 10
            marker = 'o';
        elseif STD == 20
            marker = '^';
        elseif STD == 30
            marker = 'square';
        end
        
        scatter(freq, u, 'filled', marker, 'DisplayName', caze(1:end-10))
    end
    hold off
    legend('Interpreter', 'none', 'Location', 'southeast')
    
end


%%

% Edge Detection Constants
std_tol            = 2;
background_removal = 1000;
canny_lower        = 0.1;
canny_upper        = 0.4;
grad_tol           = 0.4;
nan_dist           = 5;

% Filtering Parameters
left_bound_value   = -115;
right_bound_value  = 115;
phase_tolerance    = 5;
wave_buffer        = 1;
cmap_perc          = 0.05;

% Constants Structure
constants.background_removal = background_removal;
constants.std_tol            = std_tol;
constants.canny_lower        = canny_lower;
constants.canny_upper        = canny_upper;
constants.grad_tol           = grad_tol;
constants.nan_dist           = nan_dist;

constants.left_bound_value   = left_bound_value;
constants.right_bound_value  = right_bound_value;
constants.phase_tolerance    = phase_tolerance; 
constants.wave_buffer        = wave_buffer;

% Data paths
main_path     = '/Volumes/Zein_PIV_2/';
processing    = '/PIV_MP(2x24x24_50%ov)/PostProc';
project_names = {'Weitemyer_Alternating', 'Weitemyer_Parallel'};


for i = 1:length(project_names)
    project_path = strcat(main_path, project_names{i});
    results_path = strcat(project_path, '_results/means/');
    
    cases = dir(results_path);
    cases = {cases(3:end).name};
    
    
    for j = 1:length(cases)
        caze = cases{j};
        
        STD  = str2double(caze(strfind(caze, 'STD') + 3: strfind(caze, 'STD') + 4));
        freq = str2double(caze(strfind(caze, 'FQ') + 2: strfind(caze, 'FQ') + 3));
        
        data = load(strcat(results_path, caze));
        data = data.output;
        ti_u   = sqrt(data.ensemble.uu) ./ data.ensemble.u;
        ti_v   = sqrt(data.ensemble.vv) ./ data.ensemble.v;
        
        ti_v_mean = mean(ti_v, 'all', 'omitnan');
        ti_v_std  = std(ti_v, 0, 'all', 'omitnan');
        std_tol   = 2;
        
        
        ti_v(ti_v < (ti_v_mean - std_tol * ti_v_std)) = nan;
        ti_v(ti_v > (ti_v_mean + std_tol * ti_v_std)) = nan;
        
        ti_u(1:20, :) = nan;
        
        
        figure()
%         subplot(1,2,1)
        contourf(data.X, data.Y, ti_u, 500, 'LineStyle', 'none')
        axis equal
        colorbar()
        colormap parula
%         xlim([left_bound_value, right_bound_value])
        xlabel('x [mm]')
        ylabel('z [mm]')
        T = title(strcat(project_names{i}, {' '}, caze(1:end-10), {' '}, 'streamwise Ti'));
        set(T, 'Interpreter', 'none')
        
%         subplot(1,2,2)
%         contourf(data.X, data.Y, ti_v, 500, 'LineStyle', 'none')
%         axis equal
%         colorbar()
%         colormap parula
%         ylim([min(unique(data.Y)), 205])
%         xlim([left_bound_value, right_bound_value])
% %         yline(0, 'Color', 'red', 'LineWidth', 4)
%         xlabel('x [mm]')
%         ylabel('z [mm]')
%         T = title(strcat(project_names{i}, {' '}, caze(1:end-10), {' '}, 'vertical Ti'));
%         set(T, 'Interpreter', 'none')
%         
    end

end









