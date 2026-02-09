%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/WindWave_Functions/');
experiment_log  = readcell('Offshore_Inflow_Log.xlsx');
wave_parameters = readcell('Offshore_Waves.xlsx');

s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";

top_bound_value   = 205;       % relative to Y centered at still water
left_bound_value  = -121;      % relative to X centered at DaVis default
right_bound_value = 115;       % relative to X centered at DaVis default
range = abs(left_bound_value) + abs(right_bound_value);

%% Import Data

% Cases with waves
WTs = {'4', '6', '8'};
WVs = {'A', 'B', 'C', 'D'};

BL_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/boundary_layer';
LL_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/log_law/';
MN_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means/';
figure_path = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures';

% Load cases with waves
for s = 1:length(WTs)
    speed = WTs{s};
    for i = 1:length(WVs)
        % Case
        caze = strcat('WT', num2str(speed), '_WV', WVs{i}, '_AG0');
        fprintf('Loading Case: %s...\n', caze)
    
        % Paths
        BL_path = fullfile(BL_folder, strcat(caze, '_BL.mat'));
        LL_path = fullfile(LL_folder, strcat(caze, '_LL.mat'));
        MN_path = fullfile(MN_folder, strcat(caze, '_MEANS.mat'));
        
        % Save LL to Structure
        LL_temp = load(LL_path);
        LL_temp = LL_temp.output;
        LL_Waves.(caze) = LL_temp;

        % Save BL to Structure
        BL_temp = load(BL_path);
        BL_temp = BL_temp.output;
        BL_Waves.(caze) = BL_temp;
       
        % Save MN to Structure
        means_temp = load(MN_path);
        means_temp = means_temp.output;
        MN_Waves.(caze) = means_temp;
    end
    fprintf('\n')
end

% Load cases without waves
for s = 1:length(WTs)
    % Case
    no_wave_case = strcat('WT', num2str(WTs{s}), '_WV0_AGP');
    fprintf('Loading Case: %s...\n', no_wave_case)
    
    % Paths
    BL_path = fullfile(BL_folder, strcat(no_wave_case, '_BL.mat'));
    LL_path = fullfile(LL_folder, strcat(no_wave_case, '_LL.mat'));
    MN_path = fullfile(MN_folder, strcat(no_wave_case, '_MEANS.mat'));
    
    % Save LL to Structure
    LL_temp = load(LL_path);
    LL_temp = LL_temp.output;
    LL_No_Waves.(no_wave_case) = LL_temp;

    % Save BL to Structure
    BL_temp = load(BL_path);
    BL_temp = BL_temp.output;
    BL_No_Waves.(no_wave_case) = BL_temp;
    
    % Save MN to Structure
    means_temp = load(MN_path);
    means_temp = means_temp.output;
    MN_No_Waves.(no_wave_case) = means_temp;
end

% Get Case Names
wave_cases    = fields(MN_Waves);
no_wave_cases = fields(MN_No_Waves);

%% Get AVG Freestream from all cases
clc;
for s = 1:length(WTs)
    
    c = 1;
    all_u_infs = zeros(1,5);
    WT = WTs{s};
    % Waves
    for w = 1:length(WVs)
        WV = WVs{w};
        caze = strcat('WT', WT, '_WV', WV, '_AG0');
        u_inf = BL_Waves.(caze).u_inf;
        all_u_infs(c) = u_inf;
        c = c + 1;
    end

    % No Waves
    caze = strcat('WT', WT, '_WV0_AGP');
    u_inf = BL_No_Waves.(caze).u_inf;
    all_u_infs(5) = u_inf;
    avg_u_inf = mean(all_u_infs, 'all');
    fprintf("WT %s = %2.2f m/s\n", WT, avg_u_inf);
    u_infs.(strcat("WT", WT)) = avg_u_inf;
end

%% Cf Phase Subplots

WT = 4;
gauss_window = 1;
WV_colors = cool(length(WVs));

ax = figure('name', 'Cf Subplots');
t = tiledlayout(2,2);
sgtitle(strcat('WT', num2str(WT)))
for i = 1:4
    h(i) = nexttile();
    for w = 1:length(WVs)
        caze = strcat('WT', num2str(WT), '_WV', WVs{w}, '_AG0');

        if contains(caze, 'WVA') == 1
            color = WV_colors(1,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
            label = 'WVA';
        elseif contains(caze, 'WVB') == 1
            color = WV_colors(2,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
            label = 'WVB';
        elseif contains(caze, 'WVC') == 1
            color = WV_colors(3,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
            label = 'WVC';
        elseif contains(caze, 'WVD') == 1
            color = WV_colors(4,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
            label = 'WVD';
        end
    
        % Data
        u_inf = BL_Waves.(caze).u_inf;
        uv    = MN_Waves.(caze).phase(i).uv;
        x     = BL_Waves.(caze).x;
    
        % Compute u*
        u_star_profile = sqrt(max(-1 * uv, [], 1, 'omitnan'));
        u_star_profile = smoothdata(u_star_profile, 'gaussian', gauss_window);

        % Compute Cf
        cf = 2 * (u_star_profile / u_inf).^2;

        % Waves
        centered_normalized_x = (x - (range/2)) / wave_length;
        reference_wave = MN_Waves.(caze).phase(i).reference_wave / wave_amplitude;
        
        % Scale Wave
        wave_scale = 0.005;
        wave_offset = 0.008;
        reference_wave = reference_wave * wave_scale;
        reference_wave = reference_wave + wave_offset;
         
        % Plot
        hold on
        a = plot(centered_normalized_x, reference_wave, 'color', 'black', 'linewidth', 2, 'HandleVisibility', 'off');
        b(w) = plot(centered_normalized_x, cf, 'color', color, 'linewidth', 2);
        a.Color(4) = 1;
        uistack(a, 'bottom')
        xline(-0.5)
        xline(-0.25, 'linestyle', '--')
        xline(0)
        xline(0.25, 'linestyle', '--')
        xline(0.5)
        % Plot Vertical Line At Max Cf
        % [max_cf, max_cf_idx] = max(cf);
        % xline(centered_normalized_x(max_cf_idx), 'linestyle', '--', 'color', color)
        hold off
        tit = strcat('Phase', {' '}, num2str(i));
        title(tit{1})
        xlabel('$x / \lambda$', 'Interpreter', 'latex')
        ylabel('$C_f$', 'Interpreter', 'latex')
        
    end
end

leg = legend(b, 'WVA', 'WVB', 'WVC', 'WVD', 'orientation', 'horizontal');
leg.Layout.Tile = 'north';
linkaxes(h, 'xy')


%% Cf vs Wave Steepness

phase = 3;
marker_size = 75;
gauss_window = 10;

wind_speed_colors = turbo(3);
ax = figure('name', 'Cf vs Steepness');
hold on
for i = 1:length(WTs)
    
    steepnesses = zeros(1,4);
    max_cfs = zeros(1,4);
    color = wind_speed_colors(i,:);

    for w = 1:length(WVs)
        caze = strcat('WT', num2str(WTs{i}), '_WV', WVs{w}, '_AG0');

        if contains(caze, 'WVA') == 1
            
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
            
        elseif contains(caze, 'WVB') == 1
            
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
            
        elseif contains(caze, 'WVC') == 1
            
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
            
        elseif contains(caze, 'WVD') == 1
            
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
            
        end

        
        % Steepness
        ak = 2 * pi * (wave_amplitude/wave_length);
        steepnesses(1,w) = ak;
    
        % Data
        u_inf = BL_Waves.(caze).u_inf;
        uv    = MN_Waves.(caze).phase(phase).uv;
        x     = BL_Waves.(caze).x;
        centered_normalized_x = ((x - (range/2)) / wave_length) - ((i - 1)/4);
        reference_wave = MN_Waves.(caze).phase(i).reference_wave / wave_amplitude;
    
        % Compute u*
        u_star_profile = sqrt(max(-1 * uv, [], 1, 'omitnan'));
        u_star_profile = smoothdata(u_star_profile, 'gaussian', gauss_window);
    
        % Compute Cf
        cf = 2 * (u_star_profile / u_inf).^2;
        cf = smoothdata(cf, 'gaussian', gauss_window);
       
        cf(((x - (range/2)) / wave_length) < -0.5) = nan;
        cf(((x - (range/2)) / wave_length) > 0.5) = nan;

        max_cf = max(cf, [], 'omitnan');
        max_cfs(1,w) = max_cf;

        scatter(ak, max_cf, marker_size, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', color, 'HandleVisibility', 'off')
    end

    [B,I] = sort(steepnesses);
    label = sprintf('u_{\\infty} = %1.2f m/s', u_infs.(strcat("WT", WTs{i})));
    plot(B, max_cfs(I), linewidth = 2, color = color, displayname = label)

end
hold off
grid on
xlim([0.17, 0.31])
xlabel('$ak$', 'interpreter', 'Latex', 'FontSize', 16)
ylabel('$\max \left(C_f \right)$', 'Interpreter', 'Latex', 'FontSize', 16)
legend('Location', 'NorthWest')

% exportgraphics(ax, fullfile(figure_path, 'Cf_Steepness.png'), "Resolution", 300);


%% Cf vs Steepness All Phases


marker_size = 75;
gauss_window = 10;

wind_speed_colors = turbo(3);
ax = figure('name', 'Cf vs Steepness', 'position', [200,300,1000,500]);
t = tiledlayout(1,4, 'padding', 'tight');

clear h
for phase = 1:4
    h(phase) = nexttile();
    hold on
    for i = 1:length(WTs)
        
        steepnesses = zeros(1,4);
        max_cfs = zeros(1,4);
        color = wind_speed_colors(i,:);
    
        for w = 1:length(WVs)
            caze = strcat('WT', num2str(WTs{i}), '_WV', WVs{w}, '_AG0');
    
            if contains(caze, 'WVA') == 1
                
                wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
                wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
                
            elseif contains(caze, 'WVB') == 1
                
                wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
                wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
                
            elseif contains(caze, 'WVC') == 1
                
                wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
                wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
                
            elseif contains(caze, 'WVD') == 1
                
                wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
                wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
                
            end
    
            
            % Steepness
            ak = 2 * pi * (wave_amplitude/wave_length);
            steepnesses(1,w) = ak;
        
            % Data
            u_inf = BL_Waves.(caze).u_inf;
            uv    = MN_Waves.(caze).phase(phase).uv;
            x     = BL_Waves.(caze).x;
            centered_normalized_x = ((x - (range/2)) / wave_length) - ((i - 1)/4);
            reference_wave = MN_Waves.(caze).phase(i).reference_wave / wave_amplitude;
        
            % Compute u*
            u_star_profile = sqrt(max(-1 * uv, [], 1, 'omitnan'));
            u_star_profile = smoothdata(u_star_profile, 'gaussian', gauss_window);
        
            % Compute Cf
            cf = 2 * (u_star_profile / u_inf).^2;
            cf(((x - (range/2)) / wave_length) < -0.5) = nan;
            cf(((x - (range/2)) / wave_length) > 0.5) = nan;
            max_cf = max(cf, [], 'omitnan');
            max_cfs(1,w) = max_cf;
    
            scatter(ak, max_cf, marker_size, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', color, 'HandleVisibility', 'off')
        end
    
        [B,I] = sort(steepnesses);
        label = sprintf('u_{\\infty} = %1.2f m/s', u_infs.(strcat("WT", WTs{i})));
        plot(B, max_cfs(I), linewidth = 2, color = color, displayname = label)
    
    end
    hold off
    grid on
    xlim([0.17, 0.31])
    ylim([0.007, 0.03])
    xlabel('$ak$', 'interpreter', 'Latex', 'FontSize', 16)
    ylabel('$\max \left(C_f \right)$', 'Interpreter', 'Latex', 'FontSize', 16)
    legend('Location', 'NorthWest')
    title(sprintf("Phase %1.0f", phase))
end

linkaxes(h, 'xy')

%% Stitched Cf Across All Waves

WT = 4;
gauss_window = 1;
WV_colors = cool(length(WVs));

ax = figure('name', 'Cf stitched');

hold on
for i = 1:4
    for w = 1:length(WVs)
        caze = strcat('WT', num2str(WT), '_WV', WVs{w}, '_AG0');

        if contains(caze, 'WVA') == 1
            color = WV_colors(1,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
            label = 'WVA';
        elseif contains(caze, 'WVB') == 1
            color = WV_colors(2,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
            label = 'WVB';
        elseif contains(caze, 'WVC') == 1
            color = WV_colors(3,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
            label = 'WVC';
        elseif contains(caze, 'WVD') == 1
            color = WV_colors(4,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
            label = 'WVD';
        end
    
        % Data
        u_inf = BL_Waves.(caze).u_inf;
        uv    = MN_Waves.(caze).phase(i).uv;
        x     = BL_Waves.(caze).x;
    
        % Compute u*
        u_star_profile = sqrt(max(-1 * uv, [], 1, 'omitnan'));
        u_star_profile = smoothdata(u_star_profile, 'gaussian', gauss_window);

        % Compute Cf
        cf = 2 * (u_star_profile / u_inf).^2;
        

        % Waves
        centered_normalized_x = ((x - (range/2)) / wave_length) - ((i - 1)/4);
        reference_wave = MN_Waves.(caze).phase(i).reference_wave / wave_amplitude;
        
        
        % Scale Wave
        wave_scale = 0.005;
        wave_offset = 0.008;
        reference_wave = reference_wave * wave_scale;
        reference_wave = reference_wave + wave_offset;
        % 
        % Plot
        a = plot(centered_normalized_x, reference_wave, 'color', 'black', 'linewidth', 2, 'HandleVisibility', 'off');
        b(w) = plot(centered_normalized_x, cf, 'color', color, 'linewidth', 2);
        a.Color(4) = 1;
        uistack(a, 'bottom')

        % Plot Vertical Line At Max Cf
        % [max_cf, max_cf_idx] = max(cf);
        % xline(centered_normalized_x(max_cf_idx), 'linestyle', '--', 'color', color)
        
    end
end
xline(-1, 'color', 'black', 'linestyle', '--', 'linewidth', 2)
xline(-0.5, 'color', 'black', 'linestyle', '--', 'linewidth', 2)
xline(0, 'color', 'black', 'linestyle', '--', 'linewidth', 2)
hold off

legend(b, 'WVA', 'WVB', 'WVC', 'WVD');
xlabel('$x / \lambda$', 'Interpreter', 'latex')
ylabel('$C_f$', 'Interpreter', 'latex')
% xlim([-1.85 1.15])
xlim([-1.1,0.1])
ylim([0.0020 0.0220])
title(strcat('WT', num2str(WT), ' All Waves: Stitched Phases'));

%% Re Phase Subplots

WT = 4;
gauss_window = 1;
WV_colors = cool(length(WVs));
kinematic_viscosity = 1.48E-5;

ax = figure('name', 'Re subplot');
t = tiledlayout(2,2);
sgtitle(strcat('WT', num2str(WT)))
for i = 1:4
    disp(i)
    h(i) = nexttile();
    for w = 1:length(WVs)
        caze = strcat('WT', num2str(WT), '_WV', WVs{w}, '_AG0');

        if contains(caze, 'WVA') == 1
            color = WV_colors(1,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
            label = 'WVA';
        elseif contains(caze, 'WVB') == 1
            color = WV_colors(2,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
            label = 'WVB';
        elseif contains(caze, 'WVC') == 1
            color = WV_colors(3,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
            label = 'WVC';
        elseif contains(caze, 'WVD') == 1
            color = WV_colors(4,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
            label = 'WVD';
        end
    
        % Data
        u_inf = BL_Waves.(caze).u_inf;
        uv    = MN_Waves.(caze).phase(i).uv;
        theta = BL_Waves.(caze).phase(i).momentum * 1E-3;
        x     = BL_Waves.(caze).x;
    
        % Compute Re
        Re = (u_inf * theta) / kinematic_viscosity;
        Re = smoothdata(Re, 'gaussian', gauss_window);

        % Waves
        centered_normalized_x = (x - (range/2)) / wave_length;
        reference_wave = MN_Waves.(caze).phase(i).reference_wave / wave_amplitude;
        
        % Scale Wave
        wave_scale = 800;
        wave_offset = 1900;
        reference_wave = reference_wave * wave_scale;
        reference_wave = reference_wave + wave_offset;
         
        % Plot
        hold on
        a = plot(centered_normalized_x, reference_wave, 'color', 'black', 'linewidth', 2, 'HandleVisibility', 'off');
        b(w) = plot(centered_normalized_x, Re, 'color', color, 'linewidth', 2);
        a.Color(4) = 1;
        uistack(a, 'bottom')
        % Plot Vertical Line At Max Re
        [max_Re, max_Re_idx] = max(Re);
        xline(centered_normalized_x(max_Re_idx), 'linestyle', '--', 'color', color)
        hold off
        tit = strcat('Phase', {' '}, num2str(i));
        title(tit{1})
        xlabel('$x / \lambda$', 'Interpreter', 'latex')
        ylabel('$Re_{\theta}$', 'Interpreter', 'latex')
        
    end
end

leg = legend(b, 'WVA', 'WVB', 'WVC', 'WVD', 'orientation', 'horizontal');
leg.Layout.Tile = 'north';
linkaxes(h, 'xy')


%% Stitched Re_Theta Across All Waves

WT = 4;
gauss_window = 1;
WV_colors = cool(length(WVs));
kinematic_viscosity = 1.48E-5;

ax = figure('name', 'Re stitched');

hold on
for i = 2:2:4
    for w = 1:length(WVs)
        caze = strcat('WT', num2str(WT), '_WV', WVs{w}, '_AG0');

        if contains(caze, 'WVA') == 1
            color = WV_colors(1,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
            label = 'WVA';
        elseif contains(caze, 'WVB') == 1
            color = WV_colors(2,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
            label = 'WVB';
        elseif contains(caze, 'WVC') == 1
            color = WV_colors(3,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
            label = 'WVC';
        elseif contains(caze, 'WVD') == 1
            color = WV_colors(4,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
            label = 'WVD';
        end
    
        % Data
        u_inf = BL_Waves.(caze).u_inf;
        uv    = MN_Waves.(caze).phase(i).uv;
        theta = BL_Waves.(caze).phase(i).momentum * 1E-3;
        x     = BL_Waves.(caze).x;
    
        % Compute u*
        u_star_profile = sqrt(max(-1 * uv, [], 1, 'omitnan'));
        u_star_profile = smoothdata(u_star_profile, 'gaussian', gauss_window);

        % Compute Re
        Re = (u_inf * theta) / kinematic_viscosity;
        Re = smoothdata(Re, 'gaussian', gauss_window);

        % Waves
        centered_normalized_x = ((x - (range/2)) / wave_length) - ((i - 1)/4);
        reference_wave = MN_Waves.(caze).phase(i).reference_wave / wave_amplitude;

        Re(((x - (range/2)) / wave_length) < -0.25) = nan;
        Re(((x - (range/2)) / wave_length) > 0.25) = nan;

        
        % Scale Wave
        wave_scale = 800;
        wave_offset = 1900;
        reference_wave = reference_wave * wave_scale;
        reference_wave = reference_wave + wave_offset;
       

        % Plot
        a = plot(centered_normalized_x, reference_wave, 'color', 'black', 'linewidth', 2, 'HandleVisibility', 'off');
        b(w) = plot(centered_normalized_x, Re, 'color', color, 'linewidth', 2);
        a.Color(4) = 1;
        uistack(a, 'bottom')
        % Plot Vertical Line At Max Cf
        % [max_Re, max_Re_idx] = max(Re);
        % xline(centered_normalized_x(max_Re_idx), 'linestyle', '--', 'color', color)
        
    end
end
xline(-1, 'color', 'black', 'linestyle', '--', 'linewidth', 2)
xline(-0.5, 'color', 'black', 'linestyle', '--', 'linewidth', 2)
xline(0, 'color', 'black', 'linestyle', '--', 'linewidth', 2)
hold off

legend(b, 'WVA', 'WVB', 'WVC', 'WVD');
xlabel('$x / \lambda$', 'Interpreter', 'latex')
ylabel('$Re_{\theta}$', 'Interpreter', 'latex')
% xlim([-1.85 1.15])
xlim([-1.1 0.1])
title(strcat('WT', num2str(WT), ' All Waves: Stitched Phases'));


%% CHECK CROPPING

% Crop should be 0.25
crop_bound   = 0.25;
gauss_window = 5;
line_width   = 2;
kinematic_viscosity = 1.48E-5;

ax = figure('name', 'Cf vs Re Subplot', 'Position', [100,200,1200,400]);
t = tiledlayout(1,4,'TileSpacing','tight','Padding','compact');
sgtitle('$Cf$ vs $Re_{\theta}$ All Wind Speeds and Waves', 'interpreter', 'latex')

clc; close all; clear r h
for w = 1:4
    h(w) = nexttile();
    WV = WVs{w};
    disp(WV)

    hold on
    for s = 1:length(WTs)
        WT = WTs{s};

        for phase = 2:2:4
            caze = strcat('WT', WT, '_WV', WV, '_AG0');
        
            % Wave Parameters
            if contains(caze, 'WVA') == 1
                color = WV_colors(1,:);
                wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
                wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
                label = 'WVA';
            elseif contains(caze, 'WVB') == 1
                color = WV_colors(2,:);
                wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
                wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
                label = 'WVB';
            elseif contains(caze, 'WVC') == 1
                color = WV_colors(3,:);
                wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
                wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
                label = 'WVC';
            elseif contains(caze, 'WVD') == 1
                color = WV_colors(4,:);
                wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
                wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
                label = 'WVD';
            end

            % Data
            u_inf = BL_Waves.(caze).u_inf;
            uv    = MN_Waves.(caze).phase(phase).uv;
            theta = BL_Waves.(caze).phase(phase).momentum * 1E-3;
            x     = BL_Waves.(caze).x;
        
            % Compute u*
            u_star_profile = sqrt(max(-1 * uv, [], 1, 'omitnan'));
        
            % Compute Cf
            cf = 2 * (u_star_profile / u_inf).^2;
        
            % Compute Re
            Re_theta = (u_inf * theta) / kinematic_viscosity;
        
            % Smooth
            cf = smoothdata(cf, 'gaussian', gauss_window);
            Re_theta = smoothdata(Re_theta, 'gaussian', gauss_window);
        
            % Crop
            centered_normalized_x = ((x - (range/2)) / wave_length);
            cf(centered_normalized_x < -crop_bound) = nan;
            cf(centered_normalized_x > crop_bound) = nan;
            Re_theta(centered_normalized_x < -crop_bound) = nan;
            Re_theta(centered_normalized_x > crop_bound) = nan;
        
            if phase == 2
                temp = Re_theta(~isnan(Re_theta));
                Phase2_Re_theta_end_value = temp(end);

                temp = cf(~isnan(cf));
                Phase2_cf_end_value = temp(end);
                
            elseif phase == 4
                temp = Re_theta(~isnan(Re_theta));
                Phase4_Re_theta_start_value = temp(1);

                temp = cf(~isnan(cf));
                Phase4_cf_value = temp(1);

                Re_theta_offset = Phase2_Re_theta_end_value - Phase4_Re_theta_start_value;
                cf_offset = Phase2_cf_end_value - Phase4_cf_value;

                Re_theta = Re_theta + Re_theta_offset;
                cf = cf + cf_offset;
            end

            if phase == 2
                color = 'red';
            else
                color = 'black';
            end

            % Mask
            mask      = ~isnan(cf);
            cf_no_nan = cf(mask);
            Re_no_nan = Re_theta(mask);
            
            if phase == 2
                centered_normalized_x = centered_normalized_x - 0.25;
            elseif phase == 4
                centered_normalized_x = centered_normalized_x + 0.25;
            end
    
            % Plot
            plot(Re_theta, cf, 'color', color)
            % plot(centered_normalized_x, Re_theta, 'color', color)
            % xline(-crop_bound)
            % xline(crop_bound)

            title(strcat('WV', WV))
            xlabel('$Re_{\theta}$', 'Interpreter', 'latex')
            if w == 1
                ylabel('$C_f$', 'Interpreter', 'latex')
            end


        end
        
    end
    hold off
end

linkaxes(h, 'xy')


%% Cf vs Re_theta Phases 2 + 4: Subplot

% Crop should be 0.25
crop_bound   = 0.25;
gauss_window = 5;
line_width   = 2;
kinematic_viscosity = 1.48E-5;
FS = 16;

ax = figure('name', 'Cf vs Re Subplot', 'Position', [100,200,1200,400]);
t = tiledlayout(1,4,'TileSpacing','tight','Padding','compact');
% sgtitle('$Cf$ vs $Re_{\theta}$ All Wind Speeds and Waves', 'interpreter', 'latex')

clc;
h = 0;
r = 0;
for w = 1:4
    h(w) = nexttile();
    WV = WVs{w};
    disp(WV)

    for s = 1:length(WTs)
        WT = WTs{s};

        for phase = 2:2:4
            caze = strcat('WT', WT, '_WV', WV, '_AG0');
        
            % Wave Parameters
            if contains(caze, 'WVA') == 1
                color = WV_colors(1,:);
                wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
                wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
                label = 'WVA';
            elseif contains(caze, 'WVB') == 1
                color = WV_colors(2,:);
                wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
                wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
                label = 'WVB';
            elseif contains(caze, 'WVC') == 1
                color = WV_colors(3,:);
                wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
                wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
                label = 'WVC';
            elseif contains(caze, 'WVD') == 1
                color = WV_colors(4,:);
                wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
                wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
                label = 'WVD';
            end

            % Data
            u_inf = BL_Waves.(caze).u_inf;
            uv    = MN_Waves.(caze).phase(phase).uv;
            theta = BL_Waves.(caze).phase(phase).momentum * 1E-3;
            x     = BL_Waves.(caze).x;
        
            % Compute u*
            u_star_profile = sqrt(max(-1 * uv, [], 1, 'omitnan'));
        
            % Compute Cf
            cf = 2 * (u_star_profile / u_inf).^2;
        
            % Compute Re
            Re_theta = (u_inf * theta) / kinematic_viscosity;
        
            % Smooth
            cf = smoothdata(cf, 'gaussian', gauss_window);
            Re_theta = smoothdata(Re_theta, 'gaussian', gauss_window);
        
            % Crop
            centered_normalized_x = ((x - (range/2)) / wave_length);
            cf(centered_normalized_x < -crop_bound) = nan;
            cf(centered_normalized_x > crop_bound) = nan;
            Re_theta(centered_normalized_x < -crop_bound) = nan;
            Re_theta(centered_normalized_x > crop_bound) = nan;
        
            % Path Colors
            num_values = sum(~isnan(cf));

            % Choose Colormap Based on Wind Speed
            if WT == '4'
                temp_color = summer(2 * num_values);
            elseif WT == '6'
                temp_color = autumn(2 * num_values);
            elseif WT == '8'
                temp_color = winter(2 * num_values);
            end
        

            % Mask
            mask      = ~isnan(cf);
            cf_no_nan = cf(mask);
            Re_no_nan = Re_theta(mask);
    
            % Plot
            if phase == 2

                % Fit a Line to Bottom Part of Hystersis
                c   = polyfit(Re_no_nan, cf_no_nan, 1);
                fit = polyval(c, Re_no_nan);

                % Plot Line with Different Colors
                hold on
                for i = 1:length(cf_no_nan)-1
                    if i == 1
                        r(s) = plot(Re_no_nan(i:i+1), cf_no_nan(i:i+1), 'color', temp_color(i,:),'LineWidth', line_width);
                    else
                        plot(Re_no_nan(i:i+1), cf_no_nan(i:i+1), 'color', temp_color(i,:),'LineWidth', line_width);
                    end
                end
                % plot(Re_no_nan, fit, 'color', 'black', 'linestyle', '--', 'linewidth', 1)
                hold off

            elseif phase == 4
                % Plot Line with Different Colors
                hold on
                for i = 1:length(cf_no_nan)-1
                    plot(Re_no_nan(i:i+1), cf_no_nan(i:i+1), 'color', temp_color(i + num_values,:),'LineWidth', line_width)
                end
                hold off

            end

            title(strcat('WV', WV))
            xlabel('$Re_{\theta}$', 'Interpreter', 'latex', 'fontsize', FS)
            if w == 1
                ylabel('$C_f$', 'Interpreter', 'latex', 'fontsize', FS)
            end

        end
    end
end

linkaxes(h, 'xy')
leg = legend(r, 'WT4', 'WT6', 'WT8', 'orientation', 'horizontal');
leg.Layout.Tile = 'north';

% exportgraphics(ax, fullfile(figure_path, "Hysteresis.png"), resolution=300)


%% Cf / Re_theta vs x / wave length: single plot

% Crop should be 0.25
gauss_window = 10;
line_width   = 2;
kinematic_viscosity = 1.48E-5;
phase = 4;

clc;
clear r h
ax = figure('name', 'Cf vs Re Subplot');
hold on
for w = 2:4
    WV = WVs{w};
    disp(WV)

    for s = 1:length(WTs)
        WT = WTs{s};

        caze = strcat('WT', WT, '_WV', WV, '_AG0');
    
        % Wave Parameters
        if contains(caze, 'WVA') == 1
            color = WV_colors(1,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
            label = 'WVA';
        elseif contains(caze, 'WVB') == 1
            color = WV_colors(2,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
            label = 'WVB';
        elseif contains(caze, 'WVC') == 1
            color = WV_colors(3,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
            label = 'WVC';
        elseif contains(caze, 'WVD') == 1
            color = WV_colors(4,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
            label = 'WVD';
        end

        % Data
        u_inf = BL_Waves.(caze).u_inf;
        uv    = MN_Waves.(caze).phase(phase).uv;
        theta = BL_Waves.(caze).phase(phase).momentum * 1E-3;
        x     = BL_Waves.(caze).x;
        x     = (x - range/2)/wave_length;
    
        % Compute u*
        u_star_profile = sqrt(max(-1 * uv, [], 1, 'omitnan'));
    
        % Compute Cf
        cf = 2 * (u_star_profile / u_inf).^2;
    
        % Compute Re
        Re_theta = (u_inf * theta) / kinematic_viscosity;
    
        % Smooth
        cf = smoothdata(cf, 'gaussian', gauss_window);
        Re_theta = smoothdata(Re_theta, 'gaussian', gauss_window);
    

        % Linestyle Based on wind speed
        if contains(caze, 'WT4') == 1
            ls = '-';
        elseif contains(caze, 'WT6') == 1
            ls = '--';
        elseif contains(caze, 'WT8') == 1
            ls = ':';
        end
    
        % Plot
        plot(x, cf ./ Re_theta, 'LineWidth', line_width, 'color', color, 'linestyle', ls)
   
    end
end
hold off
ylabel('$Cf / Re_{\theta}$', 'interpreter', 'latex')
xlabel('$x / \lambda$', 'interpreter', 'latex')
xlim([-0.25, 0.25])


%% Wave-Specific Profiles: Shear Stress


% Crop should be 0.25
crop_bound   = 0.25;
gauss_window = 10;
line_width   = 2;
kinematic_viscosity = 1.48E-5;
FS = 16;

clear h
ax = figure('name', 'Profiles Along Wave', 'Position', [100,200,1200,400]);
t = tiledlayout(1,3,'TileSpacing','tight','Padding','compact');

% Profiles At Trough
h(1) = nexttile();
ylim([-15,200])
title("Trough")
xlabel("$-\overline{u'v'} / {u_{\infty}}^2$", 'Interpreter', 'Latex', 'fontsize', FS)
ylabel("y [mm]", 'Interpreter', 'latex', 'fontsize', FS)
hold on
for w = 1:4
    WV = WVs{w};
    for s = 1:length(WTs)
        WT = WTs{s};
        caze = strcat('WT', WT, '_WV', WV, '_AG0');

        % Wave Parameters
        if contains(caze, 'WVA') == 1
            color = WV_colors(1,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
            label = 'WVA';
        elseif contains(caze, 'WVB') == 1
            color = WV_colors(2,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
            label = 'WVB';
        elseif contains(caze, 'WVC') == 1
            color = WV_colors(3,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
            label = 'WVC';
        elseif contains(caze, 'WVD') == 1
            color = WV_colors(4,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
            label = 'WVD';
        end

        % Data
        u_inf = BL_Waves.(caze).u_inf;
        uv    = MN_Waves.(caze).phase(3).uv;
        theta = BL_Waves.(caze).phase(3).momentum * 1E-3;
        x     = BL_Waves.(caze).x;
        Y     = MN_Waves.(caze).Y;
        u_inf = u_infs.(strcat("WT", WT));

        uv_profile = -uv(:, range/2) / (u_inf^2);
        uv_profile = smoothdata(uv_profile, 'gaussian', gauss_window);

        plot(uv_profile, Y(:, range/2), 'color', color, 'linewidth', line_width)
    end
end
hold off


% Profiles At Peak
h(2) = nexttile();
ylim([-15,200])
title("Peak")
xlabel("$-\overline{u'v'} / {u_{\infty}}^2$", 'Interpreter', 'Latex', 'fontsize', FS)
hold on
for w = 1:4
    WV = WVs{w};
    for s = 1:length(WTs)
        WT = WTs{s};
        caze = strcat('WT', WT, '_WV', WV, '_AG0');

        % Wave Parameters
        if contains(caze, 'WVA') == 1
            color = WV_colors(1,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
            label = 'WVA';
        elseif contains(caze, 'WVB') == 1
            color = WV_colors(2,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
            label = 'WVB';
        elseif contains(caze, 'WVC') == 1
            color = WV_colors(3,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
            label = 'WVC';
        elseif contains(caze, 'WVD') == 1
            color = WV_colors(4,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
            label = 'WVD';
        end

        % Data
        u_inf = BL_Waves.(caze).u_inf;
        uv    = MN_Waves.(caze).phase(1).uv;
        theta = BL_Waves.(caze).phase(1).momentum * 1E-3;
        x     = BL_Waves.(caze).x;
        Y     = MN_Waves.(caze).Y;
        u_inf = u_infs.(strcat("WT", WT));

        uv_profile = -uv(:, range/2) / (u_inf^2);
        uv_profile = smoothdata(uv_profile, 'gaussian', gauss_window);

        plot(uv_profile, Y(:, range/2), 'color', color, 'linewidth', line_width)
    end
end
hold off


% Profiles At Max Cf
h(3) = nexttile();
ylim([-15,200])
title("Max Skin Friction")
xlabel("$-\overline{u'v'} / {u_{\infty}}^2$", 'Interpreter', 'Latex', 'fontsize', FS)
hold on
for w = 1:4
    WV = WVs{w};
    disp(WV)
    for s = 1:length(WTs)
        WT = WTs{s};
        caze = strcat('WT', WT, '_WV', WV, '_AG0');
    
        % Wave Parameters
        if contains(caze, 'WVA') == 1
            color = WV_colors(1,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
            label = 'WVA';
        elseif contains(caze, 'WVB') == 1
            color = WV_colors(2,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
            label = 'WVB';
        elseif contains(caze, 'WVC') == 1
            color = WV_colors(3,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
            label = 'WVC';
        elseif contains(caze, 'WVD') == 1
            color = WV_colors(4,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
            label = 'WVD';
        end

        % Data
        u_inf = BL_Waves.(caze).u_inf;
        uv    = MN_Waves.(caze).phase(phase).uv;
        theta = BL_Waves.(caze).phase(phase).momentum * 1E-3;
        x     = BL_Waves.(caze).x;
        Y     = MN_Waves.(caze).Y;
        u_inf = u_infs.(strcat("WT", WT));
    
        % Compute u*
        u_star_profile = sqrt(max(-1 * uv, [], 1, 'omitnan'));
    
        % Compute Cf
        cf = 2 * (u_star_profile / u_inf).^2;
   
        % Smooth
        cf = smoothdata(cf, 'gaussian', gauss_window);

        % Crop
        centered_normalized_x = ((x - (range/2)) / wave_length);
        cf(centered_normalized_x < -crop_bound) = nan;
        cf(centered_normalized_x > crop_bound) = nan;
        
        % Find Maximum Cf
        [Max_Cf, I] = max(cf, [], 'all');
        
        uv_profile = -uv(:, I) / (u_inf^2);
        uv_profile = smoothdata(uv_profile, 'gaussian', gauss_window);

        plot(uv_profile, Y(:, I), 'color', color, 'linewidth', line_width)

    end
end
hold off
linkaxes(h, 'xy')

% exportgraphics(ax, fullfile(figure_path, "Hysteresis.png"), resolution=300)



%% Wave-Specific Profiles: Log Law


% Crop should be 0.25
crop_bound   = 0.25;
gauss_window = 10;
line_width   = 2;
kinematic_viscosity = 1.48E-5;
FS = 16;

clear h
ax = figure('name', 'Profiles Along Wave', 'Position', [100,200,1200,400]);
t = tiledlayout(1,3,'TileSpacing','tight','Padding','compact');

% Profiles At Trough
h(1) = nexttile();
title("Trough")
ylabel('$u^+$', 'Interpreter', 'Latex')
xlabel('$y^+$', 'Interpreter', 'Latex')
hold on
for w = 1:4
    WV = WVs{w};
    for s = 1:length(WTs)
        WT = WTs{s};
        caze = strcat('WT', WT, '_WV', WV, '_AG0');

        % Wave Parameters
        if contains(caze, 'WVA') == 1
            color = WV_colors(1,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
            label = 'WVA';
        elseif contains(caze, 'WVB') == 1
            color = WV_colors(2,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
            label = 'WVB';
        elseif contains(caze, 'WVC') == 1
            color = WV_colors(3,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
            label = 'WVC';
        elseif contains(caze, 'WVD') == 1
            color = WV_colors(4,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
            label = 'WVD';
        end

        % Data
        u     = MN_Waves.(caze).phase(3).u;
        Y     = MN_Waves.(caze).Y * 1E-3;
        u_inf = u_infs.(strcat("WT", WT));

        u_star = LL_Waves.(caze).phase(3).u_star;
        u_profile = u(:, range/2);
        u_plus = u_profile / u_star;
        y_plus = (Y(:,range/2) * u_star) / kinematic_viscosity;

        plot(y_plus, u_plus, 'color', color, 'linewidth', line_width)
    end
end
hold off
set(gca,'xscale','log')


% Profiles At Peak
h(2) = nexttile();
title("Peak")
xlabel('$y^+$', 'Interpreter', 'Latex')
hold on
for w = 1:4
    WV = WVs{w};
    for s = 1:length(WTs)
        WT = WTs{s};
        caze = strcat('WT', WT, '_WV', WV, '_AG0');

        % Wave Parameters
        if contains(caze, 'WVA') == 1
            color = WV_colors(1,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
            label = 'WVA';
        elseif contains(caze, 'WVB') == 1
            color = WV_colors(2,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
            label = 'WVB';
        elseif contains(caze, 'WVC') == 1
            color = WV_colors(3,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
            label = 'WVC';
        elseif contains(caze, 'WVD') == 1
            color = WV_colors(4,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
            label = 'WVD';
        end

        % Data
        u     = MN_Waves.(caze).phase(1).u;
        theta = BL_Waves.(caze).phase(1).momentum * 1E-3;
        x     = BL_Waves.(caze).x;
        Y     = MN_Waves.(caze).Y * 1E-3;
        u_inf = u_infs.(strcat("WT", WT));

        u_star = LL_Waves.(caze).phase(1).u_star;
        u_profile = u(:, range/2);
        u_plus = u_profile / u_star;
        y_plus = (Y(:,range/2) * u_star) / kinematic_viscosity;

        plot(y_plus, u_plus, 'color', color, 'linewidth', line_width)
    end
end
hold off
set(gca,'xscale','log')


% Profiles At Max Cf
h(3) = nexttile();
title("Max Skin Friction")
xlabel("$-\overline{u'v'} / {u_{\infty}}^2$", 'Interpreter', 'Latex', 'fontsize', FS)
hold on
for w = 1:4
    WV = WVs{w};
    disp(WV)
    for s = 1:length(WTs)
        WT = WTs{s};
        caze = strcat('WT', WT, '_WV', WV, '_AG0');
    
        % Wave Parameters
        if contains(caze, 'WVA') == 1
            color = WV_colors(1,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'A') == 1), 3};
            label = 'WVA';
        elseif contains(caze, 'WVB') == 1
            color = WV_colors(2,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'B') == 1), 3};
            label = 'WVB';
        elseif contains(caze, 'WVC') == 1
            color = WV_colors(3,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'C') == 1), 3};
            label = 'WVC';
        elseif contains(caze, 'WVD') == 1
            color = WV_colors(4,:);
            wave_length = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 2};
            wave_amplitude = wave_parameters{find(strcmp(wave_parameters, 'D') == 1), 3};
            label = 'WVD';
        end

        % Data
        u     = MN_Waves.(caze).phase(phase).u;
        theta = BL_Waves.(caze).phase(phase).momentum * 1E-3;
        x     = BL_Waves.(caze).x;
        Y     = MN_Waves.(caze).Y * 1E-3;
        u_inf = u_infs.(strcat("WT", WT));
    
        % Compute u*
        u_star_profile = sqrt(max(-1 * uv, [], 1, 'omitnan'));
    
        % Compute Cf
        cf = 2 * (u_star_profile / u_inf).^2;
   
        % Smooth
        cf = smoothdata(cf, 'gaussian', gauss_window);

        % Crop
        centered_normalized_x = ((x - (range/2)) / wave_length);
        cf(centered_normalized_x < -crop_bound) = nan;
        cf(centered_normalized_x > crop_bound) = nan;
        
        % Find Maximum Cf
        [Max_Cf, I] = max(cf, [], 'all');
        
        u_star = u_star_profile(I);
        u_profile = u(:, I);
        u_plus = u_profile / u_star;
        y_plus = (Y(:,range/2) * u_star) / kinematic_viscosity;

        plot(y_plus, u_plus, 'color', color, 'linewidth', line_width)

    end
end
hold off
set(gca,'xscale','log')
linkaxes(h, 'xy')

% exportgraphics(ax, fullfile(figure_path, "Hysteresis.png"), resolution=300)










