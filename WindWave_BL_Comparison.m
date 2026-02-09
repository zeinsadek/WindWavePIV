%% 2D Offshore PIV Processing: Zein Sadek, 1/2023

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Functions/');
experiment_log  = readcell('Offshore_Inflow_Log.xlsx');
wave_parameters = readcell('Offshore_Waves.xlsx');


top_bound_value   = 205;       % relative to Y centered at still water
left_bound_value  = -121;      % relative to X centered at DaVis default
right_bound_value = 115;       % relative to X centered at DaVis default
range = abs(left_bound_value) + abs(right_bound_value);

%% Import Data

% Cases with waves
WTs = {'4', '6', '8'};
WVs = {'A', 'B', 'C', 'D'};

BL_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/boundary_layer/';
MN_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means/';
path      = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/paper_figures';

% Load cases with waves
for s = 1:length(WTs)
    speed = WTs{s};
    for i = 1:length(WVs)
        % Case
        caze = strcat('WT', num2str(speed), '_WV', WVs{i}, '_AG0');
        fprintf('Loading Case: %s...\n', caze)
    
        % Paths
        BL_path = fullfile(BL_folder, strcat(caze, '_BL.mat'));
        MN_path = fullfile(MN_folder, strcat(caze, '_MEANS.mat'));
        
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
    MN_path = fullfile(MN_folder, strcat(no_wave_case, '_MEANS.mat'));
    
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

%% Castillo Correlation

WT_colors = turbo(length(WTs));
WV_colors = cool(length(wave_cases) + 1);
markers   = {'o', 'square', 'diamond', '^'};
alpha = 0.1;
size  = 10;

figure('Name', 'Castillo Correlation')
hold on
for i = 1:length(wave_cases)
    
    % marker = markers{i};
    color = WV_colors(i,:);   
    caze  = wave_cases{i};
    disp(caze)
    
    % Ensemble
    thickness    = BL_Waves.(caze).ensemble.thickness;
    displacement = BL_Waves.(caze).ensemble.displacement;
    momentum     = BL_Waves.(caze).ensemble.momentum;
    G = displacement ./ thickness;
    H = displacement ./ momentum;
    scatter(G, H, size, 'o', 'filled', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha)

    % Phase
    for j = 1:4
        thickness    = BL_Waves.(caze).phase(j).thickness;
        displacement = BL_Waves.(caze).phase(j).displacement;
        momentum     = BL_Waves.(caze).phase(j).momentum;
        G = displacement ./ thickness;
        H = displacement ./ momentum;
        scatter(G, H, size, 'o', 'filled', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha)
    end
end

for i = 1:length(no_wave_cases)
    % marker = markers{i};
    color  = WV_colors(i,:);   
    caze   = no_wave_cases{i};
    disp(caze)
    
    % Ensemble
    thickness    = BL_No_Waves.(caze).ensemble.thickness;
    displacement = BL_No_Waves.(caze).ensemble.displacement;
    momentum     = BL_No_Waves.(caze).ensemble.momentum;
    G = displacement ./ thickness;
    H = displacement ./ momentum;
    scatter(G, H, size, 'o', 'filled', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', alpha)
end
hold off

title('Castillo Correltaion')
xlabel('$\delta^* / \delta$', 'Interpreter', 'latex')
ylabel('$H$', 'Interpreter', 'latex')
grid on
xlim([0, 0.6])
ylim([1,4])

%% Boundary Layer Thicknesses

gauss_window = 20;
phase_colors = turbo(4);
line_width   = 1.5;

% ax = figure('name', 'Thicknesses', 'Position', [50,200,1500,400]);
ax = figure('name', 'Thicknesses', 'Position', [200,200,1000,700]);
t = tiledlayout(3,4,'TileSpacing', 'tight', 'Padding', 'tight');
hold on

c = 1;
for s = 1:length(WTs)
    WT = WTs{s};
    for w = 1:length(WVs)
        h(c) = nexttile();
        WV = WVs{w};

        
        % Wave Data 
        caze = strcat('WT', num2str(WT), '_WV', WVs{w}, '_AG0');
        x    = BL_Waves.(caze).x;

        % No Wave Data
        no_wave_caze      = strcat('WT', num2str(WT), '_WV0_AGP');
        no_wave_thickness = BL_No_Waves.(no_wave_caze).ensemble.thickness;
        no_wave_thickness = smoothdata(no_wave_thickness, 'gaussian', gauss_window);

        hold on
        % Phase Thicknesses
        for i = 1:4
            thickness = BL_Waves.(caze).phase(i).thickness;
            thickness = smoothdata(thickness, 'gaussian', gauss_window);
            r(i) = plot(x, thickness, 'Linewidth', line_width, 'color', phase_colors(i,:));
        end

        % Ensemble Thickness
        thickness = BL_Waves.(caze).ensemble.thickness;
        thickness = smoothdata(thickness, 'gaussian', gauss_window);
        r(5) = plot(x, thickness, 'linewidth', line_width, 'color', 'black');

        % Ensemble No Waves Thickness
        r(6) = plot(x, no_wave_thickness, 'LineStyle', '--', 'linewidth', line_width, 'color', 'black');

        hold off
        % axis equal
        xlim([min(x) max(x)])   
        ylim([60,130])

        % Wave Titles
        if WT == '4'
            title(strcat('WV', WV))
        end

        % X axis titles
        if WT == '8'
            xlabel('x [mm]')
        end

        % Y axis titles
        if WV == 'A'
            ylabel('$\delta$ [mm]', 'interpreter', 'latex')
        end
        c = c + 1;

    end
    
end
linkaxes(h, 'xy')
leg = legend(r, 'Phase 1', 'Phase 2', 'Phase 3', 'Phase 4', 'Ensemble', 'Ensemble No Waves', 'orientation', 'horizontal');
leg.Layout.Tile = 'north';

% exportgraphics(ax, fullfile(path, strcat('Thicknesses.png')), 'resolution', 200);


%% Velocity Profiles

lw = 1; 
wave_colors = cool(4);
speed_styles = {'-', '--', '-.'};

% Ensemble
ax = figure('Name', 'Ensemble u profiles');
hold on

% With Waves
for s = 1:length(WTs)
    WT = WTs{s};

    if WT == "4"
        style = speed_styles{1};
    elseif WT == "6"
        style = speed_styles{2};
    elseif WT == "8"
        style = speed_styles{3};
    end


    for w = 1:length(WVs)
        WV = WVs{w};

        if WV == "A"
            color = wave_colors(1,:);
        elseif WV == "B"
            color = wave_colors(2,:);
        elseif WV == "C"
            color = wave_colors(3,:);
        elseif WV == "D"
            color = wave_colors(4,:);
        end

        % Wave Data 
        caze = strcat('WT', num2str(WT), '_WV', WV, '_AG0');
        X    = MN_Waves.(caze).X;
        Y    = MN_Waves.(caze).Y;
        y    = flipud(unique(Y));
        u    = MN_Waves.(caze).ensemble.u;
        u_inf = BL_Waves.(caze).u_inf;

        % Slice
        slice = u(:, round(range/2));
        slice = slice / u_inf;

        plot(slice, y, color=color, LineStyle=style, LineWidth=lw)
    end
end

% No Waves
for s = 1:length(WTs)
    WT = WTs{s};

    if WT == "4"
        style = speed_styles{1};
    elseif WT == "6"
        style = speed_styles{2};
    elseif WT == "8"
        style = speed_styles{3};
    end

    % Wave Data 
    caze = strcat('WT', num2str(WT), '_WV0_AGP');
    X    = MN_No_Waves.(caze).X;
    Y    = MN_No_Waves.(caze).Y;
    y    = flipud(unique(Y));
    u    = MN_No_Waves.(caze).ensemble.u;
    u_inf = BL_No_Waves.(caze).u_inf;

    % Slice
    slice = u(:, round(range/2));
    slice = slice / u_inf;

    plot(slice, y, color='black', LineStyle=style, LineWidth=lw)

end
hold off
ylim([0, 200])
xlim([0.55,1.025])

%% Boundary Layer Thicknesses

gauss_window = 20;
phase_colors = turbo(4);
line_width   = 1.5;

% ax = figure('name', 'Thicknesses', 'Position', [50,200,1500,400]);
ax = figure('name', 'Thicknesses', 'Position', [200,200,1000,700]);
t = tiledlayout(3,4,'TileSpacing', 'tight', 'Padding', 'tight');
hold on

c = 1;
for s = 1:length(WTs)
    WT = WTs{s};
    for w = 1:length(WVs)
        h(c) = nexttile();
        WV = WVs{w};

        
        % Wave Data 
        caze = strcat('WT', num2str(WT), '_WV', WVs{w}, '_AG0');
        x    = BL_Waves.(caze).x;

        % No Wave Data
        no_wave_caze      = strcat('WT', num2str(WT), '_WV0_AGP');
        no_wave_thickness = BL_No_Waves.(no_wave_caze).ensemble.thickness;
        no_wave_thickness = smoothdata(no_wave_thickness, 'gaussian', gauss_window);

        hold on
        % Phase Thicknesses
        for i = 1:4
            thickness = BL_Waves.(caze).phase(i).thickness;
            thickness = smoothdata(thickness, 'gaussian', gauss_window);
            r(i) = plot(x, thickness, 'Linewidth', line_width, 'color', phase_colors(i,:));
        end

        % Ensemble Thickness
        thickness = BL_Waves.(caze).ensemble.thickness;
        thickness = smoothdata(thickness, 'gaussian', gauss_window);
        r(5) = plot(x, thickness, 'linewidth', line_width, 'color', 'black');

        % Ensemble No Waves Thickness
        r(6) = plot(x, no_wave_thickness, 'LineStyle', '--', 'linewidth', line_width, 'color', 'black');

        hold off
        % axis equal
        xlim([min(x) max(x)])   
        ylim([60,130])

        % Wave Titles
        if WT == '4'
            title(strcat('WV', WV))
        end

        % X axis titles
        if WT == '8'
            xlabel('x [mm]')
        end

        % Y axis titles
        if WV == 'A'
            ylabel('$\delta$ [mm]', 'interpreter', 'latex')
        end
        c = c + 1;

    end
    
end
linkaxes(h, 'xy')
leg = legend(r, 'Phase 1', 'Phase 2', 'Phase 3', 'Phase 4', 'Ensemble', 'Ensemble No Waves', 'orientation', 'horizontal');
leg.Layout.Tile = 'north';

% exportgraphics(ax, fullfile(path, strcat('Thicknesses.png')), 'resolution', 200);



%% Ensemble Shape Factors

clc;
gauss_window = 15;
WT_colors = turbo(length(WTs));
WV_colors = cool(length(WVs));
markers   = {'o', 'square', 'diamond', '^'};
linewidth = 1.5;

figure('Name', 'Ensemble Shape Factors')
hold on
for i = 1:length(wave_cases)
    
    % color = WV_colors(i,:);   
    caze  = wave_cases{i};
    disp(caze)

    % Color based on wave
    if contains(caze, 'WVA') == 1
        color = WV_colors(1,:);
    elseif contains(caze, 'WVB') == 1
        color = WV_colors(2,:);
    elseif contains(caze, 'WVC') == 1
        color = WV_colors(3,:);
    elseif contains(caze, 'WVD') == 1
        color = WV_colors(4,:);
    end
    
    % Linestyle Based on wind speed
    if contains(caze, 'WT4') == 1
        ls = '-';
    elseif contains(caze, 'WT6') == 1
        ls = '--';
    elseif contains(caze, 'WT8') == 1
        ls = ':';
    end

    % Ensemble
    x = BL_Waves.(caze).x;
    thickness    = BL_Waves.(caze).ensemble.thickness;
    displacement = BL_Waves.(caze).ensemble.displacement;
    momentum     = BL_Waves.(caze).ensemble.momentum;
    H = displacement ./ momentum;
    smooth_H = smoothdata(H, 'gaussian', gauss_window);
    plot(x, smooth_H, 'color', color, 'linewidth', linewidth, 'linestyle', ls)

end

for i = 1:length(no_wave_cases)
    caze = no_wave_cases{i};
    disp(caze)

    % Linestyle Based on wind speed
    if contains(caze, 'WT4') == 1
        ls = '-';
    elseif contains(caze, 'WT6') == 1
        ls = '--';
    elseif contains(caze, 'WT8') == 1
        ls = ':';
    end


    % Ensemble
    thickness    = BL_No_Waves.(caze).ensemble.thickness;
    displacement = BL_No_Waves.(caze).ensemble.displacement;
    momentum     = BL_No_Waves.(caze).ensemble.momentum;
    H = displacement ./ momentum;
    plot(x, H, 'color', 'black', 'linewidth', linewidth, 'linestyle', ls)
end
hold off

title('Ensemble Average Shape Factors')
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$H$', 'Interpreter', 'latex')
grid on
xlim([0, max(x)])
ylim([1,1.5])


%% Phase Average Shape Factors [x normalized by wavelength and centered]

clc;
WV_colors = cool(length(WVs));
gauss_window = 3;
linewidth = 2;
FS = 16;

ax = figure('Name', 'Phase Shape Factors Centered Normalized', 'Position', [100,100,800,400]);
t = tiledlayout(2,1,'padding', 'tight');

clear h r
for j = 1:2:3
    h(j) = nexttile();
    hold on
    for i = 1:length(wave_cases)
         
        caze  = wave_cases{i};
        disp(caze)
    
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
        
        % Linestyle Based on wind speed
        if contains(caze, 'WT4') == 1
            ls = '-';
        elseif contains(caze, 'WT6') == 1
            ls = '--';
        elseif contains(caze, 'WT8') == 1
            ls = '-.';
        end

        % Phase Average
        x = BL_Waves.(caze).x;
        displacement = BL_Waves.(caze).phase(j).displacement;
        momentum     = BL_Waves.(caze).phase(j).momentum;

        H = displacement ./ momentum;
        smooth_H = smoothdata(H, 'gaussian', gauss_window);
        r(i) = plot((x - range/2)/wave_length, smooth_H, 'color', color, 'linewidth', linewidth, 'linestyle', ls);

        % Reference Wave
        if contains(caze, 'WVD') == 1
            wave_scale = 0.17;
            wave_offset = 1.32;
            ref_wave = MN_Waves.(caze).phase(j).reference_wave / wave_amplitude;
            ref_wave = (wave_scale * ref_wave) + wave_offset;
            p = plot((x - range/2)/wave_length, ref_wave, 'color', 'black', 'linewidth', 4);
            p.Color(4) = 0.2;
            uistack(p, 'bottom')
        end

    end
    hold off
    tit = strcat('Phase', {' '}, num2str(j));
    title(tit{1})
    ylabel('$H$', 'Interpreter', 'latex', 'fontsize', FS)
    grid on
    xlim([-1, 1])
end
xlabel('$x / \lambda$', 'Interpreter', 'latex', 'fontsize', FS)
linkaxes(h, 'xy')
leg = legend(r, 'WVA', 'WVB', 'WVC', 'WVD', 'orientation', 'horizontal');
leg.Layout.Tile = 'north';

exportgraphics(ax, fullfile(path, strcat('Phase_Shape.png')), 'resolution', 200);

%% Stitched Momentum Phase Average 

WT = 4;
gauss_window = 1;
WV_colors = cool(length(WVs));
kinematic_viscosity = 1.48E-5;

ax = figure('name', 'momentum stitched');

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
        x        = BL_Waves.(caze).x;
        u_inf    = BL_Waves.(caze).u_inf;
        momentum = BL_Waves.(caze).phase(i).momentum;
        momentum = smoothdata(momentum, 'gaussian', gauss_window);
        Re_momentum = (u_inf * momentum * 1E-3) / kinematic_viscosity;
    
        % Waves
        centered_normalized_x = ((x - (range/2)) / wave_length) - ((i - 1)/4);
        reference_wave = MN_Waves.(caze).phase(i).reference_wave / wave_amplitude;
        
        % Scale Wave
        wave_scale = 600;
        wave_offset = 1800;
        reference_wave = reference_wave * wave_scale;
        reference_wave = reference_wave + wave_offset;
     
        % Plot
        a = plot(centered_normalized_x, reference_wave, 'color', 'black', 'linewidth', 2, 'HandleVisibility', 'off');
        b(w) = plot(centered_normalized_x, Re_momentum, 'color', color, 'linewidth', 2);
        a.Color(4) = 1;
        uistack(a, 'bottom')
        
    end
end
xline(-1, 'color', 'black', 'linestyle', '--', 'linewidth', 2)
xline(-0.5, 'color', 'black', 'linestyle', '--', 'linewidth', 2)
xline(0, 'color', 'black', 'linestyle', '--', 'linewidth', 2)
hold off

xlim([-1.1,0.1])
legend(b, 'WVA', 'WVB', 'WVC', 'WVD');
xlabel('$x / \lambda$', 'Interpreter', 'latex')
ylabel('$Re_{\theta}$', 'Interpreter', 'latex')
title(strcat('WT', num2str(WT), ' All Waves: Stitched Phases'));


%% Ensemble Displacement and Momentum Thicknesses

clc;
gauss_window = 15;
WT_colors = turbo(length(WTs));
WV_colors = cool(length(WVs));
markers   = {'o', 'square', 'diamond', '^'};
linewidth = 1.5;

ax = figure('Name', 'Ensemble displacment and momentum');
t = tiledlayout(2,1,'TileSpacing','compact','Padding', 'tight');
clear h
for j = 1:2
    f(j) = nexttile();
    hold on
    c = 1;
    for i = 1:length(wave_cases)
        
        % color = WV_colors(i,:);   
        caze  = wave_cases{i};
        disp(caze)
    
        % Color based on wave
        if contains(caze, 'WVA') == 1
            color = WV_colors(1,:);
        elseif contains(caze, 'WVB') == 1
            color = WV_colors(2,:);
        elseif contains(caze, 'WVC') == 1
            color = WV_colors(3,:);
        elseif contains(caze, 'WVD') == 1
            color = WV_colors(4,:);
        end
        
        % Linestyle Based on wind speed
        if contains(caze, 'WT4') == 1
            ls = '-';
        elseif contains(caze, 'WT6') == 1
            ls = '--';
        elseif contains(caze, 'WT8') == 1
            ls = ':';
        end
    
        % Ensemble
        x = BL_Waves.(caze).x;
        thickness    = BL_Waves.(caze).ensemble.thickness;
        displacement = BL_Waves.(caze).ensemble.displacement;
        momentum     = BL_Waves.(caze).ensemble.momentum;

        smooth_displacement = smoothdata(displacement, 'gaussian', gauss_window);
        smooth_momentum     = smoothdata(momentum, 'gaussian', gauss_window);
        
        if j == 1
            plot(x, smooth_displacement, 'color', color, 'linewidth', linewidth, 'linestyle', ls)
        elseif j == 2
            plot(x, smooth_momentum, 'color', color, 'linewidth', linewidth, 'linestyle', ls)
        end
    end
    
    c = 1;
    for i = 1:length(no_wave_cases)
        caze = no_wave_cases{i};
        disp(caze)
    
        % Linestyle Based on wind speed
        if contains(caze, 'WT4') == 1
            ls = '-';
        elseif contains(caze, 'WT6') == 1
            ls = '--';
        elseif contains(caze, 'WT8') == 1
            ls = ':';
        end
    
    
        % Ensemble
        thickness    = BL_No_Waves.(caze).ensemble.thickness;
        displacement = BL_No_Waves.(caze).ensemble.displacement;
        momentum     = BL_No_Waves.(caze).ensemble.momentum;

        smooth_displacement = smoothdata(displacement, 'gaussian', gauss_window);
        smooth_momentum     = smoothdata(momentum, 'gaussian', gauss_window);
        
        if j == 1
            plot(x, smooth_displacement, 'color', color, 'linewidth', linewidth, 'linestyle', ls)

        elseif j == 2
            plot(x, smooth_momentum, 'color', color, 'linewidth', linewidth, 'linestyle', ls)
        end
        c = c + 1;
    end
    hold off
    xlim([0, max(x)])
    
    if j == 1
        title('Displacement Thickness')
        ylabel('$\delta^* [mm]$', 'Interpreter', 'latex', 'fontsize', 16)
        % legend(h, 'WVA', 'WVB', 'WVC', 'WVD')
    elseif j ==2
        title('Momentum Thickness')
        xlabel('$x [mm]$', 'Interpreter', 'latex', 'fontsize', 16)
        ylabel('$\theta [mm]$', 'Interpreter', 'latex', 'fontsize', 16)
    end
end
linkaxes(f,'y')




%% Phase Average Displacement, Momentum [x normalized by wavelength and centered]

clc;
WV_colors    = cool(length(WVs));
gauss_window = 3;
linewidth    = 2;
phase        = 1;
FS = 16;

ax = figure('Name', 'Phase Displacement and Momentum', 'Position', [100,100,800,400]);
t = tiledlayout(2,1);

clear h r
for j = 1:2
    h(j) = nexttile();
    hold on
    for i = 1:length(wave_cases)
         
        caze  = wave_cases{i};
        disp(caze)
    
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
        
        % Linestyle Based on wind speed
        if contains(caze, 'WT4') == 1
            ls = '-';
        elseif contains(caze, 'WT6') == 1
            ls = '--';
        elseif contains(caze, 'WT8') == 1
            ls = '-.';
        end

        % Phase Average
        x = BL_Waves.(caze).x;
        centered_normalized_x = (x - range/2)/wave_length;

        displacement = BL_Waves.(caze).phase(phase).displacement;
        momentum     = BL_Waves.(caze).phase(phase).momentum;

        H = displacement ./ momentum;
        smooth_H = smoothdata(H, 'gaussian', gauss_window);
        smooth_displacement = smoothdata(displacement, 'gaussian', gauss_window);
        smooth_momentum = smoothdata(momentum, 'gaussian', gauss_window);

        if j == 1
            data = smooth_displacement;
            tit = 'Displacement Thickness';
            ylabel('$\delta^*$ [mm]', 'Interpreter', 'latex', 'fontsize', FS)
            wave_scale = 9;
            wave_offset = 18;

        elseif j == 2
            data = smooth_momentum;
            tit = 'Momentum Thickness';
            ylabel('$\theta$ [mm]', 'Interpreter', 'latex', 'fontsize', FS)
            wave_scale = 5;
            wave_offset = 14;
            
        elseif j == 3
            data = smooth_H;
            tit = 'Shape Factor';
            ylabel('$H$', 'Interpreter', 'latex')
            wave_scale = 0.17;
            wave_offset = 1.32;
        end

        r(i) = plot(centered_normalized_x, data, 'color', color, 'linewidth', linewidth, 'linestyle', ls);
        title(tit)

        % Reference Wave
        if contains(caze, 'WVD') == 1 
            ref_wave = MN_Waves.(caze).phase(phase).reference_wave / wave_amplitude;
            ref_wave = (wave_scale * ref_wave) + wave_offset;
            p = plot((x - range/2)/wave_length, ref_wave, 'color', 'black', 'linewidth', 4);
            p.Color(4) = 0.2;
            uistack(p, 'bottom')
        end

    end
    hold off
    grid on
    xlim([-1, 1])
end
xlabel('$x / \lambda$', 'Interpreter', 'latex', 'fontsize', FS)
linkaxes(h, 'xy')
ylim([0, 35])
leg = legend(r, 'WVA', 'WVB', 'WVC', 'WVD', 'orientation', 'horizontal');
leg.Layout.Tile = 'north';



exportgraphics(ax, fullfile(path, strcat("Phase_", num2str(phase), "_Displacement_Momentum.png")), resolution=300)












