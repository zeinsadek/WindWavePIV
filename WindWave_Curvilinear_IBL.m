%%% Curvilinear Integrated Boundary Layer Equation

addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
means_path = fullfile(project_path, 'means');
curvilinear_path = fullfile(project_path, 'curvilinear');

% Case
caze = 'WT8_WVD_AG0';

% Load data
curvilinear = load(fullfile(curvilinear_path, strcat(caze, '_CURVILINEAR.mat')));
curvilinear = curvilinear.output;

% Wave Parameters
wave_type       = caze(strfind(caze, 'WV') + 2);
wave_parameters = readcell("Offshore_Waves.xlsx");
wavelength      = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 2};
amplitude       = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 3};

clear means_path curvilinear_path project_path wave_parameters wave_type

%% Compute Curvilinear IBL

% Cropping bounds
left = 5E-3;
right = 230E-3;
buffer = 3E-3;
top = 200E-3;

% Kinematic Viscosity
nu = 1.48E-5;

for phase = 1:4

    % Monitor
    fprintf("Phase %1.0f\n\n", phase)

    % X (xi) in mm
    X = curvilinear.X * 1E-3;
    x = X(1,:);
    dx = mean(diff(X(1,:)));
    vertical_lines = curvilinear.phase(phase).vertical_lines * 1E-3;
    
    % Y (zeta) in mm
    Y = curvilinear.Y * 1E-3;
    y = Y(:,1);
    dy = mean(diff(Y(:,1)));
    horizontal_lines = curvilinear.phase(phase).horizontal_lines * 1E-3;
    
    % Waves
    wave_profile = curvilinear.phase(phase).wave_profile * 1E-3;
    max_wave_profile = curvilinear.phase(phase).max_wave_profile * 1E-3;
    bottom = max_wave_profile + buffer;
    
    % Velocities
    u = curvilinear.phase(phase).u;
    v = curvilinear.phase(phase).v;
    
    % Stresses
    uu = curvilinear.phase(phase).uu;
    vv = curvilinear.phase(phase).vv;
    uv = curvilinear.phase(phase).uv;
    
    % dx
    [dxi, ~] = gradient(vertical_lines, 1, 1);
    
    % dy
    [~, dzeta] = gradient(horizontal_lines, 1, 1);


    %%% Term 1: Viscous Stress
    [~, du_dy] = gradient(u, 1, 1);
    du_dzeta = du_dy ./ dzeta;
    term1 = nu * du_dzeta;
    term1 = curvilinearCrop(term1, vertical_lines, horizontal_lines, left, right, bottom, top);
    
    % Save
    IBL.phase(phase).term1 = term1;
    clear du_dy du_dzeta


    %%% Term 2: Shear Stress
    term2 = -uv;
    term2 = curvilinearCrop(term2, vertical_lines, horizontal_lines, left, right, bottom, top);
    
    % Save
    IBL.phase(phase).term2 = term2;


    %%% Term 3: Mean Momentum Flux
    [du2_dx, ~] = gradient(u.^2, 1, 1);
    du2_dxi = du2_dx ./ dxi;
    int_du2_dxi = curvilinearIntegrate(du2_dxi, horizontal_lines);
    
    [du_dx, ~] = gradient(u, 1, 1);
    du_dxi = du_dx ./ dxi;
    int_du_dxi = curvilinearIntegrate(du_dxi, horizontal_lines);
    u_int_du_dxi = u .* int_du_dxi;
    
    term3 = -1 * int_du2_dxi + u_int_du_dxi;
    term3 = curvilinearCrop(term3, vertical_lines, horizontal_lines, left, right, bottom, top);
    
    % Save
    IBL.phase(phase).term3 = term3;
    clear du2_dx du2_dxi int_du2_dxi du_dx du_dxi int_du_dxi u_int_du_dxi

    
    %%% Term 4: Turbulent Momentum Flux
    [duu_dx, ~] = gradient(uu, 1, 1);
    duu_dxi = duu_dx ./ dxi;
    int_duu_dxi = curvilinearIntegrate(duu_dxi, horizontal_lines);
    
    term4 = -1 * int_duu_dxi;
    term4 = curvilinearCrop(term4, vertical_lines, horizontal_lines, left, right, bottom, top);

    % Save
    IBL.phase(phase).term4 = term4;
    clear duu_dx duu_dxi int_duu_dxi


    %%% Term 5: Pressure Correction
    [dvv_dx, ~] = gradient(vv, 1, 1);
    dvv_dxi = dvv_dx ./ dxi;
    int_dvv_dxi = curvilinearIntegrate(dvv_dxi, horizontal_lines);
    
    term5 = int_dvv_dxi;
    term5 = curvilinearCrop(term5, vertical_lines, horizontal_lines, left, right, bottom, top);

    % Save
    IBL.phase(phase).term5 = term5;
    clear dvv_dx dvv_dxi int_dvv_dxi

end


%% Plotting

levels = 100;
lw = 3;

for phase = 1:4

    % Monitor
    fprintf("Phase %1.0f\n\n", phase)

    % X (xi) in mm
    X = curvilinear.X * 1E-3;
    x = X(1,:);
    dx = mean(diff(X(1,:)));
    vertical_lines = curvilinear.phase(phase).vertical_lines * 1E-3;
    
    % Y (zeta) in mm
    Y = curvilinear.Y * 1E-3;
    y = Y(:,1);
    dy = mean(diff(Y(:,1)));
    horizontal_lines = curvilinear.phase(phase).horizontal_lines * 1E-3;
    
    % Waves
    wave_profile = curvilinear.phase(phase).wave_profile * 1E-3;
    max_wave_profile = curvilinear.phase(phase).max_wave_profile * 1E-3;
    
    terms = {'term1', 'term2', 'term3', 'term4', 'term5'};
    figure()
    sgtitle(sprintf('IBL %s: Phase %1.0f', caze, phase))
    tiledlayout(1, length(terms))
    for t = 1:length(terms)
        term = terms{t};
        h(t) = nexttile;
        hold on
        contourf(vertical_lines, horizontal_lines, IBL.phase(phase).(term), levels, 'linestyle', 'none')
        plot(x, wave_profile, 'color', 'black', 'linewidth', lw)
        plot(x, max_wave_profile, 'color', 'red', 'linewidth', lw)
        hold off
        axis equal
        title(term)
        colorbar
    end
    linkaxes(h, 'xy')
end



%% Functions

% Curvilienar Integration
function output = curvilinearIntegrate(quantity, horizontal_lines)
    output = nan(size(quantity));
    for i = 1:size(quantity, 2)
        zeta_column = horizontal_lines(:,i) * 1E-3;
        data_column = quantity(:,i);
        nan_mask = ~isnan(zeta_column) & ~isnan(data_column);
        
        if sum(nan_mask) > 1
            output(nan_mask, i) = flipud(cumtrapz(flipud(zeta_column(nan_mask)), flipud(data_column(nan_mask))));
        end
    end
end


% Crop Curvilinear Data
function output = curvilinearCrop(quantity, vertical_lines, horizontal_lines, left, right, bottom, top)
    output = quantity;
    output(vertical_lines < left) = nan;
    output(vertical_lines > right) = nan;
    output(horizontal_lines < bottom) = nan;
    output(horizontal_lines > top) = nan;
end