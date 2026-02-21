%%% Curvilinear Coordinates from ChatGPT

addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')

% Paths
clc; clear; close all;
project_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV";
means_path = fullfile(project_path, 'means');
curvilinear_path = fullfile(project_path, 'curvilinear_new');

% Case
caze = 'WT8_WVD_AG0';

% Figures
figure_path = fullfile(project_path, "figures", caze, "curvilinear_new");

if ~exist(figure_path, 'dir')
    mkdir(figure_path);
end

% Load data
means = load(fullfile(means_path, strcat(caze, '_MEANS.mat')));
means = means.output;
wave_type = caze(strfind(caze, 'WV') + 2);

% Wave Parameters
wave_parameters = readcell("Offshore_Waves.xlsx");
wavelength      = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 2};
amplitude       = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 3};
phase_offset    = [0, wavelength/4, wavelength/2, 3*wavelength/4];

%%% BIG CORRECTION
% wavenumber      = 1 / wavelength;
wavenumber = (2 * pi) / wavelength;

% OG coordinates
X  = means.X;
Y  = means.Y;
x  = X(1,:);
y  = Y(:,1);
dx = mean(diff(x));
dy = mean(diff(y));
components = {'u', 'v', 'uu', 'vv', 'uv'};

% Cartesian mesh for convenience
output.X = X;
output.Y = Y;

for p = 1:4

    % Loop through phases
    phase = p;
    fprintf("Curvilinear phase %1.0f\n", phase);

    % Centering wave profile around data
    x_shift = mean(x);
    max_wave_profile = means.phase(phase).max_wave_profile;
    top_crop = 205;
    
    wave_profile = amplitude * cos(2 * pi * (x - x_shift - phase_offset(phase)) / wavelength);
    wave_orthogonal = (2 * pi / wavelength) * amplitude * sin(2 * pi * (x - x_shift - phase_offset(phase)) / wavelength);
    
    clear x_shift wave_parameters means_path wave_type
    
    % Save wave profiles
    output.phase(phase).max_wave_profile = max_wave_profile;
    output.phase(phase).wave_profile = wave_profile;


    for c = 1:length(components)
        component = components{c};
    
        % Load component
        tmp = means.phase(phase).(component);
    
        % Shift RHS
        LHS_values = tmp(:,86);
        RHS_values = tmp(:,87);
        differences = RHS_values - LHS_values;
        differences(isnan(differences)) = 0;
        corrected = tmp;
        corrected(:,87:end) = corrected(:,87:end) - differences;
    
        % Interpolate across strip
        [~, left] = min(abs(x - 115));
        [~, right] = min(abs(x - 120));
        corrected(:, left:right) = nan;
        corrected = fillmissing(corrected, 'makima', 2);
        corrected(Y < max_wave_profile) = nan;

        % Make sure original nans are in place
        corrected(isnan(tmp)) = nan;
    
        % Resave data
        means_corrected.phase(phase).(component) = corrected;
        clear corrected
    end


    % Velocity/Stress data
    u  = means_corrected.phase(phase).u;
    v  = means_corrected.phase(phase).v;
    uu = means_corrected.phase(phase).uu;
    vv = means_corrected.phase(phase).vv;
    uv = means_corrected.phase(phase).uv;
    
    clear tmp LHS_values RHS_values differences left right
    clear c component levels

    % Build Extended Mesh 
    % Output size will match input PIV grid
    [rows, columns]  = size(X);
    horizontal_lines = nan(rows, columns);
    vertical_lines   = nan(rows, columns);
    
    % Horizontal (Zeta) Lines
    for r = 1:rows
        y0 = y(r);
        horizontal_lines(r, :) = y0 + wave_profile .* exp(-wavenumber * y0);
    end
    
    % Build Vertical (Xi) Lines 
    for c = 1:columns
        x0 = x(c);
        % Zero to wave surface
        y_zeroed = y - wave_profile(c);
        vertical_lines(:, c) = (y_zeroed .* wave_orthogonal(c) .* exp(-wavenumber * y_zeroed)) + x0;
    end
    
    % Crop below wave profile
    vertical_lines(horizontal_lines < wave_profile) = nan;
    horizontal_lines(horizontal_lines < wave_profile) = nan;
    
    clear r c x0 y0 y_zeroed wave_orthogonal

    % Save curvilinear mesh
    output.phase(phase).vertical_lines = vertical_lines;
    output.phase(phase).horizontal_lines = horizontal_lines;


    % Interpolate onto Curvilinear Coordinates 
    % Fill NaNs in original data just in case
    u_filled = u;
    v_filled = v;
    u_filled(isnan(u_filled)) = 0;
    v_filled(isnan(v_filled)) = 0;
    
    % Interpolate
    u_interp = interp2(X, Y, u_filled, vertical_lines, horizontal_lines, 'cubic', nan);
    v_interp = interp2(X, Y, v_filled, vertical_lines, horizontal_lines, 'cubic', nan);
    
    u_interp(u_interp == 0) = nan;
    v_interp(v_interp == 0) = nan;
    
    % Crop wave
    u_interp(horizontal_lines < max_wave_profile) = nan;
    v_interp(horizontal_lines < max_wave_profile) = nan;
    
    % Crop top
    u_interp(horizontal_lines > top_crop) = nan;
    v_interp(horizontal_lines > top_crop) = nan;
    
    clear u_filled v_filled 


    % Computing wave tangent angles
    % Compute slope of tangent lines across curvilinear space
    [dz_dx, dz_dy] = gradient(horizontal_lines, dx, dy);
    dz_dx(horizontal_lines < wave_profile) = nan;
    dz_dy(horizontal_lines < wave_profile) = nan;
    
    angles = atan2(dz_dx, dz_dy);
    angles(horizontal_lines < wave_profile) = nan;
    
    clear dz_dx dz_dy

    % Save angles
    output.phase(phase).angles = angles;


    % Projecting velocities from cartesian onto curvilinear
    % PASSIVE TRANSFORMATION
    
    % Tangent to ζ line (ξ direction)
    % ξ points "along the wave"
    xi_hat_x = cos(angles);    
    xi_hat_y = sin(angles);
    
    % Normal to ζ line (ζ direction, ~wave-normal)
    % ζ is perpendicular to ξ
    zeta_hat_x = -sin(angles);   
    zeta_hat_y =  cos(angles);
    
    % Project Cartesian velocities onto curvilinear axes
    u_tangent = u_interp .* xi_hat_x + v_interp .* xi_hat_y;
    u_normal  = -(u_interp .* zeta_hat_x + v_interp .* zeta_hat_y);
    
    % Clean up top of data
    u_tangent(horizontal_lines > top_crop) = nan;
    u_normal(horizontal_lines > top_crop) = nan;

    % Save curvilinear velocities
    output.phase(phase).u = u_tangent;
    output.phase(phase).v = u_normal;



    % Interpolating stresses
    % Fill NaNs in original data
    uu(isnan(uu)) = 0;
    vv(isnan(vv)) = 0;
    uv(isnan(uv)) = 0;
    
    % Interpolate
    uu_interp = interp2(X, Y, uu, vertical_lines, horizontal_lines, 'cubic', nan);
    vv_interp = interp2(X, Y, vv, vertical_lines, horizontal_lines, 'cubic', nan);
    uv_interp = interp2(X, Y, uv, vertical_lines, horizontal_lines, 'cubic', nan);
    
    % Return nans
    uu_interp(uu_interp == 0) = nan;
    vv_interp(vv_interp == 0) = nan;
    uv_interp(vv_interp == 0) = nan;
    
    % Crop below wave
    uu_interp(horizontal_lines < max_wave_profile) = nan;
    vv_interp(horizontal_lines < max_wave_profile) = nan;
    uv_interp(horizontal_lines < max_wave_profile) = nan;



    % Projecting stresses
    % Projected normal stresses
    tau_xixi   = uu_interp .* xi_hat_x.^2 + 2 * uv_interp .* xi_hat_x .* xi_hat_y + vv_interp .* xi_hat_y.^2;
    tau_zetazeta = uu_interp .* zeta_hat_x.^2 + 2 * uv_interp .* zeta_hat_x .* zeta_hat_y + vv_interp .* zeta_hat_y.^2;
    
    % Projected shear stress (symmetric in turbulence, so either direction is fine)
    tau_xizeta = uu_interp .* xi_hat_x .* zeta_hat_x + uv_interp .* (xi_hat_x .* zeta_hat_y + xi_hat_y .* zeta_hat_x) + vv_interp .* xi_hat_y .* zeta_hat_y;
    
    % Crop below wave
    tau_xixi(horizontal_lines < max_wave_profile) = nan;
    tau_zetazeta(horizontal_lines < max_wave_profile) = nan;
    tau_xizeta(horizontal_lines < max_wave_profile) = nan;
    
    % Crop top
    tau_xixi(horizontal_lines > top_crop) = nan;
    tau_zetazeta(horizontal_lines > top_crop) = nan;
    tau_xizeta(horizontal_lines > top_crop) = nan;

    % Save curvilinear stresses
    output.phase(phase).uu = tau_xixi;
    output.phase(phase).vv = tau_zetazeta;
    output.phase(phase).uv = tau_xizeta;

    clear tau_xixi tau_zetazeta tau_xizeta uu_interp vv_interp uv_interp
end

clear angles dx dy max_wave_profile means means_corrected p phase phase_offset top_crop
clear u v uu vv uv u_interp u_normal u_tangent v_interp wave_profile rows columns x y X Y
clear xi_hat_x xi_hat_y zeta_hat_x zeta_hat_y amplitude wavelength wavenumber vertical_lines horizontal_lines


%% Plot to check

curvilinear_titles = {'Wave Tangent Velocity',...
                      'Wave Normal Velocity', ...
                      'Wave Tangent, Normal Stress', ...
                      'Wave Normal, Normal Stress', ...
                      'Curvilinear Shear Stress'};

lw = 2;
spacing = 10;
clc; close all
for c = 1:length(components)
    component = components{c};
    cmin = min([output.phase(:).(component)], [], 'all');
    cmax = max([output.phase(:).(component)], [], 'all');

    ax = figure();
    tiledlayout(2,2)
    sgtitle(sprintf('%s: %s', caze, curvilinear_titles{c}), 'interpreter', 'none')
    for p = 1:4
        h(p) = nexttile;
        hold on
        contourf(output.phase(p).vertical_lines, output.phase(p).horizontal_lines, output.phase(p).(component), 100, 'linestyle', 'none')
        plot(output.X(1,:), output.phase(p).max_wave_profile, 'color', 'red', 'linewidth', lw)
        plot(output.X(1,:), output.phase(p).wave_profile, 'color', 'black', 'linewidth', lw)

        % Grid lines
        [row, col] = size(output.phase(p).vertical_lines);
        for i = 1:spacing:row
            P = plot(output.X(1,:), output.phase(p).horizontal_lines(i,:), 'color', 'black');
            P.Color(4) = 0.25;
        end

        for i = 1:spacing:col
            P = plot(output.phase(p).vertical_lines(:,i), output.Y(:,1) + output.phase(p).wave_profile(i), 'color', 'black');
            P.Color(4) = 0.25;
        end

        hold off
        axis equal
        colorbar()
        clim([cmin, cmax])
        title(sprintf("Phase %1.0f", p))
    end
    linkaxes(h,'xy')

    % Save figure
    figure_name = strcat(caze, "_Curvilinear_", component, ".png");
    exportgraphics(ax, fullfile(figure_path, figure_name), 'resolution', 300)
end

close all
clear curvilinear_titles col p h lw spacing component cmin cmax row col p i

% Save to matfile
fprintf("Saving Matfile\n")
filename = strcat(caze, "_CURVILINEAR.mat");
save(fullfile(curvilinear_path, filename), 'output')
fprintf("Matfile Saved!\n\n")

