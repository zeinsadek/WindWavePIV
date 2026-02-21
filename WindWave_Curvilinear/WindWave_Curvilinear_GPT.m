%%% Curvilinear Coordinates from ChatGPT

addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')

% Load data
clc; clear; close all;
means_path = "/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means";
caze = 'WT4_WVA_AG0';

means = load(fullfile(means_path, strcat(caze, '_MEANS.mat')));
means = means.output;
wave_type = caze(strfind(caze, 'WV') + 2);
phase = 2;

% Wave Parameters
wave_parameters = readcell("Offshore_Waves.xlsx");
wavelength      = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 2};
amplitude       = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 3};
phase_offset    = [0, wavelength/4, wavelength/2, 3*wavelength/4];

%%% BIG CORRECTION
% wavenumber      = 1 / wavelength;
wavenumber      = (2 * pi) / wavelength;

% OG coordinates
X  = means.X;
Y  = means.Y;
x  = X(1,:);
y  = Y(:,1);
dx = mean(diff(x));
dy = mean(diff(y));

% Centering wave profile around data
x_shift = mean(x);
max_wave_profile = means.phase(phase).max_wave_profile;
top_crop = 205;

wave_profile = amplitude * cos(2 * pi * (x - x_shift - phase_offset(phase)) / wavelength);
wave_orthogonal = (2 * pi / wavelength) * amplitude * sin(2 * pi * (x - x_shift - phase_offset(phase)) / wavelength);

clear x_shift wave_parameters means_path wave_type

%% Repair cat-scratch

components = {'u', 'v', 'uu', 'vv', 'uv'};

clc;
for c = 1:length(components)
    component = components{c};

    % Monitor
    fprintf("%s\n", component)

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
    % clear corrected

    fprintf("\n")
end

levels = 10;
% Contours
for c = 1:length(components)
    component = components{c};
    figure('color', 'white')
    title(component)
    hold on
    contourf(X, Y, means_corrected.phase(phase).(component), levels, 'linestyle', 'none')
    plot(x, max_wave_profile, 'linewidth', 3, 'Color', 'blue')
    plot(x, wave_profile, 'linewidth', 3)
    hold off
    colorbar
    axis equal
end


% Velocity/Stress data
u  = means_corrected.phase(phase).u;
v  = means_corrected.phase(phase).v;
uu = means_corrected.phase(phase).uu;
vv = means_corrected.phase(phase).vv;
uv = means_corrected.phase(phase).uv;

%%
clear means components tmp LHS_values RHS_values differences left right
clear c component levels

%% Build Extended Mesh 

clc;
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


%% Plot to Check 

pad = 40;
sz = 30;
spacing = 5;
colors = parula(numel(vertical_lines(1:spacing:end, 1:spacing:end)));

figure('color', 'white')
hold on

% Plot curves
for i = 1:spacing:rows
    plot(x, horizontal_lines(i,:), 'color', 'black')
end
for i = 1:spacing:columns
    plot(vertical_lines(:,i), y + wave_profile(i), 'color', 'black')
end

% Plot scatter of points that we will interpolate at
c = 1;
for i = 1:spacing:rows
    for j = 1:spacing:columns
        scatter(vertical_lines(i,j), horizontal_lines(i,j), sz, 'filled', 'MarkerFaceColor', colors(c,:))
        c = c + 1;
    end
end

plot(x, max_wave_profile, 'linewidth', 3, 'Color', 'blue')
plot(x, wave_profile, 'linewidth', 3)
hold off
axis equal
xline(min(x))
xline(max(x))
yline(min(y))
yline(max(y))
xlim([min(x) - pad, max(x) + pad])
ylim([min(y) - pad, max(y) + pad])

clear c i j pad sz colors spacing


%% Interpolate onto Curvilinear Coordinates 

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

% Plot Result 
% Plotted on Curvilinear Grid
figure('color', 'white')
tiledlayout(2,2)

nexttile
hold on
contourf(vertical_lines, horizontal_lines, u_interp, 500, 'linestyle', 'none')
plot(x, max_wave_profile, 'linewidth', 3)
plot(x, wave_profile, 'linewidth', 3)
hold off
colorbar()
title('U')
axis equal
xlabel('x [mm]')
ylabel('y [mm]')
title('u Velocity Field in Curvilinear Coordinates')

nexttile
hold on
contourf(vertical_lines, horizontal_lines, v_interp, 500, 'linestyle', 'none')
plot(x, max_wave_profile, 'linewidth', 3)
plot(x, wave_profile, 'linewidth', 3)
hold off
colorbar()
title('V')
axis equal
xlabel('x [mm]')
ylabel('y [mm]')
title('v Velocity Field in Curvilinear Coordinates')


% Plotted on Cartesian Grid
nexttile
contourf(X, Y, u_interp, 500, 'linestyle', 'none')
colorbar()
title('U')
axis equal
xlabel('x [mm]')
ylabel('y [mm]')
title('u Velocity Field in Cartesian Coordinates')

nexttile
contourf(X, Y, v_interp, 500, 'linestyle', 'none')
colorbar()
title('V')
axis equal
xlabel('x [mm]')
ylabel('y [mm]')
title('v Velocity Field in Cartesian Coordinates')


%% Computing wave tangent angles

% Compute slope of tangent lines across curvilinear space
[dz_dx, dz_dy] = gradient(horizontal_lines, dx, dy);
dz_dx(horizontal_lines < wave_profile) = nan;
dz_dy(horizontal_lines < wave_profile) = nan;

angles = atan2(dz_dx, dz_dy);
angles(horizontal_lines < wave_profile) = nan;

clear dz_dx dz_dy

figure('color', 'white')
hold on
contourf(vertical_lines, horizontal_lines, angles * (180/pi), 20)
plot(x, wave_profile, 'color', 'red', 'linewidth', 3)
for i = 1:10:150
    plot(x, horizontal_lines(i,:), 'color', 'black', 'linewidth', 2)
end
hold off
axis equal
xlabel('x [mm]')
ylabel('y [mm]')
colormap('coolwarm')
colorbar
clim([-10, 10])
title('Wave Tangent Angles [degrees]')


%% Projecting velocities from cartesian onto curvilinear
% PASSIVE TRANSFORMATION

% Compute unit vectors
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

% Plot
close all
figure('color', 'white')
tiledlayout(1,2)
nexttile
hold on
contourf(vertical_lines, horizontal_lines, u_tangent, 500, 'linestyle', 'none')
plot(x, max_wave_profile, 'LineWidth', 3, 'color', 'red')
plot(x, wave_profile, 'LineWidth', 3, 'color', 'black')
hold off
axis equal
title('Velocity Along ξ (Wave-following)')
colorbar

nexttile
hold on
contourf(vertical_lines, horizontal_lines, u_normal, 500, 'linestyle', 'none')
plot(x, max_wave_profile, 'LineWidth', 3, 'color', 'red')
plot(x, wave_profile, 'LineWidth', 3, 'color', 'black')
hold off
axis equal
title('Velocity Along ζ (Wave-normal)')
colorbar


step = 8;
scale = 0.5;
Xq = vertical_lines(1:step:end, 1:step:end);
Yq = horizontal_lines(1:step:end, 1:step:end);
theta = angles(1:step:end, 1:step:end);

figure;
hold on
contourf(vertical_lines, horizontal_lines, u_normal, 50, 'LineColor', 'none');
plot(x, max_wave_profile, 'LineWidth', 3, 'color', 'red')
plot(x, wave_profile, 'LineWidth', 3, 'color', 'black')
colormap(parula);
colorbar;

% Plot curves
for i = 1:step:rows
    plot(x, horizontal_lines(i,:), 'color', 'black')
end
for i = 1:step:columns
    plot(vertical_lines(:,i), y + wave_profile(i), 'color', 'black')
end

% Plot ξ̂ as red quivers
quiver(Xq, Yq, xi_hat_x(1:step:end, 1:step:end), xi_hat_y(1:step:end, 1:step:end), scale, 'r', 'LineWidth', 1.5);

% Plot ζ̂ as blue quivers
quiver(Xq, Yq, zeta_hat_x(1:step:end, 1:step:end), zeta_hat_y(1:step:end, 1:step:end), scale, 'b', 'LineWidth', 1.5);

title('Curvilinear Unit Vectors ξ̂ and ζ̂');
axis equal tight


%% Interpolating stresses

% Interpolate into new grid
% Fill NaNs in original data just in case
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


% Plot to check
figure('color', 'white')
tiledlayout(1,3)

nexttile
hold on
contourf(vertical_lines, horizontal_lines, uu_interp, 100, 'linestyle', 'none');
plot(x, max_wave_profile, 'LineWidth', 3, 'color', 'red')
plot(x, wave_profile, 'LineWidth', 3, 'color', 'black')
hold off
axis equal
colorbar()
title('Interpolated cartesian uu stress')

nexttile
hold on
contourf(vertical_lines, horizontal_lines, vv_interp, 100, 'linestyle', 'none');
plot(x, max_wave_profile, 'LineWidth', 3, 'color', 'red')
plot(x, wave_profile, 'LineWidth', 3, 'color', 'black')
hold off
axis equal
colorbar()
title('Interpolated cartesian vv stress')

nexttile
hold on
contourf(vertical_lines, horizontal_lines, uv_interp, 100, 'linestyle', 'none');
plot(x, max_wave_profile, 'LineWidth', 3, 'color', 'red')
plot(x, wave_profile, 'LineWidth', 3, 'color', 'black')
hold off
axis equal
colorbar()
title('Interpolated cartesian uv stress')


%% Projecting stresses

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

% Plot to show
figure('color', 'white')
tiledlayout(1,3)

nexttile
hold on
contourf(vertical_lines, horizontal_lines, tau_xixi, 100, 'linestyle', 'none');
plot(x, max_wave_profile, 'LineWidth', 3, 'color', 'red')
plot(x, wave_profile, 'LineWidth', 3, 'color', 'black')
hold off
axis equal
colorbar()
title('Wave Tangent, Normal Stress')

nexttile
hold on
contourf(vertical_lines, horizontal_lines, tau_zetazeta, 100, 'linestyle', 'none');
plot(x, max_wave_profile, 'LineWidth', 3, 'color', 'red')
plot(x, wave_profile, 'LineWidth', 3, 'color', 'black')
hold off
axis equal
colorbar()
title('Wave Normal, Normal Stress')


nexttile
hold on
contourf(vertical_lines, horizontal_lines, tau_xizeta, 100, 'linestyle', 'none');
plot(x, max_wave_profile, 'LineWidth', 3, 'color', 'red')
plot(x, wave_profile, 'LineWidth', 3, 'color', 'black')
hold off
axis equal
colorbar()
title('Curvilinear Shear Stress')

%%

lw = 2;
spacing = 10;

figure('color', 'white')
hold on
contourf(vertical_lines, horizontal_lines, tau_xizeta, 100, 'linestyle', 'none')
plot(x, max_wave_profile, 'LineWidth', 3, 'color', 'red')
plot(x, wave_profile, 'LineWidth', 3, 'color', 'blue')
for i = 1:spacing:rows
    plot(x, horizontal_lines(i,:), 'color', 'black', 'linewidth', lw)
end
for i = 1:spacing:columns
    plot(vertical_lines(:,i), y + wave_profile(i), 'color', 'black', 'linewidth', lw)
end
hold off
axis equal
xlim([min(x), max(x)])
ylim([min(y), top_crop])


