%%% Curvilinear Coordinates
% Zein Sadek
% 11 / 2023

clc; clear; close all;

%%% Wave Parameters
% Amplitde
a = 10;

% Wavelength
l = 200;
k = 1/l;
steepness = (2*a)/l;
fprintf('Wave Steepness: %3.3f\n', steepness)

% Cartesian Coordinates
x = 0:1:500;
z = rot90(-10:1:1000);
[X, Z] = meshgrid(x, z);

% Wave profile / Normal wave profile
wave = a * cos((2 * pi * k) * x);
wave_ortho = ((2 * pi) / l) * a * sin((2 * pi * k) * x);

% Mask Z meshgrid
Z(Z < wave) = nan;

%%% Plotting Parameters
zeta_step = 10;
xi_step = zeta_step;

linewidth = 0.5;
zeta_color = 'red';
xi_color = 'blue';
normal_scale = 10;
normal_color = 'green';

figure()
hold on

% Wave surface
plot(x, wave, 'LineWidth', linewidth, 'color', 'black')

% Zeta Lines (Horizontal Lines)
for i = zeta_step:zeta_step:length(z(z > 0))
    zeta = i;
    plot(x, zeta + wave * exp(-k * zeta), 'LineWidth', linewidth, 'Color', zeta_color)
end

% Xi Lines (Vertical Lines)
for i = 1:xi_step:length(x)
    % Zero Z so that wave surface is zero
    zero_slice = Z(:,i) - min(Z(:,i), [], 'all', 'omitnan');

    % Lines of constant xi
    plot((zero_slice .* wave_ortho(i) .* exp(-k * zero_slice)) + x(i) * ones(1, length(zero_slice)), Z(:,i), 'color', xi_color, 'LineWidth', linewidth)

    % Wave normal vectors
    plot([x(i), x(i) + (normal_scale * wave_ortho(i))], [wave(i), wave(i) + (normal_scale * 1)], 'linewidth', linewidth, 'Color', normal_color)
end

hold off
axis equal
ylim([-20, max(z)])
xlim([min(x), max(x)])
xlabel('x [mm]')
ylabel('y [mm]')


