%% compare_stokes_phases.m
%  Reads digitized eta_max data at four wave phases and overlays
%  the corresponding phase of a 3rd-order Stokes wave for comparison.

clear; close all; clc;

%% ============== FILE / PHASE SETUP ==============
files = {
    'phase1.csv', ...   % phi = 0
    'phase2.csv', ...   % phi = pi/2
    'phase3.csv', ...   % phi = pi
    'phase4.csv'  ...   % phi = 3*pi/2
};

% Stokes phase that matches each dataset
% (phi=0 and phi=pi are swapped to align crests/troughs with data)
phi_stokes = [pi+pi/4, pi/2+pi/4, 0+pi/4, 3*pi/2+pi/4];

% Subplot titles match the original figure labels
phase_labels = {'\phi = 0', '\phi = \pi/2', ...
                '\phi = \pi', '\phi = 3\pi/2'};

%% ============== USER-ADJUSTABLE WAVE PARAMETERS ==============
lambda   = 189.6/1000;             % Wavelength [m]
a0       = 9.2/1000;               % Wave amplitude [m]
d_lambda = 1.27;                   % Depth-to-wavelength ratio
k        = 2*pi/lambda;            % Wavenumber [rad/m]
d        = d_lambda * lambda;
kd       = k * d;
epsilon  = a0 * k;                 % Steepness
chi      = 1 / tanh(kd);           % 1/tanh(kd)

fprintf('lambda = %.1f mm,  a0 = %.1f mm,  kd = %.3f,  eps = %.4f,  chi = %.4f\n', ...
    lambda*1e3, a0*1e3, kd, epsilon, chi);

%% ============== READ DATA & PLOT ==============
figure('Position', [100 100 1200 800]);

z0_avg = 0;
rms_error_avg = 0;

for idx = 1:4
    subplot(2,2,idx);

    % --- Read digitized data (units: mm) ---
    data        = readmatrix(files{idx});
    x_mm        = data(:,1);
    eta_data_mm = data(:,2);

    % --- Compute 3rd-order Stokes wave at the DATA x-points ---
    x_m   = x_mm * 1e-3;                       % convert to [m]
    phi   = phi_stokes(idx);
    theta = k * x_m - phi;

    % First-order: O(1)
    term1  = a0 * cos(theta);
    % Second-order: O(epsilon)
    term2  = epsilon * (a0/4) * chi * (3*chi^2 - 1) * cos(2*theta);
    % Third-order: O(epsilon^2)
    term3a = -epsilon^2 * a0 * (3/8) * (chi^4 - 3*chi^2 + 3) * cos(theta);
    term3b =  epsilon^2 * a0 * (3/64) * (8*chi^6 + (chi^2 - 1)^2) * cos(3*theta);

    eta_stokes_mm = (term1 + term2 + term3a + term3b) * 1e3;  % [mm]

    % --- Plot ---
    plot(x_mm, eta_data_mm, 'o', 'Color', [0.83 0.63 0.09], ...
        'MarkerSize', 4, 'LineWidth', 1.2, 'DisplayName', '\eta_{max} (data)');
    hold on;
    plot(x_mm, eta_stokes_mm, 'b-', 'LineWidth', 2, ...
        'DisplayName', '3rd-order Stokes');
    hold off;

    xlabel('x [mm]');
    ylabel('\eta [mm]');
    title(phase_labels{idx}, 'FontSize', 14);
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    xlim([0 240]);

    % --- RMS metrics (both evaluated at the same x-points) ---
    rms_data   = sqrt(mean(eta_data_mm.^2));
    rms_stokes = sqrt(mean(eta_stokes_mm.^2));

    eta_diff = eta_data_mm-eta_stokes_mm;
    rms_error = sqrt(mean(eta_diff.^2) - mean(eta_diff).^2);
    % rms_error  = sqrt(mean((eta_data_mm - eta_stokes_mm).^2));

    z0 = rms_error * exp(-0.4*8.5) / 1000;
    fprintf('Phase %d: RMS_diff = %.4e mm, z0 = %.4e m\n', idx, rms_error, z0);

    z0_avg = z0_avg + z0;
    rms_error_avg = rms_error_avg  +rms_error;
end

sgtitle('Digitized \eta_{max}  vs  3rd-Order Stokes Wave at Each Phase', 'FontSize', 16);

z0_avg = z0_avg / 4;
rms_error_avg = (rms_error_avg/4);

rms_OG = 7.8e-5*1000;   %mm
z0_OG  = rms_OG * exp(-0.4*8.5)/1000;

fprintf('\n Fig4: RMS = %.4e mm, z0 = %.4e m\n', rms_error_avg, z0_avg);
fprintf('OG: RMS = %.4e mm, z0 = %.4e m\n', rms_OG, z0_OG);


