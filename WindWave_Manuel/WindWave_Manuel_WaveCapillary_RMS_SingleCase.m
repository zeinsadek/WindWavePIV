% Try to seperate capillary action on the surface of the forced waves and
% determine the RMS of there surface ripples

clc; close all; clear

% Get case to process
wind_speed = 'WT8';
wave_type = 'A';
caze = [wind_speed, '_WV', wave_type, '_AG0'];

% Load all useful data
waves  = load(strcat('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/wave/', caze, '_WAVE.mat'));
means  = load(strcat('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means/', caze, '_MEANS.mat'));
phases = load(strcat('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/phase/', caze, '_PHASE.mat'));

% Rename structure
waves  = waves.output;
means  = means.output;
phases = phases.output;
clc; fprintf('Data loaded\n\n')

% Get wave parameters
wave_parameters = readcell("/Users/zeinsadek/Desktop/Experiments/PIV/Processing/WindWave/Offshore_Waves.xlsx");
wavelength      = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 2}; %#ok<*FNDSB>
amplitude       = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 3};
frequency       = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 4};
phase_offset    = [0, wavelength/4, wavelength/2, 3*wavelength/4];
wavenumber      = (2 * pi) / wavelength;

% Get freestream
freestreams.('WT4') = 2.4181;
freestreams.('WT6') = 3.8709;
freestreams.('WT8') = 5.4289;

freestream = freestreams.(wind_speed);

clear wave_type


%% Load relevent data

% Coordinates
x = waves.x;

% Bin instantaneous snapshots by phase
clc;
for phase = 1:4

    % Indicies of images (out of the 6000 frames) that belong to 'phase'
    indicies = find(phases.phase_average_idx == phase);
    fprintf('Phase %1.0f: %4.0f Images\n\n', phase, length(indicies))

    % Save wave profiles
    PhaseBinnedWaves(phase).waves = waves.wave_profiles(indicies,:); %#ok<*SAGROW>

    % Save fitted phase shifts
    PhaseBinnedWaves(phase).shifts = waves.fitted_phases(indicies);

    % Save reference wave profiles
    PhaseBinnedWaves(phase).reference = imresize(means.phase(phase).reference_wave.', [1, length(x)]);

end

clear phases wave_phases_indicies phase indicies

% Compute phase-averaged wave surface
for phase = 1:4
    PhaseBinnedWaves(phase).phase_averaged_waves = mean(PhaseBinnedWaves(phase).waves, 1, 'omitnan');
end

clear phase


%% Plot one profile with corresponding stokes reference and compute RMS

% Which phase and which image
phase = 1;
image = 1;

% Line width
lw = 2;

% High-pass filter cutoff wavelength [mm]
% cutoff [mm] — set below forced wavelength
lambda_c = 30;                



% Phase averages surface profile
phase_average_wave = PhaseBinnedWaves(phase).phase_averaged_waves;

% Measured profile and its fitted phase shift
eta_measured = PhaseBinnedWaves(phase).waves(image, :);
phi_snap = PhaseBinnedWaves(phase).shifts(image);   % fitted phase [mm]

% Reference profile and its nominal phase shift (phase-bin center)
phi_ref = (phase - 1) * wavelength/4;

% Shift reference onto the snapshot's phase
delta_phi    = phi_snap - phi_ref;
eta_ref_shifted = interp1(x, phase_average_wave, x - delta_phi, 'spline', nan);

% Find valid (non-NaN) region after shift
valid = ~isnan(eta_ref_shifted);
x_crop          = x(valid);
eta_meas_crop   = eta_measured(valid);
eta_ref_crop    = eta_ref_shifted(valid);

% Residual on cropped domain
residual = eta_meas_crop - eta_ref_crop;



% High-pass the residual in wavenumber space
dx       = mean(diff(x));
N        = length(residual);
k_fft    = [0:floor(N/2), -floor(N/2)+1:-1] / (N * dx);  % two-sided
R_fft    = fft(residual);
R_fft(abs(1./abs(k_fft + eps)) > lambda_c) = 0;
residual_hp = real(ifft(R_fft));



% Plots
figure('color', 'white')
tile = tiledlayout(2,1,'tilespacing','tight');
sgtitle(caze, 'interpreter', 'none')

% Profiles
h(1) = nexttile;
hold on
plot(x, eta_measured, 'linewidth', 1, 'linewidth', lw, 'DisplayName', '$\eta_{i}$')
plot(x, phase_average_wave, 'linewidth', lw, 'DisplayName', '$\langle \eta \rangle$')
plot(x_crop, eta_ref_crop, 'linewidth', lw, 'DisplayName', '$\langle \eta \rangle + \phi$')
hold off
axis equal
ylabel('$y$ [mm]', 'interpreter', 'latex')
legend('interpreter', 'latex', 'location', 'eastoutside', 'box', 'off', 'fontsize', 14)

% Residuals
h(2) = nexttile;
hold on
plot(x_crop, residual, 'linewidth', lw, 'displayname', '$\eta_{i} - \langle \eta \rangle$')
plot(x_crop, residual_hp, 'linewidth', lw, 'displayname', 'High-Pass: $\eta_{i} - \langle \eta \rangle$')
hold off
legend('interpreter', 'latex', 'location', 'eastoutside', 'box', 'off', 'fontsize', 14)
ylim([-5,5])

axis equal
linkaxes(h, 'x')
xlim([-120, 120])
ylim([-20, 20])
xlabel('$x$ [mm]', 'interpreter', 'latex')
ylabel('$y$ [mm]', 'interpreter', 'latex')


% Print our RMS from different differences
clc;
fprintf('Full Difference RMS: %1.3f mm\n', rms(residual))
fprintf('High-Pass Difference RMS: %1.3f mm\n', rms(residual_hp))



%% Loop through phase and instantaneous images to compute average RMS


% High-pass filter cutoff wavelength [mm]
% cutoff [mm] — set below forced wavelength
lambda_c = 30;                

% FFT sizes
dx = mean(diff(x));

% Save values
all_residuals = [];
hp_residuals = [];

clc;
% Loop through phases
for phase = 1:4

    fprintf('Phase %1.0f\n', phase)
    % Phase averages surface profile
    phase_average_wave = PhaseBinnedWaves(phase).phase_averaged_waves;

    % Reference profile and its nominal phase shift (phase-bin center)
    phi_ref = (phase - 1) * wavelength/4;

    % Loop through images
    for image = 1:size(PhaseBinnedWaves(phase).waves, 1)

        % Measured profile and its fitted phase shift
        eta_measured = PhaseBinnedWaves(phase).waves(image, :);
        phi_snap = PhaseBinnedWaves(phase).shifts(image);

        % Shift reference onto the snapshot's phase
        delta_phi    = phi_snap - phi_ref;
        eta_ref_shifted = interp1(x, phase_average_wave, x - delta_phi, 'spline', nan);
        
        % Find valid (non-NaN) region after shift
        valid = ~isnan(eta_ref_shifted);
        x_crop          = x(valid);
        eta_meas_crop   = eta_measured(valid);
        eta_ref_crop    = eta_ref_shifted(valid);
        
        % Residual on cropped domain
        residual = eta_meas_crop - eta_ref_crop;



        % High-pass the residual in wavenumber space
        N_crop  = length(residual);
        k_fft   = [0:floor(N_crop/2), -floor(N_crop/2)+1:-1] / (N_crop * dx);
        R_fft    = fft(residual);
        R_fft(abs(1./abs(k_fft + eps)) > lambda_c) = 0;
        residual_hp = real(ifft(R_fft));



        % Save
        all_residuals = [all_residuals, residual];
        hp_residuals = [hp_residuals, residual_hp];

    end
end
clc;

% Compute RMS from residuals across all images
all_RMS = rms(all_residuals, 'omitnan');
hp_RMS = rms(hp_residuals, 'omitnan');

% Print our RMS from different differences
clc;
fprintf('Full Difference RMS: %1.3f mm\n', all_RMS)
fprintf('High-Pass Difference RMS: %1.3f mm\n\n', hp_RMS)







%% Stokes shii
% --- Compute 3rd-order Stokes wave at the DATA x-points ---
% % First-order: O(1)
% term1  = a0 * cos(theta);
% 
% % Second-order: O(epsilon)
% term2  = epsilon * (a0/4) * chi * (3*chi^2 - 1) * cos(2*theta);
% 
% % Third-order: O(epsilon^2)
% term3a = -epsilon^2 * a0 * (3/8) * (chi^4 - 3*chi^2 + 3) * cos(theta);
% term3b =  epsilon^2 * a0 * (3/64) * (8*chi^6 + (chi^2 - 1)^2) * cos(3*theta);
% 
% % Third-order stokes wave
% eta_stokes_mm = (term1 + term2 + term3a + term3b);
