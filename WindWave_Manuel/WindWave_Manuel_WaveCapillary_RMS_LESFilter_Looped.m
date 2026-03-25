% Wrap up the curvilinear-grid Cartesian velocities instantaneous 
% and phase averages for Cerrina

clc; close all; clear

% Freestream velocities
freestreams.('WT4') = 2.4181;
freestreams.('WT6') = 3.8709;
freestreams.('WT8') = 5.4289;

% Approximate wavelengths in mm for labeling plots
wavelengths.A = '410';
wavelengths.B = '313';
wavelengths.C = '189';
wavelengths.D = '124';

steepnesses.A = '0.180';
steepnesses.B = '0.211';
steepnesses.C = '0.305';
steepnesses.D = '0.267';

 
% Cutoff frequency specified by JHU LES
LES_cutoff.A = 38.5125;
LES_cutoff.B = 29.3719;
LES_cutoff.C = 35.5688;
LES_cutoff.D = 34.9594;


% Loop through cases
wave_types = {'A', 'B', 'C', 'D'};
wind_speeds = {'WT4', 'WT6', 'WT8'};

for w = 1:length(wave_types)
    wave_type = wave_types{w};

    % High-pass filter cutoff wavelength [mm]
    % cutoff [mm] — set below forced wavelength
    lambda_c = LES_cutoff.(wave_type);  

    for s = 1:length(wind_speeds)
        wind_speed = wind_speeds{s};
        caze = [wind_speed, '_WV', wave_type, '_AG0'];
        disp(caze)

        % Get freestream
        freestream = freestreams.(wind_speed);

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



        %%% Load relevent data
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
        
        
        
        %%% Loop through phase and instantaneous images to compute average RMS
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

        

        % Save to output file
        output.(caze).parameters.wavelength_mm = wavelength;
        output.(caze).parameters.amplitude_mm = amplitude;
        output.(caze).parameters.steepness = (2 * pi * amplitude) / wavelength;
        output.(caze).parameters.frequency_hz = frequency;
        output.(caze).parameters.freestream_ms = freestream;
        output.(caze).parameters.high_pass_cutoff_wavelength_mm = lambda_c;
        output.(caze).parameters.wave_speed_ms = frequency * wavelength * 1E-3;
        output.(caze).parameters.wave_age = (frequency * wavelength * 1E-3) / freestream;

        % RMS in mm
        output.(caze).RMS_mm.all = all_RMS;
        output.(caze).RMS_mm.high_pass = hp_RMS;

    end
end



%% Plot to see if things change with wavelength

wind_speed_markers = {'^', 'square', 'o'};
wave_colors = {'#FE6202', '#DC2680', '#775EEF', '#648FFF'};

figure('color', 'white')
hold on
for w = 1:length(wave_types)
    wave_type = wave_types{w};

    for s = 1:length(wind_speeds)
        wind_speed = wind_speeds{s};
        caze = [wind_speed, '_WV', wave_type, '_AG0'];
        disp(caze)

        frequency = output.(caze).parameters.frequency_hz;
        wavelength = output.(caze).parameters.wavelength_mm;
        freestream = output.(caze).parameters.freestream_ms;
        steepness = output.(caze).parameters.steepness;

        label = sprintf('$\\lambda_{%s}, \\hspace{1mm} ak_{%s}$', wavelengths.(wave_type), steepnesses.(wave_type));
        scatter((frequency * wavelength * 1E-3) / freestream, output.(caze).RMS_mm.high_pass, ...
                50, wind_speed_markers{s}, 'filled', 'MarkerFaceColor', wave_colors{w}, ...
                'HandleVisibility', 'off')

    end

    hLeg = plot(nan, nan, 'o', ...
        'MarkerFaceColor',  wave_colors{w}, ...
        'MarkerEdgeColor','none', ...
        'MarkerSize', 4, ...        % <-- CONTROL LEGEND SIZE HERE
        'LineWidth', 1, ...
        'LineStyle','none', ...
        'DisplayName',label);
end
% Add legend for markers
plot(nan, nan, 'color', 'white', 'displayname', '')

% Dummy lines to get wind-speed legend
for s = 1:length(wind_speeds)
    wind_speed = wind_speeds{s};

    if ismember(wind_speed(end), {'4'})
        u_inf = 2.4181;
    elseif ismember(wind_speed(end), {'6'})
        u_inf = 3.8709;
    elseif ismember(wind_speed(end), {'8'})
        u_inf = 5.4289;
    end

    hLeg = plot(nan, nan, wind_speed_markers{s}, ...
    'MarkerFaceColor','black', ...
    'MarkerEdgeColor','none', ...
    'MarkerSize', 4, ...        % <-- CONTROL LEGEND SIZE HERE
    'LineWidth', 1, ...
    'LineStyle','none', ...
    'DisplayName', sprintf('$u_{\\infty} = %1.2f$ m/s', u_inf));

end

leg = legend('interpreter', 'latex', 'location', 'eastoutside', 'box', 'off', 'fontsize', 10);
leg.IconColumnWidth = 19;

leg.ItemTokenSize(1) = 10;

hold off
axis square
xlim([0.05, 0.35])
xticks(0.1:0.1:0.3)
xlabel('$c \mathbin{/} u_{\infty}$', 'interpreter', 'latex', 'fontsize', 12)
ylabel('$RMS(\eta)$ [mm]', 'interpreter', 'latex', 'fontsize', 12)
hold off

%% Save to matfile

save_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/manuel';
file_name = 'WindWave_ForcedWaveCapillaryRMS_LESCutoff.mat';
fprintf('Saving matfile...\n')
pause(3)
save(fullfile(save_folder, file_name), 'output')
fprintf('Matfile saved!\n\n')


