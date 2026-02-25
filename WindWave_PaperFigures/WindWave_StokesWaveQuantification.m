% Stokes wave quantification

clc; close all; clear
caze = 'WT6_WVB_AG0';
means = load(strcat('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/means/', caze, '_MEANS.mat'));
waves = load(strcat('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/wave/', caze, '_WAVE.mat'));
phases = load(strcat('/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/phase/', caze, '_PHASE.mat'));

means  = means.output;
waves = waves.output;
phases = phases.output;

wave_type = caze(7);
wave_parameters = readcell("Offshore_Waves.xlsx");
wavelength      = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 2};
amplitude       = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 3};
wavenumber      = (2 * pi) / wavelength;

clc; fprintf('Data loaded\n\n')

%% Save a subset of frames and waves

num_frames = 6000;

clear output
X = means.X;
x = X(1,:);
xc = x - mean(x, 'all');

for f = 1:num_frames
    output.waves(f, :) = imresize(waves.wave_profiles(f,:), [1, size(means.X, 2)]);
    output.phase_fits(f) = waves.fitted_phases(f);
    output.phases(f) = phases.phase_average_idx(f);
    clear f
end

clear num_frames



%% 'Phase average' the wave profiles

phase_shifts = [0, wavelength/4, wavelength/2, 3*wavelength/4];
single_cosine = @(b,x) amplitude * cos(wavenumber * (x - b(1)));

phase = 4;
phase_shift = phase_shifts(phase);
phase_wave_profiles = output.waves((output.phases == phase),:);
averaged_wave_profile = mean(phase_wave_profiles, 1, 'omitnan');

% figure('color', 'white')
% hold on
% % Individual wave profiles
% for i = 1:10:size(phase_wave_profiles, 1)
%     P = plot(xc, phase_wave_profiles(i,:), 'linewidth', 0.5, 'Color', 'black');
%     P.Color(:,4) = 0.25;
% end
% % Averaged wave profile
% plot(xc, averaged_wave_profile, 'LineWidth', 2, 'Color', 'blue')
% % Cosine
% plot(xc, single_cosine(phase_shift, xc), 'linewidth', 2, 'color', 'red');
% hold off
% 
% axis equal
% xlim([-120, 120])
% ylim([-50, 50])
% xticks(-100:50:100)
% yticks(-50:25:50)
% xlabel('x [mm]')
% ylabel('y [mm]')
% 
% 
% figure('color', 'white')
% plot(xc, averaged_wave_profile - single_cosine(phase_shift, xc))
% axis equal


% %% GPT Linear least-squares
% clc;
% y = averaged_wave_profile.';
% xc = xc(:);
% 
% % Build design matrix for cos(k(x-phi)), cos(2k(x-phi)), cos(3k(x-phi))
% X = [ ...
%     cos(1*wavenumber*(xc - phase_shift)), ...
%     cos(2*wavenumber*(xc - phase_shift)), ...
%     cos(3*wavenumber*(xc - phase_shift))  ];
% 
% % Solve LS: minimize ||X*b - y||_2
% b_hat = X \ y;
% 
% % Fitted stokes series (your function, now using best-fit amplitudes)
% stokes = @(b,xq) (b(1) * cos(1*wavenumber*(xq - phase_shift))) + ...
%                  (b(2) * cos(2*wavenumber*(xq - phase_shift))) + ...
%                  (b(3) * cos(3*wavenumber*(xq - phase_shift)));
% 
% y_fit = stokes(b_hat, x);
% 
% % Diagnostics
% res = y - y_fit;
% rmse = sqrt(mean(res.^2));
% fprintf('b1=%.2fmm, b2=%.2fmm, b3=%.2fmm \n\n', b_hat(1), b_hat(2), b_hat(3));
% % fprintf('RMSE = %.4f (units of y)\n', rmse);
% fprintf('A2/A1 = %.4f\nA3/A1 = %.3f\n\n', b_hat(2)/b_hat(1), b_hat(3)/b_hat(1));
% fprintf('sqrt{A2^2 + A3^2}/A1 = %.3f\n\n', sqrt(b_hat(2)^2 + b_hat(3)^2) / b_hat(1))
% 
% 
% figure('color','white'); 
% title(sprintf('Phase = %1.0f', phase))
% % tiledlayout(2,1)
% % h(1) = nexttile;
% hold on
% plot(xc, y, 'Linewidth', 2, 'color', 'black')
% plot(xc, single_cosine(phase_shift, xc), 'linewidth', 2, 'color', 'red', 'linestyle', '--');
% plot(xc, stokes(b_hat, xc), 'linewidth', 2, 'color', 'blue', 'linestyle', ':');
% hold off
% axis equal
% xlim([-120, 120]); ylim([-50, 50])
% xticks(-100:50:100); yticks(-50:25:50)
% xlabel('x [mm]'); ylabel('y [mm]')
% legend('Phase Averaged Wave Profile','Cosine Wave','LS 1–3 Harmonics','location','north', 'box', 'off')


%% Looped over all 4 phases

figure('color','white'); 
tiledlayout(4,1)

clc;
for phase = 1:4
    phase_shift = phase_shifts(phase);
    phase_wave_profiles = output.waves((output.phases == phase),:);
    averaged_wave_profile = mean(phase_wave_profiles, 1, 'omitnan');

    y = averaged_wave_profile.';
    xc = xc(:);
    
    % Build design matrix for cos(k(x-phi)), cos(2k(x-phi)), cos(3k(x-phi))
    X = [ ...
        cos(1*wavenumber*(xc - phase_shift)), ...
        cos(2*wavenumber*(xc - phase_shift)), ...
        cos(3*wavenumber*(xc - phase_shift))  ];
    
    % Solve LS: minimize ||X*b - y||_2
    b_hat = X \ y;
    
    % Fitted stokes series (your function, now using best-fit amplitudes)
    stokes = @(b,xq) (b(1) * cos(1*wavenumber*(xq - phase_shift))) + ...
                     (b(2) * cos(2*wavenumber*(xq - phase_shift))) + ...
                     (b(3) * cos(3*wavenumber*(xq - phase_shift)));
    
    y_fit = stokes(b_hat, xc);
    
    % Diagnostics
    % rmse = sqrt(mean(res.^2));
    nonlinearity = sqrt(b_hat(2)^2 + b_hat(3)^2) / b_hat(1);
    tmp(phase) = nonlinearity;

    % fprintf('b1=%.2fmm, b2=%.2fmm, b3=%.2fmm \n\n', b_hat(1), b_hat(2), b_hat(3));
    % % fprintf('RMSE = %.4f (units of y)\n', rmse);
    % fprintf('A2/A1 = %.4f\nA3/A1 = %.3f\n\n', b_hat(2)/b_hat(1), b_hat(3)/b_hat(1));
    % fprintf('sqrt{A2^2 + A3^2}/A1 = %.3f\n\n', nonlinearity)

    
    first_order_res = y - single_cosine(phase_shift, xc);
    stokes_res = y - y_fit;

    % Normalized first-order RMS
    first_order_RMS = rms(first_order_res);
    stokes_RMS = rms(stokes_res);
    wave_profile_rms = rms(y);

    normalized_first_order_rms = first_order_RMS / wave_profile_rms;
    normalized_stokes_RMS = stokes_RMS / wave_profile_rms;

    normalized_first_order_RMS_per_phase(phase) = normalized_first_order_rms;
    normalized_stokes_RMS_per_phase(phase) = normalized_stokes_RMS;

    % Compute R^2 per phase first
    R2_first_order(phase) = 1 - (first_order_RMS / wave_profile_rms);
    R2_stokes(phase) = 1 - (stokes_RMS / wave_profile_rms);


    h(phase) = nexttile;
    title(sprintf('Phase = %1.0f', phase))
    hold on
    plot(xc, y, 'Linewidth', 2, 'color', 'black')
    plot(xc, single_cosine(phase_shift, xc), 'linewidth', 2, 'color', 'red', 'linestyle', '--');
    plot(xc, stokes(b_hat, xc), 'linewidth', 2, 'color', 'blue', 'linestyle', ':');
    hold off
    axis equal
    xlim([-120, 120]); ylim([-25, 25])
    xticks(-100:50:100); yticks(-50:25:50)
    xlabel('x [mm]'); ylabel('y [mm]')
    if phase == 1
        legend('Phase Averaged Wave Profile','Cosine Wave','LS 1–3 Harmonics','location','northeastoutside', 'box', 'off')
    end
end

avg_nonlinearity = mean(tmp, 'all');
avg_first_order_RMS = mean(normalized_first_order_RMS_per_phase, 'all');
avg_stokes_RMS = mean(normalized_stokes_RMS_per_phase, 'all');
% fprintf('Avg nonlinearity = %1.4f\n\n', avg_nonlinearity)


first_order_RS = 1 - avg_first_order_RMS;
stokes_RS = 1 - avg_stokes_RMS;

fprintf('%s\n\n', caze)
fprintf('Average RMS First:\n')
fprintf('First Order R^2 = %1.3f\n', first_order_RS)
fprintf('Stokes R^2 = %1.3f\n\n', stokes_RS)


avg_R2_first_order = mean(R2_first_order, 'all');
avg_R2_stokes = mean(R2_stokes, 'all');

fprintf('Average R2:\n')
fprintf('First Order R^2 = %1.3f\n', avg_R2_first_order)
fprintf('Stokes R^2 = %1.3f\n', avg_R2_stokes)
