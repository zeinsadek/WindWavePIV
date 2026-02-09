%% Wave Detection for Cropping PIV Images
% Zein Sadek, 1/23

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Functions/');

% Constantss
background_removal = 1000;
canny_lower        = 0.1;
canny_upper        = 0.4;
grad_tol           = 0.4;
nan_dist           = 5;
left_bound_value   = -120;
right_bound_value  = 120;
phase_tolerance    = 5;        % Units of [mm]

% Image Path
raw_path           = '/Volumes/Zein_PIV/Offshore_Inflow/wave_test/Correction';
image_name         = dir([raw_path, '/*.im7']);
D                  = length(image_name);

% Define Image Dimensions
raw                = readimx([raw_path, '\', image_name(1).name]);
raw_image          = raw.Frames{1,1}.Components{1,1}.Planes{1,1};
nf                 = size(raw_image);
x                  = raw.Frames{1,1}.Scales.X.Slope.*linspace(1, nf(1), nf(1)).* raw.Frames{1,1}.Grids.X + raw.Frames{1,1}.Scales.X.Offset;
y                  = raw.Frames{1,1}.Scales.Y.Slope.*linspace(1, nf(2), nf(2)).* raw.Frames{1,1}.Grids.Y + raw.Frames{1,1}.Scales.Y.Offset;
[X, Y]             = meshgrid(x, y);
[~, left_bound]    = min(abs(x - left_bound_value));
[~, right_bound]   = min(abs(x - right_bound_value));
x                  = x(left_bound:right_bound);

% Reference Wave
wave_amplitude     = 3.351;
wave_length        = 125.871;
% NEED TO TAKE STILL WATER IMAGE TO MEASURE
vertical_offset    = -103;
phase_offset       = [0, wave_length/4, wave_length/2, 3*wave_length/4];
cos_fit            = @(b, v) wave_amplitude * cos(2 * pi * (v - b(1)) / wave_length) + vertical_offset;

% Saves
phase_average_idx    = zeros(1,D);
wave_profiles        = zeros(D,length(x));
wave_profile_lengths = zeros(1,D);


%%% Bad Frame was 1691

fprintf('<phaseaverage> PROGRESS: \n');
for frame_number = 1:10
    
    % Print Progress.
    progressbarText(frame_number/D);
            
    % Load Images
    raw       = readimx([raw_path, '\', image_name(frame_number).name]);
    raw_image = raw.Frames{1,1}.Components{1,1}.Planes{1,1};
    nf        = size(raw_image);
    x         = raw.Frames{1,1}.Scales.X.Slope.*linspace(1, nf(1), nf(1)).* raw.Frames{1,1}.Grids.X + raw.Frames{1,1}.Scales.X.Offset;
    y         = raw.Frames{1,1}.Scales.Y.Slope.*linspace(1, nf(2), nf(2)).* raw.Frames{1,1}.Grids.Y + raw.Frames{1,1}.Scales.Y.Offset;
    [X, Y]    = meshgrid(x, y);
    raw_image = fliplr(rot90(raw_image, -1));

    % Edge Detection
    raw_image(raw_image < background_removal) = 0;
    wave_edge = edge(raw_image, 'Canny', [canny_lower, canny_upper]);
    wave_profile = zeros(1, nf(1));
    
    for i = 1:nf(1)
        vertical_slice = wave_edge(:, i);
        [ones_r, ~]    = find(vertical_slice == 1);
        size_ones      = size(ones_r);
        len            = size_ones(1);
        if len == 0
            wave_profile(1, i) = nan;
        else
            wave_profile(1, i) = y(min(ones_r));
        end
    end
    
    % Crop slightly to avoid nans at ends
    x            = x(left_bound:right_bound);
    wave_profile = wave_profile(left_bound:right_bound);

    % Replace kinks with nan values
    spikes = find(abs(gradient(wave_profile)) >= grad_tol);
    wave_profile(spikes) = nan;
    
    
    
    % Remove values between nans
    for i = 1:length(spikes)
        first_nan = spikes(i);
        for j = 1:nan_dist
            second_nan = first_nan + j;
            if second_nan > length(x)
                second_nan = length(x);
            end
            if isnan(wave_profile(second_nan))
                wave_profile(first_nan:second_nan) = nan;
            end
        end
    end
    
    
    
    
    % Interpolate over nans
    wave_profile(isnan(wave_profile)) = interp1(x(~isnan(wave_profile)), wave_profile(~isnan(wave_profile)), x(isnan(wave_profile)));
    
    figure()
%     hold on
    plot(x, wave_profile)
    
%     plot(x, wave_profile)
%     hold off
   
    % Phase Fit
    fit_nan_crop = ~isnan(wave_profile);
    fcn          = @(b) sum((cos_fit(b, x(fit_nan_crop)) - wave_profile(fit_nan_crop)).^2);                             
    fitted_phase = fminsearch(fcn, 0);
    
    % Avoid extraneous nans in profiles
    if sum(~isnan(wave_profile)) < length(x)
        wave_profiles(frame_number,:)  = cos_fit(fitted_phase, x);
    else    
        wave_profiles(frame_number,:)  = wave_profile;
    end
    wave_profile_lengths(frame_number) = sum(~isnan(wave_profiles(frame_number,:)));
    
    % Convert negative phases to positive
    if fitted_phase < -1 * phase_tolerance
        fitted_phase = wave_length + fitted_phase;
    else
    end
    
    % Bin frames to phases
    for phase = 1:length(phase_offset)
        % Set tolerance bounds    
        lower_phase = phase_offset(phase) - phase_tolerance;
        upper_phase = phase_offset(phase) + phase_tolerance;
        
        % Evaluate 
        if fitted_phase >= lower_phase && fitted_phase <= upper_phase
            phase_average_idx(frame_number) = phase;
        else
        end
    end
%     
    figure()
    hold on
    plot(x, wave_profile)
    plot(x, cos_fit(fitted_phase, x))
    hold off
    title(frame_number)
    axis equal
    xlim([left_bound_value, right_bound_value])
    
    pause(2)
    
    disp(fitted_phase)
    close all

end


%%

% figure(3)
% clc;
% fprintf('\nTolerance Range: %4.2f mm\n', 2 * phase_tolerance)
% 
% for phase = 1:length(phase_offset)
%     reference_profile = wave_amplitude * cos(2 * pi * (x - phase_offset(phase)) / wave_length) + vertical_offset;
%     idx = find(phase_average_idx == phase);
%     
%     fprintf('\nPhase %2u: %3f Images\n', phase, length(idx))
%     
%     subplot(4,1,phase)
%     hold on 
%     for f = 1:length(idx)
%         frame      = idx(f);
%         p          = plot(x, wave_profiles(frame, :), 'Color', 'blue', 'HandleVisibility','off');
%         p.Color(4) = 0.25;
%     end
%     plot(x, reference_profile, 'Color', 'red', 'LineWidth', 2, 'DisplayName', 'Reference')
%     plot(x, mean(wave_profiles(idx, :)), 'Color', 'green', 'LineWidth', 2, 'DisplayName', 'Average')
%     plot(x, max(wave_profiles(idx,:)), 'Color', 'magenta', 'LineWidth', 2, 'DisplayName', 'Max')
%     hold off
%     axis equal
%     legend()
%     xlim([min(x), max(x)])
%     ylim([min(y), -70])
%     title(phase)
% end
% 







