%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Functions/');
experiment_log  = readcell('Offshore_Inflow_Log.xlsx');
wave_parameters = readcell('Offshore_Waves.xlsx');

% Data paths
project_path   = '/Volumes/WT4_WVD_AG0/WT4_WVD_AG0';
processing     = '/PIV_MP(2x24x24_50%ov)/PostProc';
recording_name = 'WT4_WVD_AG0';
still_name     = strcat(experiment_log{find(strcmp(experiment_log, recording_name) == 1), 2}, '_waterlevel');
% inpt_name      = recording_name;
inpt_name      = strcat(recording_name, '_4');

% Experiment Specifics
tunnel_freq    = recording_name(strfind(recording_name, 'WT') + 2);
wave_type      = recording_name(strfind(recording_name, 'WV') + 2);
% active_grid    = recording_name(strfind(recording_name, 'AG') + 2: strfind(recording_name, 'AG') + 3);

% Image paths
still_path     = strcat('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Offshore/Water_Level/', still_name, '/Correction');
raw_path       = strcat(project_path, '/', recording_name, '/Correction');
piv_path       = strcat(project_path, '/', recording_name, processing);

% Save paths
results_path   = '/Users/zeinsadek/Desktop/Experiments/Offshore/Hopkins/Initial_Results/';
mtlb_file      = strcat(results_path, 'data'   , '/', inpt_name, '_DATA.mat');
mean_file      = strcat(results_path, 'means'  , '/', inpt_name, '_MEANS.mat');
phase_file     = strcat(results_path, 'phase'  , '/', inpt_name, '_PHASE.mat');
figure_file    = strcat(results_path, 'figures', '/', recording_name);

% Edge Detection Constants
std_tol            = 2.5;
background_removal = 1000;
canny_lower        = 0.1;
canny_upper        = 0.4;
grad_tol           = 0.4;
nan_dist           = 5;

% Filtering Parameters
left_bound_value   = -115;
right_bound_value  = 115;
phase_tolerance    = 2;
wave_buffer        = 1;
cmap_perc          = 0.05;

% Wave Parameters
if wave_type == '0'
    wave_length    = 0;
    wave_amplitude = 0;
else
    wave_length    = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 2};
    wave_amplitude = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 3};
    wave_frequency = wave_parameters{find(strcmp(wave_parameters, wave_type) == 1), 4};
end


phase_offset       = [0, wave_length/4, wave_length/2, 3*wave_length/4];
% cos_fit            = @(b, v) wave_amplitude * cos(2 * pi * (v - b(1)) / wave_length) + vertical_offset;

% Constants Structure
constants.background_removal = background_removal;
constants.std_tol            = std_tol;
constants.canny_lower        = canny_lower;
constants.canny_upper        = canny_upper;
constants.grad_tol           = grad_tol;
constants.nan_dist           = nan_dist;

constants.left_bound_value   = left_bound_value;
constants.right_bound_value  = right_bound_value;
constants.phase_tolerance    = phase_tolerance; 
constants.phase_offset       = phase_offset;

constants.wave_type          = wave_type;
constants.wave_amplitude     = wave_amplitude;
constants.wave_length        = wave_length;
constants.vertical_offset    = waterlevel(still_path, piv_path, still_name, constants); 
constants.wave_buffer        = wave_buffer;

vertical_offset    = constants.vertical_offset;
cos_fit            = @(b, v) wave_amplitude * cos(2 * pi * (v - b(1)) / wave_length) + vertical_offset;


% Image Path
image_name         = dir([raw_path, '/*.im7']);
D                  = length(image_name);

% Define Image Dimensions
raw                = readimx([raw_path, '\', image_name(1).name]);
raw_image          = raw.Frames{1,1}.Components{1,1}.Planes{1,1};
nf                 = size(raw_image);
x                  = raw.Frames{1,1}.Scales.X.Slope.*linspace(1, nf(1), nf(1)).* raw.Frames{1,1}.Grids.X + raw.Frames{1,1}.Scales.X.Offset;
y                  = raw.Frames{1,1}.Scales.Y.Slope.*linspace(1, nf(2), nf(2)).* raw.Frames{1,1}.Grids.Y + raw.Frames{1,1}.Scales.Y.Offset;
[~, left_bound]    = min(abs(x - left_bound_value));
[~, right_bound]   = min(abs(x - right_bound_value));
x                  = x(left_bound:right_bound);

%%% Check how many planes needed
x_range            = abs(left_bound_value) + abs(right_bound_value);
num_planes         = ceil(wave_length / x_range);
planes_x           = zeros(num_planes, length(x));
extended_x         = linspace(num_planes * left_bound_value, num_planes * right_bound_value, num_planes * length(x));

for i = 1:num_planes
    planes_x(i, :) = extended_x(1 + (length(x) * (i - 1)):length(x) * i);
end

% Saves
phase_average_idx       = zeros(1,D);
phase_average_plane_idx = zeros(1,D);
wave_profiles           = zeros(D,length(x));

fprintf('<phaseaverage> PROGRESS: \n');
D = 500;
for frame_number = 1:D

    % Print Progress.
    progressbarText(frame_number/D);

    % Load Images
    raw       = readimx([raw_path, '\', image_name(frame_number).name]);
    raw_image = raw.Frames{1,1}.Components{1,1}.Planes{1,1};
    nf        = size(raw_image);
    x         = raw.Frames{1,1}.Scales.X.Slope.*linspace(1, nf(1), nf(1)).* raw.Frames{1,1}.Grids.X + raw.Frames{1,1}.Scales.X.Offset;
    y         = raw.Frames{1,1}.Scales.Y.Slope.*linspace(1, nf(2), nf(2)).* raw.Frames{1,1}.Grids.Y + raw.Frames{1,1}.Scales.Y.Offset;
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

    % Crop to avoid nans at ends
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

    % Phase Fit
    for plane = 1:num_planes
        planar_x = planes_x(plane, :);
        fit_nan_crop = ~isnan(wave_profile);
        fcn          = @(b) sum((cos_fit(b, planar_x(fit_nan_crop)) - wave_profile(fit_nan_crop)).^2);                             
        fitted_phase = fminsearch(fcn, 0);
        
%         
%         figure(frame_number)
%         hold on
%         plot(planar_x, wave_profile)
%         plot(planar_x, cos_fit(fitted_phase, planar_x))
%         hold off
%         pause(2)

        % Avoid extraneous nans in profiles
        if sum(~isnan(wave_profile)) < length(x)
            wave_profiles(frame_number, :) = cos_fit(fitted_phase, planar_x);
        else    
            wave_profiles(frame_number, :) = wave_profile;
        end

        % Convert negative phases to positive
        if fitted_phase < -1 * phase_tolerance
            fitted_phase = wave_length + fitted_phase;
        end

%         disp(fitted_phase)
        
        % Bin frames to phases
        for phase = 1:length(phase_offset)   
            lower_phase = phase_offset(phase) - phase_tolerance;
            upper_phase = phase_offset(phase) + phase_tolerance;
            if fitted_phase >= lower_phase && fitted_phase <= upper_phase
                phase_average_idx(frame_number)       = phase;
                phase_average_plane_idx(frame_number) = plane;
            end
        end
    end

%%%%% This Part Works
%             % Phase Fit
%             fit_nan_crop = ~isnan(wave_profile);
%             fcn          = @(b) sum((cos_fit(b, x(fit_nan_crop)) - wave_profile(fit_nan_crop)).^2);                             
%             fitted_phase = fminsearch(fcn, 0);
% 
%             % Avoid extraneous nans in profiles
%             if sum(~isnan(wave_profile)) < length(x)
%                 wave_profiles(frame_number, :) = cos_fit(fitted_phase, x);
%             else    
%                 wave_profiles(frame_number, :) = wave_profile;
%             end
% 
%             % Convert negative phases to positive
%             if fitted_phase < -1 * phase_tolerance
%                 fitted_phase = wave_length + fitted_phase;
%             else
%             end
% 
%             % Bin frames to phases
%             for phase = 1:length(phase_offset)   
%                 lower_phase = phase_offset(phase) - phase_tolerance;
%                 upper_phase = phase_offset(phase) + phase_tolerance;
% 
%                 if fitted_phase >= lower_phase && fitted_phase <= upper_phase
%                     phase_average_idx(frame_number) = phase;
%                 else
%                 end
%             end
%         end
%%%%%
end
 
%%

%%% Plot
fprintf('\nTolerance Range: %4.2f mm\n', 2 * phase_tolerance)

figure()
hold on 
for phase = 1:length(phase_offset)
% for phase = 1:1
    reference_profile = wave_amplitude * cos(2 * pi * (extended_x - phase_offset(phase)) / wave_length);
    phase_idx = find(phase_average_idx == phase);

    subplot(length(phase_offset), 1, phase)
    
    ref = plot(extended_x, reference_profile, 'Color', 'red', 'LineWidth', 2);
    
    
    for plane = 1:num_planes
        planar_x = planes_x(plane, :);

        % Of the phase sorted images identify which plane they belong to
        idx = find(phase_average_plane_idx(phase_idx) == plane);
        fprintf('\nPhase %.f Plane %.f: %.f Images\n', phase, plane, length(idx))

        hold on
        for f = 1:length(idx)
            frame      = phase_idx(f);
            p          = plot(planar_x, wave_profiles(frame, :) - (vertical_offset - 5), 'Color', 'black', 'HandleVisibility','off');
            p.Color(4) = 0.25;
        end
        hold off
%                 plot(planar_x, mean(wave_profiles(idx, :)), 'Color', 'green', 'LineWidth', 2, 'DisplayName', 'Average')
%                 plot(planar_x, max(wave_profiles(idx,:)), 'Color', 'magenta', 'LineWidth', 2, 'DisplayName', 'Max') 

    title(strcat(num2str(length(idx)), {' '}, 'Images'))
    uistack(ref, 'top')
    end

    hold off
    axis equal
    legend()
    xlim([min(x), max(x)])
%     ylim([min(y), -70])
    ylabel('z [mm]')
%     title(phase)
end




















