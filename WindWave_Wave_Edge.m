%% Wave Detection for Cropping PIV Images
% Zein Sadek, 1/23

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Functions/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Functions/Inpaint_nans/Inpaint_nans/')
addpath('/Users/zeinsadek/audio-orchestrator-ffmpeg/bin/')

% Edge Detection Constantss
background_removal = 1000;
canny_lower        = 0.1;
canny_upper        = 0.4;
grad_tol           = 2;
nan_dist           = 3;

% Image Path
raw_path         = '/Volumes/WT4_WVD_AG0/WT4_WVD_AG0/WT4_WVD_AG0/Correction';
piv_path         = '/Volumes/WT4_WVD_AG0/WT4_WVD_AG0/WT4_WVD_AG0/PIV_MP(2x24x24_50%ov)/PostProc';
raw_image_name   = dir([raw_path, '/*.im7']);
piv_image_name   = dir([piv_path, '/*.vc7']);
D                = length(raw_image_name);

% Detect bad crops by wave_profile length
wave_profile_length = zeros(1,D);

% Peaks 
max_wave_heights    = zeros(1,D);

% Save wave profiles
wave_profiles = zeros(D,171);
std_tol = 4;

video_name = 'WT4_WVD_AG0_streamwise_inst_movie_inpaint_nans.avi';
video_path = strcat('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Offshore/', video_name);


tst = [];
writerObj = VideoWriter(video_path, 'Uncompressed AVI');
writerObj.FrameRate = 3;


%%% Bad Frame was 1691

D = 20;
for frame_number = 1:D
    
    %%% Load Images
    raw = readimx([raw_path, '\', raw_image_name(frame_number).name]);
    piv = readimx([piv_path, '\', piv_image_name(frame_number).name]);
    % Raw Image
    raw_image = raw.Frames{1,1}.Components{1,1}.Planes{1,1};
    % PIV
    names       = piv.Frames{1,1}.ComponentNames;        
    U0_index    = find(strcmp(names, 'U0'));
    V0_index    = find(strcmp(names, 'V0'));
    UF          = piv.Frames{1,1}.Components{U0_index,1}.Scale.Slope.*piv.Frames{1,1}.Components{U0_index,1}.Planes{1,1} + piv.Frames{1,1}.Components{U0_index,1}.Scale.Offset;
    VF          = piv.Frames{1,1}.Components{V0_index,1}.Scale.Slope.*piv.Frames{1,1}.Components{V0_index,1}.Planes{1,1} + piv.Frames{1,1}.Components{V0_index,1}.Scale.Offset;
 
    
    % PIV Image Coordinates
    nf_piv = size(UF);
    x_piv  = piv.Frames{1,1}.Scales.X.Slope.*linspace(1, nf_piv(1), nf_piv(1)).*piv.Frames{1,1}.Grids.X + piv.Frames{1,1}.Scales.X.Offset;
    y_piv  = piv.Frames{1,1}.Scales.Y.Slope.*linspace(1, nf_piv(2), nf_piv(2)).*piv.Frames{1,1}.Grids.Y + piv.Frames{1,1}.Scales.Y.Offset;
    [X_piv, Y_piv] = meshgrid(x_piv, y_piv);
    % Align Images: Flip 'n Dip
    raw_image = fliplr(rot90(raw_image, -1));
    UF        = -1 * fliplr(rot90(UF, -1));
    VF        = fliplr(rot90(VF, -1));
    
    
    %%% Resize Raw Image to Fit into PIV Frame
    raw_image = imresize(raw_image, size(UF));
    raw_image(raw_image < background_removal) = 0;
    wave_edge = edge(raw_image, 'Canny', [canny_lower, canny_upper]);
    
    
    %%% Picking top most edge as free surface
    wave_profile = zeros(1, nf_piv(1));
    
    for i = 1:nf_piv(1)
        vertical_slice   = wave_edge(:, i);
        [ones_r, ~] = find(vertical_slice == 1);
        size_ones = size(ones_r);
        len       = size_ones(1);
        if len == 0
            wave_profile(1, i) = nan;
        else
        % If multiple edges are found take value closest to the top
            wave_profile(1, i) = y_piv(min(ones_r));
        end
    end
    
    
    %%% Interpolate over any weird kinks and nan values to have continuous curve
    spikes = find(abs(gradient(wave_profile)) >= grad_tol);
    wave_profile(spikes) = nan;
    
    % Remove values between nans if they are within a specified distance
    for i = 1:length(spikes)
        first_nan = spikes(i);
        for j = 1:nan_dist
            second_nan = first_nan + j;
            if isnan(wave_profile(second_nan))
                wave_profile(first_nan:second_nan) = nan;
            else
                continue
            end
        end
    end
    
    
    % Interpolate over nans
    wave_profile(isnan(wave_profile)) = interp1(x_piv(~isnan(wave_profile)), wave_profile(~isnan(wave_profile)), x_piv(isnan(wave_profile)));
    
    
    % Save Max wave height and profile
    wave_profiles(frame_number,:)  = wave_profile;
    max_wave_heights(frame_number) = max(wave_profile);
    
    % Length of non-nan values in each wave profile
    % Pre-emptive, bad frame detection
    wave_profile_length = sum(~isnan(wave_profile));
    
    frame_UF_mean = mean(UF, 'all', 'omitnan');
    frame_UF_stdv = std(UF, 0, 'all', 'omitnan');

    frame_VF_mean = mean(VF, 'all', 'omitnan');
    frame_VF_stdv = std(VF, 0, 'all', 'omitnan');

    UF(UF < (frame_UF_mean - std_tol * frame_UF_stdv)) = nan;
    UF(UF > (frame_UF_mean + std_tol * frame_UF_stdv)) = nan;

    VF(VF < (frame_VF_mean - std_tol * frame_VF_stdv)) = nan;
    VF(VF > (frame_VF_mean + std_tol * frame_VF_stdv)) = nan;
    
   %%% INPAINT NANS
   UF = inpaint_nans(im2gray(double(UF)));
   VF = inpaint_nans(im2gray(double(VF)));

    UF(Y_piv < wave_profile) = nan;
    VF(Y_piv < wave_profile) = nan;
   
    %%% Plot/Video
    fig = figure(1);
    contourf(X_piv, Y_piv, UF, 200, 'lineStyle', 'none');
    axis equal
    colormap parula
    ax = gca;
    set(ax,'color','w');
    ax.CLim = [0, 2.25];
%     ax.CLim = [-0.5, 0.5];
    title(strcat('Frame Number', {' '}, num2str(frame_number)))
    hold on
    plot(x_piv, wave_profile, 'Color', 'red', 'LineWidth', 2)
    hold off
    xlim([-110,110])
    ylim([-120,110])
    xlabel('x [mm]')
    ylabel('y [mm]')
    colorbar()
    
    drawnow;
    set(fig,'renderer','zbuffer')
    F = getframe(fig);
    tst = [tst F];
    clc;

    
    
    
end

% open the video writer
open(writerObj);
for i=1:length(tst)
    frame = tst(i) ;    
    writeVideo(writerObj, frame);
end
close(writerObj);
%%

pathVideoMP4 = regexprep(video_path,'\.avi','.mp4'); % generate mp4 filename
[~,~] = system(sprintf('ffmpeg -i %s -y -an -c:v libx264 -crf 0 -preset slow %s',video_path,pathVideoMP4)); % for this to work, you should have installed ffmpeg and have it available on PATH

% [~,~] = system(sprintf('ffmpeg -i %s %s',video_path,pathVideoMP4))


%%
% hold on
% for i = 1:D
%     figure(i)
%     plot(x_piv, wave_profiles(i,:))
% end
% hold off
% axis equal
% xlim([-150,150])
% ylim([-150,150])
% % yline(max(max_wave_heights))

%%

%%% Data paths
project_path   = '/Volumes/Zein_PIV/Offshore_Inflow/';
processing     = '/PIV_MP(2x24x24_50%ov)';
recording_name = 'WT4_WV0_AG0';
still_name     = '2_4_2023_waterlevel';
inpt_name      = 'iurfvno';

% Image paths
wave_type      = recording_name(strfind(recording_name, 'WV') + 2);
still_path     = strcat(project_path, '/', still_name, '/Correction');
raw_path       = strcat(project_path, '/', recording_name, '/Correction');
piv_path       = strcat(project_path, '/', recording_name, processing);

% Save paths
results_path   = '/Volumes/Zein_PIV/results/';
mtlb_file      = strcat(results_path, 'data', '/', inpt_name, '_DATA.mat');
mean_file      = strcat(results_path, 'means', '/', inpt_name, '_MEANS.mat');
phase_file     = strcat(results_path, 'phase', '/', inpt_name, '_PHASE.mat');

% Edge detection constants
std_tol            = 2;
background_removal = 1000;
canny_lower        = 0.1;
canny_upper        = 0.4;
grad_tol           = 0.4;
nan_dist           = 5;

% Wave Parameters
left_bound_value   = -115;
right_bound_value  = 115;
phase_tolerance    = 10;
wave_amplitude     = 11.78;
wave_length        = 410;
phase_offset       = [0, wave_length/4, wave_length/2, 3*wave_length/4];
cos_fit            = @(b, v) wave_amplitude * cos(2 * pi * (v - b(1)) / wave_length) + vertical_offset;


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
% constants.vertical_offset    = waterlevel(still_path, piv_path, constants);

background_removal = constants.background_removal;
canny_lower        = constants.canny_lower;
canny_upper        = constants.canny_upper;
grad_tol           = constants.grad_tol;
nan_dist           = constants.nan_dist;
left_bound_value   = constants.left_bound_value;
right_bound_value  = constants.right_bound_value;


% Image Path
image_name         = dir([still_path, '/*.im7']);
piv_image_name     = dir([piv_path, '/*.vc7']);
D                  = length(image_name);
clc;
fprintf('<waterlevel> %.f images used \n', D);

%     % Define Image Dimensions
%     raw                = readimx([still_water_path, '\', image_name(1).name]);
%     raw_image          = raw.Frames{1,1}.Components{1,1}.Planes{1,1};
%     nf                 = size(raw_image);
%     x                  = raw.Frames{1,1}.Scales.X.Slope.*linspace(1, nf(1), nf(1)).* raw.Frames{1,1}.Grids.X + raw.Frames{1,1}.Scales.X.Offset;
%     y                  = raw.Frames{1,1}.Scales.Y.Slope.*linspace(1, nf(2), nf(2)).* raw.Frames{1,1}.Grids.Y + raw.Frames{1,1}.Scales.Y.Offset;
%     [~, left_bound]    = min(abs(x - left_bound_value));
%     [~, right_bound]   = min(abs(x - right_bound_value));
%     x                  = x(left_bound:right_bound);
%     wave_profiles      = zeros(D, length(x));

% Define Image Dimensions Based on PIV Frame
data             = readimx([piv_path, '\', piv_image_name(1).name]);
names            = data.Frames{1,1}.ComponentNames;        
U0_index         = find(strcmp(names, 'U0'));
UF               = data.Frames{1,1}.Components{U0_index,1}.Scale.Slope.*data.Frames{1,1}.Components{U0_index,1}.Planes{1,1} + data.Frames{1,1}.Components{U0_index,1}.Scale.Offset;
nf               = size(UF);

x                = data.Frames{1,1}.Scales.X.Slope.*linspace(1, nf(1), nf(1)).* data.Frames{1,1}.Grids.X + data.Frames{1,1}.Scales.X.Offset;
y                = data.Frames{1,1}.Scales.Y.Slope.*linspace(1, nf(2), nf(2)).* data.Frames{1,1}.Grids.Y + data.Frames{1,1}.Scales.Y.Offset;
[~, left_bound]  = min(abs(x - left_bound_value));
[~, right_bound] = min(abs(x - right_bound_value));
x                = x(left_bound:right_bound);
wave_profiles    = zeros(D, length(x));

figure(1)
hold on
for frame_number = 1:D

    % Load Images
    raw       = readimx([still_path, '\', image_name(frame_number).name]);
    raw_image = raw.Frames{1,1}.Components{1,1}.Planes{1,1};
%             nf        = size(raw_image);
%             x         = raw.Frames{1,1}.Scales.X.Slope.*linspace(1, nf(1), nf(1)).* raw.Frames{1,1}.Grids.X + raw.Frames{1,1}.Scales.X.Offset;
%             y         = raw.Frames{1,1}.Scales.Y.Slope.*linspace(1, nf(2), nf(2)).* raw.Frames{1,1}.Grids.Y + raw.Frames{1,1}.Scales.Y.Offset;
    raw_image = imresize(raw_image, size(UF));
    raw_image = fliplr(rot90(raw_image, -1));

    % Edge Detection
    raw_image(raw_image < background_removal) = 0;
    wave_edge = edge(raw_image, 'Canny', [canny_lower, canny_upper]);
    wave_profile = zeros(1, length(x));

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

%             x            = x(left_bound:right_bound);
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
    % wave_profile(isnan(wave_profile)) = interp1(x(~isnan(wave_profile)), wave_profile(~isnan(wave_profile)), x(isnan(wave_profile)));

    % Plot
    plot(x, wave_profile)

    % Save Profiles
    wave_profiles(frame_number,:) = wave_profile;
end

hold off
axis equal
xlim([min(x), max(x)])
xlabel('x [mm]')
ylabel('z [mm]')
title('Still Surface Profile')

water_line = mean(mean(wave_profiles, 1, 'omitnan'), 'omitnan');

ylim([water_line - 10, water_line + 10])

fprintf('\n<waterlevel> Water Line at %3.2f mm\n', water_line);
































