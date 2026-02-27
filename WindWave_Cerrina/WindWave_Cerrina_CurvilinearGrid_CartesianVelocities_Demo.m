% Demo code opening the curvilinear-projected instantaneous and phase averages for
% Montreal Polytechnique


clc; close all; clear

% Get case to process
wind_speed = 'WT6';
wave_type = 'C';
caze = [wind_speed, '_WV', wave_type, '_AG0'];

matfile_folder = '/Users/zeinsadek/Desktop/Experiments/Offshore/wind_wave_PIV/curvilinear_cerrina';
matfile_name   = [caze, '_CartesianInstantaneous_curvilinearProjected.mat'];

% Load data
output = load(fullfile(matfile_folder, matfile_name));
output = output.output;


%% Breakdown constants

% Wave constants
wavelength = output.constants.wavelength_mm;
amplitude  = output.constants.amplitude_mm;
frequency  = output.constants.frequency_hz;
steepness  = output.constants.steepness;
wave_speed = output.constants.wave_speed;

% Freestream velocity
u_inf = output.constants.freestream;


%% Break down curvililinear grid

curvilinear_grid = output.curvilinear_grid;

% curvilinear_grid(phase).vertical_lines ~ curvilinear 'x'
% curvilinear_grid(phase).horizontal_lines ~ curvilinear 'y'

% curvilinear_grid(phase).wave_profile ~ reference wave profiles
% curvilinear_grid(phase).max_wave_profile ~ max wave profile per phase bin


%% Load Cartesian coordinates

X = output.Cartesian_grid.X;
Y = output.Cartesian_grid.Y;
x = output.Cartesian_grid.x;
y = output.Cartesian_grid.y;


%% Pull out curvilinear means/instantaneous

Cartesian_means = output.Cartesian_means;
Cartesian_instantaneous = output.Cartesian_instantaneous;


%% Plot curvilinear means

% Components to plot
components = {'u', 'v', 'uu', 'vv', 'uv'};

clc; close all
% Loop through components
for c = 1:length(components)
    component = components{c};

    % Figure
    figure('color', 'white')
    tiledlayout(1,4)
    sgtitle(component)

    % Loop through different phases
    for phase = 1:4
        h(phase) = nexttile;
        hold on

        % Plot quantity
        contourf(curvilinear_grid(phase).vertical_lines, ...
                 curvilinear_grid(phase).horizontal_lines, ...
                 Cartesian_means(phase).(component), ...
                 100, 'linestyle', 'none')

        % Plot reference wave profile
        plot(x, curvilinear_grid(phase).wave_profile, 'linewidth', 2, 'color', 'black')

        % Plot max wave profile
        plot(x, curvilinear_grid(phase).max_wave_profile, 'linewidth', 2, 'color', 'red')
        hold off

        title(sprintf('Phase %1.0f', phase))
        colorbar()
        axis equal
        xlim([0, 237])
        ylim([-20, 200])
    end
end

clear phase h c component 


%% Plot instantaneous curvilinear snapshot, and compare the two different crop methods

% The instantaneous images that belong to each phase bin have also been
% projected into wave-tangent (u) and wave-normal components (v). I have
% two different methods for cropping them here, as seen in the means
% contours the misalignment between each snpshot cause a region of
% unresoled data near the surface. 

% The 'instantaneous_crop' field crops each snapshot only to the
% instantaneous wave profile. The 'maximum_crop' removed this unresolved
% region from each snapshpt ahead of time. Any data that lies below the
% reference wave profile is also remvoed from both versions.

phase = 2;
frame = 1;

% Figure
clc; close all
figure('color', 'white')
tiledlayout(2,2)

fprintf('Phase %1.0f: %4.0f Images\n', phase, Cartesian_instantaneous(phase).num_images)

% Instantaneous cropping
% Wave-tangent velocity
h(1) = nexttile;
hold on
contourf(curvilinear_grid(phase).vertical_lines, ...
         curvilinear_grid(phase).horizontal_lines, ...
         Cartesian_instantaneous(phase).instantaneous_crop.u(:,:,frame), ...
         100, 'linestyle', 'none')

% Instantaneous wave profile
plot(x, Cartesian_instantaneous(phase).waves(frame,:), 'linewidth', 2, 'color', 'red')

% Reference wave profile
plot(x, curvilinear_grid(phase).wave_profile, 'linewidth', 2, 'color', 'black')
hold off
axis equal
colorbar()
title(sprintf('$u$ velocity:\nInstantaneous surface cropping'), 'interpreter', 'latex')


% Wave-normal velocity
h(2) = nexttile;
hold on
contourf(curvilinear_grid(phase).vertical_lines, ...
         curvilinear_grid(phase).horizontal_lines, ...
         Cartesian_instantaneous(phase).instantaneous_crop.v(:,:,frame), ...
         100, 'linestyle', 'none')

% Instantaneous wave profile
plot(x, Cartesian_instantaneous(phase).waves(frame,:), 'linewidth', 2, 'color', 'red')

% Reference wave profile
plot(x, curvilinear_grid(phase).wave_profile, 'linewidth', 2, 'color', 'black')
hold off
axis equal
colorbar()
title(sprintf('$v$ velocity:\nInstantaneous surface cropping'), 'interpreter', 'latex')



% Instantaneous cropping
% Wave-tangent velocity
h(3) = nexttile;
hold on
contourf(curvilinear_grid(phase).vertical_lines, ...
         curvilinear_grid(phase).horizontal_lines, ...
         Cartesian_instantaneous(phase).maximum_crop.u(:,:,frame), ...
         100, 'linestyle', 'none')

% Instantaneous wave profile
plot(x, Cartesian_instantaneous(phase).waves(frame,:), 'linewidth', 2, 'color', 'red')

% Reference wave profile
plot(x, curvilinear_grid(phase).wave_profile, 'linewidth', 2, 'color', 'black')

% Max wave profile
plot(x, curvilinear_grid(phase).max_wave_profile, 'linewidth', 2, 'color', 'blue')
hold off
axis equal
colorbar()
title(sprintf('$u$ velocity:\nMaximum surface cropping'), 'interpreter', 'latex')


% Wave-normal velocity
h(4) = nexttile;
hold on
contourf(curvilinear_grid(phase).vertical_lines, ...
         curvilinear_grid(phase).horizontal_lines, ...
         Cartesian_instantaneous(phase).maximum_crop.v(:,:,frame), ...
         100, 'linestyle', 'none')

% Instantaneous wave profile
plot(x, Cartesian_instantaneous(phase).waves(frame,:), 'linewidth', 2, 'color', 'red')

% Reference wave profile
plot(x, curvilinear_grid(phase).wave_profile, 'linewidth', 2, 'color', 'black')

% Max wave profile
plot(x, curvilinear_grid(phase).max_wave_profile, 'linewidth', 2, 'color', 'blue')
hold off
axis equal
colorbar()
title(sprintf('$v$ velocity:\nMaximum surface cropping'), 'interpreter', 'latex')


linkaxes(h, 'xy')
xlim([0, 237])
ylim([-20, 200])

clear h phase frame 
